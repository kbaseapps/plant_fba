import logging
import os
import uuid
import sys
import re
import json

#For Bokeh figures
from bokeh.io import output_file, save
from bokeh.layouts import grid

#For stats
import pandas as pd
import scipy as sp
import scipy.linalg
import scipy.stats
import numpy as np

from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.KBaseReportClient import KBaseReport

from plant_fba.core.fetch_plantseed_impl import FetchPlantSEEDImpl
from plant_fba.core.generate_figure_impl import GenerateFigureImpl

class IntegrateAppImpl:
    @staticmethod
    def _validate_params(params, required, optional=set()):
        """Validates that required parameters are present. Warns if unexpected parameters appear"""
        required = set(required)
        optional = set(optional)
        pkeys = set(params)
        if required - pkeys:
            raise ValueError("Required keys {} not in supplied parameters"
                             .format(", ".join(required - pkeys)))
        defined_param = required | optional
        for param in params:
            if param not in defined_param:
                logging.warning("Unexpected parameter {} supplied".format(param))

    def _build_figure(self,file_path,figure_matrix):
        
        # Make figure matrix html file and embed
        file_name='integrated_scatterplot_output.html'
        figure_html_path = os.path.join(file_path,file_name)
        output_file(figure_html_path)
        save(grid(figure_matrix))
        
        return file_name

    def _build_table(self, table_dict, stats_df):

        html_lines=list()        
        html_lines.append('<table class="table table-bordered table-striped">')

        header_list = ["Enzymes","Compartments","Reactions","EC numbers","Subsystems"]+self.conditions_ids+["Mahalanobis distance","p-value"]

        html_lines.append('<thead>')            
        internal_header_line = "</td><td>".join(header_list)
        html_lines.append('<tr><td>'+internal_header_line+'</td></tr>')
        html_lines.append('</thead>')

        html_lines.append("<tbody>")
        print_row = True
        for complex_row in sorted(table_dict.keys()):
            print_row = True
            cpts = ", ".join(sorted(list(table_dict[complex_row])))

            ecs = []
            subsystems = []
            reactions = []
            conditions = []
            mahal_list=[]
            pvalue_list=[]
            mahalanobis_dist="0.00"
            pvalue="0.00"
            for cpt in table_dict[complex_row]:
                for rxn in table_dict[complex_row][cpt]:
                    
                    if(rxn not in reactions):
                        reactions.append(rxn)

                    if(len(conditions)==0):
                        conditions = table_dict[complex_row][cpt][rxn]

                    if(rxn in self.reactions_data):
                        for ss in self.reactions_data[rxn]['subsystems']:
                            ss=ss.replace("_"," ")
                            ss=ss.replace(" in plants","")
                            if(ss not in subsystems):
                                subsystems.append(ss)

                        for ec in self.reactions_data[rxn]['ecs']:
                            if(ec not in ecs):
                                ecs.append(ec)

                    str_md = "0.00"
                    str_pv = "0.00"
                    if(rxn+'_'+cpt not in stats_df.index):
                        print("MISSING REACTION: ",complex_row,rxn+"_"+cpt)
                        print_row = False
                    else:
                        str_md = "{0:.2f}".format(stats_df.loc[rxn+'_'+cpt]['mahalanobis'])
                        str_pv = "{0:.2f}".format(stats_df.loc[rxn+'_'+cpt]['pvalue'])
                        if(str_pv == "0.00"):
                            str_pv = "{0:.2e}".format(stats_df.loc[rxn+'_'+cpt]['pvalue'])
                        if(mahalanobis_dist != "0.00" and str_md != mahalanobis_dist):
                            print("WARNING: CHANGING STATS FOR SAME PROTEIN COMPLEXES\n")
                            print("===================================================\n\n")
                            print(complex_row,cpts,rxn,conditions,stats_df.loc[rxn+'_'+cpt]['mahalanobis'],mahalanobis_dist,"\n")
                            print("===================================================\n\n")

                    mahalanobis_dist=str_md
                    pvalue=str_pv

            reactions= ", ".join(sorted(reactions))
            subsystems=", ".join(sorted(subsystems))
            ecs=", ".join(sorted(ecs))

            conditions_strings = list()
            for i in range(len(conditions)):
                conditions[i][0] = "{0:.2f}".format(conditions[i][0])
                conditions_strings.append(" | ".join(conditions[i]))

            # some complexes may have zero features predicted
            if(print_row is True):
                html_lines.append("<tr>")
                internal_row_line = "</td><td>".join([complex_row,cpts,reactions,ecs,subsystems]+conditions_strings+[mahalanobis_dist,pvalue])
                html_lines.append("<td>"+internal_row_line+"</td>")
                html_lines.append("</tr>")

        html_lines.append("</tbody>")
        html_lines.append("</table>")

        return "\n".join(html_lines)

    def _build_report(self, figure_matrix, table_dict, stats_df, saved_object_list, workspace_name):
        """
        _generate_report: generate summary report
        """

        # Make report directory and copy over files
        report_file_path = os.path.join(self.scratch, self.report_uuid)
        os.mkdir(report_file_path)

        table_html_string = self._build_table(table_dict, stats_df)

        if(len(self.conditions_ids)>1):
            figure_html_file = self._build_figure(report_file_path,figure_matrix)
            output_html_files = self._generate_report_html(report_file_path,
                                                           figure_html_file=figure_html_file,
                                                           table_string=table_html_string)
        else:
            output_html_files = self._generate_report_html(report_file_path,
                                                           table_string=table_html_string)

        report_params = { 'direct_html_link_index' : 0, #Use to refer to index of 'html_links'
                          'workspace_name' : workspace_name,
                          'report_object_name' : 'plant_fba_' + self.report_uuid,
                          'objects_created' : saved_object_list,
                          'html_links' : output_html_files}

        output = self.kbr.create_extended_report(report_params)

        return {'report_name': output['name'], 'report_ref': output['ref']}

    def _generate_report_html(self, file_path, figure_html_file = None, table_string = None):
        """
            _generate_report: generates the HTML for the upload report
        """
        html_report_list = list()

        ##############################################################
        # Write table html file
        ##############################################################
        # Read in template html
        with open(os.path.join('/kb/module/data','app_report_templates','integrate_abundances_report_tables_template.html')) as report_template_file:
            report_template_string = report_template_file.read()

        # Generate and Insert html title
        title_string = "-".join([self.input_params['input_expression_matrix']]+self.conditions_ids)
        report_template_string = report_template_string.replace('*TITLE*', title_string)

        # Insert html table
        table_report_string = report_template_string.replace('*TABLES*', table_string)

        # Write html file
        table_html_file="integrated_table_output.html"
        with open(os.path.join(file_path,table_html_file),'w') as table_file:
            table_file.write(table_report_string)

        ##############################################################
        # Write summary index.html file
        ##############################################################
        # Begin composing html
        html_lines = list()
        html_lines.append('<h3 style="text-align: center">Integrate Abundances with Metabolism Report</h3>')
        html_lines.append("<p>The \"Integrate Abundances with Metabolism\" app has finished running.</br>")
        html_lines.append("The app integrated the values from the <b>"+self.input_params['input_expression_matrix']+"</b> ExpressionMatrix")
        html_lines.append(" with the <b>"+self.input_params['input_fbamodel']+"</b> FBAModel</br>")
        html_lines.append("Specifically, the app integrated the values from these chosen conditions in the ExpressionMatrix: <b>"+"</b>, <b>".join(self.conditions_ids)+"</b></br>")
        html_lines.append("The results of the integration are stored in the <b>"+self.input_params['output_reaction_matrix']+"</b> ReactionMatrix.</p><br/>")
        html_lines.append('The results of the integration are also tabulated in this <a href="'+table_html_file+'" target="_blank">Table</a></br>')

        if(len(self.conditions_ids)>1):
            html_lines.append('The results of the integration can be also be visualized in these <a href="'+figure_html_file+'" target="_blank">Scatterplots</a>')

        # Read in template html
        with open(os.path.join('/kb/module/data','app_report_templates','integrate_abundances_report_template.html')) as report_template_file:
            report_template_string = report_template_file.read()

        # Insert html
        summary_report_string = report_template_string.replace('*TEXT*', "\n".join(html_lines))

        summary_html_file="index.html"
        with open(os.path.join(file_path,summary_html_file),'w') as index_file:
            index_file.write(summary_report_string)

        ##############################################################
        # Upload files and compose html report object
        ##############################################################
        # Cache it in shock as an archive
        upload_info = self.dfu.file_to_shock({'file_path': file_path,
                                              'pack': 'zip'})

        # HTML Link objects
        html_link = dict()
        # Index
        # html_link = {'shock_id' : upload_info['shock_id'],
        #              'name' : summary_html_file,
        #              'label' : 'HTML report for integrate_abundances_with_metabolism app',
        #              'description' : 'HTML report for integrate_abundances_with_metabolism app'}
        # html_report_list.append(html_link)

        if(len(self.conditions_ids)>1):
            # Figures
            html_link = {'shock_id' : upload_info['shock_id'],
                         'name' : figure_html_file,
                         'label' : 'Scatterplot figures generated by Integrate Abundances with Metabolism app',
                         'description' : 'Scatterplot figures generated by Integrate Abundances with Metabolism app'}
            html_report_list.append(html_link)

        # Table
        html_link = {'shock_id' : upload_info['shock_id'],
                     'name' : table_html_file,
                     'label' : 'HTML table generated by Integrate Abundances with Metabolism app',
                     'description' : 'HTML table generated by Integrate Abundances with Metabolism app'}
        html_report_list.append(html_link)

        print("REPORT:",html_report_list)
        return html_report_list
        
    def _load_fbamodel(self, model_ref):
        
        model_obj = self.dfu.get_objects({'object_refs': [model_ref]})['data'][0]
        print("Number of reactions: "+str(len(model_obj['data']['modelreactions'])))

        model_reaction_lookup_dict = dict()
        for index in range(len(model_obj['data']['modelreactions'])):
            model_reaction_lookup_dict[model_obj['data']['modelreactions'][index]['id']]=index

        return [model_obj,model_reaction_lookup_dict]

    def _load_expression_matrix(self, expdata_ref):

        expdata_obj = self.dfu.get_objects({'object_refs': [expdata_ref]})['data'][0]
        conditions_ids = expdata_obj['data']['data']['col_ids']
        features_ids = expdata_obj['data']['data']['row_ids']

        feature_lookup_dict=dict()
        for index in range(len(features_ids)):
            feature_lookup_dict[features_ids[index]]=index

        condition_lookup_dict=dict()
        for index in range(len(conditions_ids)):
            condition_lookup_dict[conditions_ids[index]]=index

        if(len(self.conditions_ids) == 0):
            self.conditions_ids = conditions_ids

        return [expdata_obj,features_ids,feature_lookup_dict,condition_lookup_dict]
        
    def _compile_genome_scores(self, data, conditions_indices):

        Feature_Comparison_Dict=dict()
        for feature_index in range(len(data)):

            scores_dict=dict()
            for condition in self.conditions_ids:
                condition_index = conditions_indices[condition]

                #Retrieve value from 2D matrix
                score = data[feature_index][condition_index]
        
                #Force into string for easier comparison
                str_score = "{0:.2f}".format(score)

                if(str_score == "0.00"):
                    continue

                scores_dict[condition]=score

            #Here we skip features where there aren't enough scores (should be same number of conditions)
            if(len(scores_dict) < len(self.conditions_ids)):
                continue

            for condition in scores_dict:

                if(condition not in Feature_Comparison_Dict):
                    Feature_Comparison_Dict[condition]=list()
                        
                Feature_Comparison_Dict[condition].append(scores_dict[condition])
            
        return Feature_Comparison_Dict

    def _compile_model_scores_percentiles(self, data):

        # I want to compute percentile rank for each feature under each condition
        # The Conditions_Score_Dicts variable is used to "bin" identical scores
        # (to two decimal points, can be changed)

        # First, we iterate through the conditions for computing percentile rank
        # for each condition
        model_conditions_score_lists=dict()
        model_conditions_score_pct_dicts=dict()
        for condition_index in range(len(self.conditions_ids)):
        
            # For each condition, we "bin" the scores

            score_reaction_dict=dict()
            score_reaction_list=list()
            # The counting of features is done independently because we skip scores of zero
            # (which this affect how percentile rank distributes)
            n_ftrs=0 
            for reaction_index in range(len(data)):
            
                # Retrieve value from 2D matrix
                score = data[reaction_index][condition_index]

                # Many reactions are not assigned a score, and instead have a default tiny score
                if(score == float(-sys.maxsize-1)):
                    continue
            
                # Force into string for easier comparison
                str_score = "{0:.2f}".format(score)

                # I skip the relatively large number of reactions that have a value of zero
                # to prevent the computation of the percentile rank skewing towards zero
                if(str_score == "0.00"):
                    continue
        
                n_ftrs+=1
                if(str_score not in score_reaction_dict):
                    score_reaction_dict[str_score]=list()
                score_reaction_dict[str_score].append(reaction_index)
                score_reaction_list.append(float(str_score))
    
            condition = self.conditions_ids[condition_index]
            model_conditions_score_lists[condition]=score_reaction_list

            # Then for each condition, we use the binned scores to compute
            # percentile rank
            if(condition not in model_conditions_score_pct_dicts):
                model_conditions_score_pct_dicts[condition]=dict()

            sorted_scores = sorted(score_reaction_dict.keys(),key=float)
            less_than_score_ftrs_count=0
            for score_index in range(len(sorted_scores)):

                n_score_ftrs = len(score_reaction_dict[sorted_scores[score_index]])
                half_n_score_ftrs = float(n_score_ftrs) * 0.5
                cumulative_n_score_ftrs = float(less_than_score_ftrs_count) + half_n_score_ftrs
                percentile_rank = cumulative_n_score_ftrs / float(n_ftrs)

                less_than_score_ftrs_count += len(score_reaction_dict[sorted_scores[score_index]])
                model_conditions_score_pct_dicts[condition][sorted_scores[score_index]]=percentile_rank

        # This next part of the code is to re-iterate through the data and to compose the dicts
        # that become ColumnDataStores, and also with default values

        # The reaction_percentile_comparison_dict is for the reaction percentile plot
        reaction_percentile_comparison_dict=dict()
        if('All' not in reaction_percentile_comparison_dict):
            reaction_percentile_comparison_dict['All']=dict()

        # The reaction_score_comparison_dict works for the genome features plot
        reaction_score_comparison_dict=dict()

        for reaction_index in range(len(data)):
            
            scores_dict=dict()
            for condition_index in range(len(self.conditions_ids)):
                
                #Retrieve value from 2D matrix
                score = data[reaction_index][condition_index]
                
                #Many reactions are not assigned a score, and instead a default tiny score
                if(score == float(-sys.maxsize-1)):
                    continue
                
                condition = self.conditions_ids[condition_index]
                scores_dict[condition]=score
                
            #Here we skip reactions where there aren't enough scores (should be same number of conditions)
            if(len(scores_dict) < len(self.conditions_ids)):
                continue
            
            for condition in scores_dict:

                # Collect reaction scores
                if(condition not in reaction_score_comparison_dict):
                    reaction_score_comparison_dict[condition]=list()
                reaction_score_comparison_dict[condition].append(scores_dict[condition])

                # Collect reaction percentiles
                if(condition not in reaction_percentile_comparison_dict['All']):
                    reaction_percentile_comparison_dict['All'][condition]=list()

                #Force into string for easier comparison
                str_score = "{0:.2f}".format(scores_dict[condition])

                #We skip zero scores when computing the percentiles
                #So we have to check for them here
                condition_pct = 0.00
                if(str_score != '0.00'):
                    condition_pct = model_conditions_score_pct_dicts[condition][str_score]
                reaction_percentile_comparison_dict['All'][condition].append(condition_pct)

                if('reactions' not in reaction_percentile_comparison_dict['All']):
                    reaction_percentile_comparison_dict['All']['reactions']=list()
                if(self.reactions_ids[reaction_index] not in \
                       reaction_percentile_comparison_dict['All']['reactions']):
                    reaction_percentile_comparison_dict['All']['reactions'].append(self.reactions_ids[reaction_index])
    
                base_rxn = self.reactions_ids[reaction_index].split('_')[0]
                for ss in self.reactions_data[base_rxn]['subsystems']:
                    if(ss not in reaction_percentile_comparison_dict):
                        reaction_percentile_comparison_dict[ss]=dict()
                    if(condition not in reaction_percentile_comparison_dict[ss]):
                        reaction_percentile_comparison_dict[ss][condition]=list()
                    reaction_percentile_comparison_dict[ss][condition].append(condition_pct)

                    if('reactions' not in reaction_percentile_comparison_dict[ss]):
                        reaction_percentile_comparison_dict[ss]['reactions']=list()
                    if(self.reactions_ids[reaction_index] not in \
                           reaction_percentile_comparison_dict[ss]['reactions']):
                        reaction_percentile_comparison_dict[ss]['reactions'].append(self.reactions_ids[reaction_index])

            self.mh_reactions_ids.append(self.reactions_ids[reaction_index])

        # We set the default values here at the end of the loop because we don't know 
        # how many reactions there will be for each category
        for category in reaction_percentile_comparison_dict:
            for key in ['color','size','tooltip','fill_alpha']:
                reaction_percentile_comparison_dict[category][key]=list()

            for index in range(len(reaction_percentile_comparison_dict[category][self.conditions_ids[0]])):

                reaction_percentile_comparison_dict[category]['fill_alpha'].append(1.0)
                
                # format string of subsystems for tooltip
                rxn = reaction_percentile_comparison_dict[category]['reactions'][index]
                base_rxn = rxn.split('_')[0]
                ss_string = ", ".join(self.reactions_data[base_rxn]['subsystems'])
                reaction_percentile_comparison_dict[category]['tooltip'].append(rxn+", "+ss_string)

                if(category == 'All'):

                    reaction_percentile_comparison_dict[category]['color'].append('black')
                    reaction_percentile_comparison_dict[category]['size'].append(6)

                else:

                    reaction_percentile_comparison_dict[category]['color'].append('red')
                    reaction_percentile_comparison_dict[category]['size'].append(8)
                    
        return [reaction_score_comparison_dict, reaction_percentile_comparison_dict]

    def _compile_mahalanobis_dist_pvalue(self, data, threshold):

        data_df = pd.DataFrame(data,columns=self.conditions_ids,index=self.mh_reactions_ids)

        # I don't know the math well enough to follow what's going on, but I used 
        # the recipe described here:
        # https://www.machinelearningplus.com/statistics/mahalanobis-distance/

        # Covariance matrix via numpy
        cov_mat = np.cov(data_df.values.T)

        # Inverse covariance matrix via scipy.linalg
        # It won't accept a 1x1 matrix hence the if/else
        if(len(self.conditions_ids)>1):
            inv_cov_mat = sp.linalg.inv(cov_mat)
        else:
            inv_cov_mat = 1/cov_mat

        # two terms required, second using dot product
        data_minus_mean = data_df - np.mean(data_df)
        left_term = np.dot(data_minus_mean,inv_cov_mat)

        # dot product
        mahalanobis = np.dot(left_term, data_minus_mean.T)
        data_df['mahalanobis'] = mahalanobis.diagonal()

        # chi-squared p-values with one degree of freedom (two sets of variables)
        data_df['pvalue'] = 1-sp.stats.chi2.cdf(data_df['mahalanobis'], 1)

        # find the outliers below a given threshold, i.e. p < 0.01
        outliers=data_df.loc[data_df.pvalue < threshold]
        # this is used when you want to just plot the p-values alone
        data_df.index.name = 'reactions'
        outliers.index.name = 'reactions'

        #Need to return the mapping between reactions and the p-values
        return [data_df, outliers]

    def _integrate_abundances(self, model_obj, feature_lookup_dict, expdata_obj):

        reaction_values_matrix=list()
        reactions_ids=list()
        minmax_expscore_dict=dict()
        model_complexes_dict=dict()
        fh = open(self.scratch+'/output.txt','w')
        fh2 = open(self.scratch+'/rxn01486.txt','w')
        print_data=False
        for mdlrxn in range(len(model_obj['data']['modelreactions'])):
            mdlrxn_obj=model_obj['data']['modelreactions'][mdlrxn]
            reactions_ids.append(mdlrxn_obj['id'])
            [base_rxn,cpt_id]=mdlrxn_obj['id'].split('_')

#            if(base_rxn == 'rxn01486' or base_rxn == 'rxn37610'):
#                print_data=True

            rxndata_row=list()
            for experiment in range(len(self.conditions_ids)):
                if(self.conditions_ids[experiment] not in minmax_expscore_dict):
                    minmax_expscore_dict[self.conditions_ids[experiment]]={'max':-sys.maxsize-1,'min':sys.maxsize}

                # Maximal gene expression for a reaction
                reaction_score=['nan',""]
                prots_str_list=list()
                for prt in mdlrxn_obj['modelReactionProteins']:

                    # Minimal gene expression for a complex
                    complex_score=['nan',""]
                    subs_str_list=list()
                    for sbnt in prt['modelReactionProteinSubunits']:

                        # Maximal gene expression for a subunit
                        subunit_score=['nan',""]
                        ftrs_str_list=list()
                        for feature in sbnt['feature_refs']:
                            feature=feature.split('/')[-1]
                            ftrs_str_list.append(feature)

                            ftr_score = expdata_obj['data']['data']['values'][feature_lookup_dict[feature]][experiment]

                            if(print_data is True):
                                fh2.write(mdlrxn_obj['id']+':'+feature+':'+str(ftr_score)+'\n')

                            if(ftr_score < minmax_expscore_dict[self.conditions_ids[experiment]]['min']):
                                minmax_expscore_dict[self.conditions_ids[experiment]]['min'] = ftr_score

                            if(ftr_score > minmax_expscore_dict[self.conditions_ids[experiment]]['max']):
                                minmax_expscore_dict[self.conditions_ids[experiment]]['max'] = ftr_score

                            # Maximal gene expression for a subunit
                            if(subunit_score[0] == 'nan' or subunit_score[0] < ftr_score):
                                subunit_score = [ftr_score, feature]
                        
                        if(print_data is True):
                            fh2.write(subunit_score,'\n')

                        ftr_str = "("+", ".join(ftrs_str_list)+")"
                        subs_str_list.append(ftr_str)

                        # Minimal gene expression for a complex
                        if(subunit_score[0] != 'nan'):
                            if(complex_score[0] == 'nan' or complex_score[0] > subunit_score[0]):
                                complex_score[0] = subunit_score[0]
                                complex_score[1] = subunit_score[1]

                    if(print_data is True):
                        fh2.write(complex_score,'\n')

                    sub_str = "["+", ".join(subs_str_list)+"]"
                    prots_str_list.append(sub_str)

                    # Maximal gene expression for a reaction
                    if(complex_score[0] != 'nan'):
                        if(reaction_score[0] == 'nan' or reaction_score[0] < complex_score[0]):
                            reaction_score[0] = complex_score[0]
                            reaction_score[1] = complex_score[1]
            
                if(reaction_score[0] == 'nan'):
                    reaction_score[0] = float(-sys.maxsize-1)

                if(print_data is True):
                    fh2.write(self.conditions_ids[experiment]+':'+str(reaction_score[0])+'('+reaction_score[1]+')\n')

                #Putting together dict for table
                proteins_string=', '.join(prots_str_list)
                if(len(prots_str_list)>0 and proteins_string != "[]" and proteins_string != "[()]"):
                    if(proteins_string not in model_complexes_dict):
                        model_complexes_dict[proteins_string]=dict()
                    if(cpt_id not in model_complexes_dict[proteins_string]):
                        model_complexes_dict[proteins_string][cpt_id]=dict()
                    if(base_rxn not in model_complexes_dict[proteins_string][cpt_id]):
                        model_complexes_dict[proteins_string][cpt_id][base_rxn]=list()
                    fh.write('\t'.join([self.conditions_ids[experiment],proteins_string,cpt_id,base_rxn,str(reaction_score[0]),reaction_score[1],'\n']))
                    model_complexes_dict[proteins_string][cpt_id][base_rxn].append(reaction_score)

                rxndata_row.append(reaction_score[0])

            print_data=False

            reaction_values_matrix.append(rxndata_row)

        fh.close()

        self.reactions_ids=reactions_ids
        return (reaction_values_matrix, model_complexes_dict)

    def __init__(self, config, ctx, input_params):
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        self.dfu = DataFileUtil(self.callback_url)
        self.kbr = KBaseReport(self.callback_url)

        self.scratch = config['scratch']
        self.report_uuid = str(uuid.uuid4())
        
        # There is a bug in the UI that won't let me collect a
        # a clean list of conditions, so I have to parse them
        # from a comma-separated string
        if("input_columns" in input_params and input_params["input_columns"] != ""):
            conditions = list()
            for condition in input_params["input_columns"].split(','):
                conditions.append(condition)
            input_params["input_columns"]=conditions

        self.input_params = input_params

        # set in _load_expression_matrix()
        self.conditions_ids = list()
        
        # this is an optional parameter, but restricts the 
        # number of chosen columns in the matrix
        if('input_columns' in input_params and len(input_params['input_columns'])>0):
            self.conditions_ids = input_params['input_columns']

        # set in _integrate_abundances()
        self.reactions_ids = list()

        # set in _compile_model_scores_percentiles
        self.mh_reactions_ids = list()

        with open(os.path.join("/kb/module/PlantSEED",
                               "Data/PlantSEED_v3",
                               "PlantSEED_Roles.json")) as plsd_fh:
            PS_Roles = json.load(plsd_fh)

        plantseed = FetchPlantSEEDImpl()
        self.reactions_data = plantseed.fetch_reactions(PS_Roles)

    def integrate_abundances_with_metabolism(self):

        self._validate_params(self.input_params, 
                              {'input_ws',
                               'input_fbamodel',
                               'input_expression_matrix',
                               'output_reaction_matrix'}, 
                              {'input_columns'})

        ##############################################################
        # Load model and expression objects
        ##############################################################
        model_ref = self.input_params['input_ws']+'/'+self.input_params['input_fbamodel']
        [model_obj,reaction_index] = self._load_fbamodel(model_ref)

        # The columns / conditions_ids are set in this function if not set via user parameter
        expression_ref = self.input_params['input_ws']+'/'+self.input_params['input_expression_matrix']
        [expdata_obj,features_ids,feature_index,condition_index] = self._load_expression_matrix(expression_ref)

        ##############################################################
        # Extract expression abundances for use in first scatter plot
        ##############################################################
        feature_comparison_dict = self._compile_genome_scores(expdata_obj['data']['data']['values'], condition_index)
 
        ####################################################################
        # Actually integrate abundances and build new ReactionMatrix object
        ####################################################################
        (reaction_values_matrix, model_complexes_dict) = self._integrate_abundances(model_obj, feature_index, expdata_obj)
        
        rxndata_obj = {'row_ids':self.reactions_ids,
                       'col_ids':self.conditions_ids,
                       'values':reaction_values_matrix}

        ##########################################################################################
        # Extract / organize reaction expression scores for use in first and second scatter plot
        ##########################################################################################
        [reaction_scores_dict, reaction_percentiles_dict] = self._compile_model_scores_percentiles(reaction_values_matrix)

        #############################################################################################################
        # Multi-variate mahalanobis distances computed along with outliers depending on chi-squared p-value of 0.01
        #############################################################################################################
        [mahal_dist_df,outliers] = self._compile_mahalanobis_dist_pvalue(reaction_percentiles_dict['All'],0.01)

        ##############################################################
        # Figure generator
        ##############################################################
        subsystem_select_list=["None"]
        for category in sorted(list(reaction_percentiles_dict.keys())):
            if(category == 'All'):
                continue
            subsystem_select_list.append(category)

            for rxn_idx in range(len(reaction_percentiles_dict[category]['reactions'])):
                rxn = reaction_percentiles_dict[category]['reactions'][rxn_idx]
                pval = mahal_dist_df.loc[rxn]['pvalue']
                # reaction_percentiles_dict[category]['fill_alpha'][rxn_idx] = 1-pval

        figure_generator = GenerateFigureImpl()
        figure_grid = figure_generator.generate_figure(self.conditions_ids, category_select=subsystem_select_list,
                                                       genome_features=feature_comparison_dict,
                                                       reaction_scores=reaction_scores_dict,
                                                       reaction_percentiles=reaction_percentiles_dict)

        ##############################################################
        # Finishing and Saving ReactionMatrix
        ##############################################################
        ReactionMatrix_obj = {'type':'KBaseMatrices.ReactionMatrix','name':self.input_params['output_reaction_matrix'],
                              'data':{'scale':'raw','description':'reaction expression score',
                                      'fbamodel_ref':model_ref,
                                      'expression_ref':expression_ref,
                                      'data':rxndata_obj}}

        ws_id = self.dfu.ws_name_to_id(self.input_params['input_ws'])
        saved_matrix_dict = self.dfu.save_objects({'id':ws_id,'objects':[ReactionMatrix_obj]})[0]
        saved_matrix_ref = "{}/{}/{}".format(saved_matrix_dict[6],saved_matrix_dict[0],saved_matrix_dict[4])
        saved_matrix_desc = "Reaction matrix: "+self.input_params['output_reaction_matrix']

        #####################################################################
        # Building the report with figures, tables, and saved_objects (to be improved)
        # We pass in a dict where each key is a row for the table
        #####################################################################

        output_object_files = list()
        output_object_files.append({'ref':saved_matrix_ref,'description':saved_matrix_desc})

        return self._build_report(figure_grid, model_complexes_dict, 
                                  mahal_dist_df, output_object_files, 
                                  self.input_params['input_ws'])
