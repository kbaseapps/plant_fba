import logging
import os
import uuid
import sys
import re

#For Bokeh figures
from bokeh.io import output_file, save, export_png
from bokeh.plotting import figure
from bokeh.layouts import grid
from bokeh.models import Slope

#For stats
import pandas as pd
import scipy as sp
import scipy.linalg
import scipy.stats
import numpy as np

from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.KBaseReportClient import KBaseReport

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

        [rxns_subsystems,rxns_roles,rxns_ecs] = self._load_subsystems()

        html_lines=list()        
        html_lines.append('<table class="table table-bordered table-striped">')

        header_list = ["Complex","Compartments","Reactions","EC numbers","Subsystems"]+self.conditions_ids+["Mahalanobis distance","p-value"]

        html_lines.append('<thead>')            
        internal_header_line = "</td><td>".join(header_list)
        html_lines.append('<tr><td>'+internal_header_line+'</td></tr>')
        html_lines.append('</thead>')

        html_lines.append("<tbody>")
        for complex_row in sorted(table_dict.keys()):
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

                    if(rxn in rxns_subsystems):
                        for ss in rxns_subsystems[rxn]:
                            if(ss not in subsystems):
                                subsystems.append(ss)

                    if(rxn in rxns_ecs):
                        for ec in rxns_ecs[rxn]:
                            if(ec not in ecs):
                                ecs.append(ec)

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

            for i in range(len(conditions)):
                conditions[i] = "{0:.2f}".format(conditions[i])

            html_lines.append("<tr>")
            internal_row_line = "</td><td>".join([complex_row,cpts,reactions,ecs,subsystems]+conditions+[mahalanobis_dist,pvalue])
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

        figure_html_file = self._build_figure(report_file_path,figure_matrix)
        table_html_string = self._build_table(table_dict, stats_df)

        output_html_files = self._generate_report_html(report_file_path,figure_html_file,table_html_string)

        report_params = { 'direct_html_link_index' : 0, #Use to refer to index of 'html_links'
                          'workspace_name' : workspace_name,
                          'report_object_name' : 'plant_fba_' + self.report_uuid,
                          'objects_created' : saved_object_list,
                          'html_links' : output_html_files}

        output = self.kbr.create_extended_report(report_params)

        return {'report_name': output['name'], 'report_ref': output['ref']}

    def _generate_report_html(self, file_path, figure_html_file, table_html_string):
        """
            _generate_report: generates the HTML for the upload report
        """
        html_report_list = list()

        ##############################################################
        # Write figures html file (redundant)
        ##############################################################
        #with open(os.path.join('/kb/module/data','app_report_templates','integrate_abundances_report_figures_template.html')) as report_template_file:
        #    report_template_string = report_template_file.read()

        # embed html of figures
        # figure_html_string="<iframe src=\""+figure_html_file+"\" width='auto' height='auto'</iframe>" #width=\"100%\" height=\"100%\"></iframe>"
        # index_report_string = report_template_string.replace('*FIGURES*',figure_html_string)

        #Write html file
        #with open(os.path.join(file_path,"index.html"),'w') as index_file:
        #    index_file.write(index_report_string)

        ##############################################################
        # Write table html file
        ##############################################################
        # Read in template html
        with open(os.path.join('/kb/module/data','app_report_templates','integrate_abundances_report_tables_template.html')) as report_template_file:
            report_template_string = report_template_file.read()

        # Insert html table
        table_report_string = report_template_string.replace('*TABLES*', table_html_string)

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
        html_lines.append("<p>The \"Integrate Abundances with Metabolism\" app has finished running:</br>")
        html_lines.append("The app integrated the gene abundances from the "+self.input_params['input_expression_matrix']+" ExpressionMatrix with the ")
        html_lines.append(self.input_params['input_fbamodel']+" FBAModel, resulting in the "+self.input_params['output_reaction_matrix']+" ReactionMatrix.</p><br/>")

        html_lines.append('<p>The data within the output '+self.input_params['output_reaction_matrix']+' ReactionMatrix is also presented in two forms:</p>')
        html_lines.append('<list><ul><a href="'+table_html_file+'" target="_blank">Table</a></ul></list>')
        html_lines.append('<list><ul><a href="'+figure_html_file+'" target="_blank">Scatterplots</a></ul></list>')

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
    
    def _build_scatterplot(self,data,condition_1,condition_2,title="Default Title"):

        od_fig = figure()
        od_fig.xaxis.axis_label = self.conditions_ids[condition_1]
        od_fig.yaxis.axis_label = self.conditions_ids[condition_2]
        
        df = pd.DataFrame(data,columns=self.conditions_ids)
        od_fig.title.text = title
        od_fig.circle(x=self.conditions_ids[condition_1], y=self.conditions_ids[condition_2], source=df, color="black")

        slope_line = Slope(gradient=1, y_intercept=0, line_color="red")
        od_fig.add_layout(slope_line)

        return od_fig

    def _add_to_scatterplot(self,figure,data,condition_1,condition_2,color="black"):

        df = pd.DataFrame(data,columns=self.conditions_ids)
        figure.circle(x=self.conditions_ids[condition_1], y=self.conditions_ids[condition_2], source=df, fill_color=color, line_color="black", size=6)
        
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

        self.conditions_ids = conditions_ids
        return [expdata_obj,features_ids,feature_lookup_dict]
        
    def _compile_genome_scores(self, data):

        Feature_Comparison_Dict=dict()
        for feature_index in range(len(data)):

            scores_dict=dict()
            for condition_index in range(len(self.conditions_ids)):

                #Retrieve value from 2D matrix
                score = data[feature_index][condition_index]
        
                #Force into string for easier comparison
                str_score = "{0:.2f}".format(score)

                if(str_score == "0.00"):
                    continue

                condition = self.conditions_ids[condition_index]
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

        #I want to compute percentile rank for each feature under each condition
        #The Conditions_Score_Dicts variable is used to "bin" identical scores
        #(to two decimal points, can be changed)
        #For computing percentile rank
        
        model_conditions_score_dicts=dict()
        model_conditions_score_lists=dict()
        for condition_index in range(len(self.conditions_ids)):
        
            score_reaction_dict=dict()
            score_reaction_list=list()
            n_ftrs=0 #This is accumulated independently because we skip scores of zero (as this affects how percentile rank distributes)
            for reaction_index in range(len(data)):
            
                #Retrieve value from 2D matrix
                score = data[reaction_index][condition_index]

                #Many reactions are not assigned a score, and instead a default tiny score
                if(score == float(-sys.maxsize-1)):
                    continue
            
                #Force into string for easier comparison
                str_score = "{0:.2f}".format(score)

                if(str_score == "0.00"):
                    continue
        
                n_ftrs+=1
                if(str_score not in score_reaction_dict):
                    score_reaction_dict[str_score]=list()
                score_reaction_dict[str_score].append(reaction_index)
                score_reaction_list.append(float(str_score))
    
            condition = self.conditions_ids[condition_index]
            model_conditions_score_lists[condition]=score_reaction_list
            if(condition not in model_conditions_score_dicts):
                model_conditions_score_dicts[condition]=dict()

            #Start calculating percentile rank
            sorted_scores = sorted(score_reaction_dict.keys(),key=float)
            less_than_score_ftrs_count=0
            for score_index in range(len(sorted_scores)):

                n_score_ftrs = len(score_reaction_dict[sorted_scores[score_index]])
                half_n_score_ftrs = float(n_score_ftrs) * 0.5
                cumulative_n_score_ftrs = float(less_than_score_ftrs_count) + half_n_score_ftrs
                percentile_rank = cumulative_n_score_ftrs / float(n_ftrs)

                less_than_score_ftrs_count += len(score_reaction_dict[sorted_scores[score_index]])
                model_conditions_score_dicts[condition][sorted_scores[score_index]]=percentile_rank

        reaction_score_comparison_dict=dict()
        reaction_percentile_comparison_dict=dict()
        reactions_rows=list()
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
                
                if(condition not in reaction_percentile_comparison_dict):
                    reaction_percentile_comparison_dict[condition]=list()

                #Force into string for easier comparison
                str_score = "{0:.2f}".format(scores_dict[condition])

                #We skip zero scores when computing the percentiles
                #So we have to check for them here
                condition_pct = 0.00
                if(str_score != '0.00'):
                    condition_pct = model_conditions_score_dicts[condition][str_score]
                reaction_percentile_comparison_dict[condition].append(condition_pct)
    
            reactions_rows.append(self.reactions_ids[reaction_index])

        reaction_score_comparison_df = pd.DataFrame(reaction_score_comparison_dict,columns=self.conditions_ids,index=reactions_rows)
        reaction_percentile_comparison_df = pd.DataFrame(reaction_percentile_comparison_dict,columns=self.conditions_ids,index=reactions_rows)

        return [reaction_score_comparison_df, reaction_percentile_comparison_df]

    def _compile_mahalanobis_dist_pvalue(self, data, threshold):

        # I don't know the math well enough to follow what's going on, but I used 
        # the recipe described here:
        # https://www.machinelearningplus.com/statistics/mahalanobis-distance/

        # Covariance matrix via numpy
        cov_mat = np.cov(data.values.T)

        # Inverse covariance matrix via scipy.linalg
        inv_cov_mat = sp.linalg.inv(cov_mat)
        
        # two terms required, second using dot product
        data_minus_mean = data - np.mean(data)
        left_term = np.dot(data_minus_mean,inv_cov_mat)

        # dot product
        mahalanobis = np.dot(left_term, data_minus_mean.T)
        data['mahalanobis'] = mahalanobis.diagonal()

        # chi-squared p-values with one degree of freedom (two sets of variables)
        data['pvalue'] = 1-sp.stats.chi2.cdf(data['mahalanobis'], 1)

        # find the outliers below a given threshold, i.e. p < 0.01
        outliers=data.loc[data.pvalue < threshold]
        # this is used when you want to just plot the p-values alone
        data.index.name = 'reactions'
        outliers.index.name = 'reactions'

        #Need to return the mapping between reactions and the p-values
        return [data, outliers]
        
    def _load_subsystems(self):

        # Collect Core Subsystems and Reactions (From PlantSEED Reference Genome)
        Core_SS_Classes=('Central Carbon','Amino acids','Nucleic acids')
        Reactions_Subsystems=dict()
        Reactions_Roles=dict()

        genome_ref = 'PlantSEED_v2/PlantSEED_Arabidopsis'
        genome_obj = self.dfu.get_objects({'object_refs': [genome_ref]})['data'][0]
        for ftr in genome_obj['data']['features']:
            if(len(ftr['functional_descriptions'])>0):
                for fd in ftr['functional_descriptions']:
                    (mclass,subsystem,reaction_str)=fd.split('::')
                    subsystem=subsystem.replace(" in plants","")
                    # if(mclass not in Core_SS_Classes):
                    #    continue
                    for rxn in reaction_str.split('|'):
                        if(rxn not in Reactions_Subsystems):
                            Reactions_Subsystems[rxn]=dict()
                        Reactions_Subsystems[rxn][subsystem]=1

                        if(rxn not in Reactions_Roles):
                            Reactions_Roles[rxn]=dict()

                        #Split out comments
                        Function_Comments = ftr['functions'][0].split("#")
                        for i in range(len(Function_Comments)):
                            Function_Comments[i]=Function_Comments[i].strip()

                        Function = Function_Comments.pop(0)
                        roles = re.split("\s*;\s+|\s+[\@\/]\s+", Function)
                        for role in roles:
                            Reactions_Roles[rxn][role]=1

        Reactions_ECs=dict()
        for rxn in Reactions_Roles:
            roles = list(Reactions_Roles[rxn].keys())
            for role in roles:
                if('EC' in role):
                    match = re.search(r"\d+\.[\d-]+\.[\d-]+\.[\d-]+", role)
                    if(match is not None):
                        if(rxn not in Reactions_ECs):
                            Reactions_ECs[rxn]=dict()
                        Reactions_ECs[rxn][match.group(0)]=1

        print("Collected "+str(len(Reactions_Subsystems))+" core reactions")
        return [Reactions_Subsystems,Reactions_Roles,Reactions_ECs]

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
                reaction_score='nan'
                prt_list=list()
                for prt in mdlrxn_obj['modelReactionProteins']:

                    # Minimal gene expression for a complex
                    complex_score='nan'
                    prt_sbnt_list=list()
                    for sbnt in prt['modelReactionProteinSubunits']:

                        # Maximal gene expression for a subunit
                        subunit_score='nan'
                        sbnt_ftr_list=list()
                        for feature in sbnt['feature_refs']:
                            feature=feature.split('/')[-1]
                            sbnt_ftr_list.append(feature)

                            ftr_score = expdata_obj['data']['data']['values'][feature_lookup_dict[feature]][experiment]

                            if(print_data is True):
                                fh2.write(mdlrxn_obj['id']+':'+feature+':'+str(ftr_score)+'\n')

                            if(ftr_score < minmax_expscore_dict[self.conditions_ids[experiment]]['min']):
                                minmax_expscore_dict[self.conditions_ids[experiment]]['min'] = ftr_score

                            if(ftr_score > minmax_expscore_dict[self.conditions_ids[experiment]]['max']):
                                minmax_expscore_dict[self.conditions_ids[experiment]]['max'] = ftr_score

                            # Maximal gene expression for a subunit
                            if(subunit_score == 'nan' or subunit_score < ftr_score):
                                subunit_score = ftr_score
                        
                        if(print_data is True):
                            fh2.write(str(subunit_score)+'\n')

                        if(len(sbnt_ftr_list)>0):
                            prt_sbnt_list.append(', '.join(sbnt_ftr_list))

                        # Minimal gene expression for a complex
                        if(subunit_score != 'nan'):
                            if(complex_score == 'nan' or complex_score > subunit_score):
                                complex_score = subunit_score

                    if(print_data is True):
                        fh2.write(str(complex_score)+'\n')

                    if(len(prt_sbnt_list)>0):
                        prt_list.append('; '.join(prt_sbnt_list))

                    # Maximal gene expression for a reaction
                    if(complex_score != 'nan'):
                        if(reaction_score == 'nan' or reaction_score < complex_score):
                            reaction_score = complex_score
            
                if(reaction_score == 'nan'):
                    reaction_score = float(-sys.maxsize-1)

                if(print_data is True):
                    fh2.write(self.conditions_ids[experiment]+':'+str(reaction_score)+'\n')

                #Putting together dict for table
                if(len(prt_list)>0):
                    proteins_string=' | '.join(prt_list)
                    if(proteins_string not in model_complexes_dict):
                        model_complexes_dict[proteins_string]=dict()
                    if(cpt_id not in model_complexes_dict[proteins_string]):
                        model_complexes_dict[proteins_string][cpt_id]=dict()
                    if(base_rxn not in model_complexes_dict[proteins_string][cpt_id]):
                        model_complexes_dict[proteins_string][cpt_id][base_rxn]=list()
                    fh.write('\t'.join([self.conditions_ids[experiment],proteins_string,cpt_id,base_rxn,str(reaction_score),'\n']))
                    model_complexes_dict[proteins_string][cpt_id][base_rxn].append(reaction_score)

                rxndata_row.append(reaction_score)

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
        self.input_params = input_params
        
        # set in _load_expression_matrix()
        self.conditions_ids = list()
        # set in _integrate_abundances()
        self.reactions_ids = list()


    def integrate_abundances_with_metabolism(self):

        self._validate_params(self.input_params, {'input_ws',
                                                  'input_fbamodel',
                                                  'input_expression_matrix',
                                                  'output_reaction_matrix'}, {})

        ##############################################################
        # Load model and expression objects
        ##############################################################
        model_ref = self.input_params['input_ws']+'/'+self.input_params['input_fbamodel']
        [model_obj,reaction_index] = self._load_fbamodel(model_ref)

        expression_ref = self.input_params['input_ws']+'/'+self.input_params['input_expression_matrix']
        [expdata_obj,features_ids,feature_index] = self._load_expression_matrix(expression_ref)

        # Matrix of figures to be saved under one `grid` command
        figure_matrix= list()

        ##############################################################
        # Extract expression abundances for use in first scatter plot
        ##############################################################
        feature_comparison_dict = self._compile_genome_scores(expdata_obj['data']['data']['values'])
 

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
        [reaction_score_comparison_df, reaction_percentile_comparison_df] = self._compile_model_scores_percentiles(reaction_values_matrix)

        #############################################################################################################
        # Multi-variate mahalanobis distances computed along with outliers depending on chi-squared p-value of 0.01
        #############################################################################################################
        [mahal_dist_df,outliers] = self._compile_mahalanobis_dist_pvalue(reaction_percentile_comparison_df,0.01)

        ##############################################################
        # Iterate through every unique pair of experiments/columns
        ##############################################################
        for first_condition_index in range(len(self.conditions_ids)):
            for second_condition_index in range(len(self.conditions_ids)):
                if(first_condition_index <= second_condition_index):
                    continue
  
                # Row of figures: at this point, its just two, each with two datasets
                figure_array = list()

                # Raw transcript abundance for genome in first figure of array
                feature_figure = self._build_scatterplot(feature_comparison_dict,
                                                         first_condition_index, second_condition_index,
                                                         title="Genome Features Expression Abundances")

                # Raw reaction expression scores mapped over genome abundances in first figure
                self._add_to_scatterplot(feature_figure, reaction_score_comparison_df, 
                                         first_condition_index, second_condition_index,
                                         color="lightgreen")

                # Add first figure of row
                figure_array.append(feature_figure)

                # Second scatterplot built for normalized expression scores
                reaction_figure = self._build_scatterplot(reaction_percentile_comparison_df,
                                                          first_condition_index, second_condition_index,
                                                          title="Model Reactions Percentile Rank (p<0.01)")

                # Add outliers (data matches, but its a different color)
                self._add_to_scatterplot(reaction_figure, outliers, 
                                         first_condition_index, second_condition_index,
                                         color="red")

                # Add second figure of row
                figure_array.append(reaction_figure)

                # Add row to matrix of figures
                figure_matrix.append(figure_array)

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

        ###############################################################
        # Compiling the list of complexes for use in the report table
        ###############################################################
        # complexes_dict = self._compile_model_complexes(model_obj)

        #####################################################################
        # Building the report with figures, tables, and saved_objects (to be improved)
        # We pass in a dict where each key is a row for the table
        #####################################################################

        output_object_files = list()
        output_object_files.append({'ref':saved_matrix_ref,'description':saved_matrix_desc})

        return self._build_report(figure_matrix, model_complexes_dict, mahal_dist_df, output_object_files, self.input_params['input_ws'])
