import logging
import os
import uuid
import sys

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

    def _build_report(self, figure_matrix, workspace_name, saved_obj_ref, saved_obj_desc):
        """
        _generate_report: generate summary report with counts
        """
        output_html_files = self._generate_report_html(figure_matrix)

        output_object_files = list()
        output_object_files.append({'ref':saved_obj_ref,'description':saved_obj_desc})

        report_params = { 'direct_html_link_index' : 0, #Use to refer to index of 'html_links'
                          'workspace_name' : workspace_name,
                          'report_object_name' : 'plant_fba_' + self.report_uuid,
                          'objects_created' : output_object_files,
                          'html_links' : output_html_files}

        output = self.kbr.create_extended_report(report_params)

        return {'report_name': output['name'], 'report_ref': output['ref']}

    def _generate_report_html(self, figure_matrix):
        """
            _generate_report: generates the HTML for the upload report
        """
        html_report_list = list()

        # Make report directory and copy over files
        report_file_path = os.path.join(self.scratch, self.report_uuid)
        os.mkdir(report_file_path)

        #Begin composing html
        html_lines = list()

        # Make figure matrix html file and embed
        file_name='integrated_scatterplot_output.html'
        figure_html_path = os.path.join(report_file_path,file_name)
        output_file(figure_html_path)
        save(grid(figure_matrix))

        # build html of figures
        html_lines.append("<iframe src=\""+file_name+"\" width=\"50%\" height=\"50%\"></iframe>")

        # Build HTML tables for results
        #table_lines = []
        #table_lines.append(f'<h3 style="text-align: center">Object Counts</h3>')
        #table_lines.append('<table class="table table-bordered table-striped">')
        #header = "</td><td>".join(['Amplicon', 'Taxon', 'Total'] + self.object_categories)
        #table_lines.append(f'\t<thead><tr><td>{header}</td></tr></thead>')
        #table_lines.append('\t<tbody>')
        #for taxon in taxon_list:
        #    tax_counts = [object_counts[taxon['id']].get(ws_type, 0)
        #                  for ws_type in self.object_categories]
        #    row = [taxon['id'], taxon['name'], str(sum(tax_counts))]
        #    row += [str(x) for x in tax_counts]
        #    line = "</td><td>".join(row)
        #    table_lines.append(f'\t\t<tr><td>{line}</td></tr>')
        #table_lines.append('\t</tbody>')
        #table_lines.append('</table>\n')

        # Fill in template HTML
        #with open(os.path.join(os.path.dirname(__file__), 'table_report_template.html')
        #          ) as report_template_file:
        #    report_template = report_template_file.read() \
        #        .replace('*TABLES*', "\n".join(table_lines))


        #with open(os.path.join('/kb/module/data','app_report_templates','integrate_abundances_report_template.html')) as report_template_file:
        #report_string = report_template_file.read().replace('*TABLES*', html_string)

        report_string = "\n".join(html_lines)

        #Write html file
        with open(os.path.join(report_file_path,"index.html"),'w') as index_file:
            index_file.write(report_string)

        #Cache it in shock as an archive
        upload_info = self.dfu.file_to_shock({'file_path': report_file_path,
                                              'pack': 'zip'})

        #Html Link object
        html_link = {'shock_id' : upload_info['shock_id'],
                     'name' : 'index.html',
                     'label' : 'HTML report for integrate_abundances_with_metabolism app',
                     'description' : 'HTML report for integrate_abundances_with_metabolism app'}

        html_report_list.append({'path': report_file_path,
                                 'name': 'index.html',
                                 'description': 'HTML report for integrate_abundances_with_metabolism app'})

        return html_report_list
    
    def _build_scatterplot(self,data,title="Default Title"):

        od_fig = figure()
        od_fig.xaxis.axis_label = self.conditions_ids[0]
        od_fig.yaxis.axis_label = self.conditions_ids[1]
        
        df = pd.DataFrame(data,columns=self.conditions_ids)
        od_fig.title.text = title
        od_fig.circle(x=self.conditions_ids[0], y=self.conditions_ids[1], source=df, color="black")

        slope_line = Slope(gradient=1, y_intercept=0, line_color="red")
        od_fig.add_layout(slope_line)

        return od_fig

    def _add_to_scatterplot(self,figure,data,color="black"):

        df = pd.DataFrame(data,columns=self.conditions_ids)
        figure.circle(x=self.conditions_ids[0], y=self.conditions_ids[1], source=df, fill_color=color, line_color="black", size=6)
        
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
        return [expdata_obj,conditions_ids,features_ids,feature_lookup_dict]
        
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

        return [reaction_score_comparison_dict, reaction_percentile_comparison_dict]

    def _compile_mahalanobis_dist_pvalue(self, data, threshold):

        # I don't know the math well enough to follow what's going on, but I used 
        # the recipe described here:
        # https://www.machinelearningplus.com/statistics/mahalanobis-distance/

        df = pd.DataFrame(data,columns=self.conditions_ids)

        # Covariance matrix via numpy
        cov_mat = np.cov(df.values.T)

        # Inverse covariance matrix via scipy.linalg
        inv_cov_mat = sp.linalg.inv(cov_mat)
        
        # two terms required, second using dot product
        data_minus_mean = df - np.mean(df)
        left_term = np.dot(data_minus_mean,inv_cov_mat)

        # dot product
        mahalanobis = np.dot(left_term, data_minus_mean.T)
        df['mahalanobis'] = mahalanobis.diagonal()

        # chi-squared p-values with one degree of freedom (two sets of variables)
        df['p_value'] = 1-sp.stats.chi2.cdf(df['mahalanobis'], 1)

        # find the outliers below a given threshold, i.e. p < 0.01
        outliers=df.loc[df.p_value < threshold]
        # this is used when you want to just plot the p-values alone
        outliers.index.name = 'reactions'

        #Need to return the mapping between reactions and the p-values
        return [df, outliers]
        
    def _load_subsystems(self):

        # Collect Core Subsystems and Reactions (From PlantSEED Reference Genome)
        # Core_SS_Classes=('Central Carbon','Amino acids','Nucleic acids')
        # Reactions_Subsystems=dict()

        # genome_ref = 'PlantSEED_v2/PlantSEED_Arabidopsis'
        # genome_obj = self.dfu.get_objects({'object_refs': [genome_ref]})['data'][0]
        # for ftr in genome_obj['data']['features']:
        #     if(len(ftr['functional_descriptions'])>0):
        #         for fd in ftr['functional_descriptions']:
        #             (mclass,subsystem,reaction_str)=fd.split('::')
        #             if(mclass not in Core_SS_Classes):
        #                 continue
        #             for rxn in reaction_str.split('|'):
        #                 if(rxn not in Reactions_Subsystems):
        #                     Reactions_Subsystems[rxn]=dict()
        #                 Reactions_Subsystems[rxn][subsystem]=1
        # print("Collected "+str(len(Reactions_Subsystems))+" core reactions")
        pass
    
    def _integrate_abundances(self, model_obj, conditions_ids, feature_lookup_dict, expdata_obj):

        reaction_values_matrix=list()
        reactions_ids=list()
        minmax_expscore_dict=dict()
        reaction_classification_dict=dict()
        for mdlrxn in range(len(model_obj['data']['modelreactions'])):
            mdlrxn_obj=model_obj['data']['modelreactions'][mdlrxn]
            reactions_ids.append(mdlrxn_obj['id'])

            #Reactions in one compartment vs reactions in multiple compartments
            #Reactions with multi-subunit enzymes vs Reactions with individual enzymes
            [base_rxn,cpt_id]=mdlrxn_obj['id'].split('_')
            if(base_rxn not in reaction_classification_dict):
                reaction_classification_dict[base_rxn]={'cpts':{},'single':0,'multiple':0}

            for prt in mdlrxn_obj['modelReactionProteins']:
                if(len(prt['modelReactionProteinSubunits'])==1):
                    reaction_classification_dict[base_rxn]['single']+=1
                if(len(prt['modelReactionProteinSubunits'])>1):
                    reaction_classification_dict[base_rxn]['multiple']+=1

            rxndata_row=list()
            for experiment in range(len(conditions_ids)):
                if(conditions_ids[experiment] not in minmax_expscore_dict):
                    minmax_expscore_dict[conditions_ids[experiment]]={'max':-sys.maxsize-1,'min':sys.maxsize}

                # Maximal gene expression for a reaction
                reaction_score='nan'
                for prt in mdlrxn_obj['modelReactionProteins']:

                    # Minimal gene expression for a complex
                    complex_score='nan'
                    for sbnt in prt['modelReactionProteinSubunits']:

                        # Maximal gene expression for a subunit
                        subunit_score='nan'
                        for feature in sbnt['feature_refs']:
                            feature=feature.split('/')[-1]
                            ftr_score = expdata_obj['data']['data']['values'][feature_lookup_dict[feature]][experiment]

                            if(ftr_score < minmax_expscore_dict[conditions_ids[experiment]]['min']):
                                minmax_expscore_dict[conditions_ids[experiment]]['min'] = ftr_score

                            if(ftr_score > minmax_expscore_dict[conditions_ids[experiment]]['max']):
                                minmax_expscore_dict[conditions_ids[experiment]]['max'] = ftr_score

                            # Maximal gene expression for a subunit
                            if(subunit_score == 'nan' or subunit_score < ftr_score):
                                subunit_score = ftr_score

                        # Minimal gene expression for a complex
                        if(subunit_score != 'nan'):
                            if(complex_score == 'nan' or complex_score > subunit_score):
                                complex_score = subunit_score

                    # Maximal gene expression for a reaction
                    if(complex_score != 'nan'):
                        if(reaction_score == 'nan' or reaction_score < complex_score):
                            reaction_score = complex_score
            
                if(reaction_score == 'nan'):
                    reaction_score = float(-sys.maxsize-1)
                else:
                    if(cpt_id not in reaction_classification_dict[base_rxn]['cpts']):
                        reaction_classification_dict[base_rxn]['cpts'][cpt_id]=[]
                    reaction_classification_dict[base_rxn]['cpts'][cpt_id].append({conditions_ids[experiment]:reaction_score})

                rxndata_row.append(reaction_score)
            reaction_values_matrix.append(rxndata_row)

        return (reaction_values_matrix, reactions_ids)

    def __init__(self, config, ctx):
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        self.scratch = config['scratch']
        self.dfu = DataFileUtil(self.callback_url)
        self.kbr = KBaseReport(self.callback_url)
        self.report_uuid = str(uuid.uuid4())

    def integrate_abundances_with_metabolism(self, params):

        self._validate_params(params, {'input_ws',
                                       'input_fbamodel',
                                       'input_expression_matrix',
                                       'output_reaction_matrix'}, {})

        model_ref = params['input_ws']+'/'+params['input_fbamodel']
        [model_obj,reaction_index] = self._load_fbamodel(model_ref)

        expression_ref = params['input_ws']+'/'+params['input_expression_matrix']
        [expdata_obj,conditions_ids,features_ids,feature_index] = self._load_expression_matrix(expression_ref)

        #Matrix of Figures
        figure_matrix=list()
        figure_array = list()
    
        feature_comparison_dict = self._compile_genome_scores(expdata_obj['data']['data']['values'])
        feature_figure = self._build_scatterplot(feature_comparison_dict,title="Genome Features Expression Abundances")
        figure_array.append(feature_figure)

        rxndata_obj = {'row_ids':[],'col_ids':[],'values':[]}
        rxndata_obj['col_ids']=conditions_ids

        (reaction_values_matrix, reactions_ids) = self._integrate_abundances(model_obj,conditions_ids,feature_index,expdata_obj)
        [reaction_score_comparison_dict, reaction_percentile_comparison_dict] = self._compile_model_scores_percentiles(reaction_values_matrix)

        self._add_to_scatterplot(feature_figure, reaction_score_comparison_dict, color="lightgreen")

        reaction_figure = self._build_scatterplot(reaction_percentile_comparison_dict,title="Model Reactions Percentile Rank (p<0.01)")
        figure_array.append(reaction_figure)

        [data_frame,outliers] = self._compile_mahalanobis_dist_pvalue(reaction_percentile_comparison_dict,0.01)

        self._add_to_scatterplot(reaction_figure, outliers, color="red")

        figure_matrix.append(figure_array)

        rxndata_obj['row_ids']=reactions_ids
        rxndata_obj['values']=reaction_values_matrix

        ReactionMatrix_obj = {'type':'KBaseMatrices.ReactionMatrix','name':params['output_reaction_matrix'],
                              'data':{'scale':'raw','description':'reaction expression score',
                                      'fbamodel_ref':model_ref,
                                      'expression_ref':expression_ref,
                                      'data':rxndata_obj}}

        ws_id = self.dfu.ws_name_to_id(params['input_ws'])
        saved_matrix_dict = self.dfu.save_objects({'id':ws_id,'objects':[ReactionMatrix_obj]})[0]
        saved_matrix_ref = "{}/{}/{}".format(saved_matrix_dict[6],saved_matrix_dict[0],saved_matrix_dict[4])
        saved_matrix_desc = "Reaction matrix: "+params['output_reaction_matrix']

        return self._build_report(figure_matrix, params['input_ws'], saved_matrix_ref, saved_matrix_desc)

"""
#A = single compartment, single subunits
#B = multiple compartments, single subunits
#C = single compartment, multiple subunits
#D = multiple compartments, multple subunits
row_counts={'A':0,'B':0,'C':0,'D':0}
for rxn in sorted(reaction_classification_dict):
rc_dict = reaction_classification_dict[rxn]
if(len(rc_dict['cpts'])==1 and rc_dict['single']>0):
row_counts['A']+=1
if(len(rc_dict['cpts'])>1 and rc_dict['single']>0):
row_counts['B']+=1
if(len(rc_dict['cpts'])==1 and rc_dict['multiple']>0):
row_counts['C']+=1
if(len(rc_dict['cpts'])>1 and rc_dict['multiple']>0):
row_counts['D']+=1
                                             
#Compose report string
html_lines.append('<h3 style="text-align: center">Integrate Abundances with Metabolism Report</h3>')
html_lines.append("<p>The \"Integrate Abundances with Metabolism\" app has finished running:</br>")
html_lines.append("The app integrated the gene abundances from the "+input_params['input_expression_matrix']+" ExpressionMatrix with the ")
html_lines.append(input_params['input_fbamodel']+" FBAModel, resulting in the "+input_params['output_reaction_matrix']+" ReactionMatrix.</p>")

html_lines.append('<table class="table table-bordered table-striped">')

#Reaction Compartment Number of Enzymatic Subunits Experiments
internal_header_line = "</td><td>".join(['Reaction','Compartment','# Subunits'] + conditions_ids)
html_lines.append('<thead>')
html_lines.append('<tr><td>'+internal_header_line+'</td></tr>')
html_lines.append('</thead>')

html_lines.append('<tbody>')
for rxn in sorted(reaction_classification_dict):
rc_dict = reaction_classification_dict[rxn]
for cpt in rc_dict['cpts']:
scores = list()
for exp in rc_dict['cpts'][cpt]:
scores.append(str(list(exp.values())[0]))
row = '</td><td>'.join([rxn,cpt,'0'] + scores)
html_lines.append('<tr><td>'+row+'</td></tr>')
html_lines.append('</tbody>')
html_lines.append('</table>')
html_string="\n".join(html_lines)
"""
