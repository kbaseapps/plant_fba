# -*- coding: utf-8 -*-
#BEGIN_HEADER
import os
import sys
from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.DataFileUtilClient import DataFileUtil
#END_HEADER


class plant_fba:
    '''
    Module Name:
    plant_fba

    Module Description:
    A KBase module: plant_fba
    '''

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "0.0.1"
    GIT_URL = "git@github.com:kbaseapps/plant_fba.git"
    GIT_COMMIT_HASH = "ba8d952cf73ca319caa9d8dad6f713b5cded098e"

    #BEGIN_CLASS_HEADER
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        self.shared_folder = config['scratch']
        self.dfu = DataFileUtil(self.callback_url)
        #END_CONSTRUCTOR
        pass


    def integrate_abundances_with_metabolism(self, ctx, input_params):
        """
        :param input_params: instance of type "IntegrateAbundancesParams" ->
           structure: parameter "input_ws" of String, parameter
           "input_expression_matrix" of String, parameter "input_fbamodel" of
           String, parameter "output_reaction_matrix" of String
        :returns: instance of type "ReportResults" -> structure: parameter
           "report_name" of String, parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: output_report
        #BEGIN integrate_abundances_with_metabolism

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

        model_ref = input_params['input_ws']+'/'+input_params['input_fbamodel']
        model_obj = self.dfu.get_objects({'object_refs': [model_ref]})['data'][0]
        print("Number of reactions: "+str(len(model_obj['data']['modelreactions'])))

        # Fetch expression data
        expdata_ref = input_params['input_ws']+'/'+input_params['input_expression_matrix']
        expdata_obj = self.dfu.get_objects({'object_refs': [expdata_ref]})['data'][0]
        exp_ids = expdata_obj['data']['data']['col_ids']
        ftr_ids = expdata_obj['data']['data']['row_ids']
        feature_lookup_dict=dict()
        for index in range(len(ftr_ids)):
            feature_lookup_dict[ftr_ids[index]]=index

        # Initialize reaction data
        rxndata_obj = {'row_ids':[],'col_ids':[],'values':[]}
        rxndata_obj['col_ids']=exp_ids
        for mdlrxn in range(len(model_obj['data']['modelreactions'])):
            mdlrxn_obj=model_obj['data']['modelreactions'][mdlrxn]
            rxndata_obj['row_ids'].append(mdlrxn_obj['id'])

        minmax_expscore_dict=dict()
        for mdlrxn in range(len(model_obj['data']['modelreactions'])):
            mdlrxn_obj=model_obj['data']['modelreactions'][mdlrxn]
    
            rxndata_row=list()
            for experiment in range(len(exp_ids)):
                if(exp_ids[experiment] not in minmax_expscore_dict):
                    minmax_expscore_dict[exp_ids[experiment]]={'max':-sys.maxsize-1,'min':sys.maxsize}

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
                            
                            if(ftr_score < minmax_expscore_dict[exp_ids[experiment]]['min']):
                                minmax_expscore_dict[exp_ids[experiment]]['min'] = ftr_score

                            if(ftr_score > minmax_expscore_dict[exp_ids[experiment]]['max']):
                                minmax_expscore_dict[exp_ids[experiment]]['max'] = ftr_score

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

                rxndata_row.append(reaction_score)
                rxndata_obj['values'].append(rxndata_row)

        rxnvalue_matrix = {'type':'KBaseFeatureValues.ExpressionMatrix','name':input_params['output_reaction_matrix'],
                           'data':{'scale':'1.0','type':'reaction expression score','data':rxndata_obj}}

        ws_id = self.dfu.ws_name_to_id(input_params['input_ws'])
        self.dfu.save_objects({'id':ws_id,'objects':[rxnvalue_matrix]})

        output_report=dict()

        #END integrate_abundances_with_metabolism

        # At some point might do deeper type checking...
        if not isinstance(output_report, dict):
            raise ValueError('Method integrate_abundances_with_metabolism return value ' +
                             'output_report is not type dict as required.')
        # return the results
        return [output_report]
    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK",
                     'message': "",
                     'version': self.VERSION,
                     'git_url': self.GIT_URL,
                     'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
