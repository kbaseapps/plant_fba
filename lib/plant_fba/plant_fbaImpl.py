# -*- coding: utf-8 -*-
#BEGIN_HEADER

import os
import sys
import re
import copy

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
    GIT_COMMIT_HASH = "0b8a7fdc6d1b98240c8e9aa1d8d81b492111d503"

    #BEGIN_CLASS_HEADER

    def convert_search_role(self,role):

        searchrole = role

        #Remove spaces
        searchrole = searchrole.strip()
        searchrole = searchrole.replace(' ','')

        #Make all lowercase
        searchrole = searchrole.lower()

        #Remove EC and parentheses
        searchrole = re.sub(r'\(ec[\d-]+\.[\d-]\.[\d-]\.[\d-]\)', '', searchrole)

        return searchrole

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

    def reconstruct_plant_metabolism(self, ctx, input_params):
        """
        :param input_params: instance of type "ReconstructMetabolismParams"
           -> structure: parameter "input_ws" of String, parameter
           "input_genome" of String, parameter "output_ws" of String,
           parameter "output_fbamodel" of String, parameter "template" of
           String, parameter "template_ws" of String
        :returns: instance of type "ReportResults" -> structure: parameter
           "report_name" of String, parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: output_report
        #BEGIN reconstruct_plant_metabolism

        #Compile biochemistry information
        abbrev_cpt_dict=dict()
        cpt_name_dict=dict()
        with open('/kb/module/data/compartments.txt') as fh:
            for line in fh.readlines():
                line=line.strip('\r\n')
                array=line.split('\t')
                abbrev_cpt_dict[array[3]]=array[0]
                cpt_name_dict[array[0]]=array[2]

        rxn_name_dict=dict()
        with open('/kb/module/data/reactions.txt') as fh:
            for line in fh.readlines():
                line=line.strip('\r\n')
                array=line.split('\t')
                if(array[1] == ""):
                    array[1] = array[0]
                rxn_name_dict[array[0]]=array[1]

        cpd_props_dict=dict()
        with open('/kb/module/data/compounds.txt') as fh:
            for line in fh.readlines():
                line=line.strip('\r\n')
                array=line.split('\t')
                if(array[2] == ''):
                    array[2] = array[0]
                cpd_props_dict[array[0]]={'name':array[2],
                                          'formula':array[3],
                                          'charge':array[4]}

        #Retrieve Template, and compile indexes of roles and complexes
        template_ref = input_params['template_ws']+'/'+input_params['template']
        template_obj = self.dfu.get_objects({'object_refs': [template_ref]})['data'][0]

        searchroles_dict = dict()
        roles_dict = dict()
        for role in template_obj['data']['roles']:
            searchrole = self.convert_search_role(role['name'])
            searchroles_dict[searchrole]=role['id']
            roles_dict[role['id']]=role['name']

        complex_dict = dict()
        for cpx in template_obj['data']['complexes']:
            complex_dict[cpx['id']]=cpx

        #Retrieve Genome annotation as dict
        role_cpt_ftr_dict=dict()
        genome_ref = input_params['input_ws']+'/'+input_params['input_genome']
        genome_obj = self.dfu.get_objects({'object_refs': [genome_ref]})['data'][0]
        for feature in genome_obj['data']['features']:
            if(len(feature['functions'])>0):
                for function_comment in feature['functions']:

                    #Split for comments and retrieve compartments
                    function_cpt_list = function_comment.split("#")
                    for i in range(len(function_cpt_list)):
                        function_cpt_list[i]=function_cpt_list[i].strip()

                    function = function_cpt_list.pop(0)
                    roles = re.split("\s*;\s+|\s+[\@\/]\s+", function)
                    for role in roles:
                        
                        searchrole = self.convert_search_role(role)
                        if(searchrole not in searchroles_dict):
                            continue

                        role_id = searchroles_dict[searchrole]

                        if(role_id not in role_cpt_ftr_dict):
                            role_cpt_ftr_dict[role_id]=dict()

                        for cpt in function_cpt_list:
                            abbrev_cpt=cpt
                            if(cpt not in abbrev_cpt_dict):
                                print("No compartmental abbreviation found for "+cpt)
                            else:
                                abbrev_cpt = abbrev_cpt_dict[cpt]

                            if(abbrev_cpt not in role_cpt_ftr_dict[role_id]):
                                role_cpt_ftr_dict[role_id][abbrev_cpt]=dict()

                            role_cpt_ftr_dict[role_id][abbrev_cpt][feature['id']]=1

        #Default dictionaries for objects needed for a model reaction
        default_mdlcpt_dict = { 'id': 'u0', 'label': 'unknown',
                                'pH': 7, 'potential': 0, 'compartmentIndex': 0,
                                'compartment_ref': '~//' }

        default_mdlcpd_dict = { 'id': '', 'charge': 0, 'formula': '',
                                'name': '', 'compound_ref': '',
                                'modelcompartment_ref': '~/modelcompartments/id/u0' }

        default_mdlrxn_dict = { 'id': '', 'direction': '', 'protons': 0,
                                'name': '', 'reaction_ref': '', 'probability': 0,
                                'modelcompartment_ref': '',
                                'modelReactionReagents': [], 'modelReactionProteins': [] }

        #Lookup dictionaries for compartments and compounds, to avoid duplicating them
        mdlcpts_dict = dict()
        mdlcpds_dict = dict()

        #Create New, but Empty Plant Reconstruction
        new_model_obj = { 'id' : input_params['output_fbamodel'], 'type' : "GenomeScale", 'source' : "KBase", 'source_id' : "PlantSEED_v2", 
                          'template_ref' : template_ref, 'genome_ref' : genome_ref, 'name' : input_params['output_fbamodel'],
                          'modelreactions' : [], 'modelcompounds' : [], 'modelcompartments' : [], 'biomasses' : [],
                          'gapgens' : [], 'gapfillings' : [] }
        
        for template_rxn in template_obj['data']['reactions']:
            if(template_rxn['type'] == 'gapfilling'):
                continue

            template_rxn_cpt = template_rxn['templatecompartment_ref'].split('/')[-1]

            proteins_list = list()
            #complex_ref and source are optional fields
            default_protein_dict = {'note':template_rxn['type'], 'complex_ref':'', 'modelReactionProteinSubunits':[]}
            for cpx_ref in template_rxn['templatecomplex_refs']:
                cpx_id = cpx_ref.split('/')[-1]
                model_complex_ref = "~/template/complexes/id/"+cpx_id

                new_protein_dict = copy.deepcopy(default_protein_dict)
                new_protein_dict['complex_ref']=model_complex_ref

                complex_present=False
                subunits_list = list()
                default_subunit_dict = {'role':'','triggering':0,'optionalSubunit':0,'note':'','feature_refs':[]}
                matched_role_dict = dict()

                for cpxrole in complex_dict[cpx_id]['complexroles']:
                    role_id = cpxrole['templaterole_ref'].split('/')[-1]

                    if(role_id in role_cpt_ftr_dict):

                        for role_cpt in role_cpt_ftr_dict[role_id]:
                            role_cpt_present=False
                            if(template_rxn_cpt == role_cpt and cpxrole['triggering'] == 1):
                                complex_present=True
                                role_cpt_present=True

                            if(role_cpt_present == True):
                                new_subunit_dict = copy.deepcopy(default_subunit_dict)
                                new_subunit_dict['triggering'] = cpxrole['triggering']
                                new_subunit_dict['optionalSubunit'] = cpxrole['optional_role']
                                new_subunit_dict['role'] = roles_dict[role_id]
                                
                                for ftr in role_cpt_ftr_dict[role_id][role_cpt]:
                                    feature_ref = "~/genome/features/id/"+ftr
                                    new_subunit_dict['feature_refs'].append(feature_ref)
                                
                                matched_role_dict[role_id]=1
                                subunits_list.append(new_subunit_dict)

                if(complex_present == True):
                    for cpxrole in complex_dict[cpx_id]['complexroles']:
                        role_id = cpxrole['templaterole_ref'].split('/')[-1]
                        
                        if(role_id not in matched_role_dict):
                            new_subunit_dict = copy.deepcopy(default_subunit_dict)
                            new_subunit_dict['triggering'] = cpxrole['triggering']
                            new_subunit_dict['optionalSubunit'] = cpxrole['optional_role']
                            new_subunit_dict['note'] = "Complex-based-gapfilling"
                            subunits_list.append(new_subunit_dict)

                if(len(subunits_list)>0):
                    new_protein_dict['modelReactionProteinSubunits']=subunits_list

                proteins_list.append(new_protein_dict)

            #This is important, we need to use role-based annotation to determine whether
            #a reaction should even be added to the model
            if(template_rxn['type'] == 'conditional' and len(proteins_list)==0):
                continue

            #If the check passes, then, here, we instatiate the actual reaction that goes into the model
            new_mdlrxn_id = template_rxn['id']+'0'
            new_mdlcpt_id = template_rxn_cpt+'0'
            base_rxn_id = template_rxn['id'].split('_')[0]

            new_mdlrxn_dict = copy.deepcopy(default_mdlrxn_dict)
            new_mdlrxn_dict['id'] = new_mdlrxn_id
            new_mdlrxn_dict['name'] = rxn_name_dict[base_rxn_id]
            new_mdlrxn_dict['direction'] = template_rxn['direction']
            new_mdlrxn_dict['reaction_ref']='~/template/reactions/id/'+template_rxn['id']
            new_mdlrxn_dict['modelcompartment_ref']='~/modelcompartments/id/'+new_mdlcpt_id

            #Here we check and instantiate a new modelcompartment
            if(new_mdlcpt_id not in mdlcpts_dict):
                new_mdlcpt_dict = copy.deepcopy(default_mdlcpt_dict)
                new_mdlcpt_dict['id']=new_mdlcpt_id
                new_mdlcpt_dict['label']=cpt_name_dict[template_rxn_cpt]
                new_mdlcpt_dict['compartment_ref']='~/template/compartments/id/'+template_rxn_cpt
                mdlcpts_dict[new_mdlcpt_id]=new_mdlcpt_dict

            #Add Proteins as previously determined
            new_mdlrxn_dict['modelReactionProteins']=proteins_list

            #Add Reagents
            for template_rgt in template_rxn['templateReactionReagents']:
                template_rgt_cpd_cpt_id = template_rgt['templatecompcompound_ref'].split('/')[-1]
                (template_rgt_cpd,template_rgt_cpt)=template_rgt_cpd_cpt_id.split('_')

                #Check and add new model compartment 
                new_mdlcpt_id = template_rgt_cpt+'0'
                if(new_mdlcpt_id not in mdlcpts_dict):
                    new_mdlcpt_dict = copy.deepcopy(default_mdlcpt_dict)
                    new_mdlcpt_dict['id']=new_mdlcpt_id
                    new_mdlcpt_dict['label']=cpt_name_dict[template_rgt_cpt]
                    new_mdlcpt_dict['compartment_ref']='~/template/compartments/id/'+template_rgt_cpt
                    mdlcpts_dict[new_mdlcpt_id]=new_mdlcpt_dict
               
                #Add new model compounds
                new_mdlcpd_id = template_rgt_cpd_cpt_id+'0'
                base_cpd_id = template_rgt_cpd_cpt_id.split('_')[0]

                if(new_mdlcpd_id not in mdlcpds_dict):
                    new_mdlcpd_dict = copy.deepcopy(default_mdlcpd_dict)
                    new_mdlcpd_dict['id']=new_mdlcpd_id
                    new_mdlcpd_dict['name'] = cpd_props_dict[base_cpd_id]['name']
                    new_mdlcpd_dict['charge'] = float(cpd_props_dict[base_cpd_id]['charge'])
                    new_mdlcpd_dict['formula'] = cpd_props_dict[base_cpd_id]['formula']
                    new_mdlcpd_dict['compound_ref']='~/template/compounds/id/'+template_rgt_cpd
                    new_mdlcpd_dict['modelcompartment_ref']='~/modelcompartments/id/'+new_mdlcpt_id
                    mdlcpds_dict[new_mdlcpd_id]=new_mdlcpd_dict

                new_rgt_dict = {'coefficient' : template_rgt['coefficient'],
                                'modelcompound_ref' : '~/modelcompounds/id/'+new_mdlcpd_id}

                new_mdlrxn_dict['modelReactionReagents'].append(new_rgt_dict)

            new_model_obj['modelreactions'].append(new_mdlrxn_dict)

        #Having populated with list of reactions and biomass (to come), then add all compartments and compounds
        for cpt_id in mdlcpts_dict:
            new_model_obj['modelcompartments'].append(mdlcpts_dict[cpt_id])

        #Last, but key modelcompound is the biomass, need to add it explicitly
        biocpd_id = "cpd11416"
        mdlbiocpd_dict = copy.deepcopy(default_mdlcpd_dict)
        mdlbiocpd_dict['id'] = biocpd_id+'_c0'
        mdlbiocpd_dict['name'] = 'Biomass'
        mdlbiocpd_dict['compound_ref'] = "~/template/compounds/id/"+biocpd_id
        mdlbiocpd_dict['modelcompartment_ref'] = "~/modelcompartments/id/c0"
        mdlcpds_dict[mdlbiocpd_dict['id']] = mdlbiocpd_dict

        for cpd_id in mdlcpds_dict:
            new_model_obj['modelcompounds'].append(mdlcpds_dict[cpd_id])

        default_biomass_dict = { 'id': 'bio1', 'name': 'Plant leaf biomass', 'other': 1,
                                 'dna': 0, 'rna': 0, 'protein': 0, 'cellwall': 0,
                                 'lipid': 0, 'cofactor': 0, 'energy': 0, 'biomasscompounds': [] }

        default_biocpd_dict = { 'modelcompound_ref' : '', 'coefficient' : 0 }

        for template_biomass in template_obj['data']['biomasses']:
            new_template_biomass = copy.deepcopy(default_biomass_dict)
            new_template_biomass['id'] = template_biomass['id']
            new_template_biomass['name'] = template_biomass['name']

            for entry in ['dna','rna','protein','cellwall','lipid','cofactor','energy','other']:
                new_template_biomass[entry] = template_biomass[entry]

            for template_cpd in template_biomass['templateBiomassComponents']:
                new_biocpd_dict = copy.deepcopy(default_biocpd_dict)
                mdlcpd_id = template_cpd['templatecompcompound_ref'].split('/')[-1]+'0'
                if(mdlcpd_id not in mdlcpds_dict):
                    print("Missing: ",template_cpd)
                    continue
                new_biocpd_dict['modelcompound_ref'] = '~/modelcompounds/id/'+mdlcpd_id
                new_biocpd_dict['coefficient'] = template_cpd['coefficient']
                new_template_biomass['biomasscompounds'].append(new_biocpd_dict)
        
            new_model_obj['biomasses'].append(new_template_biomass)

        print("Saving metabolic reconstruction")
        model_ws_object = {'type' : 'KBaseFBA.FBAModel', 'name' : input_params['output_fbamodel'],
                           'data' : new_model_obj }

        ws_id = self.dfu.ws_name_to_id(input_params['input_ws'])
        self.dfu.save_objects({'id':ws_id,'objects':[model_ws_object]})

        output_report = dict()

        #END reconstruct_plant_metabolism

        # At some point might do deeper type checking...
        if not isinstance(output_report, dict):
            raise ValueError('Method reconstruct_plant_metabolism return value ' +
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
