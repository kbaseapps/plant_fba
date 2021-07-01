from urllib.request import urlopen
import time
import json
import os
import re

PS_url = 'https://raw.githubusercontent.com/ModelSEED/PlantSEED/'
PS_tag = 'v3.0.0'
  
def _load_reactions():

    reactions_data = dict()

    # Load these directly from PlantSEED_Roles.json
    PS_Roles = json.load(urlopen(PS_url+PS_tag+'/Data/PlantSEED_v3/PlantSEED_Roles.json'))

    for entry in PS_Roles:
        if(entry['include'] is False):
#            print(entry['role'])
            continue

        main_class_ss = list()
        main_class = list()
        for metabolic_class in entry['classes']:
            if(len(entry['classes'][metabolic_class].keys())>0):
                main_class.append(metabolic_class)

            for ss in entry['classes'][metabolic_class].keys():
                if(ss not in main_class_ss):
                    main_class_ss.append(ss)

        for rxn in entry['reactions']:
            if(rxn not in reactions_data):
                reactions_data[rxn]={'ecs':[],
                                     'roles':[],
                                     'classes':[],
                                     'subsystems':[],
                                     'compartments':[]}

            if(entry['role'] not in reactions_data[rxn]['roles']):
                reactions_data[rxn]['roles'].append(entry['role'])

            for mclass in main_class:
                if(mclass not in reactions_data[rxn]['classes']):
                    reactions_data[rxn]['classes'].append(mclass)

            for subsystem in main_class_ss:
                if(subsystem not in reactions_data[rxn]['subsystems']):
                    reactions_data[rxn]['subsystems'].append(subsystem)

            for cpt in entry['localization']:
                if(cpt not in reactions_data[rxn]['compartments']):
                    reactions_data[rxn]['compartments'].append(cpt)

    for rxn in reactions_data:
        for role in reactions_data[rxn]['roles']:
            if('EC' in role):
                match = re.search(r"\d+\.[\d-]+\.[\d-]+\.[\d-]+", role)
                if(match is not None):
                    if(match.group(0) not in reactions_data[rxn]['ecs']):
                        reactions_data[rxn]['ecs'].append(match.group(0))

    print("Collected "+str(len(reactions_data))+" core reactions")
    return reactions_data

# These should be retrieved from the Template data
Template_Compartment_Mapping={'c':'cytosol', 'g':'golgi', 'w':'cellwall',
                              'n':'nucleus', 'r':'endoplasm',
                              'v':'vacuole', 'cv':'vacuole',
                              'd':'plastid', 'cd':'plastid',
                              'm':'mitochondria','cm':'mitochondria',
                              'mj':'mitointer', 'ce':'extracellular',
                              'x':'peroxisome', 'cx':'peroxisome',
                              'e':'extracellular','de':'plastid'}

class GenerateTableImpl:

    def log(self,message, prefix_newline=False):
        time_str = time.strftime('%Y-%m-%d %H:%M:%S', time.gmtime(time.time()))
        print(('\n' if prefix_newline else '') + time_str + ': ' + message)

    def __init__(self):
        pass

    def generate_table(self, complexes=None):
        reactions_data = _load_reactions()

        html_lines=list()        
        html_lines.append('<table class="table table-bordered table-striped">')

        header_list = ["Reactions","Complexes","Roles","Compartments","EC numbers","Subsystems","Classes"]

        html_lines.append('<thead>')            
        internal_header_line = "</td><td>".join(header_list)
        html_lines.append('<tr><td>'+internal_header_line+'</td></tr>')
        html_lines.append('</thead>')

        html_lines.append("<tbody>")

        for reaction in sorted(reactions_data.keys()):

            compartmentalized_complexes=list()

            #######################################################################
            # The code below is taken from the code for generating a ModelTemplate
            # It identifies the single letter reaction id for transport reactions
            # Accordingly, the compartments should all be sorted
            # So a compartment index of 0 matches the first position in the compartment list
            # The order is curated in the PlantSEED database

            compartments = list()
            for compartment in reactions_data[reaction]['compartments']:
                compartments.append(Template_Compartment_Mapping[compartment])

                reaction_cpt_id = compartment

                # If its a transporter, need to update the reaction compartment id
                if(len(compartment)==2):

                    # The rule is that it is always the non-cytosolic compartment
                    if('c' in compartment):
                        for cpt in compartment:
                            if(cpt != 'c'):
                                reaction_cpt_id = cpt

                    # With two main exceptions:
                    # 1) whether its an extracellular transporter
                    if('e' in compartment):
                        for cpt in compartment:
                            if(cpt != 'e'):
                                reaction_cpt_id = cpt
                                
                    # 2) whether its an intraorganellar transporter
                    if('j' in compartment):
                        reaction_cpt_id = 'j'

                model_reaction_id = reaction+"_"+reaction_cpt_id+"0"
                if(model_reaction_id in complexes and complexes[model_reaction_id] not in compartmentalized_complexes):
                    compartmentalized_complexes.append(complexes[model_reaction_id])

            complexes_str="; ".join(compartmentalized_complexes)
            compartments="; ".join(sorted(compartments))

            #######################################################################

            roles="; ".join(sorted(reactions_data[reaction]['roles']))
            ecs="; ".join(sorted(reactions_data[reaction]['ecs']))
            classes="; ".join(sorted(reactions_data[reaction]['classes']))

            subsystems="; ".join(sorted(reactions_data[reaction]['subsystems']))
            subsystems = subsystems.replace("_in_plants","")
            subsystems = subsystems.replace("_"," ")

            html_lines.append("<tr>")
            internal_row_line = "</td><td>".join([reaction,complexes_str,roles,compartments,ecs,subsystems,classes])
            html_lines.append("<td>"+internal_row_line+"</td>")
            html_lines.append("</tr>")

        html_lines.append("</tbody>")
        html_lines.append("</table>")

        return "\n".join(html_lines)

def main():
    table = GenerateTableImpl()
    table_html_string = table.generate_table()

    with open(os.path.join('app_report_templates','annotation_report_tables_template.html')) as report_template_file:
        report_template_string = report_template_file.read()

    # Insert html table
    table_report_string = report_template_string.replace('*TABLES*', table_html_string)

    table_html_file="annotation_table_output.html"
    with open(os.path.join('app_report_templates',table_html_file),'w') as table_file:
        table_file.write(table_report_string)

    table.log(message="Table written to file: "+table_html_file)

if(__name__ == "__main__"):
    main()
