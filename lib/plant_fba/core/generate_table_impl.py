from urllib.request import urlopen
import time
import json
import os
import re

from plant_fba.core.fetch_plantseed_impl import FetchPlantSEEDImpl

PS_url = 'https://raw.githubusercontent.com/ModelSEED/PlantSEED/'
PS_tag = 'v3.0.0'
  
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

    def generate_table(self, reactions, complexes=None):

        html_lines=list()        
        html_lines.append('<table class="table table-bordered table-striped">')

        header_list = ["Reaction","Compartment","Roles","Enzymes","EC numbers","Subsystems","Classes"]

        html_lines.append('<thead>')            
        internal_header_line = "</td><td>".join(header_list)
        html_lines.append('<tr><td>'+internal_header_line+'</td></tr>')
        html_lines.append('</thead>')

        html_lines.append("<tbody>")

        #######################################################################
        # Each row in the table is unique to the reaction/compartment combination
        for reaction in sorted(reactions.keys()):
            for compartment in reactions[reaction]['compartments']:
                compartment_str = Template_Compartment_Mapping[compartment]

                compartmentalized_complexes=list()

                #######################################################################
                # The code below is taken from the code for generating a ModelTemplate
                # It identifies the single letter reaction id for transport reactions
                # Accordingly, the compartments should all be sorted
                # So a compartment index of 0 matches the first position in the compartment list
                # The order is curated in the PlantSEED database

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

                #######################################################################

                model_reaction_id = reaction+"_"+reaction_cpt_id+"0"
                if(complexes is not None and \
                       model_reaction_id in complexes and \
                       complexes[model_reaction_id] not in compartmentalized_complexes):
                    compartmentalized_complexes.append(complexes[model_reaction_id])
            
                complexes_str="; ".join(compartmentalized_complexes)
                roles_str="; ".join(sorted(reactions[reaction]['roles']))
                ecs_str="; ".join(sorted(reactions[reaction]['ecs']))
                classes_str="; ".join(sorted(reactions[reaction]['classes']))

                subsystems_str ="; ".join(sorted(reactions[reaction]['subsystems']))
                subsystems_str = subsystems_str.replace("_in_plants","")
                subsystems_str = subsystems_str.replace("_"," ")

                html_lines.append("<tr>")
                internal_row_line = "</td><td>".join([reaction,compartment_str,roles_str,
                                                      complexes_str,ecs_str,subsystems_str,classes_str])
                html_lines.append("<td>"+internal_row_line+"</td>")
                html_lines.append("</tr>")

        html_lines.append("</tbody>")
        html_lines.append("</table>")

        return "\n".join(html_lines)

def main():

    plantseed = FetchPlantSEEDImpl()
    reactions_data = plantseed.fetch_reactions()

    table = GenerateTableImpl()
    table_html_string = table.generate_table(reactions_data)

    with open(os.path.join('../../../data','app_report_templates',
                           'integrate_abundances_report_tables_template.html')) as report_template_file:
        report_template_string = report_template_file.read()

    # Insert html table
    table_report_string = report_template_string.replace('*TABLES*', table_html_string)

    table_html_file="annotation_table_output.html"
    with open(os.path.join('../../../data','app_report_templates',
                           table_html_file),'w') as table_file:
        table_file.write(table_report_string)

    table.log(message="Table written to file: "+table_html_file)

if(__name__ == "__main__"):
    main()
