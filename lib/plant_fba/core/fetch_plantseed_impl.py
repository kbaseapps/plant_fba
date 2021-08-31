from urllib.request import urlopen
import time
import json
import os
import re

PS_url = 'https://raw.githubusercontent.com/ModelSEED/PlantSEED/'
PS_tag = 'kbase_release'

class FetchPlantSEEDImpl:

    def fetch_reactions(self, PS_Roles):

        reactions_data = dict()

        for entry in PS_Roles:
            if(entry['include'] is False):
                # print(entry['role'])
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

    def log(self,message, prefix_newline=False):
        time_str = time.strftime('%Y-%m-%d %H:%M:%S', time.gmtime(time.time()))
        print(('\n' if prefix_newline else '') + time_str + ': ' + message)

    def __init__(self):
        pass

def main():

    # Load these directly from PlantSEED_Roles.json
    PS_Roles = json.load(urlopen(PS_url+PS_tag+'/Data/PlantSEED_v3/PlantSEED_Roles.json'))

    plantseed = FetchPlantSEEDImpl()
    plantseed_data = plantseed.fetch_reactions(PS_Roles)
    print(plantseed_data)

if(__name__ == "__main__"):
    main()
