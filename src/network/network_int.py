'''MAIN SCRIPT
input: confident interface interaction between proteins predictd by alphabridge
ouput: json file containing the 3 different networks:
    1) network that shows confident contact interface between protein (not merged)
    2) network that shows confident contact interface between protein (merged)
    3) network that shows a general overview of the confident physical interaction between proteins'''

#import required packages

import src.network.complex_not_merged as no_merg
import src.network.complex_merged as merg
import src.network.protein_network as protein_net
import sys 
import os
import pandas as pd
import json



#CONVERTING THE JSON FILE (output of alphabridge) IN A DF

class INTERACTIVE_NETWORK:
    
    def __init__(self, alphabridge_dict):
        self.alphabridge_dict = alphabridge_dict
    
    def from_json_to_df(self, threshold):
    
        alphabridge_dict = self.alphabridge_dict
        all_data = []
        # Process JSON data and store in a list
        #from auth to label
        auth2label = {}
        chains = alphabridge_dict['structure'][0]['chains']
        for entity_type, rec_list in chains.items():
            for rec in rec_list:
                #print(rec)
                #print(rec['auth_asym_id'], rec['label_asym_id'])
                auth_asym_id = rec['auth_asym_id']
                label_asym_id = rec['label_asym_id']
                if not auth_asym_id in auth2label:
                    auth2label[auth_asym_id] = str()
                auth2label[auth_asym_id] = label_asym_id
        #from label to auth
        label2auth = {value:key for key, value in auth2label.items()}
        rows = []
        if "interactions" in alphabridge_dict:
            for interaction in alphabridge_dict["interactions"]:
                if interaction.get("cut-off") == threshold:
                    for interface in interaction["interfaces"]:
                        interface_name = interface["interface_name"]
                        for link in interface["links"]:
                            row = {
                                #folder_name": folder_name,
                                "interface": interface_name,
                                "interface_id": interface["interface_id"],
                                "prot_1": link["first"]["asym_id"],
                                "start_1": link["first"]["link_range"]["start"],
                                "end_1": link["first"]["link_range"]["end"],
                                "prot_2": link["second"]["asym_id"],
                                "start_2": link["second"]["link_range"]["start"],
                                "end_2": link["second"]["link_range"]["end"],
                            }
                            rows.append(row)

        # Add rows to the main list
        all_data.extend(rows)
        #create the df for storing the score
        pairwise_interaction = alphabridge_dict['structure'][0]['pairwise_interaction']
        #print(pairwise_interaction)
        df_pairwise_interaction = pd.DataFrame(pairwise_interaction)
        #print(df_pairwise_interaction)


        # Convert the list to a DataFrame, THIS IS THE ONE THAT IS NEEDED FOR BOTH
        df = pd.DataFrame(all_data)
        #print(df)


        # form auth to label for whatever is not a protein
        #for column in ['prot_1', 'prot_2']:
        #    df[column] = df[column].replace(auth2label)
        #print(df)
        return df, label2auth, df_pairwise_interaction, auth2label
    
    def get_network_info(self):
        
        all_threshold = {}
        threshold_list = [0.5, 0.75, 0.9]
        
        for threshold in threshold_list:
            
            df, label2auth, df_pairwise_interaction, auth2label = self.from_json_to_df(threshold)
            
            j_not_merged = no_merg.get_protein_network_no_merging(df,label2auth)
            j_proteins = protein_net.get_protein_network(df,label2auth,df_pairwise_interaction,threshold,auth2label)
            
            # Combine them into a single dictionary
            combined_networks = {"cut-off": threshold,
                "network_not_merged": j_not_merged,
            # "network_merged": j_merged,
                "protein_network": j_proteins 
            }
            #combine everything
            all_threshold[f"network at {threshold}"] = combined_networks
        
        return all_threshold


'''MAIN FUNCTION FOR OBTANING BOTH NETWORK: WITHOUT MERGING (MORE NODES) AND WITH MERGIN (LESS NODES)'''


'''    in_dir = args.in_dir
    threshold = [float(t) for t in args.threshold]
    if not os.path.exists(args.o_dir):
        os.makedirs(args.o_dir)
    all_threshold = {}
    #iterate for every possible threshold
    for th in threshold:
        #print(f"th:{th}")
        # funtion that creates the df that will be used as an inpt for the next fuction
        df, label2auth, df_pairwise_interaction, auth2label = from_json_to_df(in_dir, th)
        #function --> INFORMATION FOR THE NETWORK WITH MORE NODES(NO MERGED)
        j_not_merged = no_merg.get_protein_network_no_merging(df,label2auth)
        print("FINISHED THE NOT MERGED NETWORK")
        #function --> INFORMATION FOR THE NETWORK WITH LESS NODES(MERGED)
        #j_merged = merg.get_protein_network_merging(df)
        #print("FINISHED THE MERGED NETWORK")
        #function --> INFORMATION FOR THE NETWORK WHERE NOES == PROTEIN
        j_proteins = protein_net.get_protein_network(df,label2auth,df_pairwise_interaction,th,auth2label)
        print("FINISHED THE WHOLE PROTEIN NETWORK")
        # Combine them into a single dictionary
        combined_networks = {"cut-off": th,
            "network_not_merged": j_not_merged,
        # "network_merged": j_merged,
            "protein_network": j_proteins 
        }
        #combine everything
        all_threshold[f"network at {th}"] = combined_networks

    # Save to a JSON file
    with open(f"{args.o_dir}/combined_networks.json", "w") as f:
        json.dump(all_threshold, f, indent=4)'''

