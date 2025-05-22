#script for obtaining ONLY THE INTERACTIVE NETWORK WITH MERGE
#import required packages
import pandas as pd
import igraph as ig
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
from networkx.readwrite import json_graph;

#FUCNTION USED IN BOTH NETWORK WITH MERGING NODES AND NOT MERGIN NODES

#generating eges between the node tht belog to the same protein

def generate_intraprotein_edges(protein_dict):
    edges = []
    for protein, intervals in protein_dict.items():
        protein_edges = []
        for i in range(len(intervals)):
            for j in range(i + 1, len(intervals)):
                # Create a connection between every pair of intervals
                protein_edges.append((str(intervals[i]), str(intervals[j])))
        edges.extend(protein_edges)
    return edges

#formatting labels --> display the label nodes as "protein name (form aa,to aa) (form aa,to aa)""

def formatting_labels(s):
    parts = s.split(" ")
    result = [parts[0]]  
    temp_list = []

    for part in parts[1:]:  
        temp_list.append(part[:-1])
    for i in range(len(temp_list)-1):
        if i%2 == 0:
            temp_list[i]=temp_list[i].strip("(")
            if temp_list[i] == temp_list[i+1]:
                result.append(f"({temp_list[i] })")
            else:
                result.append(f"({temp_list[i]}-{temp_list[i+1]})")
    return " ".join(result)

#FUNCTION THAT ARE NEEDED IN THE MAIN FUCNTION FOR THE NETWORK WITH MERGING NODES

#FUNCTION FOR:

#1 MERGING INTERVALS (START,END)

def defining_intervals(interval1, interval2, min_value_for_merging):
    #case 1 cheeking if they are equal
    if interval1 == interval2:
        return interval1
    start_1, end_1 = interval1
    start_2, end_2 = interval2
    overlap =min(end_1,end_2)-max(start_1,start_2)+1
    #case 2 check if one intervals contains another
    #1st example wjere the interval1 contains the interval 2
    if start_1 <=start_2 and end_1 >= end_2:
        return interval1
    #2nd example where the interval2 contains the int 1
    if start_2 <=start_1 and end_2 >= end_1:
        return interval2
    #case 3 check if the two interval overlap with a length > of the min_value_for_merging value
    if overlap > min_value_for_merging:
        if not (end_1<start_2 or end_2<start_1): 
            return (min(start_1, start_2), max(end_1, end_2)) #tuple
    else:
        return interval1, interval2
    

#2 APPLY THE defining_intervals FUNCTION IN THE SMALL DF FOR A SINGLE PROTEIN

# UPDATE THE DF WITH NEW COLUMNS CONTANING THE MERGED INTERVALS


def get_merging_intervals(df, min_value_for_merging):
    for i in range(len(df['protein'])):
        interval_1 = (df.loc[i, 'start_merged'], df.loc[i, 'end_merged'])
        
        for j in range(len(df['protein'])): 
            if i == j: 
                continue
            
            interval_2 = (df.loc[j, 'start_merged'], df.loc[j, 'end_merged'])
            new_intervals = defining_intervals(interval_1, interval_2, min_value_for_merging)
            
            if isinstance(new_intervals, tuple) and all(isinstance(i, tuple) for i in new_intervals):
                pass
            else:
                df.at[i, 'start_merged'] = new_intervals[0]
                df.at[i, 'end_merged'] = new_intervals[1]
                df.at[j, 'start_merged'] = new_intervals[0]
                df.at[j, 'end_merged'] = new_intervals[1]
    
    return df



#3 MERGING INTERFACES (all the intervals that belong to the same interface are collected in a list)

#function that put  the intrval of the same interface in a list of list
def create_new_column_interface_intervals(df):
    # Group by "interface"
    grouped = df.groupby(["interface"])
    # grouped = df.groupby(["num_binary_combination", "interface"])
    # grouped = df.groupby(["numn_int", "interface"])
    
    # List for intervals for each interface
    interface_intervals = [None] * len(df)
    
    for _, group in grouped:
        # Create a list of tuples for the unique intervals
        intervals = [tuple(x) for x in group[["start_merged", "end_merged"]].drop_duplicates().values]
        
        # Assign the same intervals to every row in the group
        for idx in group.index:
            interface_intervals[idx] = intervals

    # Create a new column to store interface intervals as lists of tuples
    df["interface_intervals"] = interface_intervals

    # Create also a column for updating later the merged ones (now they have the same value!!!!)
    df["interface_intervals_merged"] = interface_intervals
    
    return df


#function that check if there are overlipping interfaces, here we updating the interface_intervals_merged column

def check_overlapping_interfaces(inte_1, inte_2):
    # check if every int in inte_1 are in inte_2
    if all(i in inte_2 for i in inte_1):
        return inte_2
    # check if every int in inte_2 are in inte_1
    if all(i in inte_1 for i in inte_2):
        return inte_1
    # check overlap
    if any(i in inte_2 for i in inte_1) or any(i in inte_1 for i in inte_2):
        #union of intervals (removing duplicates)
        union = inte_1 + [i for i in inte_2 if i not in inte_1]
        return sorted(union)
    #otherwise
    return None



'''MAIN FUNCTION '''
def get_protein_network_merging(df):
    #STEP 2 --> FROM the BIG dataframe TO create SMALL dayaframe for a sigle protein
    prot_names = set(df['prot_1']).union(set(df['prot_2']))
    df['original_index'] = df.index

    #dictionary for storing them
    protein_df_dict = {}
    
    for prot in prot_names:
        #whent the slected prtoein is in the first column
        prot1_rows = df[df['prot_1'] == prot][['interface','prot_1', 'start_1', 'end_1','prot_2','original_index']].rename(
            columns={'prot_1': 'protein', 'start_1': 'start', 'end_1': 'end', 'prot_2':'bin_comb_with'}
        )
        # when the selected protein is in the sec column
        prot2_rows = df[df['prot_2'] == prot][['interface','prot_2', 'start_2', 'end_2','prot_1','original_index']].rename(
            columns={'prot_2': 'protein', 'start_2': 'start', 'end_2': 'end','prot_1':'bin_comb_with'}
        )
        #combine the two rows information
        combined_rows = pd.concat([prot1_rows, prot2_rows], ignore_index=True)

          
        # Sort by 'original_index' to keep the rows in the correct order
        combined_rows = combined_rows.sort_values(by='original_index').reset_index(drop=True)
        #add the combined one as the value of the dictionary
        protein_df_dict[prot] = combined_rows

    #adding new columns for updating the intervals later (so we can keep track of the orginal ones)
    for prot, df in protein_df_dict.items():
            df['start_merged'] = df['start']
            df['end_merged'] = df['end']

    #print(f"protein_df_dict AFTER STEP 2 (from the big to the small df)-->{protein_df_dict}")

    # STEP 3 --> function for MERGIN INTERVALS in the SMALL dataframe

    for prot, df in protein_df_dict.items():
        df = get_merging_intervals(df,3)

    #print(f"protein_df_dict AFTER STEP 3(mergin intervals in the big dataframe)-->{protein_df_dict}")

    #STEP 4 --> all INTERVALS that belong to the SAME INTERFACE will be collected togheter in a list
    for prot, df in protein_df_dict.items():
        df = create_new_column_interface_intervals(df)

    #print(f"protein_df_dict AFTER STEP 4(merge intervals in interfaces intervals)-->{protein_df_dict}")

    #STEP 5 --> check if there are overlapping interfaces in the dataframe of the signle protein (here you check all togheter bc it could be)
    #that "different" interfaces in different binary interaction are actually overlapping
    #apply the function in the dataframe
    for prot, df in protein_df_dict.items():
        for i in range(len(df['protein'])):
            #print(f"i{i}")
            interval_1 = df.loc[i, 'interface_intervals_merged']
            for j in range(len(df['protein'])):
                #print(f"j:{j}")
                interval_2 = df.loc[j, 'interface_intervals_merged']
                new_interfaces = check_overlapping_interfaces(interval_1, interval_2)
                if new_interfaces == None:
                    pass
                else:
                    df.at[i,'interface_intervals_merged'] = new_interfaces
                    df.at[j,'interface_intervals_merged'] = new_interfaces

    #print(f"protein_df_dict AFTER STEP 5(check if there are overlapping interfaces)-->{protein_df_dict}")

    #STEP 6 --> obtaining again the output of alphabridge but merged

    #using the original index we have now a dataframe we the protein involoved in the interaction are one after the other ikn rw
    #ex --> first row (protein A) and second row (protein B) in the orginal df are interacting and they are in the same row
    merged_df = pd.concat(protein_df_dict.values(), axis=0).sort_values(by="original_index").reset_index(drop=True)
    #print(f"merged df:{merged_df}")

    # Creating a copy of the merged_df
    copied_df = merged_df.copy()

    # deletion of rows
    #in one df we will delete the rows from index 1 + step 2
    filtered_original = merged_df.drop(merged_df.index[1::2])  # delete 1, 3, 5, ...
    #print(f"filtered_original df:{filtered_original}")
    #in one df we will delete the rows from index 0 + step 2
    filtered_copy = copied_df.drop(copied_df.index[0::2])  # deete 0, 2, 4, ...
    #print(f"filtered_copy df:{filtered_copy}")


    # after we delete the row we can have merge the 2 df, based on the original index column
    # result --> output of AlphaBridge but with all the data processed
    final_df = pd.merge(filtered_original, filtered_copy, how="inner", on = 'original_index', suffixes=('_1', '_2'))
    final_df["interface_intervals_mergerd_1_labels"] = final_df["protein_1"]+ " " +final_df["interface_intervals_merged_1"].apply(lambda x: " ".join(map(str, x)) if isinstance(x, list) else str(x)).astype(str)
    final_df["interface_intervals_mergerd_1_labels"] = final_df["interface_intervals_mergerd_1_labels"].apply(formatting_labels)
    final_df["interface_intervals_mergerd_2_labels"] = final_df["protein_2"]+ " " +final_df["interface_intervals_merged_2"].apply(lambda x: " ".join(map(str, x)) if isinstance(x, list) else str(x)).astype(str)
    final_df["interface_intervals_mergerd_2_labels"] = final_df["interface_intervals_mergerd_2_labels"].apply(formatting_labels)


    #print(f"merged_protein_dataframe_dict  STEP 6(obtaining again the output of alphabridge but merged-->{final_df}")


    #START TO TAKE INFORMATION FOR PLOTTING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    # Initialize the dictionary to store proteins and their unique nodes
    protein_nodes_old = {}

    # Iterating through each row of the DataFrame
    
    for _, row in final_df.iterrows():
        # Extract protein IDs and their corresponding residue intervals
        prot_1_id = row['protein_1']
        prot_2_id = row['protein_2']
        prot_1_interface = row["interface_intervals_mergerd_1_labels"] 
        prot_2_interface = row["interface_intervals_mergerd_2_labels"] 

        # Initialize the protein in the dictionary if not already present
        if prot_1_id not in protein_nodes_old:
            protein_nodes_old[prot_1_id] = []
        if prot_2_id not in protein_nodes_old:
            protein_nodes_old[prot_2_id] = []

        # Add unique intervals for prot_1_id
        if prot_1_interface not in protein_nodes_old[prot_1_id]:  # Check if the interval is unique
            protein_nodes_old[prot_1_id].append(prot_1_interface)

        # Add unique intervals for prot_2_id
        if prot_2_interface not in protein_nodes_old[prot_2_id]:  # Check if the interval is unique
            protein_nodes_old[prot_2_id].append(prot_2_interface)

    #print(f'protein nodes old  {protein_nodes_old}')

    protein_nodes = dict(sorted(protein_nodes_old.items()))
    
    #print(f'protein nodes {protein_nodes}')


    # Calculate the number of unique nodes (intervals) for each protein
    nodes_for_each_protein = [len(value) for value in protein_nodes.values()]
    #print(f'nodes_for_each_protein {nodes_for_each_protein}')

    # INTERACTION BETWEEN DIFFERENT PROTEINS
    inter_protein_interactions = []

    # Iterate over each row in the dataframe
    for _, row in final_df.iterrows():
        # Extract data for prot_1 and prot_2
        #prot_1_id = row["protein_1"]
        aa_prot_1 = row["interface_intervals_mergerd_1_labels"]
        #prot_2_id = row["protein_2"]
        aa_prot_2 = row["interface_intervals_mergerd_2_labels"]

        # Create a list of interactions between protein pairs
        interaction = (aa_prot_1, aa_prot_2)
        inter_protein_interactions.append(interaction)

    # Print all protein interactions
    #print(f" different protein interactions: {inter_protein_interactions}")
    #print(len(inter_protein_interactions))

    unique_inter_protein_interactions = []   

    seen_interactions = set()  

    for interaction in inter_protein_interactions:
        sorted_interaction = tuple(sorted(interaction))  
        
        if sorted_interaction not in seen_interactions:
            unique_inter_protein_interactions.append(interaction)
            seen_interactions.add(sorted_interaction)  

    #print(unique_inter_protein_interactions)
    #print(f" unique_inter_protein_interactions: {unique_inter_protein_interactions}")

    #print(f"Unique protein interactions: {unique_inter_protein_interactions}")


    #INTERACTION WHITIN THE SAME PROTEIN
    edges_intraprpt = generate_intraprotein_edges(protein_nodes)
    #print(f'edges_intraprt {edges_intraprpt}')

    
    #PLOTTING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! THE NOT INTERACTIVE ONE
    g = ig.Graph()
    number_of_nodes = 0
    labels = []

    
    #NODES/LABELS
    #calculate the number of nodes that will be plotted
    for nodes in nodes_for_each_protein:
        number_of_nodes = number_of_nodes + nodes
        
    #labels--> each node has intervals of aa
    for key, value in protein_nodes.items():
        for interval_aa_str in value:
            labels.append(interval_aa_str)

    #print(f"labels: {labels}")
    #print(f"number of nodes:{number_of_nodes}")
    
    #adding nodes
    g.add_vertices(number_of_nodes)
    
    #assign labels to each node (the order follow the dictionary protein_nodes)
    g.vs['label'] = labels
    #for vertex in g.vs:     
    #    print(f"ID: {vertex.index}, Label: {vertex['label']}")


    #COLORS 
    # Generate as many unique colors as the number of protein groups
    num_proteins = len(nodes_for_each_protein)
    cmap = plt.get_cmap("rainbow") 
    colors = [mcolors.rgb2hex(cmap(i)) for i in np.linspace(0, 1, num_proteins)]
    
    # Assign colors to nodes based on their protein group
    node_colors = []
    for protein_index, size in enumerate(nodes_for_each_protein):
        node_colors.extend([colors[protein_index]] * size) 
        
    # Assign colors to graph nodes
    g.vs["color"] = node_colors
    #print(f"node colors{node_colors}")

    
    #EDGES
    # different interfaces
    edges_diff = []
    for source, target in unique_inter_protein_interactions:
        # Convert source and target to match the label format
        source_index = g.vs.find(label=f"{source}").index
        target_index = g.vs.find(label=f"{target}").index
        edges_diff.append((source_index, target_index))
    #print(f"edges diff:{edges_diff}")
    g.add_edges(edges_diff)
    
    # same protein edges
    edges_same = []
    for key, value in edges_intraprpt:
        # Convert source and target to match the label format
        source_index_1 = g.vs.find(label=f"{key}").index
        target_index_1 = g.vs.find(label=f"{value}").index
        edges_same.append((source_index_1, target_index_1))
    #print(f"edges sapython3 define_interfaces.py -i path-to-AlphaFold3-folderme:{edges_same}")
    g.add_edges(edges_same)
    #g.es['weight'] =[500.0] *g.ecount()
    g.es['weight'] = 500
    g.es['color'] = "black"
    g.es["width"] = 1.5
    g.es["label"] = "interface"
    #iterating all ossible edges
    for edge in edges_same:
        #for every tuple (=edge inside the list) we take the FROM node and the TO node
        source_index, target_index = edge
        #get ID f the edge
        edges_id = g.get_eid(source_index, target_index)
        
        # take the actual color of the node
        node_color = g.vs[source_index]["color"]
        g.es[edges_id]["weight"] = 3
        g.es[edges_id]["color"] = node_color  #we give for each edge the color of his node
        g.es[edges_id]["width"] = 0.5
        g.es[edges_id]["label"] = "intraface"



   #LAYAOUT

    #layout = g.layout("fruchterman_reingold")

    
   #PLOT 
    
    #fig, ax = plt.subplots(figsize=(10, 10))
    #ig.plot(g,target = ax, layout = layout, bbox =(800,800), margin = 20, vertex_color = g.vs['color'],vertex_size=10, edge_weight = g.es["weight"], vertex_label = None, color =g.es["color"], widht = g.es["widht"])


    #LEGEND 

    #PROTEINS for the legend
    #proteins_leg = []
    #for key, value in protein_nodes.items():
    #    proteins_leg.append(key)

        #creating the content of the first legend
    #nodes_patches = [mlines.Line2D([0], [0], marker='o', color='w',
    #                            markerfacecolor=colors[i], markersize=12,
    #                            label=proteins_leg[i]) for i in range(len(proteins_leg))]
    
    #creating the content of the second legend
    #edge_patches = [
    #    mpatches.Patch(color='lightgray', label='Same Protein Connections'),
    #    mpatches.Patch(color='black', label='Interface Connections')
    #]
   

    #first legend
    #first_legend = ax.legend(handles=nodes_patches, loc="upper right")
    #ax.add_artist(first_legend)
    
    #second legend
    #ax.legend(handles=edge_patches, loc="upper left", title=r"$\bf{Interfaces}$")

    #network in network x
    g_networkx = g.to_networkx()
    
    #informaton for the json file
    jobs = json_graph.node_link_data(g_networkx)

    return jobs
