#script for obtaining the protein interaction network
#EVERY NODE IS A PROTEIN
#EVERY EDGE IS A PHYSICAL INTERECTION BETWEEN THE PROTEINS

#import required packages

import igraph as ig
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
from networkx.readwrite import json_graph;
import src.network.obtaining_colors as col

#FUCNTION USED IN NETWORK WITH NOT MERGIN NODES AND THE PROTEIN==NODE NETWORK

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

#all the intervals that belong to the same interface --> togheter

def create_new_column_interface_intervals_no_merge(df):
    grouped = df.groupby(["interface"])
    
    # List for intervals for each interface (start_1, end_1)
    interface_intervals_1 = [None] * len(df)
    for _, group in grouped:
        # Convert to a list of tuples
        intervals = [tuple(x) for x in group[["start_1", "end_1"]].drop_duplicates().values]
        for idx in group.index:
            interface_intervals_1[idx] = intervals
    df["interface_intervals_1"] = interface_intervals_1

    grouped = df.groupby(["interface"])
    
    # List for intervals for each interface (start_2, end_2)
    interface_intervals_2 = [None] * len(df)
    for _, group in grouped:
        # Convert to a list of tuples
        intervals = [tuple(x) for x in group[["start_2", "end_2"]].drop_duplicates().values]
        for idx in group.index:
            interface_intervals_2[idx] = intervals
    df["interface_intervals_2"] = interface_intervals_2

    return df

'''MAIN FUNCTION!!!'''

 #STEP1--> HAVING THE DF WITH ALL POSSIBLE INFORMATION
def get_protein_network(df,label2auth,df_pairwise_interaction,threshold,auth2label):
    label_color_dict = col.get_label_color_dict(label2auth)
    df['prot_1_lab'] = df['prot_1']
    df['prot_2_lab'] = df['prot_2']
    for column in ['prot_1_lab', 'prot_2_lab']:
        df[column]= df[column].replace(label2auth)
    #print(df)
    # Initialize the dictionary to store proteins and their unique nodes
    protein_labels = {}

    # Iterating through each row of the DataFrame
    for _, row in df.iterrows():
        # Extract protein IDs and their corresponding residue intervals
        prot_1_id = row['prot_1']
        prot_2_id = row['prot_2']
        prot_1_interface = row["prot_1_lab"] 
        prot_2_interface = row["prot_2_lab"] 

        # Initialize the protein in the dictionary if not already present
        if prot_1_id not in protein_labels:
            protein_labels[prot_1_id] = []
        if prot_2_id not in protein_labels:
            protein_labels[prot_2_id] = []

        # Add unique intervals for prot_1_id
        if prot_1_interface not in protein_labels[prot_1_id]:  
            protein_labels[prot_1_id].append(prot_1_interface)

        # Add unique intervals for prot_2_id
        if prot_2_interface not in protein_labels[prot_2_id]:  
            protein_labels[prot_2_id].append(prot_2_interface)

    #print(f'Protein labels before formatting: {protein_labels}')

    # Format labels correctly
    for protein in protein_labels:
        intervals = protein_labels[protein]
        #print(protein)
        
        # take out the name of the protein in the other string (I only want it once)
        formatted_intervals = [intervals[0]]  # only the first interval keep he name
        #print(formatted_intervals)
        #for the rest of them:
        for interval in intervals[1:]:
            #replace protein == name of the proten with nothing
            formatted_intervals.append(interval.replace(f"{protein} ", "")) 
            #print(formatted_intervals)
        
        # formatted labels --> sring separated by , --> joing them in an unique one
        protein_labels[protein] = ", ".join(formatted_intervals)

    #print(f'Protein labels after formatting: {protein_labels}')

    #orering them --> for having the same color as the other network
    protein_labels_sorted = dict(sorted(protein_labels.items()))

    #print(f'Protein labels after sorting: {protein_labels_sorted}')

    # INTERACTIONS with name of the protein
    protein_pairs = []
    protein_pairs_unique = []

    for prot1, prot2 in zip(df["prot_1"], df["prot_2"]):
        pair = [prot1, prot2]
        protein_pairs.append(pair)

    for protein in protein_pairs:
        if protein not in protein_pairs_unique:
            protein_pairs_unique.append(protein)

    #print(protein_pairs)
    #print(protein_pairs_unique)

    #print("All interactions:", protein_pairs)
    #print("Unique interactions:", protein_pairs_unique)

    # INTERACTIONS with the all label of the protein
    protein_pairs_unique_labels = []

    for prot1, prot2 in protein_pairs_unique:
        labels_1 = protein_labels.get(prot1, "No Label")
        labels_2 = protein_labels.get(prot2, "No Label")
        protein_pairs_unique_labels.append([labels_1, labels_2])

    #print("Unique interactions with labels:", protein_pairs_unique_labels)
    g = ig.Graph()
    number_of_nodes = (len(protein_labels_sorted))
    labels = []
    #adding nodes
    g.add_vertices(number_of_nodes)

    #assign labels
    for key, value in protein_labels_sorted.items():
        labels.append(value)

    g.vs['label'] = labels

    #for vertex in g.vs:     
    #    print(f"ID: {vertex.index}, Label: {vertex['label']}")

    #COLORS 
    # Generate as many unique colors as the number of protein groups
    #cmap = plt.get_cmap("rainbow") 
    #colors = [mcolors.rgb2hex(cmap(i)) for i in np.linspace(0, 1, number_of_nodes)]
        
    # Assign colors to graph nodes
    colors = []

    for v in g.vs:
        label = v["label"]
        found_color = "gray"  #dafault color

        for key_letter, color in label_color_dict.items():
            if label.startswith(key_letter):
                found_color = color
                break 

        colors.append(found_color)

    g.vs["color"] = colors

    #print(f"node colors{colors}")

    #edges
    edges = []
    for source, target in protein_pairs_unique_labels: #Unique interactions with labels: [['A', 'B'], ['A', 'lig_C'], ['B', 'DNA_E'], ['B', 'DNA_F'], ['DNA_E', 'DNA_F']]
        # Convert source and target to match the label format
        source_index = g.vs.find(label=f"{source}").index #form the label give me the index 
        #print(f"source:{source_index}")
        target_index = g.vs.find(label=f"{target}").index
        #print(f"target:{target_index}")
        edges.append((source_index, target_index))

    #print(f"edges :{edges}")
    g.add_edges(edges)
    g.es['color'] = "black"
    unique_id_dict = {}
    #adding the unique identifier for the interface as the attribute type for the edges/links of the network
    #creating a dictionary for storing them
    for _, row in df.iterrows():
        from_prefix_df = row['prot_1']
        to_prefix_df = row['prot_2']
        unique_identifier = row['interface_id']

        key = (from_prefix_df, to_prefix_df)

        if key not in unique_id_dict:
            unique_id_dict[key] = []

        if unique_identifier not in unique_id_dict[key]:  # avoid duplicates
            unique_id_dict[key].append(unique_identifier)
    for edge in edges:  
        source_index, target_index = edge

        source_label_raw = g.vs[source_index]["label"]
        target_label_raw = g.vs[target_index]["label"]

        source_label = auth2label.get(source_label_raw, source_label_raw)
        target_label = auth2label.get(target_label_raw, target_label_raw)

        edge_key = (source_label, target_label)

        if edge_key in unique_id_dict:
            edge_id = g.get_eid(source_index, target_index)
            interface_ids = unique_id_dict[edge_key]

            g.es[edge_id]["interaction"] = ",".join(interface_ids)

    #print(unique_id_dict)
    #g.es["label"] = "physical interaction"
    for _, row in df_pairwise_interaction.iterrows():
        from_prefix = row['first']
        to_prefix = row['second']
        score = row['pairwise_score'] 
        #print(from_prefix, to_prefix,score)

        for edge in edges:  
            source_index, target_index = edge
        
            source_label = g.vs[source_index]["label"]
            target_label = g.vs[target_index]["label"]
        
            if from_prefix in source_label and to_prefix in target_label:
                edge_id = g.get_eid(source_index, target_index)
                if score >= threshold:
                    g.es[edge_id]["label"] = score
                    #edges_labels =g.es[edge_id]["label"]
                    #print(g.es[edge_id]["label"])

    #network in network x
    g_networkx = g.to_networkx()
    
    #informaton for the json file
    jobs = json_graph.node_link_data(g_networkx)
    for link in jobs['links']:
        link['interaction'] = link['interaction'].split(",")

    return jobs

