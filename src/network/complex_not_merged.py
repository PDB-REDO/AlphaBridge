#script for obtaining ONLY THE INTERACTIVE NETWORK WITH NO MERGE(more nodes)
#import required packages
import igraph as ig
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
from networkx.readwrite import json_graph;

#FUCNTION USED IN NETWORK WITH MERGING NODES AND NOT MERGIN NODES

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
def get_protein_network_no_merging(df,dit):

    #STEP1--> HAVING THE DF WITH ALL POSSIBLE INFORMATION
    #df_interfaces --> DF THAT IS USED FOR OBTAINING THE NO MERGED ETWORK (==MORE NODES)
    # .apply(lambda x: ...) to convert lists to space-separated strings.

    df_interfaces = create_new_column_interface_intervals_no_merge(df)
    df['prot_1_lab'] = df['prot_1']
    df['prot_2_lab'] = df['prot_2']
    for column in ['prot_1_lab', 'prot_2_lab']:
        df[column]= df[column].replace(dit)
    # Save the DataFrame with the updated protein names
    df_interfaces["interface_intervals_1_labels"] = df["prot_1_lab"]+" "+ df_interfaces["interface_intervals_1"].apply(lambda x: " ".join(map(str, x)) if isinstance(x, list) else str(x)).astype(str)
    df_interfaces["interface_intervals_1_labels"] = df_interfaces["interface_intervals_1_labels"].apply(formatting_labels)
    df_interfaces["interface_intervals_2_labels"] = df["prot_2_lab"]+ " "+ df_interfaces["interface_intervals_2"].apply(lambda x: " ".join(map(str, x)) if isinstance(x, list) else str(x)).astype(str)
    df_interfaces["interface_intervals_2_labels"] = df_interfaces["interface_intervals_2_labels"].apply(formatting_labels)

    
    #STE2 -> START TO TAKE INFORMATION FOR PLOTTING

    # Initialize the dictionary to store proteins and their unique nodes
    protein_nodes = {}

    # Iterating through each row of the DataFrame
    
    for _, row in df_interfaces.iterrows():
        # Extract protein IDs and their corresponding residue intervals
        prot_1_id = row['prot_1']
        prot_2_id = row['prot_2']
        prot_1_interface = row["interface_intervals_1_labels"] 
        prot_2_interface = row["interface_intervals_2_labels"] 

        # Initialize the protein in the dictionary if not already present
        if prot_1_id not in protein_nodes:
            protein_nodes[prot_1_id] = []
        if prot_2_id not in protein_nodes:
            protein_nodes[prot_2_id] = []

        # Add unique intervals for prot_1_id
        if prot_1_interface not in protein_nodes[prot_1_id]:  # Check if the interval is unique
            protein_nodes[prot_1_id].append(prot_1_interface)

        # Add unique intervals for prot_2_id
        if prot_2_interface not in protein_nodes[prot_2_id]:  # Check if the interval is unique
            protein_nodes[prot_2_id].append(prot_2_interface)

    #print(f'protein nodes NO SORTED {protein_nodes}')

    #sort the protein in order so the olor will be the same in the network
    protein_nodes_sorted = dict(sorted(protein_nodes.items()))
    #print(f'protein nodes sorted {protein_nodes_sorted}')

    # Calculate the number of unique nodes (intervals) for each protein
    nodes_for_each_protein = [len(value) for value in protein_nodes_sorted.values()]

    # INTERACTION BETWEEN DIFFERENT PROTEINS
    inter_protein_interactions = []

    # Iterate over each row in the dataframe
    for _, row in df_interfaces.iterrows():
        # Extract data for prot_1 and prot_2
        #prot_1_id = row["protein_1"]
        aa_prot_1 = row["interface_intervals_1_labels"]
        #prot_2_id = row["protein_2"]
        aa_prot_2 = row["interface_intervals_2_labels"]

        # Create a list of interactions between protein pairs
        interaction = (aa_prot_1, aa_prot_2)
        inter_protein_interactions.append(interaction)

    # Print all protein interactions
    #print(f"All protein interactions: {inter_protein_interactions}")

    # To get unique interactions, use a set
    unique_inter_protein_interactions = []

    for interaction in inter_protein_interactions:
        if interaction not in unique_inter_protein_interactions:
            unique_inter_protein_interactions.append(interaction)

    # Print unique protein interactions
    #print(f"Unique protein interactions: {unique_inter_protein_interactions}")


    #INTERACTION WHITIN THE SAME PROTEIN
    edges_intraprpt = generate_intraprotein_edges(protein_nodes_sorted)
    #print(f'edges_intraprt {edges_intraprpt}')

    
    #STEP 9 --> PLOTTING
    g = ig.Graph()
    number_of_nodes = 0
    labels = []

    
    #NODES/LABELS
    #calculate the number of nodes that will be plotted
    for nodes in nodes_for_each_protein:
        number_of_nodes = number_of_nodes + nodes
        
    #labels--> each node has intervals of aa
    for key, value in protein_nodes_sorted.items():
        for interval_aa_str in value:
            labels.append(interval_aa_str)

    #print(f"labels: {labels}")
    #print(f"number of nodes:{number_of_nodes}")
    
    #adding nodes
    g.add_vertices(number_of_nodes)
    
    #assign labels to each node (the order follow the dictionary protein_nodes_sorted)order
    g.vs['label'] = labels
    #for vertex in g.vs:     
     #   print(f"ID: {vertex.index}, Label: {vertex['label']}")


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
   # print(f"node colors{node_colors}")

    
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
    #g.es["label"] = "interface"
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
        #g.es[edges_id]["label"] = "intraface"


   #LAYAOUT

    #layout = g.layout("fruchterman_reingold")

    
   #PLOT 
    
    #fig, ax = plt.subplots(figsize=(10, 10))
    #ig.plot(g,target = ax, layout = layout, bbox =(800,800), margin = 20, vertex_color = g.vs['color'],vertex_size=10, edge_weight = g.es["weight"], vertex_label = None, color =g.es["color"], widht = g.es["widht"])


    #LEGEND 

    #PROTEINS for the legend
    #proteins_leg = []
    #for key, value in protein_nodes_sorted.items():
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

#g = get_protein_network_no_merging('/local_data/carmen/complexes_trial_old/complex1_pfdn1_vbp1_pfdn4_pfdn5', 0.75)
#print(g)