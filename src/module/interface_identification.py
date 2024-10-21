# -*- coding: utf-8 -*-
"""
Created on Wed Feb 7 16:53:05 2024

@author: Dan_salv
"""
import numpy as np
from itertools import chain, repeat, count, islice
from scipy import ndimage
from collections import Counter
import pandas as pd
import networkx 
from networkx.algorithms.components.connected import connected_components

from src.module.parsers import MMCIFPARSER
from src.module.alingment_utils import compare_protein_seq
from src.module.scoring import calculate_probability_contact_link, calculate_probability_contact_interface, plot_probability_histplot

from Bio.Data import IUPACData


protein_letters_1to3 = IUPACData.protein_letters_1to3

upper_protein_letters_1to3 = {k.upper():v.upper() for k,v in protein_letters_1to3.items()}

upper_protein_letters_1to3


class interface_identification():
    
    def __init__(self,
                 interacting_coevolutionary_domains, 
                 entity_region_dict, 
                 plddt_dict, 
                 rec_sequence_list, 
                 list_sequence_info, 
                 contact_matrix,
                 threshold, 
                 chain_dict, 
                 polymer_chain_dict):
        
        self.interacting_coevolutionary_domains = interacting_coevolutionary_domains
        self.entity_region_dict = entity_region_dict 
        self.plddt_dict = plddt_dict
        self.rec_sequence_list = rec_sequence_list
        self.list_sequence_info = list_sequence_info
        self.contact_matrix  = contact_matrix
        self.threshold = threshold
        self.chain_dict = chain_dict
        self.polymer_chain_dict = polymer_chain_dict
        
        
    def extract_contacts(self):
        data = []
        
        interacting_coevolutionary_domains = self.interacting_coevolutionary_domains
        contact_matrix = self.contact_matrix
        entity_region_dict = self.entity_region_dict
        
        
        cluster_group_name_list = list(interacting_coevolutionary_domains.keys())

        for cluster_group_name in cluster_group_name_list:
            
            overlap_cluster = interacting_coevolutionary_domains[cluster_group_name]['overlap_complex']
            proteins_interacting_list = list(overlap_cluster.keys())
            if len(proteins_interacting_list) >= 2:
                complex_combinations = list(unique_combinations(proteins_interacting_list, 2))


                for protA, protB in complex_combinations:
                    
                    protA_range_list = overlap_cluster[protA]
                    protB_range_list = overlap_cluster[protB]
                    
                    
                    for protA_range in protA_range_list:
                        for protB_range in protB_range_list:
                            
                            distance_submatrix = contact_matrix[protA_range[0]:protA_range[1], protB_range[0]:protB_range[1]]
                            interfaces, interface_range_list = find_interfaces(distance_submatrix, self.threshold)
                            
                            if not len(interface_range_list) == 0:
                            
                                for interface_range in interface_range_list:
                                    
                                    #interaction_dimension = interfaces[interface_range].shape
                                    
                                    protA_start, protA_end = map_residue_range(entity_region_dict[protA], protA_range, interface_range[0])
                                    protB_start, protB_end = map_residue_range(entity_region_dict[protB], protB_range, interface_range[1])
                                    
                                    data.append([cluster_group_name,protA, protA_start, protA_end, protB, protB_start, protB_end])
        
        contact_df = pd.DataFrame(data, columns=['cluster_group_name','prot_1','start_1', 'end_1' ,'prot_2','start_2', 'end_2'])
        contact_df['interfaces'] = contact_df.groupby(by=['cluster_group_name','prot_1','prot_2'], as_index=False).ngroup() +1

        contact_df['interfaces'] = contact_df['interfaces'].apply(lambda x: f'interface_{x}')

        contact_df = contact_df.reindex(
            columns = ['cluster_group_name','interfaces','prot_1', 'start_1', 'end_1', 'prot_2', 'start_2', 'end_2'])
        
        return contact_df
    
    
    def extract_interface(self):
        
        contact_df = self.extract_contacts()
        chain_dict = self.chain_dict
        entity_region_dict = self.entity_region_dict
        rec_sequence_list = self.rec_sequence_list
        contact_matrix = self.contact_matrix
        
        interaction_link_dict = {}
        for fasta_name, seq in rec_sequence_list:
            if not fasta_name in interaction_link_dict:
                interaction_link_dict[fasta_name] = []
            
        
        interface_dict = {}
        protein_interface_dict = {}
        interface_count = 0
        
        interface_group = contact_df.groupby(['interfaces'])
        
        
        for _, df_group in interface_group:
            
            ranges_prot_1 = mergeIntervals([list(range_tuple) for range_tuple in zip(df_group.start_1, df_group.end_1)])
            ranges_prot_2 = mergeIntervals([list(range_tuple) for range_tuple in zip(df_group.start_2, df_group.end_2)])
            
            interface_list, interface_link_list = merge_split_interfaces(ranges_prot_1,ranges_prot_2,df_group)
            
            proteins_involved = [df_group.prot_1.unique()[0],df_group.prot_2.unique()[0]]
            
            for interface, interface_link in zip(interface_list, interface_link_list):
                matrix_probability_interface = []
                link_id = 0
                
                interface_name = f'interface_{interface_count}'
                interface_count += 1
                
                if not interface_name in interface_dict:
                    interface_dict[interface_name] = {'prot_1':{'accesion_id':df_group.prot_1.unique()[0], 'chain' : chain_dict[df_group.prot_1.unique()[0]], 'interface_range':interface[0]},
                                                'prot_2':{'accesion_id':df_group.prot_2.unique()[0], 'chain' : chain_dict[df_group.prot_2.unique()[0]], 'interface_range':interface[1]},
                                                'links':interface_link,
                                                'interface_prob': float()}

                for link in interface_link:
                    link_id +=1
                    coord_link = map_link2coord(link, proteins_involved, entity_region_dict)
                    matrix_probability_link, link_probability = calculate_probability_contact_link(coord_link,contact_matrix)
                    matrix_probability_interface.append(matrix_probability_link)
                    
                    for residue_range, prot in zip(link,proteins_involved):
                        if not prot in interaction_link_dict:
                            interaction_link_dict[prot]= []
                        
                        link_data = (residue_range, interface, link_id, link_probability)
                        interaction_link_dict[prot].append(link_data)
                
                
                probability_contact_interface, flatten_matrix_probability_interface = calculate_probability_contact_interface(matrix_probability_interface)
                interface_dict[interface_name]['interface_prob'] = probability_contact_interface
                
                print(df_group.prot_1.unique()[0],df_group.prot_2.unique()[0])
                #plot_probability_histplot(probability_contact_interface,flatten_matrix_probability_interface)
                
                
                for index, protein_involved in enumerate(proteins_involved):
                    if not protein_involved in protein_interface_dict:
                        protein_interface_dict[protein_involved] = {}
                    if not interface_name in protein_interface_dict[protein_involved]:
                        protein_interface_dict[protein_involved][interface_name] = {'interface_range':interface[index]}
                
        return interface_dict, protein_interface_dict, interaction_link_dict
         
        
    def get_interface_info_dataframes(self, interface_dict, interaction_link_dict):
        
        rec_sequence_list = self.rec_sequence_list
        list_sequence_info = self.list_sequence_info
        chain_dict = self.chain_dict
        polymer_chain_dict = self.polymer_chain_dict
        
        plddt_dict = self.plddt_dict
        
        interface_df_per_token = get_interface_df_per_token(rec_sequence_list, list_sequence_info, chain_dict, polymer_chain_dict, interaction_link_dict, plddt_dict)
        interface_df = get_interface_df(interface_dict)
        
        return interface_df_per_token, interface_df


        

def repeat_chain(values, counts):
    return chain.from_iterable(map(repeat, values, counts))


def unique_combinations_from_value_counts(values, counts, r):
    n = len(counts)
    indices = list(islice(repeat_chain(count(), counts), r))
    if len(indices) < r:
        return
    while True:
        yield tuple(values[i] for i in indices)
        for i, j in zip(reversed(range(r)), repeat_chain(reversed(range(n)), reversed(counts))):
            if indices[i] != j:
                break
        else:
            return
        j = indices[i] + 1
        for i, j in zip(range(i, r), repeat_chain(count(j), counts[j:])):
            indices[i] = j


def unique_combinations(iterable, r):
    values, counts = zip(*Counter(iterable).items())
    return unique_combinations_from_value_counts(values, counts, r)


def find_interfaces(matrix, threshold):
    
    binary = matrix > threshold
    
    centrosymmetric_matrix = ndimage.generate_binary_structure(2,2)
        
    interfaces, number_of_interfaces = ndimage.label(binary, structure = centrosymmetric_matrix)

    interface_range_list  = ndimage.find_objects(interfaces)
    
    return interfaces, interface_range_list

def map_residue_range(protein_region, submatrix_range, iteraction_range):
    
    delta_index = submatrix_range[0] - protein_region[0]
   
    start = iteraction_range.start + delta_index + 1
    end = iteraction_range.stop  +  delta_index 
    
    return start, end


def fix_intervals(overlap_list):
    index = 0
    for i in range(1, len(overlap_list)):
        delta = overlap_list[i][0] - overlap_list[i-1][1]
        
        if delta <= 2:
            overlap_list[index][1] = max(overlap_list[index][1], overlap_list[i][1])
            
        else:
            index = index + 1
            overlap_list[index] = overlap_list[i]
    
    overlap_list = overlap_list[:index+1]
    
    return overlap_list

def mergeIntervals(arr):
 
    # Sorting based on the increasing order
    # of the start intervals
    arr.sort(key=lambda x: x[0])
 
    # Stores index of last element
    # in output array (modified arr[])
    index = 0
 
    # Traverse all input Intervals starting from
    # second interval
    for i in range(1, len(arr)):
        
        # If this is not first Interval and overlaps
        # with the previous one, Merge previous and
        # current Intervals
        if (arr[index][1] >= arr[i][0]):
            arr[index][1] = max(arr[index][1], arr[i][1])
        else:
            index = index + 1
            arr[index] = arr[i]
    
    overlap_list = fix_intervals(arr[:index+1])
    
    return overlap_list


def range_subset(range1, range2):
    """Whether range1 is a subset of range2."""
    if not range1:# empty range is subset of anything
        if range1.start == range1.stop:
            return range2.start <= range1.start <= range2.stop
        else:
            raise ValueError
        
    if not range2:
        return False  # non-empty range can't be subset of empty range
    if len(range1) > 1 and range1.step % range2.step:
        return False  # must have a single value or integer multiple step
    return range1.start in range2 and range1[-1] in range2

def common_member(a, b):
    a_set = set(a)
    b_set = set(b)
    if (a_set & b_set):
        return True
    else:
        return False
    
def to_graph(l):
    G = networkx.Graph()
    for part in l:
        # each sublist is a bunch of nodes
        G.add_nodes_from(part)
        # it also imlies a number of edges:
        G.add_edges_from(to_edges(part))
    return G

def to_edges(l):
    """ 
        treat `l` as a Graph and returns it's edges 
        to_edges(['a','b','c','d']) -> [(a,b), (b,c),(c,d)]
    """
    it = iter(l)
    last = next(it)

    for current in it:
        yield last, current
        last = current    

    
def extract_indexes_interfaces(unclustered_matrix):
    
    clusters = [np.where(row==1)[0] for row in unclustered_matrix]
    
    G = to_graph(clusters)
    index_interfaces = []
    for column_component in connected_components(G):
        row_list = set()
        for row, cluster in enumerate(clusters):
            
            if common_member(cluster, list(column_component)):
                #print(row,cluster,list(column_component))
                row_list.add(row)
        index_interfaces.append((list(row_list),list(column_component)))
    
    return index_interfaces
        
def interface2links(unclustered_matrix, index_interface, ranges_prot_1,ranges_prot_2):
    interface_links = []
    
    links_index = np.argwhere(unclustered_matrix==1)
    for link_index in links_index:
        if link_index[0] in index_interface[0] and link_index[1] in index_interface[1]:
            
            link = [ranges_prot_1[link_index[0]], ranges_prot_2[link_index[1]]]
            interface_links.append(link)
    
    
    return interface_links
        
     
def merge_split_interfaces(ranges_prot_1,ranges_prot_2,df_group):
    interface_list = []
    interface_link_list = []
    unclustered_matrix = np.zeros([len(ranges_prot_1), len(ranges_prot_2)], dtype=int)
    for i,range_1 in enumerate(ranges_prot_1):
        for j,range_2 in enumerate(ranges_prot_2):
            for start_1, end_1, start_2, end_2 in zip(df_group.start_1, df_group.end_1, df_group.start_2, df_group.end_2):
                if range_subset(range(start_1,end_1), range(range_1[0], range_1[1])) and range_subset(range(start_2,end_2), range(range_2[0], range_2[1])):
                    unclustered_matrix[i][j] = 1

    
    indexes_interfaces = extract_indexes_interfaces(unclustered_matrix)
    
    
    for index_interface in indexes_interfaces:
        
        range_interface_1 =  [ranges_prot_1[i] for i in index_interface[0]]
        range_interface_2 =  [ranges_prot_2[j] for j in index_interface[1]]
        
        interface_links = interface2links(unclustered_matrix, index_interface, ranges_prot_1,ranges_prot_2)
        
        interface_link_list.append(interface_links)
        interface_list.append([range_interface_1,range_interface_2])
        
    
    return interface_list, interface_link_list


def map_link2coord(link, binary_prot, entity_region_dict):
    
    coord_link = [list(np.array(link_terminal) + entity_region_dict[accesion_id][0] - 1) for accesion_id, link_terminal in zip(binary_prot, link)]
    return coord_link

def get_interface_and_link_per_residue(res,accesion_id,interaction_link_dict):
    link_data_list = interaction_link_dict[accesion_id]
    
    link_interface_list = []
    
    for residue_range, interface, link, link_probability in link_data_list:
        
        if residue_range[0] <= res <= residue_range[1]:
            
            link_interface = (interface, int(link),link_probability)
            
            link_interface_list.append(link_interface)
    
    if not link_interface_list:
        link_interface_list = [(None, None,None)]
           
    return link_interface_list


def get_interface_df_per_token(rec_sequence_list, list_sequence_info, chain_dict, polymer_chain_dict, interaction_link_dict, plddt_dict):
            
        list_fasta_name, list_fasta_acclen, list_fasta_centerticks, list_fasta_len = tuple(list_sequence_info)
        fasta_sequence_dict = {}
        
        for fasta_name, seq in rec_sequence_list:
            if not fasta_name in fasta_sequence_dict:
                fasta_sequence_dict[fasta_name] = ''
            fasta_sequence_dict[fasta_name] = seq
        
        res_link_data = []

        for fasta_name , prot_len in zip(list_fasta_name, list_fasta_len):
            
            for res in range(prot_len):
                
                
                if polymer_chain_dict[chain_dict[fasta_name]]['entity_type'] == 'polypeptide(L)':
                    comp_id = upper_protein_letters_1to3[fasta_sequence_dict[fasta_name][res]]
                
                else:
                
                    comp_id = fasta_sequence_dict[fasta_name][res]
                
                pldtt = plddt_dict[fasta_name][res]
                
                interface_link_list = get_interface_and_link_per_residue(res + 1,fasta_name,interaction_link_dict)
                
                for interface_nr, link_id,link_probability in interface_link_list:
                    
                    row = [fasta_name, chain_dict[fasta_name] ,res + 1, comp_id,interface_nr,link_id,pldtt,link_probability]
                
                    res_link_data.append(row)
        
        interface_df_per_token = pd.DataFrame(res_link_data, columns=['prot_name','asym_id','res_index','comp_id','interface', 'link_id' , 'plddt', 'link_probability'])
        
        return interface_df_per_token



def get_interface_df(interface_dict):

    interaction_link_data = []
    for interface in interface_dict:
        links =  interface_dict[interface]['links']
        for link in links:
            
        
            
            link_prot_1 = link[0]
            link_prot_2 =  link[1]
            
            

            accesion_id_1 = interface_dict[interface]['prot_1']['accesion_id']
            accesion_id_2 = interface_dict[interface]['prot_2']['accesion_id']
                    
                
            row = [interface, accesion_id_1,link_prot_1[0], link_prot_1[1],accesion_id_2,link_prot_2[0], link_prot_2[1]]
        
            interaction_link_data.append(row)  
    
    interface_info_df = pd.DataFrame(data=interaction_link_data, columns=['interface','prot_1','start_1','end_1','prot_2','start_2','end_2'])
    
    return interface_info_df