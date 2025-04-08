import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from itertools import chain

def calculate_probability_contact_link(coord_link, contact_probability):
    matrix_probability_link = contact_probability[np.s_[coord_link[0][0]:coord_link[0][1]+1], np.s_[coord_link[1][0]:coord_link[1][1]+1]]
    
    link_probability = np.mean([prob for prob in matrix_probability_link.flatten() if prob >= 0])
    return matrix_probability_link, link_probability

def calculate_pmc_link(coord_link, pae_plddt):
    matrix_pmc_link = pae_plddt[np.s_[coord_link[0][0]:coord_link[0][1]+1], np.s_[coord_link[1][0]:coord_link[1][1]+1]]
    
    norm_matrix_pmc_link = 1 - ( 3 * matrix_pmc_link / (2 * 100) )
    pmc_link = np.mean(norm_matrix_pmc_link.flatten())
    
    return norm_matrix_pmc_link, pmc_link
    
def calculate_probability_contact_interface(matrix_probability_interface):
    
    flatten_matrix_probability_interface = []
    filtered_flatten_matrix_probability_interface = []
    
    for matrix_probability_link in matrix_probability_interface:
        link_probability = [prob for prob in matrix_probability_link.flatten() if prob > 0]
        flatten_matrix_probability_interface += list(matrix_probability_link.flatten())
        filtered_flatten_matrix_probability_interface += link_probability
    
    quantile = np.quantile(filtered_flatten_matrix_probability_interface, 0.75)
    
    interface_probability = np.mean([prob for prob in filtered_flatten_matrix_probability_interface if prob >= quantile])
    
    contact_nr = len([prob for prob in flatten_matrix_probability_interface if prob >= quantile])
    #contact_ratio = contact_nr / len(flatten_matrix_probability_interface)
    
    return interface_probability, contact_nr , flatten_matrix_probability_interface

def calculate_pmc_interface(matrix_pmc_interface):
    
    flatten_matrix_pmc_interface = []
    
    for pmc_link in matrix_pmc_interface:
        flatten_matrix_pmc_interface += list(pmc_link.flatten())
    
    pmc_interface = np.mean(flatten_matrix_pmc_interface)
    
    return pmc_interface, flatten_matrix_pmc_interface
    
def calculate_interface_scores(interface_probability, pmc_interface, chain_pair_iptm):
    
    interface_score_iptm = interface_probability * chain_pair_iptm
    interface_score_pmc = interface_probability * pmc_interface
    
    return interface_score_iptm, interface_score_pmc

def calculate_scores_stucture(probability_structure_list, pmc_structure_list, iptm):
    
    probability_contact_structure = float()
    pmc_structure = float()
    structure_score_iptm = float()
    structure_score_pmc = float()

    if probability_structure_list and pmc_structure_list:

        flattened_probability_structure = []
        flattened_pmc_structure_list = []
        for contact_submatrix, pmc_submatrix in zip(probability_structure_list, pmc_structure_list):
            
            flattened_probability_structure += [prob for prob in contact_submatrix.flatten() if prob > 0]
            
            norm_pmc_submatrix = 1 - ( 3 * pmc_submatrix / (2 * 100) )
            flattened_pmc_structure_list += list(norm_pmc_submatrix.flatten())
        
        
        quantile = np.quantile(flattened_probability_structure, 0.75)
        #print(quantile)
        probability_contact_structure = np.mean([prob for prob in flattened_probability_structure if prob >= quantile])
        
        pmc_structure = np.mean(flattened_pmc_structure_list)
        
        structure_score_iptm = probability_contact_structure * iptm
        structure_score_pmc = probability_contact_structure * pmc_structure
        
    #plot_probability_histplot(quantile,probability_contact_structure, flattened_probability_structure, flattened_pmc_structure_list)
    
    
    return probability_contact_structure, pmc_structure, structure_score_iptm, structure_score_pmc


def plot_probability_histplot(quantile,probability_contact_structure,flattened_probability_structure, flattened_pmc_structure_list):
    
    #print(sorted(flatten_matrix_probability_interface),interface_probability, quantile)
    fig, ax = plt.subplots()
    g = sns.histplot(data=flattened_probability_structure, binwidth=0.01, cumulative = True, fill=False, stat='density')
    plt.axhline(y = 0.5, color = 'b')
    plt.axvline(x = probability_contact_structure, color = 'r')
    ax.set(xlim=(0,1))
    fig, ax = plt.subplots()
    g = sns.histplot(data=flattened_probability_structure, binwidth=0.01, cumulative = True, fill=False, stat='density')
    plt.axhline(y = 0.75, color = 'b')
    plt.axvline(x = probability_contact_structure, color = 'r')
    ax.set(xlim=(0,1))
    
    fig, ax = plt.subplots()
    g = sns.histplot(data=flattened_probability_structure, binwidth=0.01, fill=False, stat='density')
    plt.axvline(x = quantile, color = 'b')
    plt.axvline(x = probability_contact_structure, color = 'r')
    ax.set(xlim=(0,1))
    plt.show()
