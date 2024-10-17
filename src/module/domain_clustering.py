# -*- coding: utf-8 -*-
"""
Created on Wed Feb 7 16:53:05 2024

@author: Dan_salv
"""
import os
import numpy as np
import json
import igraph
import seaborn as sns
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import itertools

class domain_clustering():
    
    def __init__(self, matrix_dict, list_fasta_files, plotting: bool = True, bool_mask_clusters: bool = True, outdir: str = '' , alphafold_version: str = 'AF2'):
        
        self.matrix_dict = matrix_dict
        self.list_fasta_files = list_fasta_files
        self.plotting = plotting
        self.outdir = outdir
        self.alphafold_version = alphafold_version
        self.bool_mask_clusters = bool_mask_clusters
        
    def run_domain_clustering(self):
        matrix_dict = self.matrix_dict
        list_fasta_files = self.list_fasta_files
        
        matrix_input = matrix_dict['confidance_matrix']
        masked_confidance_matrix = matrix_dict['masked_confidance_matrix']
        masked_contact_matrix = matrix_dict['masked_contact_matrix']
        
        
        if self.alphafold_version == 'AF2':
            graph_resolution = 0.5
            matrix_cutoff = 2.6
            cmap_confidance = 'RdPu_r'
            cmap_contact = 'Blues_r'

            coevolutionary_domains = get_coevolutionary_domains(matrix_input, graph_resolution = graph_resolution, matrix_cutoff = matrix_cutoff)\
            
            interacting_coevultionary_cluster_dict, entity_region_dict = self.get_interacting_coevolutionary_domains(coevolutionary_domains)
            
            interacting_mask_cluster = self.get_interacting_mask_cluster(coevolutionary_domains,interacting_coevultionary_cluster_dict, self.bool_mask_clusters)
            
        elif self.alphafold_version == 'AF3':
            
            graph_resolution = 0.25
            matrix_cutoff = 25
            cmap_confidance =  'Blues_r'
            cmap_contact = "RdPu"
            
            coevolutionary_domains = get_coevolutionary_domains(matrix_input, graph_resolution = graph_resolution, matrix_cutoff = matrix_cutoff)
            
            interacting_coevultionary_cluster_dict, entity_region_dict = self.get_interacting_coevolutionary_domains(coevolutionary_domains)
            
            interacting_mask_cluster = self.get_interacting_mask_cluster(coevolutionary_domains,interacting_coevultionary_cluster_dict, self.bool_mask_clusters)
            
            
        if  self.plotting:
            
            plot_combination_matrix(coevolutionary_domains,masked_confidance_matrix,masked_contact_matrix,interacting_mask_cluster,list_fasta_files, self.outdir, self.alphafold_version)
            plot_separate_matrix(matrix_dict,list_fasta_files, self.outdir)
        
        return coevolutionary_domains, interacting_coevultionary_cluster_dict, entity_region_dict
    
    def get_interacting_coevolutionary_domains(self, coevolutionary_domains):
        
        list_fasta_name, list_fasta_acclen, list_fasta_centerticks, list_fasta_len = tuple(self.list_fasta_files)
        unique_coevolutionary_domains = np.unique(coevolutionary_domains)
        cluster_index_dict = {}
        interacting_coevultionary_cluster_dict = {}
        entity_region_dict = {}


        #iterate through all coevolutionary_domains found
        for cluster_name in unique_coevolutionary_domains:

            ii = np.where(coevolutionary_domains == cluster_name)[0]
            if len (ii) > 1 :
                if not cluster_name in cluster_index_dict:
                    cluster_index_dict[cluster_name] = ii
                
                #define start and end for each cluster    
                groups = (list(x) for _, x in
                        itertools.groupby(cluster_index_dict[cluster_name], lambda x, c=itertools.count(): x - next(c)))
                
                group_index_range_list = [(item[0], item[-1])[:len(item)] for item in groups]
                
                #build interacting_coevultionary_cluster_dict
                if not cluster_name in interacting_coevultionary_cluster_dict:
                    
                    cluster_group_name = f'cluster_{cluster_name}'
                    
                    interacting_coevultionary_cluster_dict[cluster_group_name] = {'range_index_list':group_index_range_list, 'overlap_complex':{}}


                #find overlapping sequences between all proteins and regions of coevolutionary_domains 
                for protein, acclen, fasta_len in zip(list_fasta_name, list_fasta_acclen, list_fasta_len):
                    
                    protein_start = acclen - fasta_len
                    protein_end = acclen - 1
                    
                    entity_region_dict[protein] = (protein_start,protein_end)
                        
                    for group_range in interacting_coevultionary_cluster_dict[cluster_group_name]['range_index_list']:
                        
                        if len(group_range) == 2: 
                            group_index_start = group_range[0]
                            group_index_end = group_range[1]
                    
                    
                            if group_index_start < protein_end and protein_start < group_index_end:
                                
                                protein_range = range(protein_start,protein_end)
                                cluster_range = range(group_index_start, group_index_end)
                                overlapping_range = (max(protein_range[0], cluster_range[0]), min(protein_range[-1], cluster_range[-1])+1)
                                
                                if not protein in interacting_coevultionary_cluster_dict[cluster_group_name]['overlap_complex']:
                                    interacting_coevultionary_cluster_dict[cluster_group_name]['overlap_complex'][protein] = []
                                
                                interacting_coevultionary_cluster_dict[cluster_group_name]['overlap_complex'][protein].append(overlapping_range)

        return interacting_coevultionary_cluster_dict, entity_region_dict
    
    def get_interacting_mask_cluster(self, clusters, interacting_coevultionary_cluster_dict ,bool_mask_clusters):
        
        if not bool_mask_clusters:
            
            interacting_mask_cluster = np.full(clusters.shape, False)
                        
            return interacting_mask_cluster
        
        else:
            
            interacting_mask_cluster = np.full(clusters.shape, True)
        
            interacting_range_list = []
            for cluster in interacting_coevultionary_cluster_dict:
                if len(interacting_coevultionary_cluster_dict[cluster]['overlap_complex']) > 1:
                    
                    for entity in interacting_coevultionary_cluster_dict[cluster]['overlap_complex']:
                        
                        for interacting_range in interacting_coevultionary_cluster_dict[cluster]['overlap_complex'][entity]:
                            
                            interacting_range_list.append(interacting_range)
        
            for index, cluster in enumerate(clusters):
                
                for interacting_range in interacting_range_list:
                    
                    if interacting_range[0] <= index <= interacting_range[1]:
                        interacting_mask_cluster[index] = False
        
        return interacting_mask_cluster
    
                    

def get_coevolutionary_domains(matrix_input, pae_power = 1, graph_resolution = 0.5, matrix_cutoff = 2.6 ):
                
        weights = 1/matrix_input**pae_power

        g = igraph.Graph()
        size = weights.shape[0]
        g.add_vertices(range(size))
        edges = np.argwhere(matrix_input < matrix_cutoff)
        sel_weights = weights[edges.T[0], edges.T[1]]
        g.add_edges(edges)
        g.es['weight']=sel_weights

        vc = g.community_leiden(weights='weight', resolution_parameter=graph_resolution/100, n_iterations=-1)
        coevolutionary_domains = np.array(vc.membership)
        
        return coevolutionary_domains


def plot_combination_matrix(coevolutionary_domains, confidance_matrix, contact_matrix, interacting_mask_cluster ,list_fasta_files, outdir, alphafold_version):
        
        labels = np.array(coevolutionary_domains)
        label_data = np.tile(labels, (2,1))
        
        mask_data = np.tile(interacting_mask_cluster, (2,1))
        
        list_fasta_name = list_fasta_files[0]
        list_fasta_acclen = list_fasta_files[1]
        list_fasta_centerticks = list_fasta_files[2]
        list_fasta_len = list_fasta_files[3]
        
        inv_fasta_acclen = [sum(list_fasta_len) - acclen for acclen in list_fasta_acclen]
        inv_fasta_acclen.insert(0, sum(list_fasta_len) -1 )

        list_fasta_len_names = [f'0 / {length}' if i != len(list_fasta_len)-1 else f'{length}' for i,length in enumerate(list_fasta_len)]
        
        fig, ax  = plt.subplots(figsize = (15,15))
        
        ax.set_xticks(list_fasta_acclen)
        ax.set_xticklabels('')            
        ax.set_xticks(list_fasta_centerticks,minor=True)
        ax.set_xticklabels(list_fasta_name, rotation=45, ha='right',va = 'center_baseline', fontsize=8,minor=True)
        
        ax.set_yticks(np.array(list_fasta_acclen))
        ax.set_yticklabels(list_fasta_len_names)
        #ax.plot([0, 1], [0, 1], transform=ax.transAxes) 
        
        divider = make_axes_locatable(ax)
        cax1 = divider.append_axes("right", size="5%", pad=0.5)
        cax2 = divider.append_axes("bottom", size="5%", pad=0.5) 
        cax3 = divider.append_axes("top", size="3%", pad=0.5) 
        
        sns.heatmap(label_data,mask= mask_data, ax=cax3 ,cbar= False ,cmap = sns.color_palette('tab20b') )
        
        for acclen in list_fasta_acclen:
            
            cax3.axvline(acclen, color = 'white', linewidth = 2)
            ax.axvline(acclen, color = 'black', linewidth = 1)
            ax.axhline(acclen , color = 'black', linewidth = 1)
        
        if alphafold_version == 'AF2':
            
            img2 = ax.imshow(contact_matrix, cmap="RdPu_r")
        
        elif alphafold_version == 'AF3':
            img2 = ax.imshow(contact_matrix, cmap="RdPu")
            
        img1 = ax.imshow(confidance_matrix, cmap="Blues_r")
        
        
        fig.colorbar(img1, orientation='vertical', cax = cax1)
        fig.colorbar(img2, orientation='horizontal', cax = cax2)
        
        cax3.set(yticklabels=[])
            
        cax3.set_xticklabels('')            
        cax3.set_xticks(list_fasta_centerticks,minor=True)
        cax3.set_xticklabels(list_fasta_name, va = 'center_baseline', fontsize=8,minor=True)
        
        cax3.tick_params(left = False, bottom =False, labelbottom=True) 
        plt.subplots_adjust(hspace=0.05)
        plt.savefig(f"{outdir}/Confidence-contact_plot.png",dpi=300)

def plot_separate_matrix(matrix_dict,list_fasta_files, outdir):
    N_COL = 3
    N_ROW = 1
    matrix_list =  ['pae','confidance_matrix','contact_matrix']
    cmap_list = ['Greens_r', 'Blues_r', "RdPu"]
    
    
    list_fasta_name = list_fasta_files[0]
    list_fasta_acclen = list_fasta_files[1]
    list_fasta_centerticks = list_fasta_files[2]
    list_fasta_len = list_fasta_files[3]
    list_fasta_len_names = [f'{length} / 0' if i != len(list_fasta_len)-1 else f'{length}' 
                for i,length in enumerate(list_fasta_len)]
    
    
    fig, axs = plt.subplots(N_ROW,N_COL,figsize=(30, 10))
    
    for n, (feature, cmap) in enumerate(zip(matrix_list, cmap_list)):
        
        matrix = matrix_dict[feature]
        n_col = n%3
        
        divider = make_axes_locatable(axs[n_col])
        cax = divider.append_axes("right", size="5%", pad=0.5)
        
        axs[n_col].set_xticks(list_fasta_acclen)
        axs[n_col].set_xticklabels('')
        axs[n_col].set_xticks(list_fasta_centerticks,minor=True)
        axs[n_col].set_xticklabels(list_fasta_name, rotation=45, ha='right',va = 'center_baseline', fontsize=8,minor=True) 
        axs[n_col].set_yticks(np.array(list_fasta_acclen)-1)
        axs[n_col].set_yticklabels(list_fasta_len_names)
        
        for i in list_fasta_acclen:

            axs[n_col].axvline(i, color = 'black', linewidth = 1)  
            axs[n_col].axhline(i, color = 'black', linewidth = 1)            


        img = axs[n_col].imshow(matrix , cmap=cmap)
        
        fig.colorbar(img, orientation='vertical', cax = cax)
    plt.subplots_adjust(top = 0.96, bottom=0.1, hspace=0.7, wspace=0.2)
    plt.show
    plt.savefig(f"{outdir}/feature_matrix.png",dpi=300)
    
    