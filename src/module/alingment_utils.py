# -*- coding: utf-8 -*-
"""
Created on Fri Nov 11 15:49:05 2022

@author: Daals
"""

from Bio import SeqIO,  Align
from Bio.Align import substitution_matrices
import os
import json
from collections import defaultdict

        
class protein_seq_alingment():
    
    def __init__(self, 
                 query_seq, ref_seq):
        
        self.query_seq = query_seq
        self.ref_seq = ref_seq
    
    def get_alignment(self):
        
        ref_seq = self.query_seq
        query_seq = self.ref_seq
        
        aligner = Align.PairwiseAligner(scoring="blastp")
        aligner.mode = 'global'
        matrix = substitution_matrices.load("BLOSUM62")
        aligner.substitution_matrix = matrix
        aligner_output = aligner.align(ref_seq, query_seq)
        alignment = aligner_output[0]
        
        
        return alignment

class compare_protein_seq:
    
    def __init__(self, 
                 structure_list, 
                 fasta_list):
        
       self.structure_list = structure_list
       self.fasta_list = fasta_list
    

    def extract_chain_dict(self):
        
        structure_list = self.structure_list
        fasta_list = self.fasta_list
        
        chain_dict = {}
        
        
        struct_dict = defaultdict(str)
        fasta_dict = defaultdict(str)

        chain_dict = {}

        for struct, fasta in zip(structure_list,fasta_list):
 
            struct_chain, struct_seq = struct
            fasta_chain, fasta_seq = fasta
            
            if not str(struct_seq) in struct_dict:
                struct_dict[str(struct_seq)] = []
            
            if not str(fasta_seq) in fasta_dict:
                fasta_dict[str(fasta_seq)] = []
            
            struct_dict[str(struct_seq)].append(struct_chain)
            fasta_dict[str(fasta_seq)].append(fasta_chain)



        key_seq_union = set().union(*[struct_dict,fasta_dict])


        for key_seq  in key_seq_union:
            
            for struct_chain, fasta_chain in zip(struct_dict[key_seq], fasta_dict[key_seq]):
                
                chain_dict[fasta_chain] = struct_chain
                
        return chain_dict

class msa_folder:
    
    def __init__(self, 
                 feature_folder, 
                 fasta_list):
        
       self.feature_folder = feature_folder
       self.fasta_list = fasta_list
    
    def extract_chain_id_map(self):
        
        feature_folder = self.feature_folder
        
        chain_id_map_path = os.path.join(feature_folder,'msas','chain_id_map.json')
        
        with open(chain_id_map_path) as json_file:
            chain_id_map = json.load(json_file)
        
        return chain_id_map
    
    def extract_msa_folder(self):
        
        feature_folder = self.feature_folder
        
        chain_id_map = self.extract_chain_id_map()
        
        parent_msa_folder = os.path.join(feature_folder,'msas')
        
        dedup_msa_folder = defaultdict(str)
        
        for record in chain_id_map:
            
            sequence = chain_id_map[record]['sequence']
            msa_folder = os.path.join(parent_msa_folder, record)
            if  os.path.exists(msa_folder):
                dedup_msa_folder[sequence] = record
        
        return dedup_msa_folder
        
        
    def map_description2msa_folder(self):
        
        dedup_msa_folder = self.extract_msa_folder()
        fasta_list = self.fasta_list
        
        fasta_msa_dict = {}
        
        for fasta_name ,fasta_seq in fasta_list:
            
            if not fasta_name in fasta_msa_dict:
                fasta_msa_dict[fasta_name] = dedup_msa_folder[str(fasta_seq)]
        
        return fasta_msa_dict

 
        
        

    
    
        
        
        
    


    
                
               




