# -*- coding: utf-8 -*-
"""
Created on Wed Feb 7 16:53:05 2024

@author: Dan_salv
"""
import os
import numpy as np
import json
from Bio import SeqIO
from pathlib import Path
import pandas as pd
import warnings
warnings.filterwarnings("ignore")
from Bio.Data import IUPACData
import errno

from src.module.parsers import PDBPARSER, MMCIFPARSER
from sklearn.metrics.pairwise import pairwise_distances
from src.module.alingment_utils import compare_protein_seq

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord



module_dir = os.path.dirname(os.path.realpath(__file__))




class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()      
    
class FEATURE_MATRIX:
    
    def __init__(self, in_dir):
        
        self.in_dir = in_dir
    
    
    def check_if_path_exist(self, filepath):
        
        if Path(filepath).exists():
            return True
        else:
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), filepath) 
           #sys.exit(f"{accesion_id}: {filename} does not exist in: {database}")
    
    
    def fasta_profiles(self, fasta_sequences):
        
        list_fasta_name = []
        list_fasta_acclen = []
        list_fasta_len = []
        list_fasta_files= []
        list_fasta_centerticks = []
        num_acc = 0
        
        for i, fasta in enumerate(fasta_sequences):
            name, sequence = fasta.id, str(fasta.seq)
            list_fasta_name.append(name)
            num_acc += len(sequence)
            if len(list_fasta_acclen) == 0 :
                center_tick = int((num_acc - 0)/2)
            else:
                center_tick = int((num_acc - list_fasta_acclen[-1])/2 + list_fasta_acclen[-1])
            list_fasta_centerticks.append(center_tick)
            list_fasta_acclen.append(num_acc)
            list_fasta_len.append(len(sequence))
        
        return [list_fasta_name, list_fasta_acclen, list_fasta_centerticks, list_fasta_len]
    
    def get_sequence_chain_tuple(self, sequence_list):
        
        sequence_chain_tuple = [(rec.id ,rec.seq) for rec in sequence_list]

        return sequence_chain_tuple
    
        
    def get_distance_matrix(self, ca_distances):
        
        distance_matrix = pairwise_distances(ca_distances,ca_distances)
        
        return distance_matrix
    
    def get_pae_plddt_matrix(self, pae, plddt):
        
        symmetric_pae = pae.copy()

        for i,column in enumerate(symmetric_pae):
            for j,row in enumerate(symmetric_pae):
                
                if symmetric_pae[i][j] < symmetric_pae[j][i]:
                    
                    symmetric_pae[i][j] = symmetric_pae[j][i]
        
        
        plddt_matrix = np.zeros((len(plddt), len(plddt))) 

        for i,column in enumerate(plddt_matrix):
            for j,row in enumerate(plddt_matrix):
                delta_index = i - j 
                if  -2 <= delta_index <= 2:
                    plddt_matrix[i][j] = 0
                else :
                    plddt_matrix[i][j] =  -1 * (((plddt[i] + plddt[j]) / 2) - 100)


        pae_plddt = symmetric_pae + plddt_matrix / 3

        pae_plddt[np.where(pae_plddt > 32)] = 32
        
        return symmetric_pae, pae_plddt, plddt_matrix
    
    def get_feature_matrix_dict(self, pae, plddt, chain_pair_iptm ,plddt_matrix, pae_plddt, symmetric_pae, contact_matrix, confidance_matrix, masked_confidance_matrix, masked_contact_matrix ):
        
        matrix_dict = {}
        
        matrix_dict['pae'] = pae
        matrix_dict['plddt'] = plddt
        matrix_dict['chain_pair_iptm'] = chain_pair_iptm
        matrix_dict['plddt_matrix'] = plddt_matrix
        matrix_dict['pae_plddt'] = pae_plddt
        matrix_dict['symmetric_pae'] = symmetric_pae
        matrix_dict['contact_matrix'] = contact_matrix
        matrix_dict['confidance_matrix'] = confidance_matrix
        matrix_dict['masked_confidance_matrix'] = masked_confidance_matrix
        matrix_dict['masked_contact_matrix'] = masked_contact_matrix
        
        return matrix_dict
    
    
    def get_scores_dict(self, scores_list, list_fasta_files):
        
        scores_dict = {}
        initial_acc = 0
        
        list_fasta_name, list_fasta_acclen, list_fasta_centerticks, list_fasta_len = tuple(list_fasta_files)

        for fasta_name, fasta_acclen in zip(list_fasta_name, list_fasta_acclen):
            
            if not fasta_name in scores_dict:
                scores_dict[fasta_name] = scores_list[initial_acc:fasta_acclen]
            
            initial_acc = fasta_acclen
        
        return scores_dict
    
    def print_matrix_dict(self, matrix_dict):
        
        excluded_keys = ['pae','plddt','symmetric_pae','pae_plddt','masked_confidance_matrix', 'masked_contact_matrix']
        tmp_matrix = {i:matrix_dict[i] for i in matrix_dict if i not in  excluded_keys}
                   
        if os.path.exists(self.in_dir):
            feature_object_path = os.path.join(self.in_dir, 'matrix_info.json')
        
            with open(feature_object_path, 'w') as f:
            
                f.write(json.dumps(tmp_matrix,
                                   cls=NumpyEncoder

                                   )
                        )

class CCM_AF3(FEATURE_MATRIX):
    
    def __init__(self, in_dir):
        
        super().__init__(in_dir)
    
    def extract_feature_filepath(self):
        
        folder_path = self.in_dir
        
        if self.check_if_path_exist(folder_path):
        
            feature_path = list(Path(folder_path).glob( "*full_data_0.json"))[0]
            structure_path = list(Path(folder_path).glob( "*model_0.cif"))[0]
            job_request_path = list(Path(folder_path).glob("*job_request*.json"))[0]
            summary_request_path = list(Path(folder_path).glob("*summary_confidences*.json"))[0]
        
            return feature_path, structure_path, job_request_path, summary_request_path
    
    def extract_sequences(self, job_request_path):

            request_file = json.load(open(job_request_path))[0] 

            chains =(i for i in "ABCDEFGHIJKLMNOPQRSTUVWXYZ")
            known_keys = ["proteinChain", 'ligand' ,"rnaSequence", "ion", "dnaSequence"]
            seq_types_symbols = {"proteinChain": "", "rnaSequence": "RNA_", "dnaSequence": "DNA_", 'ligand':'ligand_', 'ion':'ion_'}

            rec_list = []
            
            for macromolecule in request_file['sequences']:
                key = list(macromolecule.keys())[0]
                for item in range(macromolecule[key]['count']):
                    if key in ['proteinChain', 'dnaSequence', 'rnaSequence']:
                        seq = Seq(macromolecule[key]['sequence'])
                        chain_id = seq_types_symbols[key] + next(chains)
                        rec = SeqRecord(seq, id= chain_id, name = chain_id , description= chain_id)
                        rec_list.append(rec)

            return rec_list
    
    def extract_sequence_info(self):
        
        feature_path, structure_path, job_request_path, summary_request_path = self.extract_feature_filepath()
        
        structure = MMCIFPARSER(structure_path)
        
        feature_dict = read_json_file(feature_path)
        
        structure_sequence_list = structure.get_sequence_list()
        polymer_chain_dict = structure.get_polypeptide_chain_dict()
        
        tmp_rec_list = self.extract_sequences(job_request_path)
        
        rec_list = self.reorder_sequences_by_token_list(tmp_rec_list, structure_sequence_list, feature_dict)
        
        chain_tuple = self.get_sequence_chain_tuple(rec_list)
        
        list_sequence_info = self.fasta_profiles(rec_list)
        
        return list_sequence_info, chain_tuple, structure_sequence_list, polymer_chain_dict
        
    
    def reorder_sequences_by_token_list(self,tmp_rec_list, structure_sequence_list,feature_dict ):
        
        fasta_sequence_list = self.get_sequence_chain_tuple(tmp_rec_list)
        
        fasta_to_struct = compare_protein_seq(structure_sequence_list, fasta_sequence_list).extract_chain_dict()
        
        struct_to_fasta = {v: k for k, v in fasta_to_struct.items()}
        
        token_chain_ids = feature_dict['token_chain_ids']
        unique_chains = np.unique(token_chain_ids)
        
        unique_fasta_chains = [struct_to_fasta[unique_chain] for unique_chain in unique_chains if unique_chain in struct_to_fasta]
                
        reorder_fasta_sequence_list = [tuple for x in unique_fasta_chains for tuple in fasta_sequence_list if tuple[0] == x]
        
        reorder_rec_list = [ SeqRecord(reorder_tuple[1], id= reorder_tuple[0], name = reorder_tuple[0] , description= reorder_tuple[0]) for reorder_tuple in reorder_fasta_sequence_list]
        
        return reorder_rec_list
    
    def extract_plddt_per_residue(self, structure):
        structure_coordinates = structure.get_coordinates()
        polymer_chain_dict = structure.get_polypeptide_chain_dict()
        
        data = []

        for asym_id in structure_coordinates:
            for seq_id in structure_coordinates[asym_id]:
                for atom_id in structure_coordinates[asym_id][seq_id]['atom_id']:
                    
                    if polymer_chain_dict[asym_id]['entity_type'] in ['polypeptide(L)','polydeoxyribonucleotide', 'polyribonucleotide']:
                        
                        plddt = float(structure_coordinates[asym_id][seq_id]['atom_id'][atom_id]['plddt'])
                        data.append([asym_id ,seq_id,plddt])
                
        plddt_df = pd.DataFrame(data, columns=['asym_id','seq_id','plddt'])               

        plddt_per_residue_df = plddt_df.groupby(by=["asym_id", 'seq_id']).mean().reset_index()

        residue_plddts = plddt_per_residue_df['plddt'].tolist()
    
        return residue_plddts
        
    def fix_matrix_size(self, structure,feature_dict, chain_pair_iptm):
        
        pae = np.array(feature_dict['pae'])
        contact_probability = np.array(feature_dict['contact_probs'])
        
        token_chain_ids = feature_dict['token_chain_ids']
        unique_chains = np.unique(token_chain_ids)
        chain_index_dict = {}
        
        for chain in unique_chains:
            ii = [i for i, j in enumerate(token_chain_ids) if j == chain]
    
            if not chain in chain_index_dict:
                chain_index_dict[chain] = ii
            
        #get chains that are a polypeptide
        polymer_chain_dict = structure.get_polypeptide_chain_dict()
        
        polypeptide_chain_list = [chain  for chain in polymer_chain_dict if polymer_chain_dict[chain]['entity_type']  in ['polypeptide(L)','polydeoxyribonucleotide', 'polyribonucleotide']] 
        
        non_polypeptide_chain_list = list(set(unique_chains) - set(polypeptide_chain_list))
        
        remove_chain_idx = np.isin(unique_chains, non_polypeptide_chain_list)

        remove_idx_list = []

        for non_polypeptide_chain in non_polypeptide_chain_list:
        
            remove_idx_list += chain_index_dict[non_polypeptide_chain]
                
        mask = np.ones(pae.shape[0], bool)        
        mask[remove_idx_list] = 0
        
        mask_chain = np.ones(chain_pair_iptm.shape[0], bool)
        mask_chain[remove_chain_idx] = 0
          
        
        
        fix_size_pae = pae[mask,:][:,mask]
        fix_size_contact_probability = contact_probability[mask,:][:,mask]
        fix_size_chain_pair_iptm = chain_pair_iptm[mask_chain,:][:,mask_chain]
        
        return fix_size_pae, fix_size_contact_probability, fix_size_chain_pair_iptm
    
    def get_feature_info(self):                                                                                                                                                                                                                                    
        
        
        feature_path, structure_path, job_request_path, summary_request_path = self.extract_feature_filepath()
        
        structure = MMCIFPARSER(structure_path)
        
        feature_dict = read_json_file(feature_path)
        
        summary_request_dict =  read_json_file(summary_request_path)

        chain_pair_iptm = np.where(np.array(summary_request_dict['chain_pair_iptm'])==None, 0, np.array(summary_request_dict['chain_pair_iptm'])).astype(float) 

        plddt = self.extract_plddt_per_residue(structure)
        
        distance_matrix = self.get_distance_matrix(structure.get_ca_distances())
        
        pae, contact_probability, chain_pair_iptm = self.fix_matrix_size(structure, feature_dict, chain_pair_iptm)
        
        return distance_matrix, pae, contact_probability, plddt, chain_pair_iptm
    
    
    def extract_matrix_dict(self):   
        
        distance_matrix, pae, contact_probability, plddt, chain_pair_iptm = self.get_feature_info()
        
        symmetric_pae, pae_plddt, plddt_matrix = self.get_pae_plddt_matrix(pae, plddt)
        
        confidance_matrix = pae_plddt
            
        contact_matrix = contact_probability
        
        binary_contact = contact_probability > 0.5
        
        mask_upper=  np.triu(binary_contact, k=0)
        masked_contact_matrix = np.ma.array(binary_contact, mask=mask_upper)
    
        mask_lower =  np.tri(pae_plddt.shape[0], k=0)
        masked_confidance_matrix = np.ma.array(pae_plddt, mask=mask_lower)
        
        matrix_dict = self.get_feature_matrix_dict(pae, plddt, chain_pair_iptm, plddt_matrix, pae_plddt, symmetric_pae, contact_matrix, confidance_matrix, masked_confidance_matrix, masked_contact_matrix )
        
        return matrix_dict
    
    
    

def read_json_file(json_file):
    
    with open(json_file, 'r') as f:
        
        file = json.loads(f.read())
        
    return file

        
        