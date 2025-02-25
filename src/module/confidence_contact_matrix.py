# -*- coding: utf-8 -*-
"""
Created on Wed Feb 7 16:53:05 2024

@author: Dan_salv
"""
import os
import numpy as np
import json
from itertools import groupby
from Bio import SeqIO
from pathlib import Path
import pandas as pd
import warnings
warnings.filterwarnings("ignore")
from Bio.Data import IUPACData
import errno

from src.module.rec_input import RECORD_AF3, RECORD_SERVER
from src.module.parsers import PDBPARSER, MMCIFPARSER
from sklearn.metrics.pairwise import pairwise_distances

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Data import IUPACData


protein_letters_1to3 = IUPACData.protein_letters_1to3

upper_protein_letters_1to3 = {k.upper():v.upper() for k,v in protein_letters_1to3.items()}

upper_protein_letters_1to3


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
        confidence_matrix = pae_plddt.copy()
        
        confidence_matrix[np.where(confidence_matrix > 32)] = 32
        
        return symmetric_pae, pae_plddt, confidence_matrix , plddt_matrix
    
    def get_feature_matrix_dict(self, pae, plddt, iptm, chain_pair_iptm ,plddt_matrix, pae_plddt , symmetric_pae, contact_matrix, confidence_matrix, masked_confidence_matrix, masked_contact_matrix ):
        
        matrix_dict = {}
        
        matrix_dict['pae'] = pae
        matrix_dict['plddt'] = plddt
        matrix_dict['iptm'] = iptm
        matrix_dict['chain_pair_iptm'] = chain_pair_iptm
        matrix_dict['plddt_matrix'] = plddt_matrix
        matrix_dict['pae_plddt'] = pae_plddt
        matrix_dict['symmetric_pae'] = symmetric_pae
        matrix_dict['contact_matrix'] = contact_matrix
        matrix_dict['confidence_matrix'] = confidence_matrix
        matrix_dict['masked_confidence_matrix'] = masked_confidence_matrix
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
        
        excluded_keys = ['pae','plddt','symmetric_pae','pae_plddt','masked_confidence_matrix', 'masked_contact_matrix']
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
    
    def check_alphafold_dialect(self):
        
        folder_path = self.in_dir
        
        job_request_path = list(Path(folder_path).glob( "*_data.json"))
        
        if not job_request_path:
            return False
        
        else:
            request_job = read_json_file(job_request_path[0])
            
            if request_job['dialect'] == 'alphafold3':
                return True
            else:
                raise NotImplementedError('Format File not Valid')

    
    def extract_feature_filepath(self):
        
        folder_path = self.in_dir
        
        if self.check_if_path_exist(folder_path):
        
            if self.check_alphafold_dialect():
                
                feature_path = [file  for file in list(Path(folder_path).glob( "*_confidences.json")) if not'_summary_confidences' in str(file)][0]
                structure_path = list(Path(folder_path).glob( "*model.cif"))[0]
                job_request_path = list(Path(folder_path).glob("*data.json"))[0]
                summary_request_path = list(Path(folder_path).glob("*summary_confidences*.json"))[0]
            
            else:
        
                feature_path = list(Path(folder_path).glob( "*full_data_0.json"))[0]
                structure_path = list(Path(folder_path).glob( "*model_0.cif"))[0]
                job_request_path = list(Path(folder_path).glob("*job_request*.json"))[0]
                summary_request_path = list(Path(folder_path).glob("*summary_confidences*.json"))[0]
        
            return feature_path, structure_path, job_request_path, summary_request_path
    
    def extract_rec_list(self, job_request_path, structure_sequence_list, feature_dict):

        request_file = read_json_file(job_request_path)
        
        if self.check_alphafold_dialect():

            rec_list = RECORD_AF3(request_file, structure_sequence_list, feature_dict).process_record_file()
        else:
            rec_list = RECORD_SERVER(request_file, structure_sequence_list, feature_dict).process_record_file()
                       
        return rec_list
                     
    def extract_sequence_info(self):
        
        feature_path, structure_path, job_request_path, summary_request_path = self.extract_feature_filepath()
        
        structure = MMCIFPARSER(structure_path)
        
        feature_dict = read_json_file(feature_path)
        
        structure_sequence_list = structure.get_sequence_list()
        
        rec_list = self.extract_rec_list(job_request_path, structure_sequence_list, feature_dict)
        
        sequence_info_dict = self.extract_sequence_info_dict(feature_dict, rec_list)
     
        return rec_list, sequence_info_dict
    
    def extract_sequence_info_dict(self, feature_dict, rec_list):
        
        label_asym_id_list = []
        label_asym_id_acclen_list = []
        label_asym_id_centerticks_list = []
        label_asym_id_len_list = []
        num_acc = 0

        token_chain_id = feature_dict['token_chain_ids']

        for rec in rec_list:
            
            label_asym_id = rec['label_asym_id']
            
            if rec['rec_type'] == 'polymer':
                entity_length = len(rec['sequence']) 
            else:
                entity_length = token_chain_id.count(label_asym_id)
            
            label_asym_id_list.append(label_asym_id)
            
            num_acc += entity_length
            if len(label_asym_id_acclen_list) == 0 :
                center_tick = int((num_acc - 0)/2)
            else:
                center_tick = int((num_acc - label_asym_id_acclen_list[-1])/2 + label_asym_id_acclen_list[-1])
            label_asym_id_centerticks_list.append(center_tick)
            label_asym_id_acclen_list.append(num_acc)
            label_asym_id_len_list.append(entity_length)

        sequence_info_dict = {
            'label_asym_id' : label_asym_id_list,
            'acclen' : label_asym_id_acclen_list,
            'centerticks' : label_asym_id_centerticks_list,
            'length' : label_asym_id_len_list
        }
        return sequence_info_dict
        
    def extract_plddt_per_token(self, structure):
        #will need to be changed again 
        structure_coordinates = structure.get_coordinates()
        rec_list, sequence_info_dict= self.extract_sequence_info()
        
        token_plddt_list = []
        
        
        for rec in rec_list:
            label_asym_id = rec['label_asym_id']
            rec_type = rec['rec_type']
            for seq_id in structure_coordinates[label_asym_id]:
                
                plddt_atom_list = [float(atom_id['plddt']) for atom_id in structure_coordinates[label_asym_id][seq_id]['atom_id']]
                
                if rec_type == 'polymer':
                    
                    plddt_value = [np.mean(plddt_atom_list)]
                    
                else:
                    plddt_value = plddt_atom_list
                
                token_plddt_list += plddt_value
    
        return token_plddt_list
    
    def get_plddt_dict(self):
        #will need to be changed again
        feature_path, structure_path, job_request_path, summary_request_path = self.extract_feature_filepath()
        
        structure = MMCIFPARSER(structure_path) 
        
        structure_coordinates = structure.get_coordinates()
        
        rec_list, sequence_info_dict = self.extract_sequence_info()
        
        plddt_dict = {}
        
        for rec in rec_list:
            label_asym_id = rec['label_asym_id']
            rec_type = rec['rec_type']
            for seq_id in structure_coordinates[label_asym_id]:
                
                plddt_atom_list = [float(atom_id['plddt']) for atom_id in structure_coordinates[label_asym_id][seq_id]['atom_id']]
                atom_type_list = [atom_id['atom_type'] for atom_id in structure_coordinates[label_asym_id][seq_id]['atom_id']]
                
                if not label_asym_id in plddt_dict:
                    plddt_dict[label_asym_id] = {}
                if not seq_id in plddt_dict[label_asym_id]:
                    plddt_dict[label_asym_id][seq_id] = {'plddt': float(),
                                                         'atom_type_list': []}
                    
                
                if rec_type == 'polymer':
                    
                    plddt_dict[label_asym_id][seq_id]['plddt'] = np.mean(plddt_atom_list)
                    plddt_dict[label_asym_id][seq_id]['atom_type_list'] = atom_type_list
                else:
                    plddt_dict[label_asym_id][seq_id]['plddt'] = plddt_atom_list
                    plddt_dict[label_asym_id][seq_id]['atom_type_list'] = atom_type_list
                

        return plddt_dict
    
    def fix_matrix_size(self, feature_dict, rec_list):
        
        pae = np.array(feature_dict['pae'])
        contact_probability = np.array(feature_dict['contact_probs'])

        token_chain_ids = feature_dict['token_chain_ids']
        token_res_ids = feature_dict['token_res_ids']

        chain_index_dict = {}

        mask = []

        fixed_pae = pae.copy()
        fixed_contact_probability = contact_probability.copy()
            
        for rec in rec_list:
            
            label_aysm_id = rec['label_asym_id']
            
            if rec['macromolecule_type'] == 'protein' and rec['modifications']: 
                
                for modification in rec['modifications']:
                            
                    ptm_position = modification['ptmPosition']
                    
                    mask = np.array([False if token_chain == label_aysm_id and token_res == ptm_position else True 
                            for token_chain, token_res in zip(token_chain_ids,token_res_ids)])
                    
                    fixed_pae = summarize_ptm_matrix(fixed_pae, mask, ptm_position, np.min)
                    fixed_contact_probability = summarize_ptm_matrix(fixed_contact_probability, mask, ptm_position, np.max)
                    
                    token_chain_ids, token_res_ids = fix_token_lists(mask, token_chain_ids, token_res_ids)
            
        return fixed_pae, fixed_contact_probability
    
    def get_feature_info(self):                                                                                                                                                                                                                                    
        
        feature_path, structure_path, job_request_path, summary_request_path = self.extract_feature_filepath()
        
        structure = MMCIFPARSER(structure_path)
        
        structure_sequence_list = structure.get_sequence_list()
        
        feature_dict = read_json_file(feature_path)
        
        rec_list = self.extract_rec_list(job_request_path, structure_sequence_list, feature_dict)
        
        summary_request_dict =  read_json_file(summary_request_path)

        chain_pair_iptm = np.where(np.array(summary_request_dict['chain_pair_iptm'])==None, 0, np.array(summary_request_dict['chain_pair_iptm'])).astype(float) 

        iptm = summary_request_dict['iptm']
            
        plddt = self.extract_plddt_per_token(structure)
        
        distance_matrix = self.get_distance_matrix(structure.get_ca_distances())
            
        pae, contact_probability= self.fix_matrix_size(feature_dict, rec_list)
        
        return distance_matrix, pae, contact_probability, plddt, iptm , chain_pair_iptm
    
    
    def extract_matrix_dict(self):   
        
        distance_matrix, pae, contact_probability, plddt, iptm, chain_pair_iptm = self.get_feature_info()
        
        symmetric_pae, pae_plddt, confidence_matrix, plddt_matrix = self.get_pae_plddt_matrix(pae, plddt)
            
        contact_matrix = contact_probability
        
        binary_contact = contact_probability > 0.5
        
        mask_upper=  np.triu(binary_contact, k=0)
        masked_contact_matrix = np.ma.array(binary_contact, mask=mask_upper)
    
        mask_lower =  np.tri(pae_plddt.shape[0], k=0)
        masked_confidence_matrix = np.ma.array(confidence_matrix, mask=mask_lower)
        
        matrix_dict = self.get_feature_matrix_dict(pae, plddt, iptm ,chain_pair_iptm, plddt_matrix, pae_plddt, symmetric_pae, contact_matrix, confidence_matrix, masked_confidence_matrix, masked_contact_matrix )
        
        return matrix_dict
    
    def extract_chain_info_dict(self):
        
        
        rec_list, sequence_info_dict = self.extract_sequence_info()
        plddt_dict = self.get_plddt_dict()
        
        chain_info_dict = {
            'polymer' : [],
            'non_polymer' : []
        }
        
        for rec in rec_list:
            
            auth_asym_id = rec['auth_asym_id']
            label_asym_id = rec['label_asym_id']
            
            if rec['rec_type'] == 'polymer':
                rec['residues'] = []
                
                for index, residue in enumerate(rec['sequence']):
                    seq_id = index + 1
                    
                    residue_dict = {
                        "seq_id": int(),
                        "comp_id" : str(),
                        "plddt" : float()
                    }
                    
                    if rec['macromolecule_type'] == 'protein':
                        
                        if not rec['modifications']:
                            comp_id = upper_protein_letters_1to3[residue]
                        
                        else:
                            for modification in rec['modifications']:
                                
                                if modification['ptmPosition'] == seq_id:
                                  
                                    comp_id = modification['ptmType'].replace('CCD_', '') if modification['ptmType'].startswith('CCD_') else modification['ptmType']
                                
                                else:
                                    comp_id = upper_protein_letters_1to3[residue]
                                    
                    else:
                        comp_id = residue
                        
                    residue_dict['seq_id'] = seq_id
                    residue_dict['comp_id'] = comp_id
                    residue_dict['plddt'] = plddt_dict[label_asym_id][seq_id]['plddt']
                    
                    rec['residues'].append(residue_dict)
                chain_info_dict['polymer'].append(rec)
            else:
                rec['atoms'] = []   
                for index, (atom_type, atom_plddt) in enumerate(zip(plddt_dict[label_asym_id]['.']['atom_type_list'], plddt_dict[label_asym_id]['.']['plddt'])):
                    atom_id = index + 1
                    atom_dict = {
                    'atom_type': str(),
                    'atom_id': int(),
                    "plddt" : float()
                    }
                    
                    atom_dict['atom_type'] = atom_type
                    atom_dict['atom_id'] = atom_id
                    atom_dict['plddt'] = atom_plddt
                    
                    
                    rec['atoms'].append(atom_dict)
                chain_info_dict['non_polymer'].append(rec)
    
        return chain_info_dict, sequence_info_dict
              

def read_json_file(json_file):
    
    with open(json_file, 'r') as f:
        
        file = json.loads(f.read())
        
    return file

        
def summarize_ptm_matrix(matrix, mask, ptm_position, func):
    #column array
    scored_array = matrix[:, ~mask]
    #row array
    aligned_array = matrix[~mask]
    #matrix without ptm token
    no_ptm_matrix = matrix[mask,:][:,mask]
    
    #ptm expanded metric
    inside_values = [index for index, value in enumerate(mask) if not value]
    
    ptm_value = func(matrix[~mask,:][:,~mask])
    
    #get the minimum maximum value from each matrix array
    scored_values = func(scored_array, axis=1)
    aligned_values = func(aligned_array, axis=0)
    
    #remove ptm_expanded_values
    scored_values = np.delete(scored_values, inside_values)
    aligned_values = np.delete(aligned_values, inside_values)
    #remove ptm value to match dimesnsion array
    aligned_values =  np.insert(aligned_values, ptm_position -1,  ptm_value)
    #add ptm token as a residue
    added_column_matrix = np.insert(no_ptm_matrix, ptm_position -1, scored_values, axis=1)

    ptm_matrix = np.insert(added_column_matrix, ptm_position -1, aligned_values, axis=0)
    
    return ptm_matrix

def fix_token_lists(mask, token_chain_ids, token_res_ids):
    
    remove_ptm_index = [index for index, value in enumerate(mask) if not value][1:]
    
    token_chain_ids =  np.delete(token_chain_ids, remove_ptm_index)
    token_res_ids =  np.delete(token_res_ids, remove_ptm_index)
    
    return token_chain_ids, token_res_ids