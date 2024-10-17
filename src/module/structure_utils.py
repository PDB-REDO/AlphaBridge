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
import configparser

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


config = configparser.ConfigParser()

module_dir = os.path.dirname(os.path.realpath(__file__))

config.read(f'{module_dir}/config.ini')


WORKING_DIR = config['DEFAULT']['WORKING_DIR']
AF3_DIR = config['DEFAULT']['AF3_DIR']


class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()      
    
class FEATURE_MATRIX:
    
    def __init__(self, job_id, outdir:str = '', alphafold_version: str = 'AF2'):
        
        self.job_id = job_id
        self.alphafold_version = alphafold_version
        self.outdir = outdir

    def get_sequence_info_path(self):
        
        job_id = self.job_id
        
        if self.alphafold_version == 'AF2':
            
            filepath = f'/DATA/zata/projects/interactions/data/request/fasta/{job_id}.fasta'
        
        elif self.alphafold_version == 'AF3':
            
            filepath = self.get_feature_folder_path()
        
        if Path(filepath).exists():
            return filepath
        else:
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), filepath) 
           #sys.exit(f"{accesion_id}: {filename} does not exist in: {database}")
    
    def get_feature_folder_path(self):
        
        job_id = self.job_id
        
        if self.alphafold_version == 'AF2':
        
            filepath = f'/DATA/zata/projects/interactions/data/results/{job_id}'
        
        elif self.alphafold_version == 'AF3':
            
            filepath = f'{AF3_DIR}/{job_id}'
            
        if Path(filepath).exists():
            return filepath
        else:
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), filepath) 
           #sys.exit(f"{accesion_id}: {filename} does not exist in: {database}")
    
    def extract_sequences(self, sequence_info_path):

        if self.alphafold_version == 'AF2':
            
            rec_list = SeqIO.parse(open(sequence_info_path),'fasta')
        
        elif self.alphafold_version == 'AF3':


            folder_path = Path(sequence_info_path)

            job_request_pattern = "*job_request*.json"
            job_request = json.load(open(list(folder_path.glob(job_request_pattern))[0]))[0]


            chains =(i for i in "ABCDEFGHIJKLMNOPQRSTUVWXYZ")
            known_keys = ["proteinChain", 'ligand' ,"rnaSequence", "ion", "dnaSequence"]
            seq_types_symbols = {"proteinChain": "", "rnaSequence": "RNA_", "dnaSequence": "DNA_", 'ligand':'ligand_', 'ion':'ion_'}

            rec_list = []
            
            for macromolecule in job_request['sequences']:
                key = list(macromolecule.keys())[0]
                for item in range(macromolecule[key]['count']):
                    if key in ['proteinChain', 'dnaSequence', 'rnaSequence']:
                        seq = Seq(macromolecule[key]['sequence'])
                        chain_id = seq_types_symbols[key] + next(chains)
                        rec = SeqRecord(seq, id= chain_id, name = chain_id , description= chain_id)
                        rec_list.append(rec)

        return rec_list

    def get_sequence_chain_tuple(self):
        
        sequence_info_path =  self.get_sequence_info_path()

        sequence_list = self.extract_sequences(sequence_info_path)

        sequence_chain_tuple = [(rec.id ,rec.seq) for rec in sequence_list]

        return sequence_chain_tuple

    def reorder_sequences_by_token_list(self):
        
        structure_path, feature_path = self.find_feature_files()
        
        structure_sequence_list = MMCIFPARSER(structure_path).get_sequence_list()
        fasta_sequence_list = self.get_sequence_chain_tuple()
        
        fasta_to_struct = compare_protein_seq(structure_sequence_list, fasta_sequence_list).extract_chain_dict()
        
        struct_to_fasta = {v: k for k, v in fasta_to_struct.items()}
        
        feature_dict = read_json_file(feature_path)
        
        token_chain_ids = feature_dict['token_chain_ids']
        unique_chains = np.unique(token_chain_ids)
        
        unique_fasta_chains = [struct_to_fasta[unique_chain] for unique_chain in unique_chains if unique_chain in struct_to_fasta]
                
        reorder_fasta_sequence_list = [tuple for x in unique_fasta_chains for tuple in fasta_sequence_list if tuple[0] == x]
        
        reorder_rec_list = [ SeqRecord(reorder_tuple[1], id= reorder_tuple[0], name = reorder_tuple[0] , description= reorder_tuple[0]) for reorder_tuple in reorder_fasta_sequence_list]
        
        return reorder_rec_list
    
    def fasta_profiles(self):
        
        list_fasta_name = []
        list_fasta_acclen = []
        list_fasta_len = []
        list_fasta_files= []
        list_fasta_centerticks = []
        num_acc = 0
        sequence_info_path =  self.get_sequence_info_path()
        
        if self.alphafold_version == 'AF2' : 
            fasta_sequences = self.extract_sequences(sequence_info_path)
        
        
        elif self.alphafold_version == 'AF3':
            fasta_sequences = self.reorder_sequences_by_token_list()
        
        
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
        
    def find_feature_files(self):
        
        if self.alphafold_version == 'AF2':
            feature_folder = self.get_feature_folder_path()
            
            with open(os.path.join(feature_folder, 'ranking_debug.json')) as json_file:
                ranking = json.load(json_file)
            
            first_model = ranking['order'][0]
            feature_filename = f'result_{first_model}.pkl'
            
            structure_path = os.path.join(feature_folder, 'ranked_0.pdb')
            feature_path = os.path.join(feature_folder, feature_filename)
        
        elif self.alphafold_version == 'AF3':
            folder_path = Path(self.get_feature_folder_path())
            
            feature_path = list(Path(folder_path).glob( "*full_data_0.json"))[0]
            structure_path = list(Path(folder_path).glob( "*model_0.cif"))[0]
            
        return structure_path, feature_path
    
    def get_distance_matrix(self):
        
        structure_path, feature_path = self.find_feature_files()
        
        if self.alphafold_version == 'AF2':
            
            ca_distances = PDBPARSER(structure_path).get_ca_distances()
        
        elif self.alphafold_version == 'AF3':
            
            ca_distances = MMCIFPARSER(structure_path).get_ca_distances()


        distance_matrix = pairwise_distances(ca_distances,ca_distances)
        
        return distance_matrix
    
    
    def fix_matrix_size(self, structure_path,feature_dict):
        
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
        polymer_chain_dict = MMCIFPARSER(structure_path).get_polypeptide_chain_dict()
        
        polypeptide_chain_list = [chain  for chain in polymer_chain_dict if polymer_chain_dict[chain]['entity_type']  in ['polypeptide(L)','polydeoxyribonucleotide', 'polyribonucleotide']] 

        non_polypeptide_chain_list = list(set(unique_chains) - set(polypeptide_chain_list))

        remove_idx_list = []

        for non_polypeptide_chain in non_polypeptide_chain_list:
        
            remove_idx_list += chain_index_dict[non_polypeptide_chain]
                
        mask = np.ones(pae.shape[0], bool)
        mask[remove_idx_list] = 0
        
        fix_size_pae = pae[mask,:][:,mask]
        fix_size_contact_probability = contact_probability[mask,:][:,mask]
        
        return fix_size_pae, fix_size_contact_probability

    
    def extract_plddt_per_residue(self, structure_path):
        structure_coordinates = MMCIFPARSER(structure_path).get_coordinates()
        polymer_chain_dict = MMCIFPARSER(structure_path).get_polypeptide_chain_dict()
        
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
        
    def get_feature_info(self):
        structure_path, feature_path = self.find_feature_files()        
        
        if self.alphafold_version == 'AF2':

            distance_matrix = self.get_distance_matrix()
        
            feature_dict = np.load(feature_path, allow_pickle=True)
        
            pae = feature_dict['predicted_aligned_error']
            plddt = feature_dict['plddt']
            contact_probability = {}

        
        elif self.alphafold_version == 'AF3':

            folder_path = Path(self.get_feature_folder_path())
            feature_dict = read_json_file(feature_path)

            plddt = self.extract_plddt_per_residue(structure_path)
            distance_matrix = self.get_distance_matrix()
            
            pae, contact_probability = self.fix_matrix_size(structure_path, feature_dict)
        

        return distance_matrix, pae,contact_probability, plddt


    def extract_matrix_dict(self):
        
        matrix_dict = {}
        
        
        distance_matrix, pae, contact_probability, plddt = self.get_feature_info()

            
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
        
        if self.alphafold_version == 'AF2':
            
            log_distance = np.log10(distance_matrix, out=np.zeros_like(distance_matrix), where=(distance_matrix!=0))

            log_pae = np.log10(symmetric_pae, out=np.zeros_like(symmetric_pae), where=(pae!=0))
        
            confidance_matrix = log_distance + log_pae
            contact_matrix = distance_matrix

            modified_distance_matrix = distance_matrix.copy()
            modified_distance_matrix[np.where(distance_matrix > 40)] = 40
            
            mask_upper =  np.triu(modified_distance_matrix, k=0)
            masked_contact_matrix = np.ma.array(modified_distance_matrix, mask=mask_upper)
            
            mask_lower =  np.tri(pae.shape[0], k=0)
            masked_confidance_matrix = np.ma.array(pae, mask=mask_lower)
        
        elif self.alphafold_version == 'AF3':
            
            confidance_matrix = pae_plddt
            
            contact_matrix = contact_probability
            
            binary_contact = contact_probability > 0.5
            
            mask_upper=  np.triu(binary_contact, k=0)
            masked_contact_matrix = np.ma.array(binary_contact, mask=mask_upper)
        
            mask_lower =  np.tri(pae_plddt.shape[0], k=0)
            masked_confidance_matrix = np.ma.array(pae_plddt, mask=mask_lower)
        

        matrix_dict['pae'] = pae
        matrix_dict['plddt'] = plddt
        matrix_dict['plddt_matrix'] = plddt_matrix
        matrix_dict['pae_plddt'] = pae_plddt
        matrix_dict['symmetric_pae'] = symmetric_pae
        matrix_dict['contact_matrix'] = contact_matrix
        matrix_dict['confidance_matrix'] = confidance_matrix
        matrix_dict['masked_confidance_matrix'] = masked_confidance_matrix
        matrix_dict['masked_contact_matrix'] = masked_contact_matrix
        
        return matrix_dict
    
    def print_matrix_dict(self):
        
        
        matrix_dict = self.extract_matrix_dict()
        
        matrix_dict.pop('masked_confidance_matrix', None)
        matrix_dict.pop('masked_contact_matrix', None)        
            
        if os.path.exists(self.outdir):
            feature_object_path = os.path.join(self.outdir, 'matrix_info.json')
        
            with open(feature_object_path, 'w') as f:
                f.write(json.dumps(matrix_dict, indent=4,
                                   cls=NumpyEncoder
                                   )
                        )
        
    
    

def read_json_file(json_file):
    
    with open(json_file, 'r') as f:
        
        file = json.loads(f.read())
        
    return file


   

    
