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
from scipy.special import expit

from src.module.rec_input import RECORD_AF3, RECORD_SERVER
from src.module.parsers import PDBPARSER, MMCIFPARSER
from sklearn.metrics.pairwise import pairwise_distances

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Data import IUPACData


protein_letters_1to3 = IUPACData.protein_letters_1to3

upper_protein_letters_1to3 = {k.upper():v.upper() for k,v in protein_letters_1to3.items()}

upper_protein_letters_1to3['X'] = 'UNK'


module_dir = os.path.dirname(os.path.realpath(__file__))




class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()      
    
class FEATURE_MATRIX:

    MACROMOLECULE_NAME_DICT = {
        'protein': 'Protein',
        'dna': 'DNA',
        'rna': 'RNA',
        'glycan': 'Glycan',
        'ligand': 'Ligand',
        'ion': 'Ion',
    }

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
    
    def get_feature_matrix_dict(self, pae, plddt, iptm, chain_pair_iptm, plddt_matrix, pae_plddt, symmetric_pae, contact_matrix, confidence_matrix, masked_confidence_matrix, masked_contact_matrix, unfixed_pae, unfixed_contact_probability):

        matrix_dict = {}

        matrix_dict['pae'] = pae
        matrix_dict['unfixed_pae'] = unfixed_pae
        matrix_dict['unfixed_contact_probability'] = unfixed_contact_probability
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
    
    def _build_sequence_info_dict(self, chain_entity_lengths):
        """Build sequence_info_dict from a list of (label_asym_id, entity_length) pairs.

        This is the shared accumulation logic used by both CCM_AF3 and CCM_BOLTZ.
        Each subclass prepares the input list in its own way, then delegates here.
        """
        label_asym_id_list = []
        acclen_list = []
        centerticks_list = []
        length_list = []
        num_acc = 0

        for label_asym_id, entity_length in chain_entity_lengths:
            label_asym_id_list.append(label_asym_id)
            num_acc += entity_length
            if not acclen_list:
                center_tick = int(num_acc / 2)
            else:
                center_tick = int((num_acc - acclen_list[-1]) / 2 + acclen_list[-1])
            centerticks_list.append(center_tick)
            acclen_list.append(num_acc)
            length_list.append(entity_length)

        return {
            'label_asym_id': label_asym_id_list,
            'acclen': acclen_list,
            'centerticks': centerticks_list,
            'length': length_list,
        }

    def _build_plddt_dict(self, structure_coordinates, rec_list):
        """Build plddt_dict from CIF atom coordinates and rec_list.

        For polymer residues the pLDDT is the mean over all atoms; for non-polymer
        entities it is kept as a per-atom list (matches the CIF B_iso_or_equiv convention
        used by both AF3 and Boltz2).
        """
        plddt_dict = {}

        for rec in rec_list:
            label_asym_id = rec['label_asym_id']
            rec_type = rec['rec_type']
            plddt_dict.setdefault(label_asym_id, {})

            for seq_id, residue_data in structure_coordinates[label_asym_id].items():
                plddt_atom_list = [float(atom['plddt']) for atom in residue_data['atom_id']]
                atom_type_list  = [atom['atom_type']   for atom in residue_data['atom_id']]

                plddt_dict[label_asym_id][seq_id] = {
                    'plddt': np.mean(plddt_atom_list) if rec_type == 'polymer' else plddt_atom_list,
                    'atom_type_list': atom_type_list,
                }

        return plddt_dict

    def _get_residue_comp_id(self, rec, residue, _seq_id):
        """Return the 3-letter component ID for a single residue position.

        The base implementation covers standard polymer types (no PTMs).
        CCM_AF3 overrides this to handle post-translational modifications;
        _seq_id is unused here but required by that override's interface.
        """
        if rec['macromolecule_type'] == 'Protein':
            return upper_protein_letters_1to3[residue]
        return residue

    def _build_chain_info_dict(self, rec_list, plddt_dict):
        """Assemble chain_info_dict from rec_list and plddt_dict.

        Populates residue-level pLDDT and comp_id for polymer chains, and
        atom-level pLDDT for non-polymer entities. Calls _get_residue_comp_id,
        which subclasses can override (e.g. for PTM handling in CCM_AF3).
        Also computes and assigns entity_degree for all chains.
        """
        chain_info_dict = {'polymer': [], 'non_polymer': []}

        for rec in rec_list:
            label_asym_id = rec['label_asym_id']
            rec['macromolecule_type'] = self.MACROMOLECULE_NAME_DICT[rec['macromolecule_type']]

            if rec['rec_type'] == 'polymer':
                rec['residues'] = []
                for index, residue in enumerate(rec['sequence']):
                    seq_id = index + 1
                    rec['residues'].append({
                        'seq_id': seq_id,
                        'comp_id': self._get_residue_comp_id(rec, residue, seq_id),
                        'plddt': plddt_dict[label_asym_id][seq_id]['plddt'],
                    })
                chain_info_dict['polymer'].append(rec)
            else:
                rec['atoms'] = [
                    {'atom_type': atom_type, 'atom_id': idx + 1, 'plddt': atom_plddt}
                    for idx, (atom_type, atom_plddt) in enumerate(
                        zip(
                            plddt_dict[label_asym_id]['.']['atom_type_list'],
                            plddt_dict[label_asym_id]['.']['plddt'],
                        )
                    )
                ]
                chain_info_dict['non_polymer'].append(rec)

        monomer_name_len_dict, _ = get_monomer_info_dict(chain_info_dict)
        for recs in chain_info_dict.values():
            for rec in recs:
                rec['entity_degree'] = monomer_name_len_dict[rec['auth_asym_id']]['entity_degree']

        return chain_info_dict

    def extract_matrix_dict(self):
        """Build the full matrix dict.

        Calls get_feature_info(), which each subclass must implement and which
        must return (distance_matrix, pae, contact_matrix, plddt, iptm, chain_pair_iptm).
        Subclasses that perform PTM fixing (CCM_AF3) override this method to also
        supply unfixed_pae and unfixed_contact_probability; here they default to
        the same arrays since no fixing occurs.
        """
        _, pae, contact_matrix, plddt, iptm, chain_pair_iptm = self.get_feature_info()

        symmetric_pae, pae_plddt, confidence_matrix, plddt_matrix = \
            self.get_pae_plddt_matrix(pae, plddt)

        binary_contact = contact_matrix > 0.5
        mask_upper = np.triu(binary_contact, k=0)
        masked_contact_matrix = np.ma.array(binary_contact, mask=mask_upper)

        mask_lower = np.tri(pae_plddt.shape[0], k=0)
        masked_confidence_matrix = np.ma.array(confidence_matrix, mask=mask_lower)

        return self.get_feature_matrix_dict(
            pae, plddt, iptm, chain_pair_iptm,
            plddt_matrix, pae_plddt, symmetric_pae,
            contact_matrix, confidence_matrix,
            masked_confidence_matrix, masked_contact_matrix,
            unfixed_pae=pae,
            unfixed_contact_probability=contact_matrix,
        )

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
    
    def __init__(self, 
                 in_dir, 
                 sample: int = 0):
        
        super().__init__(in_dir)
        
        self.sample = sample
    
    def check_if_alphabridge_server(self):
        
        folder_path = self.in_dir
        
        server_source_path = os.path.join(folder_path, 'server_source')

        # Check if server source file exists
        if not os.path.isfile(server_source_path):
            return False
        else:
            return True
    
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
        sample = self.sample
        
        if self.check_if_path_exist(folder_path):
            
            if not self.check_if_alphabridge_server():
        
                if self.check_alphafold_dialect():
                    
                    feature_path = [file  for file in list(Path(folder_path).glob( "*_confidences.json")) if not'_summary_confidences' in str(file)][0]
                    structure_path = list(Path(folder_path).glob( "*model.cif"))[0]
                    job_request_path = list(Path(folder_path).glob("*data.json"))[0]
                    summary_request_path = list(Path(folder_path).glob("*summary_confidences*.json"))[0]
                    alphafold_dialect = 'AlphaFold_local'
                    
                else:
            
                    feature_path = list(Path(folder_path).glob( f"*full_data_{sample}.json"))[0]
                    structure_path = list(Path(folder_path).glob( f"*model_{sample}.cif"))[0]
                    job_request_path = list(Path(folder_path).glob(f"*job_request.json"))[0]
                    summary_request_path = list(Path(folder_path).glob(f"*summary_confidences_{sample}.json"))[0]
                    alphafold_dialect = 'AlphaFold_server'
            else:
                
                server_source_path = os.path.join(folder_path, 'server_source')
                
                feature_path = list(Path(folder_path).glob( f"*confidence_metrics.json"))[0]
                structure_path = list(Path(folder_path).glob( f"*structure.cif"))[0]
                job_request_path = list(Path(folder_path).glob(f"*job_request.json"))[0]
                summary_request_path = list(Path(folder_path).glob(f"*summary_metrics.json"))[0]
                
                with open(server_source_path, 'r') as file:
                    alphafold_dialect = file.read().strip()
                
            return feature_path, structure_path, job_request_path, summary_request_path, alphafold_dialect
    
    def extract_rec_list(self, job_request_path, structure_sequence_list, feature_dict, alphafold_dialect):

        request_file = read_json_file(job_request_path)
        if alphafold_dialect == 'AlphaFold_local':
            #print(structure_sequence_list)
            rec_list = RECORD_AF3(request_file, structure_sequence_list, feature_dict).process_record_file()
        elif alphafold_dialect == 'AlphaFold_server':
            rec_list = RECORD_SERVER(request_file, structure_sequence_list, feature_dict).process_record_file()
        else:
            raise NotImplementedError('Format File not Valid')
                       
        return rec_list
    
    def extract_job_id_name(self, job_request_path, alphafold_dialect):
        
        request_file = read_json_file(job_request_path)
        
        if alphafold_dialect == 'AlphaFold_local':
            job_id_name = request_file['name']
        
        elif alphafold_dialect == 'AlphaFold_server':
            job_id_name = request_file[0]['name']

        return job_id_name

    def extract_sequence_info(self):

        feature_path, structure_path, job_request_path, summary_request_path, alphafold_dialect = self.extract_feature_filepath()

        structure = MMCIFPARSER(structure_path)
        
        feature_dict = read_json_file(feature_path)
        
        structure_sequence_list = structure.get_sequence_list()
        
        rec_list = self.extract_rec_list(job_request_path, structure_sequence_list, feature_dict, alphafold_dialect)
        
        sequence_info_dict = self.extract_sequence_info_dict(feature_dict, rec_list)
     
        return rec_list, sequence_info_dict
    
    def extract_sequence_info_dict(self, feature_dict, rec_list):
        token_chain_id = feature_dict['token_chain_ids']
        chain_entity_lengths = [
            (rec['label_asym_id'],
             len(rec['sequence']) if rec['rec_type'] == 'polymer'
             else token_chain_id.count(rec['label_asym_id']))
            for rec in rec_list
        ]
        return self._build_sequence_info_dict(chain_entity_lengths)
        
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
        _, structure_path, _, _, _ = self.extract_feature_filepath()
        structure = MMCIFPARSER(structure_path)
        rec_list, _ = self.extract_sequence_info()
        return self._build_plddt_dict(structure.get_coordinates(), rec_list)
    
    def fix_matrix_size(self, feature_dict, rec_list):

        pae = np.array(feature_dict['pae'])
        contact_probability = np.array(feature_dict['contact_probs'])

        token_chain_ids = feature_dict['token_chain_ids']
        token_res_ids = feature_dict['token_res_ids']

        chain_index_dict = {}

        mask = []

        unfixed_pae = pae.copy()
        unfixed_contact_probability = contact_probability.copy()
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
            
        return fixed_pae, fixed_contact_probability, unfixed_pae, unfixed_contact_probability
    
    def get_feature_info(self):                                                                                                                                                                                                                                    
        
        feature_path, structure_path, job_request_path, summary_request_path, alphafold_dialect = self.extract_feature_filepath()
        
        structure = MMCIFPARSER(structure_path)
        
        structure_sequence_list = structure.get_sequence_list()
        
        feature_dict = read_json_file(feature_path)
        
        rec_list = self.extract_rec_list(job_request_path, structure_sequence_list, feature_dict, alphafold_dialect)
        
        summary_request_dict =  read_json_file(summary_request_path)

        chain_pair_iptm = np.where(np.array(summary_request_dict['chain_pair_iptm'])==None, 0, np.array(summary_request_dict['chain_pair_iptm'])).astype(float) 

        iptm = summary_request_dict['iptm']
            
        plddt = self.extract_plddt_per_token(structure)
        
        distance_matrix = self.get_distance_matrix(structure.get_ca_distances())
            
        pae, contact_probability, unfixed_pae, unfixed_contact_probability = self.fix_matrix_size(feature_dict, rec_list)

        return distance_matrix, pae, contact_probability, plddt, iptm, chain_pair_iptm, unfixed_pae, unfixed_contact_probability
    
    
    def _get_residue_comp_id(self, rec, residue, seq_id):
        """Override: handle PTMs in addition to standard residues."""
        if rec['macromolecule_type'] == 'Protein':
            if not rec['modifications']:
                return upper_protein_letters_1to3[residue]
            comp_id = upper_protein_letters_1to3[residue]
            for modification in rec['modifications']:
                if modification['ptmPosition'] == seq_id:
                    ptm_type = modification['ptmType']
                    comp_id = ptm_type.replace('CCD_', '') if ptm_type.startswith('CCD_') else ptm_type
                else:
                    comp_id = upper_protein_letters_1to3[residue]
            return comp_id
        return residue

    def extract_matrix_dict(self):

        distance_matrix, pae, contact_probability, plddt, iptm, chain_pair_iptm, unfixed_pae, unfixed_contact_probability = self.get_feature_info()

        symmetric_pae, pae_plddt, confidence_matrix, plddt_matrix = self.get_pae_plddt_matrix(pae, plddt)

        contact_matrix = contact_probability

        binary_contact = contact_probability > 0.5

        mask_upper = np.triu(binary_contact, k=0)
        masked_contact_matrix = np.ma.array(binary_contact, mask=mask_upper)

        mask_lower = np.tri(pae_plddt.shape[0], k=0)
        masked_confidence_matrix = np.ma.array(confidence_matrix, mask=mask_lower)

        matrix_dict = self.get_feature_matrix_dict(pae, plddt, iptm, chain_pair_iptm, plddt_matrix, pae_plddt, symmetric_pae, contact_matrix, confidence_matrix, masked_confidence_matrix, masked_contact_matrix, unfixed_pae, unfixed_contact_probability)

        return matrix_dict

    def extract_chain_info_dict(self):
        rec_list, sequence_info_dict = self.extract_sequence_info()
        plddt_dict = self.get_plddt_dict()
        chain_info_dict = self._build_chain_info_dict(rec_list, plddt_dict)
        return chain_info_dict, sequence_info_dict
              
class CCM_BOLTZ(FEATURE_MATRIX):
    """Confidence Contact Matrix extractor for Boltz2 predictions."""

    # Boltz2 mol_type integer → macromolecule_type string
    MOL_TYPE_MAP = {
        0: 'protein',
        1: 'rna',
        2: 'dna',
        3: 'ligand',
        4: 'ion',
        5: 'glycan',
    }

    def __init__(self, in_dir, sample: int = 0):
        super().__init__(in_dir)
        self.sample = sample

    # ------------------------------------------------------------------
    # Path helpers
    # ------------------------------------------------------------------

    def _get_job_name(self):
        predictions_dir = os.path.join(self.in_dir, 'predictions')
        job_names = [
            d for d in os.listdir(predictions_dir)
            if os.path.isdir(os.path.join(predictions_dir, d))
        ]
        if len(job_names) != 1:
            raise ValueError(
                f"Expected exactly one job directory in {predictions_dir}, "
                f"found: {job_names}"
            )
        return job_names[0]

    def extract_feature_filepath(self):
        in_dir = self.in_dir
        sample = self.sample

        if self.check_if_path_exist(in_dir):
            job_name = self._get_job_name()
            pred_dir = os.path.join(in_dir, 'predictions', job_name)

            structure_path  = os.path.join(pred_dir, f'{job_name}_model_{sample}.cif')
            pae_path        = os.path.join(pred_dir, f'pae_{job_name}_model_{sample}.npz')
            plddt_path      = os.path.join(pred_dir, f'plddt_{job_name}_model_{sample}.npz')
            pde_path        = os.path.join(pred_dir, f'pde_{job_name}_model_{sample}.npz')
            confidence_path = os.path.join(pred_dir, f'confidence_{job_name}_model_{sample}.json')
            records_path    = os.path.join(in_dir, 'processed', 'records', f'{job_name}.json')

            return structure_path, pae_path, plddt_path, pde_path, confidence_path, records_path, job_name

    def extract_job_id_name(self):
        _, _, _, _, _, records_path, _ = self.extract_feature_filepath()
        records_data = read_json_file(records_path)
        return records_data['id']

    # ------------------------------------------------------------------
    # Record / sequence info
    # ------------------------------------------------------------------

    def extract_rec_list_from_boltz(self, records_data, structure_sequence_list):
        """Build rec_list from boltz processed/records JSON and CIF sequence list."""
        seq_map = dict(structure_sequence_list)

        rec_list = []
        for chain in records_data['chains']:
            mol_type = chain['mol_type']
            macromolecule_type = self.MOL_TYPE_MAP.get(mol_type, 'ligand')
            chain_name = chain['chain_name']

            if macromolecule_type in ('protein', 'rna', 'dna'):
                rec = {
                    'rec_type': 'polymer',
                    'macromolecule_type': macromolecule_type,
                    'sequence': seq_map.get(chain_name, ''),
                    'entity_degree': int(),
                    'modifications': [],
                    'auth_asym_id': chain_name,
                    'label_asym_id': chain_name,
                }
            else:
                rec = {
                    'rec_type': 'non_polymer',
                    'macromolecule_type': macromolecule_type,
                    'non_poly_entity': seq_map.get(chain_name, ''),
                    'smiles': str(),
                    'entity_degree': int(),
                    'auth_asym_id': chain_name,
                    'label_asym_id': chain_name,
                }
            rec_list.append(rec)

        return rec_list

    def extract_sequence_info_dict(self, records_data, rec_list):
        chain_entity_lengths = [
            (rec['label_asym_id'], chain['num_residues'])
            for rec, chain in zip(rec_list, records_data['chains'])
        ]
        return self._build_sequence_info_dict(chain_entity_lengths)

    def extract_sequence_info(self):
        structure_path, _, _, _, _, records_path, _ = self.extract_feature_filepath()

        structure = MMCIFPARSER(structure_path)
        structure_sequence_list = structure.get_sequence_list()
        records_data = read_json_file(records_path)

        rec_list = self.extract_rec_list_from_boltz(records_data, structure_sequence_list)
        sequence_info_dict = self.extract_sequence_info_dict(records_data, rec_list)

        return rec_list, sequence_info_dict

    # ------------------------------------------------------------------
    # pLDDT dict (per-residue / per-atom for chain_info_dict)
    # ------------------------------------------------------------------

    def get_plddt_dict(self):
        structure_path, _, _, _, _, records_path, _ = self.extract_feature_filepath()
        structure = MMCIFPARSER(structure_path)
        records_data = read_json_file(records_path)
        rec_list = self.extract_rec_list_from_boltz(records_data, structure.get_sequence_list())
        return self._build_plddt_dict(structure.get_coordinates(), rec_list)

    # ------------------------------------------------------------------
    # chain_pair_iptm matrix
    # ------------------------------------------------------------------

    def get_chain_pair_iptm_matrix(self, confidence_data, num_chains):
        """Convert pair_chains_iptm dict-of-dicts to a numpy matrix."""
        pair_chains_iptm = confidence_data.get('pair_chains_iptm', {})
        matrix = np.zeros((num_chains, num_chains))
        for i in range(num_chains):
            for j in range(num_chains):
                matrix[i, j] = pair_chains_iptm.get(str(i), {}).get(str(j), 0.0)
        return matrix

    # ------------------------------------------------------------------
    # Core feature extraction
    # ------------------------------------------------------------------

    def get_feature_info(self):
        structure_path, pae_path, plddt_path, pde_path, confidence_path, records_path, _ = \
            self.extract_feature_filepath()

        structure = MMCIFPARSER(structure_path)
        records_data = read_json_file(records_path)

        pae   = np.load(pae_path)['pae']
        # Boltz plddt npz is in [0, 1]; scale to [0, 100] to match AF3 convention
        plddt = np.load(plddt_path)['plddt'] * 100.0
        pde   = np.load(pde_path)['pde']

        confidence_data = read_json_file(confidence_path)
        iptm            = confidence_data['iptm']
        num_chains      = len(records_data['chains'])
        chain_pair_iptm = self.get_chain_pair_iptm_matrix(confidence_data, num_chains)

        contact_matrix  = pde_to_contact_prob(pde).astype(np.float64)
        distance_matrix = self.get_distance_matrix(structure.get_ca_distances())

        return distance_matrix, pae, contact_matrix, plddt, iptm, chain_pair_iptm

    def extract_chain_info_dict(self):
        structure_path, _, _, _, _, records_path, _ = self.extract_feature_filepath()
        structure = MMCIFPARSER(structure_path)
        records_data = read_json_file(records_path)
        rec_list = self.extract_rec_list_from_boltz(records_data, structure.get_sequence_list())
        sequence_info_dict = self.extract_sequence_info_dict(records_data, rec_list)
        plddt_dict = self.get_plddt_dict()
        chain_info_dict = self._build_chain_info_dict(rec_list, plddt_dict)
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

def get_monomer_info_dict(chain_info_dict):
    
    chain_info_dict = chain_info_dict
    
    total_token_nr = 0
    ligand_nr = 0
    ion_nr = 0
    glycan_nr = 0
    
    total_ligand_token = 0
    total_ion_token = 0
    total_glycan_token = 0
    
    ion_degree = 5
    ligand_degree = 10
    glycan_degree = 10

    poly_bolean = False

    monomer_name_len_dict = {}
    poly_type_dict = {}
    
    for entity_type, rec_list in chain_info_dict.items():
    
                for rec in rec_list:
                
                    auth_asym_id = rec['auth_asym_id']
                    if not auth_asym_id in monomer_name_len_dict:
                        monomer_name_len_dict[auth_asym_id] = {'entity_length' : int(), 'entity_degree':int()}
                        poly_type_dict[auth_asym_id] = str()
                    
                    poly_type_dict[auth_asym_id] = rec['macromolecule_type']
    
                    if entity_type == 'polymer':
                    
                        entity_length = len(rec['residues'])
                        monomer_name_len_dict[auth_asym_id]['entity_length'] = entity_length
                    else:
                        entity_length = len(rec['atoms'])
                        #print(entity_length)
                        poly_bolean = True
                        
                        monomer_name_len_dict[auth_asym_id]['entity_length'] = None
    
                        if rec['macromolecule_type'] == 'Ligand':
                            ligand_nr +=1
                            total_ligand_token += entity_length
                            
                        elif rec['macromolecule_type'] == 'Ion':
                            ion_nr += 1
                            total_ion_token += entity_length
                        elif rec['macromolecule_type'] == 'Glycan':
                            glycan_nr += 1
                            total_glycan_token += entity_length
    
                    total_token_nr += entity_length
    

    if poly_bolean:

        non_poly_section = (ligand_degree * ligand_nr)  + (ion_degree * ion_nr) + (glycan_degree * glycan_nr)
        
        non_poly_section_cut_off = 30

        if non_poly_section > non_poly_section_cut_off:
            
            non_poly_entity_length = (non_poly_section_cut_off * total_token_nr) // 360
                
            ligand_length = non_poly_entity_length // (ligand_nr + 1/2 * ion_nr + glycan_nr)
            ion_length = ligand_length // 2
            glycan_length = ligand_length
    
        else:

            ligand_length = (ligand_degree * total_token_nr) // 360
            ion_length = (ion_degree * total_token_nr) // 360
            glycan_length = (glycan_degree * total_token_nr) // 360
            

        for auth_asym_id, type in poly_type_dict.items():
            
            if type == 'Ligand':
                monomer_name_len_dict[auth_asym_id]['entity_length'] = int(ligand_length)
            elif type == 'Ion':
                monomer_name_len_dict[auth_asym_id]['entity_length'] = int(ion_length)
            elif type == 'Glycan':
                monomer_name_len_dict[auth_asym_id]['entity_length'] = int(glycan_length)
        

    entity_sum_length = sum([value['entity_length'] for key, value in monomer_name_len_dict.items()])

    for key, value in monomer_name_len_dict.items():
        monomer_name_len_dict[key]['entity_degree'] = (value['entity_length'] * 360) / entity_sum_length
                   
    return monomer_name_len_dict, poly_type_dict   


def pde_to_contact_prob(pde, threshold=4.0, k=0.8):
    """Approximate contact probability from a PDE matrix via a sigmoid.

    P(contact | i,j) ≈ σ(−k · (PDE_ij − threshold))

    Parameters
    ----------
    pde : np.ndarray, shape (N, N)
        PDE matrix in Å. Lower values = higher confidence.
    threshold : float
        PDE value at which contact probability = 0.5. Default 4.0 Å.
    k : float
        Sigmoid steepness. Default 0.8.

    Returns
    -------
    contact_prob : np.ndarray, shape (N, N), dtype float32
    """
    pde = np.asarray(pde, dtype=np.float64)
    contact_prob = expit(-k * (pde - threshold)).astype(np.float32)
    return contact_prob
