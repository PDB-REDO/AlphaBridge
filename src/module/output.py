from itertools import chain, repeat, count, islice
from collections import Counter
import numpy as np
import math
from src.module.scoring import (
    build_chains_array, calculate_ipsae, calculate_actifptm
)

class OUTPUT:
    
    def __init__(self, structure_score_dict, interactions_list, sequence_info_dict, matrix_dict):
        self.structure_score_dict = structure_score_dict
        self.interactions_list = interactions_list
        self.sequence_info_dict = sequence_info_dict
        self.matrix_dict = matrix_dict

    def get_pairwise_scores(self):
        
        label_asym_id_list = []
        pairwise_score_dict = {}
        pairwise_combination_list = []
        
        for entity_type, rec_list in self.structure_score_dict['chains'].items():

            for rec in rec_list:
                
                label_asym_id = rec['label_asym_id']
                label_asym_id_list.append(label_asym_id)
        
        biomolecule_combination_list = list(unique_combinations(label_asym_id_list, 2))
        
        pairwise_combination_list = []
        for biomolecule_combination in biomolecule_combination_list:
            score_dict = {
                'first' : str(),
                'second' : str(),
                'pairwise_score' : float()
            }
            biomolecule_1 = biomolecule_combination[0]
            biomolecule_2 = biomolecule_combination[1]
            
            score_dict['first'] = biomolecule_1
            score_dict['second'] = biomolecule_2
            pairwise_combination_list.append(score_dict)
        
        
        for interactions in self.interactions_list:
        
            for interace in interactions['interfaces']:
                interface_chain = (interace['links'][0]['first']['asym_id'] + interace['links'][0]['second']['asym_id'])

                if not interface_chain in pairwise_score_dict:
                    pairwise_score_dict[interface_chain] = 0
                
                pairwise_score_dict[interface_chain] =  interactions['cut-off']
                
        for pairwise_combination in pairwise_combination_list:
            
            biomolecule_1 = pairwise_combination['first']
            biomolecule_2 = pairwise_combination['second']
            interface_chain = biomolecule_1 + biomolecule_2
            
            if interface_chain in pairwise_score_dict:
                pairwise_combination['pairwise_score'] = pairwise_score_dict[interface_chain]
            
        
    
        
        return pairwise_combination_list
    
    

    def get_alphabridge_dict(self):

        label_asym_id_list = [chain['label_asym_id'] for polymer_type in self.structure_score_dict['chains'] for chain in self.structure_score_dict['chains'][polymer_type]]

        pairwise_scores = self.get_pairwise_scores()

        # Build inputs for ipSAE / actifpTM
        chains_array   = build_chains_array(self.sequence_info_dict)
        chain_type_dict = _build_chain_type_dict(self.structure_score_dict['chains'])
        polymer_ids    = [rec['label_asym_id'] for rec in self.structure_score_dict['chains']['polymer']]

        pae           = np.array(self.matrix_dict['unfixed_pae'])
        contact_probs = np.array(self.matrix_dict['unfixed_contact_probability'])

        ipsae_results   = calculate_ipsae(pae, chains_array, polymer_ids, chain_type_dict)
        actifptm_results = calculate_actifptm(pae, contact_probs, chains_array, polymer_ids, chain_type_dict)

        for combo in pairwise_scores:
            c1, c2 = combo['first'], combo['second']

            # ipSAE: take max of both asymmetric directions
            ab = ipsae_results.get((c1, c2), {})
            ba = ipsae_results.get((c2, c1), {})
            for key in ('ipsae_d0chn', 'ipsae_d0dom', 'ipsae_d0res', 'iptm_d0chn'):
                combo[key] = max(ab.get(key, 0.0), ba.get(key, 0.0))

            # actifpTM: symmetric pair
            combo['actifptm'] = actifptm_results.get((c1, c2)) or actifptm_results.get((c2, c1), 0.0)

        self.structure_score_dict['pairwise_interaction'] = pairwise_scores
        self.structure_score_dict['AB_score'] = calculate_alphabridge_score(pairwise_scores, label_asym_id_list)

        alphabridge_dict = {
            "structure": [self.structure_score_dict],
            "interactions": self.interactions_list
        }

        return alphabridge_dict
        
        
        

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

def _build_chain_type_dict(chains_dict):
    """Map label_asym_id -> 'nucleic_acid' or 'protein' for d0 calculation."""
    NUCLEIC = {'DNA', 'RNA'}
    chain_type = {}
    for rec_list in chains_dict.values():
        for rec in rec_list:
            mol_type = rec['macromolecule_type']
            chain_type[rec['label_asym_id']] = 'nucleic_acid' if mol_type in NUCLEIC else 'protein'
    return chain_type


def calculate_alphabridge_score(combination_list, label_asym_id_list):
    """
    Calculate the AlphaBridge score based on the pairwise score.
    This is a placeholder for the actual scoring logic.
    """
    alphabridge_score = float()

    if not combination_list:
        return alphabridge_score

    else:
        interacting_partners = [combination for combination in combination_list if combination['pairwise_score'] > 0]

        interacting_asym_id_list = []
        for interacting_partner in interacting_partners:
            interacting_asym_id_list.append(interacting_partner['first'])
            interacting_asym_id_list.append(interacting_partner['second'])

        non_interacting_asym_id_list = list(set(label_asym_id_list) - set(interacting_asym_id_list))

        if not non_interacting_asym_id_list:

            pairwise_scores = [combination['pairwise_score'] for combination in interacting_partners]
            alphabridge_score = math.prod(pairwise_scores) ** (1/len(pairwise_scores))

    return alphabridge_score


