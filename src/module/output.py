from itertools import chain, repeat, count, islice
from collections import Counter


class output:
    
    def __init__(self, structure_score_dict, interactions_list):
        self.structure_score_dict = structure_score_dict
        self.interactions_list = interactions_list

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
        
        self.structure_score_dict['pairwise_interaction'] = self.get_pairwise_scores()
        
        alphabridge_dict = {
        "structure": [self.structure_score_dict],
        "interactions" : self.interactions_list
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