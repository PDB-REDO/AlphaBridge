import os
import numpy as np
import json

class RECORD():
    
    def __init__(self,
                 record_file,
                 structure_sequence_list,
                 feature_dict
                 ):
        
        self.record_file = record_file 
        self.structure_sequence_list = structure_sequence_list
        self.feature_dict = feature_dict
    
    def map_asym_id(self, rec_list):
        
        structure_sequence_list = self.structure_sequence_list
        
        feature_dict = self.feature_dict
        
        used_indices = set()
        token_chain_ids = np.unique(feature_dict['token_chain_ids']).tolist()
        structure_sequence_list_ordered = sorted(structure_sequence_list, key=lambda x: token_chain_ids.index(x[0]))
        seq_types_symbols = {"protein": "", "rna": "RNA_", "dna": "DNA_", 'ligand':'lig_', 'ion':'Ion_', 'glycan': 'Glycan_'}
        #print(structure_sequence_list_ordered)
        #print(rec_list)
       
        for record in rec_list:
            entity = record['sequence'] if record['rec_type'] == 'polymer' else record['non_poly_entity']
            #print(record)
            #print(entity,structure_sequence_list_ordered)
            index_to_match = next(
                (i for i, tup in enumerate(structure_sequence_list_ordered)
                if tup[1] == entity and i not in used_indices),
                None
            )
            macromolecule_type = record['macromolecule_type']

            if index_to_match is not None:
                used_indices.add(index_to_match)
            else:     
                print(f"⚠️ No match found for entity: {entity}")
            
            
            record['auth_asym_id'] = seq_types_symbols[macromolecule_type] + structure_sequence_list_ordered[index_to_match][0]
            record['label_asym_id'] = structure_sequence_list_ordered[index_to_match][0]
            
            
            
        
        sorted_rec_list = sorted(rec_list, key=lambda x: token_chain_ids.index(x['label_asym_id']))
    
        return sorted_rec_list
        
    

class RECORD_AF3(RECORD):

    def __init__(self, 
                record_file,
                structure_sequence_list,
                feature_dict):
        
        super().__init__(
                record_file,
                structure_sequence_list,
                feature_dict)
        
    
    def process_record_file(self):
        
        rec_list = []
        record_file = self.record_file 

        for macromolecule in record_file['sequences']:
            macromolecule_type  = list(macromolecule.keys())[0]
            
            asym_id_list = [macromolecule[macromolecule_type]['id']] if isinstance(macromolecule[macromolecule_type]['id'], str) else macromolecule[macromolecule_type]['id']

            for asym_id in asym_id_list:
                
                if macromolecule_type == 'protein':
                
                    record = macromolecule[macromolecule_type]
                    
                    protein_rec_info = PROTEIN(macromolecule_type, record).get_rec_info()
                    
                    rec_list += protein_rec_info
                
                elif macromolecule_type in ['dna', 'rna']:
                    
                    record = macromolecule[macromolecule_type]
                    
                    nucleotide_rec_info = NUCLEOTIDE(macromolecule_type, record).get_rec_info()
                    
                    rec_list += nucleotide_rec_info

                elif macromolecule_type  == 'ligand':
                                    
                    record = transform_rec(macromolecule[macromolecule_type], macromolecule_type, asym_id)
                    
                    non_poly_rec_info = NON_POLYMER(macromolecule_type, record).get_rec_info()
                    
                    rec_list += non_poly_rec_info
                else:
                    raise ValueError(f"Macromolecule type {macromolecule_type} is not supported")
        
        rec_list = self.map_asym_id(rec_list)
                  
        return rec_list

class RECORD_SERVER(RECORD):

    def __init__(self, 
                record_file,
                structure_sequence_list,
                feature_dict):
        
        super().__init__(
                record_file,
                structure_sequence_list,
                feature_dict)
        
    
    def process_record_file(self):
        
        rec_list = []
        record_file = self.record_file[0] 
        
        macromolecule_name_dict = {
            'proteinChain':'protein',
            'dnaSequence':'dna',
            'rnaSequence':'rna',
            'ligand':'ligand',
            'ion': 'ion'
        }

        for macromolecule in record_file['sequences']:
             
            macromolecule_type = list(macromolecule.keys())[0]
            for item in range(macromolecule[macromolecule_type]['count']):
                
                if macromolecule_type == 'proteinChain':
                
                    record = macromolecule[macromolecule_type]
                    
                    protein_rec_info = PROTEIN(macromolecule_name_dict[macromolecule_type], record).get_rec_info()
                    
                    rec_list += protein_rec_info
                
                elif macromolecule_type in ['dnaSequence', 'rnaSequence']:
                    
                    
                    record = macromolecule[macromolecule_type]
                    
                    nucleotide_rec_info = NUCLEOTIDE(macromolecule_name_dict[macromolecule_type], record).get_rec_info()
                    
                    rec_list += nucleotide_rec_info

                elif macromolecule_type in ['ligand', 'ion']:
                    
                    record = macromolecule[macromolecule_type]
                    
                    non_poly_rec_info = NON_POLYMER(macromolecule_name_dict[macromolecule_type], record).get_rec_info()
                    
                    rec_list += non_poly_rec_info
                else:
                    raise ValueError(f"Macromolecule type {macromolecule_type} is not supported")
        
        rec_list = self.map_asym_id(rec_list)
                  
        return rec_list

class POLYMER():
    
    def __init__(self, macromolecule_type):
        
        self.id = macromolecule_type

    def create_polymer_rec_dict(self, sequence):
        
        polymer_dict = {
            'rec_type' : 'polymer',
            'macromolecule_type'  : self.id,
            'sequence' : sequence,
            'entity_degree' : int(),
            'modifications' : [],
            'auth_asym_id' : str(),
            'label_asym_id' : str()
        }
        return polymer_dict
    
class PROTEIN(POLYMER):
    
    def __init__(self, macromolecule_type, record):
        
        super().__init__(macromolecule_type)
        
        self.record = record
        self.sequence =  record['sequence']
        
    def create_glycan_rec_dict(self, glycan, position):
        
        glycan_dict = {
            'rec_type' : 'non_polymer',
            'macromolecule_type'  : 'glycan',
            'non_poly_entity' : glycan,
            'entity_degree' : int(),
            'position': position,
            'auth_asym_id' : str(),
            'label_asym_id' : str()
        }
        
        return glycan_dict

     
    def extract_modification_info(self):
        
        record = self.record
    
        modification_info = {
            'ptm' : [],
            'glycans' : []
        }
        
        possible_modifications = {'modifications', 'glycans'}
        record_keys = set(record.keys())
        modification_list = list(set(possible_modifications) & set(record_keys))
        
        for modification_type in modification_list:
            
            if modification_type == 'modifications':
                
                ptm_list = record['modifications']
                
                modification_info['ptm'] = ptm_list
            
            elif modification_type == 'glycans':
                
                glycans_list = record['glycans']
                
                modification_info['glycans'] = glycans_list
            
        return modification_info
            
            
    def get_rec_info(self):
        
        sequence = self.sequence
        
        rec_info = []

        polymer_rec_dict = self.create_polymer_rec_dict(sequence)

        rec_info.append(polymer_rec_dict)
        modification_info = self.extract_modification_info()
        
        for modification, modification_list in modification_info.items():
            
            if not modification_list:
                continue
            
            else:
                
                if modification == 'ptm':
                
                    polymer_rec_dict['modifications'] = modification_list
            
                elif modification == 'glycans':
                    
                    glycan_rec_info = []
                    
                    for glycan in modification_list:
                        
                        glycan_id = glycan['residues']
                        glycan_position = glycan['position']
                        
                        glycan_rec_dict = self.create_glycan_rec_dict(glycan_id, glycan_position)
                        glycan_rec_info.append(glycan_rec_dict)
                        
                    rec_info += glycan_rec_info
                
        return rec_info

class NUCLEOTIDE(POLYMER):
    
    def __init__(self, macromolecule_type, record):
        
        super().__init__(macromolecule_type)
        
        self.record = record
        self.sequence =  record['sequence']
        
    def get_rec_info(self):
        
        sequence = self.sequence
        
        rec_info = []

        polymer_rec_dict = self.create_polymer_rec_dict(sequence)

        rec_info.append(polymer_rec_dict)
        
        if 'modifications' in self.record:
            
            if self.record['modifications']:
                raise NotImplementedError("FOUND MODIFIED NUCLEOTIDE; Modified nucleotides are not yet supported")

        return rec_info

class NON_POLYMER():

    def __init__(self, macromolecule_type, record):

        self.macromolecule_type = macromolecule_type
        self.record = record
        self.non_poly_entity = record[macromolecule_type].replace('CCD_', '') if record[macromolecule_type].startswith('CCD_') else record[macromolecule_type]  
            
    def create_rec_dict(self):
        
        non_polymer_dict = {
            'rec_type' : 'non_polymer',
            'macromolecule_type': self.macromolecule_type,
            'non_poly_entity' : self.non_poly_entity,
            'smiles' : str(),
            'entity_degree' : int(),
            'auth_asym_id' : str(),
            'label_asym_id' : str()
        }
        
        return non_polymer_dict

    def get_rec_info(self):
        rec_dict = self.create_rec_dict()
        if 'smiles' in self.record:
            rec_dict['smiles'] = self.record['smiles']
            
        rec_info = [rec_dict]
       
        return rec_info
    
def transform_rec(rec,macromolecule_type, asym_id):
    
    
    if 'smiles' in rec:
        ligand_rec = f'LIG_{asym_id}'
        rec_transformed ={
                            macromolecule_type: ligand_rec,
                            'smiles': rec['smiles'],
                            'count': 1
                        }
    elif 'ccdCodes' in rec:
        ligand_rec = f'{asym_id}'
        rec_transformed ={
                            macromolecule_type: rec['ccdCodes'][0],
                            'smiles': str(),
                            'count': 1
                        }
        
    else:
        raise NotImplemented("No support (yet) for this type of non-polymer entities in local AlphaFold3 version")
        
    return rec_transformed