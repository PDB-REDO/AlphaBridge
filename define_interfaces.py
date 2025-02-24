import os
import sys
import json
import pandas as pd
import networkx as nx
import numpy as np

from src.module.confidence_contact_matrix import CCM_AF3
from src.module.alingment_utils import compare_protein_seq
from src.module.domain_clustering import domain_clustering
from src.module.parsers import MMCIFPARSER, HSSPPARSER, alphafold_msa
#from src.module.conservation_score import CONSERVATION_SCORE
from src.module.interface_identification import interface_identification
from src.module.ribbon_diagram import RIBBON_DIAGRAM

import argparse 

working_dir = os.path.dirname(os.path.realpath(__file__))


def parse_args():
    #####################
    # START CODING HERE #
    #####################
    # Implement a simple argument parser (WITH help documentation!) that parses
    # the information needed by main() from commandline. 

    parser = argparse.ArgumentParser()
    
    if len(sys.argv)==1:
        parser.print_help()
        # parser.print_usage() # for just the usage line
        parser.exit()

    parser.add_argument('-i', dest='in_dir', default='',
                        help='path to a directory where input folder are stored')
    
    parser.add_argument('-c', dest='config_dir', default='',
                        help='path to a directory where input folder are stored')
    
    parser.add_argument('-m','--mode', dest='mode', choices=['AF3', 'AF2', 'ColabFold'] , default='AF3',
                        help='output from different AlphaFold Version. Options: AF3, AF2, ColabFold')
    
    #parser.add_argument('-t','--threshold', dest='contact_threshold' , default=0.7, type=restricted_float,
                        #help='contact threshold to detect a contact-link in the contact_proability matrix')
    
    # parser.add_argument(?)
    # parser.add_argument(?)
    # parser.add_argument(?)

    args = parser.parse_args()

    return args

def restricted_float(x):
    try:
        x = float(x)
    except ValueError:
        raise argparse.ArgumentTypeError("%r not a floating-point literal" % (x,))

    if x < 0.0 or x > 1.0:
        raise argparse.ArgumentTypeError("%r not in range [0.0, 1.0]"%(x,))
    return x

def write_dataframe(df, filename, outdir_path):
    
    filepath= os.path.join(outdir_path, filename)
    
    df.to_csv(f'{filepath}.csv', index = False)



def define_interfaces(in_dir, mode):
        
        
    outdir = os.path.join(in_dir, 'AlphaBridge')

    if not os.path.isdir(outdir):
        os.makedirs(outdir)
            

    if mode == 'AF3':
        
            FEATURE_OBJECT = CCM_AF3(in_dir)
            
            feature_path, structure_path, job_request_path, summary_request_path = FEATURE_OBJECT.extract_feature_filepath()
            chain_info_dict, sequence_info_dict = FEATURE_OBJECT.extract_chain_info_dict()
            
    else:
        raise  NotImplementedError("Output from AF2 or ColabFold not implemented yet")
    
    matrix_dict = FEATURE_OBJECT.extract_matrix_dict()
    contact_matrix =  matrix_dict['contact_matrix']
    confidance_matrix = matrix_dict['pae_plddt']
    iptm = matrix_dict['iptm']
    chain_pair_iptm_matrix = matrix_dict['chain_pair_iptm'] 


    coevolutionary_domains, coevolutionary_cluster_dict, entity_region_dict = domain_clustering(matrix_dict,
                                                                                sequence_info_dict,
                                                                                alphafold_version=mode,
                                                                                outdir = outdir, 
                                                                                plotting=True).run_domain_clustering()
    
    
    #elements = np.linspace(0.4, 1, 40).tolist()
    
    #contact_threshold_list = [round(x, 3) for x in elements]
    contact_threshold = [0.5, 0.75, 0.9]
    interactions_list = []
    
    
    for contact_threshold in contact_threshold_list:

        INTERFACE_IDENTIFICATION = interface_identification(coevolutionary_cluster_dict, 
                                                            entity_region_dict,
                                                            chain_info_dict,
                                                            sequence_info_dict,
                                                            iptm,
                                                            chain_pair_iptm_matrix,
                                                            confidance_matrix,
                                                            contact_matrix,
                                                            contact_threshold)

        interactions_dict= INTERFACE_IDENTIFICATION.extract_interfaces()
        interactions_list.append(interactions_dict)
    
        biomolecule_interface_dict= INTERFACE_IDENTIFICATION.map_info_interfaces(interactions_dict)
    
        #interface_info_df = INTERFACE_IDENTIFICATION.get_interface_info_dataframes(interactions_dict)
    
    
        ribbon_diagram = RIBBON_DIAGRAM(
                        interactions_dict,
                        biomolecule_interface_dict,
                        chain_info_dict,
                        contact_threshold,
                        outdir=outdir,
                        boolean_modified_non_poly_length= True)
        ribbon_diagram.plot_ribbon_diagram()
        


    structure_score_dict = INTERFACE_IDENTIFICATION.get_structure_score_dict(chain_info_dict)

    #structure_info_df = INTERFACE_IDENTIFICATION.get_structure_info_dataframes(structure_score_dict)

    alphabridge_dict = {
        "structure": [structure_score_dict],
        "interactions" : interactions_list
    }


    with open(f"{outdir}/alphabridge_data.json", "w") as file:
        file.write(json.dumps(alphabridge_dict, indent=4))
    
    #write_dataframe(structure_info_df, 'structure_scores', outdir )
        

def main():
    args = parse_args()
    
    in_dir = args.in_dir
    mode = args.mode
    #contact_threshold = args.contact_threshold
    
    define_interfaces(in_dir, mode)
    
    print('finished')
    
    
   
if __name__ == '__main__':
    main()