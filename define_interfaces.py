import os
import sys
import json
import pandas as pd
import networkx as nx

from src.module.confidance_contact_matrix import CCM_AF3
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
    
    parser.add_argument('-t','--threshold', dest='contact_threshold' , default=0.7, type=restricted_float,
                        help='contact threshold to detect a contact-link in the contact_proability matrix')
    
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

def define_interfaces(in_dir, mode, contact_threshold):
    
    if mode == 'AF3':
     
        FEATURE_OBJECT = CCM_AF3(in_dir)
        
        feature_path, structure_path, job_request_path, summary_request_path = FEATURE_OBJECT.extract_feature_filepath()
        
        list_sequence_info, rec_sequence_list, structure_sequence_list, polymer_chain_dict = FEATURE_OBJECT.extract_sequence_info()
        
    else:
        raise  NotImplementedError("Output from AF2 or ColabFold not implemented yet")

    chain_dict = compare_protein_seq(structure_sequence_list, rec_sequence_list).extract_chain_dict()
    
    matrix_dict = FEATURE_OBJECT.extract_matrix_dict()
    
    plddt_dict = FEATURE_OBJECT.get_scores_dict(matrix_dict['plddt'], list_sequence_info)
    
    contact_matrix =  matrix_dict['contact_matrix']
    
    
    
    FEATURE_OBJECT.print_matrix_dict(matrix_dict)
    
    coevolutionary_domains, interacting_coevolutionary_domains, entity_region_dict = domain_clustering(matrix_dict,
                                                                               list_sequence_info,
                                                                               alphafold_version=mode,
                                                                               outdir = in_dir, 
                                                                               plotting=True).run_domain_clustering()
    
    
    INTERFACE_IDENTIFICATION = interface_identification(interacting_coevolutionary_domains, 
                                                        entity_region_dict,
                                                        plddt_dict,
                                                        rec_sequence_list,
                                                        list_sequence_info,
                                                        contact_matrix,
                                                        contact_threshold,
                                                        chain_dict,
                                                        polymer_chain_dict)
    
    interface_dict, protein_interface_dict, interaction_link_dict = INTERFACE_IDENTIFICATION.extract_interface()   
    
    interface_df_per_token, interface_df = INTERFACE_IDENTIFICATION.get_interface_info_dataframes(interface_dict,interaction_link_dict)
    
    ribbon_diagram = RIBBON_DIAGRAM(list_sequence_info,
                   interface_dict,
                   protein_interface_dict,
                   plddt_dict,
                   outdir=in_dir)
    
    ribbon_diagram.plot_ribbon_diagram()
    
    
    
    
    
    return interface_df_per_token, interface_df
    
    
    

def main():
    args = parse_args()
    
    in_dir = args.in_dir
    mode = args.mode
    contact_threshold = args.contact_threshold
    
    interface_df_per_token, interface_df = define_interfaces(in_dir, mode, contact_threshold)
    
    write_dataframe(interface_df_per_token, 'interface_df_per_token', in_dir )
    write_dataframe(interface_df, 'binding_interfaces', in_dir )
    
    print('finished')
    
    
   
if __name__ == '__main__':
    main()