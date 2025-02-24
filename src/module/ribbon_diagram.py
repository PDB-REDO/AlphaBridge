# -*- coding: utf-8 -*-
"""
Created on Wed Feb 7 16:53:05 2024

@author: Dan_salv
"""
import os
from pycirclize import Circos
from math import degrees
from matplotlib.patches import Patch
from matplotlib.colors import rgb2hex
from matplotlib import colormaps
import matplotlib.pyplot as plt
import distinctipy
import time
import numpy as np


module_dir = os.path.dirname(os.path.realpath(__file__))




class RIBBON_DIAGRAM:
    
    def __init__(self, interactions_dict, biomolecule_interface_dict, chain_info_dict ,contact_threshold, conservation_dict: dict = {}, outdir: str = '', boolean_modified_non_poly_length: bool = False):
        
        self.interactions_dict = interactions_dict
        self.biomolecule_interface_dict = biomolecule_interface_dict
        self.chain_info_dict = chain_info_dict
        self.contact_threshold = contact_threshold
        self.conservation_dict = conservation_dict
        self.boolean_modified_non_poly_length = boolean_modified_non_poly_length
        self.outdir = outdir
    
    
    def get_monomer_info_dict(self):
        
        chain_info_dict = self.chain_info_dict
        
        monomer_name_len_dict = {}
        poly_type_dict = {}

        for entity_type, rec_list in chain_info_dict.items():

            for rec in rec_list:
                
                auth_asym_id = rec['auth_asym_id']
                
                if not auth_asym_id in monomer_name_len_dict:
                    monomer_name_len_dict[auth_asym_id] = int()
                    poly_type_dict[auth_asym_id] = str()
                        
                poly_type_dict[auth_asym_id] = rec['macromolecule_type']
                
                if entity_type == 'polymer':
                    
                    entity_length = len(rec['residues'])
                else:
                    entity_length = len(rec['atoms'])
                
                monomer_name_len_dict[auth_asym_id] = entity_length
                            
        return monomer_name_len_dict, poly_type_dict 
    
    def get_modified_monomer_info_dict(self):
        
        chain_info_dict = self.chain_info_dict
        
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
                            monomer_name_len_dict[auth_asym_id] = int()
                            poly_type_dict[auth_asym_id] = str()
                        
                        poly_type_dict[auth_asym_id] = rec['macromolecule_type']
        
                        if entity_type == 'polymer':
                        
                            entity_length = len(rec['residues'])
                            monomer_name_len_dict[auth_asym_id] = entity_length
                        else:
                            entity_length = len(rec['atoms'])
                            poly_bolean = True
                            
                            monomer_name_len_dict[auth_asym_id] = None
        
                            if rec['macromolecule_type'] == 'ligand':
                                ligand_nr +=1
                                total_ligand_token += entity_length
                                
                            elif rec['macromolecule_type'] == 'ion':
                                ion_nr += 1
                                total_ion_token += entity_length
                            elif rec['macromolecule_type'] == 'glycan':
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
                
                if type == 'ligand':
                    monomer_name_len_dict[auth_asym_id] = int(ligand_length)
                elif type == 'ion':
                    monomer_name_len_dict[auth_asym_id] = int(ion_length)
                elif type == 'glycan':
                    monomer_name_len_dict[auth_asym_id] = int(glycan_length)
                    
                            
        return monomer_name_len_dict, poly_type_dict   
    
    def get_plddt_dict(self):
        plddt_dict = {}
        chain_info_dict = self.chain_info_dict
        for entity_type, rec_list in chain_info_dict.items():

            for rec in rec_list:
                
                auth_asym_id = rec['auth_asym_id']
                
                if not auth_asym_id in plddt_dict:
                    plddt_dict[auth_asym_id] = []
                
                if entity_type == 'polymer':
                    plddt_list = [residue['plddt'] for residue in rec['residues']]
                else:
                    plddt_list = [atom['plddt'] for atom in rec['atoms']]
                    
                
                plddt_dict[auth_asym_id] = plddt_list
        
        return plddt_dict
        

    def create_ribbon_plot(self):
        
    
        chain_info_dict = self.chain_info_dict
        biomolecule_interface_dict = self.biomolecule_interface_dict
        interactions_dict = self.interactions_dict
       
        plddt_dict = self.get_plddt_dict()
        label2auth = get_label2auth(chain_info_dict)
        auth2label = get_auth2label(chain_info_dict)
        
         
        if not self.boolean_modified_non_poly_length:
            sectors, poly_type_dict =  self.get_monomer_info_dict()
        else:
            sectors, poly_type_dict =  self.get_modified_monomer_info_dict()
        
        interfaces_list = [interface['interface_id'] for interface in interactions_dict['interfaces']]
        
        
        #biomolecule_interface_dict = self.biomolecule_interface_dict
        #interface_dict = self.interface_dict
        ''' COLOUR SCHEMA NEEDED'''
        
        cmap_plddt = colormaps['Spectral'] 
        norm_plddt = plt.Normalize(vmin=0, vmax=100)
        sm_plddt = plt.cm.ScalarMappable(cmap=cmap_plddt, norm=norm_plddt)    
        
        cmap_conservation = colormaps['RdBu_r'] 
        norm_conservation = plt.Normalize(vmin=0, vmax=1)
        sm_conservation = plt.cm.ScalarMappable(cmap=cmap_conservation, norm=norm_conservation)    
        
        start_time = time.time()
        interface2color = get_interface2color(interfaces_list)
        
        end_time = time.time()
        elapsed_time = end_time - start_time
        #print(f"colour {elapsed_time:.6f} seconds")
        '''INITIALIZE CIRCOS PLOT'''
        

        
        circos = Circos(sectors, space= 0.75)
        for sector in circos.sectors:
            
            auth_asym_id = sector.name
            label_asym_id = auth2label[auth_asym_id]
            
            poly_type = poly_type_dict[auth_asym_id]
            
            tracks_position_list = [(75, 85), (88,93), (95, 100)]
            
            plddt_list = plddt_dict[auth_asym_id]
            
            if poly_type_dict[auth_asym_id] in ['ligand' , 'ion', 'glycan'] and self.boolean_modified_non_poly_length:
                
                plddt_value = np.mean(plddt_list)
                
                plddt_list = [plddt_value] * sectors[auth_asym_id]

            if not self.conservation_dict:
                
                conservation_list = []
                tracks_position_list = tracks_position_list[:-1]
                
            else:
                
                conservation_list = self.conservation_dict[auth_asym_id]
                
                circos.colorbar(bounds=(1.35, 0.29, 0.02, 0.5), vmin=0, vmax=1, cmap=cmap_conservation,
                colorbar_kws=dict(label="AlphaMissense"),
                tick_kws=dict(labelsize=8, labelrotation=0),)
            
            for track_position in tracks_position_list:
                #add tracks
                track = sector.add_track(track_position)
                #add track-axis, ticksa and text
                track.axis()
                    
            track.xticks_by_interval(100)
            track.text(auth_asym_id, color="black", size=9, r = tracks_position_list[-1][1] + 10)
            
            for i in range(int(sector.size)):

                plddt_color = rgb2hex(sm_plddt.to_rgba(plddt_list[i])[:3])
                plddt_color = get_colour_plddt(plddt_list[i])            
                sector.rect(start=i, end=i + 1, r_lim= tracks_position_list[1], color = plddt_color , lw=0)
                
                if conservation_list:
                    conservation_color = rgb2hex(sm_conservation.to_rgba(conservation_list[i])[:3])
                    sector.rect(start=i, end=i + 1, r_lim=tracks_position_list[2], color = conservation_color , lw=0)
                
                
            for interface_id in interfaces_list:
                if label_asym_id in biomolecule_interface_dict:
                    
                    if interface_id in biomolecule_interface_dict[label_asym_id]:
                        #check if label_asym id is a ligand or an ion, and boolean_modified_non_poly_length
                        if poly_type in ['ligand' , 'ion', 'glycan'] and self.boolean_modified_non_poly_length:
                            
                            interface_range = [1, sectors[auth_asym_id]]
                            
                            #print(interface_range)
                            
                            degree_range = [degrees(sector.x_to_rad(residue_number - 1)) for residue_number in interface_range] 
                            circos.rect(r_lim=(75, 85), deg_lim=(degree_range[0], degree_range[1]),fc=interface2color[interface_id], ec="black", lw=0.5)
                           
                        else:
                            
                            interface_range_list = biomolecule_interface_dict[label_asym_id][interface_id]
                            
                            for interface_range in interface_range_list:
                                
                                degree_range = [degrees(sector.x_to_rad(residue_number - 1)) for residue_number in interface_range] 
                                circos.rect(r_lim=(75, 85), deg_lim=(degree_range[0], degree_range[1]),fc=interface2color[interface_id], ec="black", lw=0.5)
                               

            
        
        for interface in interactions_dict['interfaces']:
    
            interface_id = interface['interface_id']
            
            biomolecule_1 = interface['links'][0]['first']['asym_id']
            biomolecule_2 = interface['links'][0]['second']['asym_id']
            
            
            color = interface2color[interface_id]
            
            
            for link in interface['links']:
                
                interaction_1 = link['first']
                interaction_2 = link['second']
                bridge_1, bridge_2 = (extract_bridge(interaction, 
                                                     label2auth,
                                                     sectors, 
                                                     poly_type_dict, 
                                                     self.boolean_modified_non_poly_length) for interaction in (interaction_1, interaction_2))
                                
                
                circos.link(bridge_1,bridge_2, 
                            color=color, alpha = 0.25)
        
        return circos
    
    
    def plot_ribbon_diagram(self):
        
        start_time = time.time()
        
        circos = self.create_ribbon_plot()
        
        end_time = time.time()
        
        elapsed_time = end_time - start_time
        #print(f"create ribbon plot {elapsed_time:.6f} seconds")
        
        outdir = self.outdir
        contact_threshold = self.contact_threshold
        filename = f'{outdir}/{contact_threshold}_ribbon_plot.png'
        
        fig = circos.plotfig()
        plddt_color_list = ['#0053d6','#65cbf3','#ffdb13', '#ff7d45']

        plddt_label = ['Very high','High','Low', 'Very Low']

        rect_handles = []
        for idx, color in enumerate(plddt_color_list):
            rect_handles.append(Patch(color=color, label=plddt_label[idx]))
        _ = circos.ax.legend(
            handles=rect_handles,
            bbox_to_anchor=(1.2, 0.55),
            loc="center",
            fontsize=8,
            title="Model confidence",
            ncol=1,
        )
        
        fig.savefig(filename)
                
        
        
        
def get_interface2color(interfaces_list):

        interface_nr = len(interfaces_list)
        #cmap = colormaps['tab20']  # matplotlib color palette name, n colors
        cmap = distinctipy.get_colors(interface_nr)    
        #color_list = [rgb2hex(cmap(i)[:3]) for i in range(cmap.N)]
        #reord_color_list = color_list[::2] + color_list[1::2]
        color_list = [rgb2hex(rgb) for rgb in cmap]
        #interface2color = {name:reord_color_list[index]  for index,name in enumerate(interfaces_list)}
        interface2color = {name:color_list[index]  for index,name in enumerate(interfaces_list)}
        return interface2color
             
        

def get_colour_plddt(plddt_value):
    
    if plddt_value < 50 :
        return '#ff7d45'
    elif 50 <= plddt_value < 70:
        return '#ffdb13'
    elif 70 <= plddt_value < 90:
        return '#65cbf3'
    else:
        return '#0053d6'
    

        
def get_label2auth(chain_info_dict):
    label2auth = {rec['label_asym_id']: rec['auth_asym_id'] for entity_type, rec_list in chain_info_dict.items() for rec in rec_list}
    return label2auth

def get_auth2label(chain_info_dict):
    
    auth2label = {rec['auth_asym_id']: rec['label_asym_id'] for entity_type, rec_list in chain_info_dict.items() for rec in rec_list}
    return auth2label

def extract_bridge(interaction, label2auth,monomer_length_dict,poly_type_dict, boolean_modified_non_poly_length):
    
    label_asym_id = interaction['asym_id'] 
    auth_asym_id = label2auth[label_asym_id] 
    
    if poly_type_dict[auth_asym_id] in ['ligand', 'ion', 'glycan'] and boolean_modified_non_poly_length:
        
        start = 0
        end = monomer_length_dict[auth_asym_id] - 1
    else:    
        start = interaction['link_range']['start'] - 1
        end = interaction['link_range']['end'] - 1
    
    
    return auth_asym_id, start, end