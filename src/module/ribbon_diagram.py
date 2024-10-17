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


module_dir = os.path.dirname(os.path.realpath(__file__))




class RIBBON_DIAGRAM:
    
    def __init__(self, list_sequence_info, interface_dict, protein_interface_dict, plddt_dict: dict = {}, conservation_dict: dict = {}, outdir: str = ''):
        self.list_sequence_info = list_sequence_info
        self.interface_dict = interface_dict
        self.protein_interface_dict = protein_interface_dict
        self.plddt_dict = plddt_dict
        self.conservation_dict = conservation_dict
        self.outdir = outdir
    
    
    def create_ribbon_plot(self):
        
        list_fasta_name, list_fasta_acclen, list_fasta_centerticks, list_fasta_len = tuple(self.list_sequence_info)
        interfaces_list = list(self.interface_dict.keys())
        
        protein_interface_dict = self.protein_interface_dict
        interface_dict = self.interface_dict
        ''' COLOUR SCHEMA NEEDED'''
        
        cmap_plddt = colormaps['Spectral'] 
        norm_plddt = plt.Normalize(vmin=0, vmax=100)
        sm_plddt = plt.cm.ScalarMappable(cmap=cmap_plddt, norm=norm_plddt)    
        
        cmap_conservation = colormaps['RdBu_r'] 
        norm_conservation = plt.Normalize(vmin=0, vmax=1)
        sm_conservation = plt.cm.ScalarMappable(cmap=cmap_conservation, norm=norm_conservation)    
        
        interface2color = get_interface2color(interfaces_list)
        
        '''INITIALIZE CIRCOS PLOT'''
        
        sectors = {fasta_name :fasta_len for fasta_name, fasta_len in zip(list_fasta_name,list_fasta_len)}
        
        circos = Circos(sectors, space = 5)
        
        for sector in circos.sectors:
            tracks_position_list = [(75, 85), (88,93), (95, 100)]
            
            plddt_list = self.plddt_dict[sector.name]
        
            if not self.conservation_dict:
                
                conservation_list = []
                tracks_position_list = tracks_position_list[:-1]
                
            else:
                
                conservation_list = self.conservation_dict[sector.name]
                
                circos.colorbar(bounds=(1.35, 0.29, 0.02, 0.5), vmin=0, vmax=1, cmap=cmap_conservation,
                colorbar_kws=dict(label="AlphaMissense"),
                tick_kws=dict(labelsize=8, labelrotation=0),)
            
            for track_position in tracks_position_list:
                #add tracks
                track = sector.add_track(track_position)
                #add track-axis, ticksa and text
                track.axis()
                    
            track.xticks_by_interval(100)
            track.text(sector.name, color="black", size=9, r = tracks_position_list[-1][1] + 10)
            
            for i in range(int(sector.size)):

                plddt_color = rgb2hex(sm_plddt.to_rgba(plddt_list[i])[:3])
                plddt_color = get_colour_plddt(plddt_list[i])            
                sector.rect(start=i, end=i + 1, r_lim= tracks_position_list[1], color = plddt_color , lw=0)
                
                if conservation_list:
                    conservation_color = rgb2hex(sm_conservation.to_rgba(conservation_list[i])[:3])
                    sector.rect(start=i, end=i + 1, r_lim=tracks_position_list[2], color = conservation_color , lw=0)
                
                
            for interface_name in interfaces_list:
                if sector.name in protein_interface_dict:
                    
                    if interface_name in protein_interface_dict[sector.name]:
                        
                        interface_range_list = protein_interface_dict[sector.name][interface_name]['interface_range']
                        
                        for interface_range in interface_range_list:
                            
                            degree_range = [degrees(sector.x_to_rad(residue_number - 1)) for residue_number in interface_range] 
                            circos.rect(r_lim=(75, 85), deg_lim=(degree_range[0], degree_range[1]),fc=interface2color[interface_name], ec="black", lw=0.5)

        for interface in interface_dict:
    
            interface_dict[interface]
            
            prot_1 = interface_dict[interface]['prot_1']['accesion_id']
            prot_2 = interface_dict[interface]['prot_2']['accesion_id']
            
            color = interface2color[interface]
            
            links =  interface_dict[interface]['links']
            
            for link in links:
                
                link_prot_1 = link[0]
                link_prot_2 =  link[1]
                
                
                circos.link((prot_1, link_prot_1[0]-1, link_prot_1[1]-1),(prot_2, link_prot_2[0]-1, link_prot_2[1]-1), 
                            color=color, alpha = 0.25)
        
        return circos
    
    
    def plot_ribbon_diagram(self):
        
        circos = self.create_ribbon_plot()
        outdir = self.outdir
        filename = f'{outdir}/ribbon_plot.png'
        
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
            title="Model confidance",
            ncol=1,
        )
        
        fig.savefig(filename)
                
        
        
        
def get_interface2color(interfaces_list):
    
        cmap = colormaps['tab20']  # matplotlib color palette name, n colors
            
        color_list = [rgb2hex(cmap(i)[:3]) for i in range(cmap.N)]

        reord_color_list = color_list[::2] + color_list[1::2]

        interface2color = {name:reord_color_list[index]  for index,name in enumerate(interfaces_list)}
        
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