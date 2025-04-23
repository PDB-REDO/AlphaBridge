import matplotlib.pyplot as plt 
import matplotlib.colors as mcolors
import numpy as np

#function for keeping the color constant 
#retun a dictionar --> protein : color
def get_label_color_dict(dit):
    label_list = []
    for aut, label in dit.items():
        label_list.append(label)
    num_proteins = len(label_list)
    cmap = plt.get_cmap("rainbow") 
    colors_list = [mcolors.rgb2hex(cmap(i)) for i in np.linspace(0, 1, num_proteins)]
    label_color_dict = dict(zip(label_list, colors_list))
    return label_color_dict