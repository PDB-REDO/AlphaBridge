import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

def calculate_probability_contact_link(coord_link, contact_probability):
    matrix_probability_link = contact_probability[np.s_[coord_link[0][0]:coord_link[0][1]+1], np.s_[coord_link[1][0]:coord_link[1][1]+1]]
    link_probability = np.mean([prob for prob in matrix_probability_link.flatten() if prob >= 0.1])
    return matrix_probability_link, link_probability
    

def calculate_probability_contact_interface(matrix_probability_interface):
    flatten_matrix_probability_interface = []
    
    for matrix_probability_link in matrix_probability_interface:
        link_probability = [prob for prob in matrix_probability_link.flatten() if prob > 0]
        flatten_matrix_probability_interface += link_probability
    quantile = np.quantile(flatten_matrix_probability_interface, 0.5)
    
    interface_probability = np.mean([prob for prob in flatten_matrix_probability_interface if prob >= quantile])
    contact_nr = len([prob for prob in flatten_matrix_probability_interface if prob >= quantile])
    contact_ratio = contact_nr / len(flatten_matrix_probability_interface)
    
    return interface_probability, flatten_matrix_probability_interface

def plot_probability_histplot(interface_probability,flatten_matrix_probability_interface):
    
    quantile = np.quantile(flatten_matrix_probability_interface, 0.5)
    #print(sorted(flatten_matrix_probability_interface),interface_probability, quantile)
    fig, ax = plt.subplots()
    g = sns.histplot(data=flatten_matrix_probability_interface, binwidth=0.01)
    plt.axvline(x = quantile, color = 'b')
    ax.set(xlim=(0,1))
    plt.show()
