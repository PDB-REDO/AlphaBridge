import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from itertools import chain

def calculate_probability_contact_link(coord_link, contact_probability):
    matrix_probability_link = contact_probability[np.s_[coord_link[0][0]:coord_link[0][1]+1], np.s_[coord_link[1][0]:coord_link[1][1]+1]]
    
    link_probability = np.mean([prob for prob in matrix_probability_link.flatten() if prob >= 0])
    return matrix_probability_link, link_probability

def calculate_pmc_link(coord_link, pae_plddt):
    matrix_pmc_link = pae_plddt[np.s_[coord_link[0][0]:coord_link[0][1]+1], np.s_[coord_link[1][0]:coord_link[1][1]+1]]
    
    norm_matrix_pmc_link = 1 - ( 3 * matrix_pmc_link / (2 * 100) )
    pmc_link = np.mean(norm_matrix_pmc_link.flatten())
    
    return norm_matrix_pmc_link, pmc_link
    
def calculate_probability_contact_interface(matrix_probability_interface):
    
    flatten_matrix_probability_interface = []
    filtered_flatten_matrix_probability_interface = []
    
    for matrix_probability_link in matrix_probability_interface:
        link_probability = [prob for prob in matrix_probability_link.flatten() if prob > 0]
        flatten_matrix_probability_interface += list(matrix_probability_link.flatten())
        filtered_flatten_matrix_probability_interface += link_probability
    
    quantile = np.quantile(filtered_flatten_matrix_probability_interface, 0.75)
    
    interface_probability = np.mean([prob for prob in filtered_flatten_matrix_probability_interface if prob >= quantile])
    
    contact_nr = len([prob for prob in flatten_matrix_probability_interface if prob >= quantile])
    #contact_ratio = contact_nr / len(flatten_matrix_probability_interface)
    
    return interface_probability, contact_nr , flatten_matrix_probability_interface

def calculate_pmc_interface(matrix_pmc_interface):
    
    flatten_matrix_pmc_interface = []
    
    for pmc_link in matrix_pmc_interface:
        flatten_matrix_pmc_interface += list(pmc_link.flatten())
    
    pmc_interface = np.mean(flatten_matrix_pmc_interface)
    
    return pmc_interface, flatten_matrix_pmc_interface
    
def calculate_interface_scores(interface_probability, pmc_interface, chain_pair_iptm):
    
    interface_score_iptm = np.sqrt(interface_probability * chain_pair_iptm)
    interface_score_pmc =  np.sqrt(interface_probability * pmc_interface)

    return interface_score_iptm, interface_score_pmc

def calculate_scores_stucture(probability_structure_list, pmc_structure_list, iptm):
    
    probability_contact_structure = float()
    pmc_structure = float()
    structure_score_iptm = float()
    structure_score_pmc = float()

    if probability_structure_list and pmc_structure_list:

        flattened_probability_structure = []
        flattened_pmc_structure_list = []
        for contact_submatrix, pmc_submatrix in zip(probability_structure_list, pmc_structure_list):
            
            flattened_probability_structure += [prob for prob in contact_submatrix.flatten() if prob > 0]
            
            norm_pmc_submatrix = 1 - ( 3 * pmc_submatrix / (2 * 100) )
            flattened_pmc_structure_list += list(norm_pmc_submatrix.flatten())
        
        
        quantile = np.quantile(flattened_probability_structure, 0.75)
        #print(quantile)
        probability_contact_structure = np.mean([prob for prob in flattened_probability_structure if prob >= quantile])
        
        pmc_structure = np.mean(flattened_pmc_structure_list)
        
        structure_score_iptm = np.sqrt(probability_contact_structure * iptm)
        structure_score_pmc = np.sqrt(probability_contact_structure * pmc_structure)
        
    #plot_probability_histplot(quantile,probability_contact_structure, flattened_probability_structure, flattened_pmc_structure_list)
    
    
    return probability_contact_structure, pmc_structure, structure_score_iptm, structure_score_pmc


def ptm_func(x, d0):
    return 1.0 / (1.0 + (x / d0) ** 2)

def calc_d0(L, pair_type='protein'):
    L = float(L)
    min_value = 2.0 if pair_type == 'nucleic_acid' else 1.0
    d0 = 1.24 * (L - 15.0) ** (1.0 / 3.0) - 1.8 if L > 27 else 1.0
    return max(min_value, d0)

def calc_d0_array(L_array, pair_type='protein'):
    L = np.maximum(np.array(L_array, dtype=float), 27.0)
    min_value = 2.0 if pair_type == 'nucleic_acid' else 1.0
    return np.maximum(min_value, 1.24 * (L - 15.0) ** (1.0 / 3.0) - 1.8)

def build_chains_array(sequence_info_dict):
    chains = []
    prev_acc = 0
    for chain_id, acc in zip(sequence_info_dict['label_asym_id'], sequence_info_dict['acclen']):
        chains.extend([chain_id] * (acc - prev_acc))
        prev_acc = acc
    return np.array(chains)

def _get_pair_type(chain1, chain2, chain_type_dict):
    t1 = chain_type_dict.get(chain1, 'protein')
    t2 = chain_type_dict.get(chain2, 'protein')
    return 'nucleic_acid' if t1 == 'nucleic_acid' or t2 == 'nucleic_acid' else 'protein'

def calculate_ipsae(pae, chains_array, chain_ids, chain_type_dict, pae_cutoff=10.0):
    """
    ipSAE: interface predicted SAE score (Dunbrack lab).
    Computes four variants per ordered chain pair (chain1 -> chain2):
      ipsae_d0chn  — d0 from combined chain length, PAE cutoff applied
      ipsae_d0dom  — d0 from interface-contributing residues only
      ipsae_d0res  — per-residue d0 based on each residue's contact count
      iptm_d0chn   — same d0 as d0chn but no PAE cutoff (all pairs)
    Returns {(chain1, chain2): {score_name: float}} for all ordered pairs.
    """
    results = {}

    for chain1 in chain_ids:
        for chain2 in chain_ids:
            if chain1 == chain2:
                continue

            idx1 = np.where(chains_array == chain1)[0]
            idx2 = np.where(chains_array == chain2)[0]
            if len(idx1) == 0 or len(idx2) == 0:
                continue

            pair_type = _get_pair_type(chain1, chain2, chain_type_dict)
            pae_sub = pae[np.ix_(idx1, idx2)]   # (n1, n2)
            valid   = pae_sub < pae_cutoff       # bool (n1, n2)

            # d0chn
            d0chn  = calc_d0(len(idx1) + len(idx2), pair_type)
            ptm_chn = ptm_func(pae_sub, d0chn)
            ipsae_d0chn = float(np.array([
                ptm_chn[k, valid[k]].mean() if valid[k].any() else 0.0
                for k in range(len(idx1))
            ]).mean())

            # d0dom
            n0dom  = int(valid.any(axis=1).sum() + valid.any(axis=0).sum())
            d0dom  = calc_d0(max(n0dom, 1), pair_type)
            ptm_dom = ptm_func(pae_sub, d0dom)
            ipsae_d0dom = float(np.array([
                ptm_dom[k, valid[k]].mean() if valid[k].any() else 0.0
                for k in range(len(idx1))
            ]).mean())

            # d0res
            n0res       = valid.sum(axis=1).astype(float)
            d0res       = calc_d0_array(n0res, pair_type)
            d0res_byres = np.zeros(len(idx1))
            for k in range(len(idx1)):
                if valid[k].any():
                    d0res_byres[k] = ptm_func(pae_sub[k], d0res[k])[valid[k]].mean()
            ipsae_d0res = float(d0res_byres.mean())

            # iptm_d0chn (no cutoff)
            iptm_d0chn = float(ptm_chn.mean())

            results[(chain1, chain2)] = {
                'ipsae_d0chn': ipsae_d0chn,
                'ipsae_d0dom': ipsae_d0dom,
                'ipsae_d0res': ipsae_d0res,
                'iptm_d0chn':  iptm_d0chn,
            }

    return results

def calculate_actifptm(pae, contact_probs, chains_array, chain_ids, chain_type_dict):
    """
    actifpTM: active-interface predicted TM-score (Steinegger/ColabFold).
    Adapted for AF3: uses exported PAE values directly instead of distogram logits.
    d0 is based on total token count (same as global ipTM normalization).
    contact_probs weights each pair by its soft contact probability.
    Returns {(chain1, chain2): float} for unique unordered pairs.
    """
    results = {}
    n_total = len(chains_array)
    inter_chain_mask = (chains_array[:, None] != chains_array[None, :]).astype(float)

    for i, chain1 in enumerate(chain_ids):
        for j in range(i + 1, len(chain_ids)):
            chain2 = chain_ids[j]

            idx1 = np.where(chains_array == chain1)[0]
            idx2 = np.where(chains_array == chain2)[0]
            if len(idx1) == 0 or len(idx2) == 0:
                results[(chain1, chain2)] = 0.0
                continue

            pair_type = _get_pair_type(chain1, chain2, chain_type_dict)
            d0 = calc_d0(max(n_total, 19), pair_type)

            # Build full-size pair weight matrix for this interface only
            pair_weights = np.zeros((n_total, n_total))
            pair_weights[np.ix_(idx1, idx2)] = contact_probs[np.ix_(idx1, idx2)]
            pair_weights[np.ix_(idx2, idx1)] = contact_probs[np.ix_(idx2, idx1)]

            total_weight = pair_weights.sum()
            if total_weight == 0.0:
                results[(chain1, chain2)] = 0.0
                continue

            ptm_matrix    = ptm_func(pae, d0) * inter_chain_mask
            residue_weights = (pair_weights.sum(axis=1) > 0).astype(float)
            residuewise   = (ptm_matrix * (pair_weights / total_weight)).sum(axis=1) * residue_weights

            results[(chain1, chain2)] = float(residuewise.max())

    return results


def plot_probability_histplot(quantile,probability_contact_structure,flattened_probability_structure, flattened_pmc_structure_list):
    
    #print(sorted(flatten_matrix_probability_interface),interface_probability, quantile)
    fig, ax = plt.subplots()
    g = sns.histplot(data=flattened_probability_structure, binwidth=0.01, cumulative = True, fill=False, stat='density')
    plt.axhline(y = 0.5, color = 'b')
    plt.axvline(x = probability_contact_structure, color = 'r')
    ax.set(xlim=(0,1))
    fig, ax = plt.subplots()
    g = sns.histplot(data=flattened_probability_structure, binwidth=0.01, cumulative = True, fill=False, stat='density')
    plt.axhline(y = 0.75, color = 'b')
    plt.axvline(x = probability_contact_structure, color = 'r')
    ax.set(xlim=(0,1))
    
    fig, ax = plt.subplots()
    g = sns.histplot(data=flattened_probability_structure, binwidth=0.01, fill=False, stat='density')
    plt.axvline(x = quantile, color = 'b')
    plt.axvline(x = probability_contact_structure, color = 'r')
    ax.set(xlim=(0,1))
    plt.show()
