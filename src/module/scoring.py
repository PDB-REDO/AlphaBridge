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

            # d0chn — aggregate = max over chain1 residues (matches Dunbrack lab implementation)
            d0chn   = calc_d0(len(idx1) + len(idx2), pair_type)
            ptm_chn = ptm_func(pae_sub, d0chn)
            byres_d0chn = np.array([
                ptm_chn[k, valid[k]].mean() if valid[k].any() else 0.0
                for k in range(len(idx1))
            ])
            ipsae_d0chn = float(byres_d0chn.max())

            # d0dom
            n0dom   = int(valid.any(axis=1).sum() + valid.any(axis=0).sum())
            d0dom   = calc_d0(max(n0dom, 1), pair_type)
            ptm_dom = ptm_func(pae_sub, d0dom)
            byres_d0dom = np.array([
                ptm_dom[k, valid[k]].mean() if valid[k].any() else 0.0
                for k in range(len(idx1))
            ])
            ipsae_d0dom = float(byres_d0dom.max())

            # d0res
            n0res       = valid.sum(axis=1).astype(float)
            d0res       = calc_d0_array(n0res, pair_type)
            byres_d0res = np.zeros(len(idx1))
            for k in range(len(idx1)):
                if valid[k].any():
                    byres_d0res[k] = ptm_func(pae_sub[k], d0res[k])[valid[k]].mean()
            ipsae_d0res = float(byres_d0res.max())

            # iptm_d0chn (no cutoff) — mean over all chain2 residues per chain1 residue, then max
            iptm_d0chn = float(ptm_chn.mean(axis=1).max())

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
    # Use PAE matrix size for d0 (total tokens including non-polymer atoms)
    n_total = pae.shape[0]

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

            n1, n2 = len(idx1), len(idx2)

            # Extract sub-matrices for this chain pair only — avoids any
            # dimension mismatch between chains_array and the full PAE
            pae_12   = pae[np.ix_(idx1, idx2)]            # (n1, n2)
            pae_21   = pae[np.ix_(idx2, idx1)]            # (n2, n1)
            cmap_12  = contact_probs[np.ix_(idx1, idx2)]  # (n1, n2)
            cmap_21  = contact_probs[np.ix_(idx2, idx1)]  # (n2, n1)

            # Combined (n1+n2) interface weight matrix, off-diagonal blocks only
            pair_weights = np.zeros((n1 + n2, n1 + n2))
            pair_weights[:n1, n1:] = cmap_12
            pair_weights[n1:, :n1] = cmap_21

            total_weight = pair_weights.sum()
            if total_weight == 0.0:
                results[(chain1, chain2)] = 0.0
                continue

            # TM scores for both directions
            ptm_combined = np.zeros((n1 + n2, n1 + n2))
            ptm_combined[:n1, n1:] = ptm_func(pae_12, d0)
            ptm_combined[n1:, :n1] = ptm_func(pae_21, d0)

            residue_weights = (pair_weights.sum(axis=1) > 0).astype(float)
            residuewise = (ptm_combined * (pair_weights / total_weight)).sum(axis=1) * residue_weights

            results[(chain1, chain2)] = float(residuewise.max())

    return results


def calculate_ipsae_complex(pae, chains_array, chain_ids, chain_type_dict, pae_cutoff=10.0):
    """
    Complex-level ipSAE: all inter-chain pairs across all polymer chains are
    considered simultaneously. For each residue i, the per-residue score is
    the mean ptm over all valid inter-chain partners from ANY other chain
    (pae < cutoff). The final score is the max over all polymer residues.

    Three d0 variants are reported:
      ipsae_d0chn  — d0 from total polymer residues in the complex
      ipsae_d0dom  — d0 from residues that have at least one valid inter-chain pair
      ipsae_d0res  — per-residue d0 from each residue's total valid inter-chain contacts
    """
    idx_list     = [np.where(chains_array == c)[0] for c in chain_ids]
    all_idx      = np.concatenate(idx_list)
    chain_labels = np.concatenate([[c] * len(ix) for c, ix in zip(chain_ids, idx_list)])
    n_total      = len(all_idx)

    if n_total == 0:
        return {'ipsae_d0chn': 0.0, 'ipsae_d0dom': 0.0,
                'ipsae_d0res': 0.0, 'iptm_d0chn': 0.0}

    pair_type = 'nucleic_acid' if any(
        chain_type_dict.get(c) == 'nucleic_acid' for c in chain_ids
    ) else 'protein'

    pae_sub = pae[np.ix_(all_idx, all_idx)]                         # (N, N)
    inter   = chain_labels[:, None] != chain_labels[None, :]        # bool (N, N)
    valid   = inter & (pae_sub < pae_cutoff)                        # bool (N, N)

    # d0chn: d0 from total polymer residues
    d0chn   = calc_d0(n_total, pair_type)
    ptm_chn = ptm_func(pae_sub, d0chn)
    byres_d0chn = np.array([
        ptm_chn[k, valid[k]].mean() if valid[k].any() else 0.0
        for k in range(n_total)
    ])
    ipsae_d0chn = float(byres_d0chn.max())

    # d0dom: d0 from residues with at least one valid inter-chain pair
    n0dom   = int(valid.any(axis=1).sum())
    d0dom   = calc_d0(max(n0dom, 1), pair_type)
    ptm_dom = ptm_func(pae_sub, d0dom)
    byres_d0dom = np.array([
        ptm_dom[k, valid[k]].mean() if valid[k].any() else 0.0
        for k in range(n_total)
    ])
    ipsae_d0dom = float(byres_d0dom.max())

    # d0res: per-residue d0 from each residue's total valid inter-chain contacts
    n0res = valid.sum(axis=1).astype(float)
    d0res = calc_d0_array(n0res, pair_type)
    byres_d0res = np.zeros(n_total)
    for k in range(n_total):
        if valid[k].any():
            byres_d0res[k] = ptm_func(pae_sub[k], d0res[k])[valid[k]].mean()
    ipsae_d0res = float(byres_d0res.max())

    # iptm_d0chn: no cutoff
    iptm_d0chn = float((ptm_chn * inter).sum(axis=1).max() / max(inter.sum(axis=1).max(), 1))

    return {
        'ipsae_d0chn': ipsae_d0chn,
        'ipsae_d0dom': ipsae_d0dom,
        'ipsae_d0res': ipsae_d0res,
        'iptm_d0chn':  iptm_d0chn,
    }


def calculate_actifptm_complex(pae, contact_probs, chains_array, chain_ids, chain_type_dict):
    """
    actifpTM for the whole complex: all inter-chain contacts across all polymer
    chains are considered simultaneously.
    Returns a single float.
    """
    idx_list     = [np.where(chains_array == c)[0] for c in chain_ids]
    all_idx      = np.concatenate(idx_list)
    chain_labels = np.concatenate([[c] * len(ix) for c, ix in zip(chain_ids, idx_list)])
    n_complex    = len(all_idx)

    if n_complex == 0:
        return 0.0

    pair_type = 'nucleic_acid' if any(
        chain_type_dict.get(c) == 'nucleic_acid' for c in chain_ids
    ) else 'protein'
    d0 = calc_d0(max(pae.shape[0], 19), pair_type)

    inter        = (chain_labels[:, None] != chain_labels[None, :]).astype(float)
    pae_sub      = pae[np.ix_(all_idx, all_idx)]
    cmap_sub     = contact_probs[np.ix_(all_idx, all_idx)]
    pair_weights = cmap_sub * inter

    total_weight = pair_weights.sum()
    if total_weight == 0.0:
        return 0.0

    ptm_sub         = ptm_func(pae_sub, d0) * inter
    residue_weights = (pair_weights.sum(axis=1) > 0).astype(float)
    residuewise     = (ptm_sub * (pair_weights / total_weight)).sum(axis=1) * residue_weights

    return float(residuewise.max())


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
