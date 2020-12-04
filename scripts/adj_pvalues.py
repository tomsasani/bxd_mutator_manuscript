import scipy.stats as ss
import numpy as np
import math

def adj_pvalues(pvals, sub_0_counts, sub_1_counts):

    sorted_pval_idx = np.unravel_index(np.argsort(pvals, axis=None), pvals.shape)
    sorted_pval_idx = np.vstack(sorted_pval_idx)

    adj_pvals = np.ones(pvals.shape)

    # procedure to adjust p-values based on kelley's MS
    for i in np.arange(sorted_pval_idx.shape[1]):
        # get the indices of this particular mutation type
        idx_x, idx_y = sorted_pval_idx[0][i], sorted_pval_idx[1][i]

        # if this is the mutation type with the lowest p-value, 
        # we just use its original p-value
        if i == 0:
            adj_pvals[idx_x, idx_y] = pvals[idx_x, idx_y]
        # otherwise, we'll do an adjusted Chi2 test where we'll take
        # S_0 (number of mutations of this type in subset 0)
        # S_1 (number of mutations of this type in subset 1)
        # P_0 (number of all other mutations in subset 0, from i + 1 up to 96)
        # P_1 (number of all other mutations in subset 1, from i + 1 up to 96)
        elif i > 0 and i < 96:
            idx_x, idx_y = sorted_pval_idx[0][i], sorted_pval_idx[1][i]
            s_0 = sub_0_counts[idx_x, idx_y]
            s_1 = sub_1_counts[idx_x, idx_y]

            remaining_idxs = np.arange(i+1, 96)

            p_0, p_1 = 0, 0
            for idx in remaining_idxs:
                orig_x, orig_y = sorted_pval_idx[0][idx], sorted_pval_idx[1][idx]
                orig_count_x = sub_0_counts[orig_x, orig_y]
                orig_count_y = sub_1_counts[orig_x, orig_y]

                p_0 += orig_count_x
                p_1 += orig_count_y

            if 0 in (p_0, p_1): new_p = 1

            else:
                _,new_p,_,_ = ss.chi2_contingency([ [s_0, p_0],
                                                    [s_1, p_1] ])
                adj_pvals[idx_x, idx_y] = new_p
        else:
            idx_x, idx_y = sorted_pval_idx[0][i], sorted_pval_idx[1][i]
            s_0 = sub_0_counts[idx_x, idx_y]
            s_1 = sub_1_counts[idx_x, idx_y]
            orig_x, orig_y = sorted_pval_idx[0][95], sorted_pval_idx[1][95]
            p_0 = sub_0_counts[orig_x, orig_y]
            p_1 = sub_1_counts[orig_x, orig_y]

            _,new_p,_,_ = ss.chi2_contingency([ [s_0, p_0],
                                                [s_1, p_1] ])
            adj_pvals[idx_x, idx_y] = new_p

    return adj_pvals
