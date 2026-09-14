import os
import glob
import pickle
import numpy as np
import pandas as pd
from scipy.stats import norm

import rpy2.robjects as ro
from rpy2.robjects import default_converter
from rpy2.robjects.conversion import localconverter
from rpy2.robjects import numpy2ri

# Read beta summary statistics into DataFrame for a given level of hierarchy
def get_inla_summary(postfile, level=1):
    readRDS = ro.r["readRDS"]
    colnames_fn = ro.r["colnames"]
    rownames_fn = ro.r["rownames"]
    names_fn = ro.r["names"]
    is_null = ro.r["is.null"]
    
    res = readRDS(postfile)
    
    res_names = list(names_fn(res))
    if "beta_effective" not in res_names:
        raise KeyError(
            f"'beta_summary' not found in {postfile}. "
            f"Available entries: {res_names}"
        )

    beta_levels = list(names_fn(res.rx2("beta_effective").rx2("summaries")))
    if "beta_l%d" % level not in beta_levels:
        raise KeyError(
            f"'beta_l{level}' not found in {postfile}."
            f"Available entries: {beta_levels}"
        )
    beta_summary = res.rx2("beta_effective").rx2("summaries").rx2("beta_l%d" % level)
    
    # Extract coefficient names and posterior means
    coef_names = list(beta_summary.rx2("coef"))
    beta_means = np.array(beta_summary.rx2("mean"))
    beta_sd = np.array(beta_summary.rx2("sd"))
    df = pd.DataFrame({'mu':beta_means, 'sd':beta_sd}, index=coef_names)

    return df

# Return probability of means of two distributions being different by a threshold.
# - Two-sided: p_large = P(|mu1-mu2\>thresh)
def p_large(mu1, sigma1, mu2, sigma2, thresh):
    mu_delta = mu1 - mu2
    var_delta = sigma1**2 + sigma2**2
    sd_delta = np.sqrt(var_delta)

    p_gt_thresh = 1 - norm.cdf(thresh, mu_delta, sd_delta)
    p_lt_neg_thresh = norm.cdf(-thresh, mu_delta, sd_delta)
    p_large = p_gt_thresh + p_lt_neg_thresh

    return p_large

# Return p_large for expression of condition 1 vs. condition 2.
# Conditions specified as 3-tuples of INTS: (level, condition, celltype)
# ...or 2-tuples, assuming level-1 conditions.
def compute_de_indices(postfile, cond1, cond2, thresh):
    # If level information is not provided, assume Level-1
    if len(cond1)==2:
        cond1 = (1,) + tuple(cond1)
    if len(cond2)==2:
        cond2 = (1,) + tuple(cond2)

    for c in [cond1, cond2]:
        df_inla = get_inla_summary(postfile, level=c[0])
        
        if 'beta_l%d_c%d_m%d' % tuple(c) not in df_inla.index:
            raise KeyError(
                f"'beta_l{c[0]}_c{c[1]}_m{c[2]}' not found in results."
            )

    mu1, sigma1 = df_inla.loc['beta_l%d_c%d_m%d' % tuple(cond1)].values
    mu2, sigma2 = df_inla.loc['beta_l%d_c%d_m%d' % tuple(cond2)].values
    delta = mu2-mu1

    return delta, p_large(mu1, sigma1, mu2, sigma2, thresh)

# Return p_large for expression of condition 1 vs. condition 2.
# Conditions specified as 3-tuples of STRINGS: (level, condition, celltype)
# ...or 2-tuples, assuming level-1 conditions.
def compute_de_names(postfile, sinfo, cond1, cond2, thresh):
    # If level information is not provided, assume Level-1
    if len(cond1)==2:
        cond1 = (1,) + tuple(cond1)
    if len(cond2)==2:
        cond2 = (1,) + tuple(cond2)

    index_tuples = []
    for c in [cond1, cond2]:
        beta_mapping = sinfo['beta_mapping']['beta_level_%d' % c[0]]
        if c[1] not in beta_mapping:
            raise KeyError(
                f"Condition '{c[1]}' not found in beta_level_{c[0]}."
                f"Available entries: {beta_mapping}"
            )
        c_ind = list(beta_mapping).index(c[1]) + 1

        annot_mapping = sinfo['annotation_mapping']
        if c[2] not in annot_mapping:
            raise KeyError(
                f"Annotation '{c[1]}' not found."
                f"Available entries: {annot_mapping}"
            )
        m_ind = list(annot_mapping).index(c[2]) + 1

        index_tuples.append((c[0], c_ind, m_ind))

    return compute_de_indices(postfile, index_tuples[0], index_tuples[1], thresh)
        