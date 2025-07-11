import os
import pandas as pd

def load_lmem_pred(clinvar, res_dir, metric, obsdata_df = None, sign_threshold = 0.2):
    """
    """
    
    # load results of uni and multivariate analysis
    res_dfs_dict = {}
    idx_cols_saved = False
    for file in os.listdir(res_dir):
        if not file.startswith('.') and metric in file:
            res_name = ".".join(file.split(".")[2:4])
            res_dfs_dict[res_name] = pd.read_csv(os.path.join(res_dir, file), sep = "\t", index_col = 0)
            if not idx_cols_saved:
                idxs = res_dfs_dict[res_name].index.to_list()
                cols = res_dfs_dict[res_name].columns.to_list()
                idx_cols_saved = True
            else:
                res_dfs_dict[res_name] = res_dfs_dict[res_name].reindex(idxs)
                res_dfs_dict[res_name] = res_dfs_dict[res_name].T.reindex(cols).T

    # combine uni and multi results: if multi was calculated and uni was significant, use multi calculation
    sign_mask = res_dfs_dict["univariateME.pvals"] < sign_threshold
    res_dfs_dict["multivariateME.coeffs"] = res_dfs_dict["multivariateME.coeffs"][sign_mask]
    res_dfs_dict["multivariateME.intercepts"] = res_dfs_dict["multivariateME.intercepts"][sign_mask]
    res_dfs_dict["multivariateME.pvals"] = res_dfs_dict["multivariateME.pvals"][sign_mask]
    res_dfs_dict["multivariateME.lowi"] = res_dfs_dict["multivariateME.lowi"][sign_mask]
    res_dfs_dict["multivariateME.highi"] = res_dfs_dict["multivariateME.highi"][sign_mask]

    coeffs_df = res_dfs_dict["multivariateME.coeffs"].combine_first(res_dfs_dict["univariateME.coeffs"])
    interc_df = res_dfs_dict["multivariateME.intercepts"].combine_first(res_dfs_dict["univariateME.intercepts"])
    pvals_df = res_dfs_dict["multivariateME.pvals"].combine_first(res_dfs_dict["univariateME.pvals"])
    lowi_df = res_dfs_dict["multivariateME.lowi"].combine_first(res_dfs_dict["univariateME.lowi"])
    highi_df = res_dfs_dict["multivariateME.highi"].combine_first(res_dfs_dict["univariateME.highi"])
    res_df = interc_df[[clinvar]].rename({clinvar: "intercept"}, axis = 1).merge(
        coeffs_df[[clinvar]].rename({clinvar: "coeff"}, axis = 1),
        right_index = True, left_index = True, how = "outer").merge(
        pvals_df[[clinvar]].rename({clinvar: "qval"}, axis = 1),
        right_index = True, left_index = True, how = "outer").merge(
        lowi_df[[clinvar]].rename({clinvar: "lowci"}, axis = 1), 
        right_index = True, left_index = True, how = "outer").merge(
        highi_df[[clinvar]].rename({clinvar: "highci"}, axis = 1), 
        right_index = True, left_index = True, how = "outer")
    res_df = res_df.reset_index(names = "gene")

    # calculate LMEM predictions if observed data is given
    if obsdata_df is not None:

        obsdata_df = obsdata_df.rename({"GENE": "gene"}, axis = 1)

        ## create age_decades for AGE variable
        if clinvar == "age_decades": 
            obsdata_df["age_decades"] = obsdata_df.apply(lambda row: row["AGE"] / 10, axis = 1)

        
        res_df = res_df.merge(obsdata_df[["gene", clinvar]], on = "gene", how = "inner").drop_duplicates()
        res_df["predicted"] = res_df["intercept"] + res_df["coeff"] * res_df[clinvar]
        res_df["predicted"] = res_df["intercept"] + res_df["coeff"] * res_df[clinvar]
        res_df["predicted"] = res_df["intercept"] + res_df["coeff"] * res_df[clinvar]

        if clinvar == "age_decades":
            res_df["AGE"] = res_df.apply(lambda row: row["age_decades"] * 10, axis = 1) 

    return res_df

    
def load_omega_res(omega_dir):
    """
    Loads omega results from the deepCSA output directory.
    Avoids loading multi results
    """

    res_dfs = []

    for file in os.listdir(omega_dir):
        if file.startswith("output_mle") and "multi" not in file:
            res_df = pd.read_csv(os.path.join(omega_dir, file), sep = "\t")
            res_dfs.append(res_df)

    omega_df = pd.concat(res_dfs, ignore_index = True)

    return omega_df
