import os
import pandas as pd

def load_lmem_pred(clinvar, res_dir, metric, obsdata_df = None, index2replace = ""):
    """
    """
    
    # load results of uni and multivariate analysis
    res_dfs_dict = {}
    for file in os.listdir(res_dir):
        if not file.startswith('.') and metric in file:
            res_name = ".".join(file.split(".")[2:4])
            res_dfs_dict[res_name] = pd.read_csv(os.path.join(res_dir, file), sep = "\t", index_col = 0)

    # combine uni and multi results: if multi was calculated, prioritaze this result
    coeffs_df = res_dfs_dict["multivariateME.coeffs"].combine_first(res_dfs_dict["univariateME.coeffs"])
    # coeffs_df.index = coeffs_df.index.str.replace(index2replace, '', regex = True) #TODO: bug index2replace
    interc_df = res_dfs_dict["multivariateME.intercepts"].combine_first(res_dfs_dict["univariateME.intercepts"])
    # interc_df.index = interc_df.index.str.replace(index2replace, '', regex = True) #TODO: bug index2replace
    pvals_df = res_dfs_dict["multivariateME.pvals"].combine_first(res_dfs_dict["univariateME.pvals"])
    # pvals_df.index = pvals_df.index.str.replace(index2replace, '', regex = True) #TODO: bug index2replace
    lowi_df = res_dfs_dict["multivariateME.lowi"].combine_first(res_dfs_dict["univariateME.lowi"])
    # lowi_df.index = lowi_df.index.str.replace(index2replace, '', regex = True) #TODO: bug index2replace
    highi_df = res_dfs_dict["multivariateME.highi"].combine_first(res_dfs_dict["univariateME.highi"])
    # highi_df.index = highi_df.index.str.replace(index2replace, '', regex = True) #TODO: bug index2replace
    res_df = interc_df[[clinvar]].merge(
        coeffs_df[[clinvar]], right_index = True, left_index = True, how = "outer").merge(
        pvals_df[[clinvar]], right_index = True, left_index = True, how = "outer").merge(
        lowi_df[[clinvar]], right_index = True, left_index = True, how = "outer").merge(
        highi_df[[clinvar]], right_index = True, left_index = True, how = "outer")
    res_df.columns = ["intercept", "coeff", "qval", "lowci", "highci"]
    res_df = res_df.reset_index(names = "gene")

    # calculate LMEM predictions if observed data is given
    if obsdata_df is not None:

        obsdata_df = obsdata_df.rename({"GENE": "gene"}, axis = 1)

        ## create age_decades for AGE variable
        if clinvar == "age_decades": #TODO: do this transformation in the nb/script to process the clinical data
            obsdata_df["age_decades"] = obsdata_df.apply(lambda row: row["AGE"] / 10, axis = 1)

        
        res_df = res_df.merge(obsdata_df[["gene", clinvar]], on = "gene", how = "inner").drop_duplicates()
        res_df["predicted"] = res_df["intercept"] + res_df["coeff"] * res_df[clinvar]
        res_df["predicted"] = res_df["intercept"] + res_df["coeff"] * res_df[clinvar]
        res_df["predicted"] = res_df["intercept"] + res_df["coeff"] * res_df[clinvar]

        if clinvar == "age_decades":
            res_df["AGE"] = res_df.apply(lambda row: row["age_decades"] * 10, axis = 1)  #TODO: do this transformation in the nb/script to process the clinical data

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
