import sys
import pandas as pd
import numpy as np
import os

sys.path.append('../..') 
from consensus_variables import * 

## -- Auxiliary -- ##
def get_broad_consequence(list_of_annotations):
    """
    Group variants into broader consequence types.
    """
        
    CONSEQUENCES_LIST = [
        'transcript_ablation',
        'splice_acceptor_variant',
        'splice_donor_variant',
        'stop_gained',
        'frameshift_variant',
        'stop_lost',
        'start_lost',
        'transcript_amplification',
        'inframe_insertion',
        'inframe_deletion',
        'missense_variant',
        'protein_altering_variant',
        'splice_region_variant',
        'splice_donor_5th_base_variant',
        'splice_donor_region_variant',
        'splice_polypyrimidine_tract_variant',
        'incomplete_terminal_codon_variant',
        'start_retained_variant',
        'stop_retained_variant',
        'synonymous_variant',
        'coding_sequence_variant',
        'mature_miRNA_variant',
        '5_prime_UTR_variant',
        '3_prime_UTR_variant',
        'non_coding_transcript_exon_variant',
        'intron_variant',
        'NMD_transcript_variant',
        'non_coding_transcript_variant',
        'upstream_gene_variant',
        'downstream_gene_variant',
        'TFBS_ablation',
        'TFBS_amplification',
        'TF_binding_site_variant',
        'regulatory_region_ablation',
        'regulatory_region_amplification',
        'feature_elongation',
        'regulatory_region_variant',
        'feature_truncation',
        'intergenic_variant'
    ]
    
    GROUPING_DICT = {
        'transcript_ablation': 'nonsense',
        'splice_acceptor_variant': 'nonsense',
        'splice_donor_variant': 'nonsense',
        'stop_gained': 'nonsense',
        'frameshift_variant': 'nonsense',
        'stop_lost': 'nonsense',
        'start_lost': 'nonsense',
        'missense_variant': 'missense',
        'inframe_insertion': 'indel',
        'inframe_deletion': 'indel',
        'splice_donor_variant': 'splicing',
        'splice_acceptor_variant': 'splicing',
        'splice_region_variant': 'splicing',
        'splice_donor_5th_base_variant': 'splicing',
        'splice_donor_region_variant': 'splicing',
        'splice_polypyrimidine_tract_variant': 'splicing',
        'synonymous_variant': 'synonymous',
        'incomplete_terminal_codon_variant': 'synonymous',
        'start_retained_variant': 'synonymous',
        'stop_retained_variant': 'synonymous',
        'protein_altering_variant' : 'protein_altering_variant',
        'transcript_amplification' : 'transcript_amplification', 
        'coding_sequence_variant': 'coding_sequence_variant', 
        'mature_miRNA_variant': 'non_coding_exon_region',
        '5_prime_UTR_variant': 'non_coding_exon_region',
        '3_prime_UTR_variant': 'non_coding_exon_region',
        'non_coding_transcript_exon_variant': 'non_coding_exon_region',
        'NMD_transcript_variant': 'non_coding_exon_region',
        'intron_variant': 'intron_variant',
        'non_coding_transcript_variant' : 'non_coding_transcript_variant',
        'upstream_gene_variant': 'non_genic_variant',
        'downstream_gene_variant': 'non_genic_variant',
        'TFBS_ablation': 'non_genic_variant',
        'TFBS_amplification': 'non_genic_variant',
        'TF_binding_site_variant': 'non_genic_variant',
        'regulatory_region_ablation': 'non_genic_variant',
        'regulatory_region_amplification': 'non_genic_variant',
        'feature_elongation': 'non_genic_variant',
        'regulatory_region_variant': 'non_genic_variant',
        'feature_truncation': 'non_genic_variant',
        'intergenic_variant': 'non_genic_variant',
        '-'  : '-'
    }
    
    consequence_rank_dict = { consequence : rank for rank, consequence in enumerate(CONSEQUENCES_LIST) }
    rank_consequence_dict = { rank : consequence for rank, consequence in enumerate(CONSEQUENCES_LIST) }
    
    list_of_single_annotations = []
    list_of_broad_annotations = []
    for x in list_of_annotations:
        all_consequences = x.split(",")
        all_consequences_ranks = map(lambda x: consequence_rank_dict[x], all_consequences)
        single_consequence = rank_consequence_dict[min(all_consequences_ranks)]
        list_of_single_annotations.append(single_consequence)
        if single_consequence in GROUPING_DICT:
            list_of_broad_annotations.append(GROUPING_DICT[single_consequence])
        else:
            list_of_broad_annotations.append(single_consequence)

    return list_of_broad_annotations

## -- Main -- ##

def load_cohort_data(clinvars_file, muts_file, mutrate_file):
    """
    """
    
    # load maf file and count number of muts per sample
    muts_df = pd.read_csv(muts_file, sep = "\t")
    print(f"Total number of mutations {len(muts_df)}")
    print(f"Types of variants included: {muts_df['TYPE'].unique()}")
    mnvs = muts_df.loc[muts_df["TYPE"] == "MNV"]["MUT_ID"].unique()
    mnv_maxlen = max([len(mut.split("_")[-1].split(">")[0]) for mut in mnvs])
    print(f"\tMNV max length: {mnv_maxlen}")
    
    muts_df["CLEAN_SAMPLE_ID"] = muts_df["SAMPLE_ID"].map(old2new_sample_names)
    nmuts_df = muts_df.groupby("CLEAN_SAMPLE_ID").size().to_frame("N_MUTS")
    print(f"Number of mutations:")
    print(f"\tRange: {nmuts_df['N_MUTS'].min()} - {nmuts_df['N_MUTS'].max()}")
    print(f"\tMedian: {nmuts_df['N_MUTS'].median()}")

    # load mutrate calculations
    ## deepCSA current config: all_types mutation rate includes whatever type is remaining in the somatic MAF
    mutrate_df = pd.read_csv(mutrate_file, sep = "\t")
    mutrate_df_PA = mutrate_df.loc[(mutrate_df["SAMPLE_ID"].str.contains("P19_")) &
                                (mutrate_df["GENE"] == "ALL_GENES") &
                                (mutrate_df["MUTTYPES"] == "all_types") & 
                                (mutrate_df["REGIONS"] == "protein_affecting")]
    mutrate_df = mutrate_df.loc[(mutrate_df["SAMPLE_ID"].str.contains("P19_")) &
                                (mutrate_df["GENE"] == "ALL_GENES") &
                                (mutrate_df["MUTTYPES"] == "all_types") & 
                                (mutrate_df["REGIONS"] == "all")]
    mutrate_df["CLEAN_SAMPLE_ID"] = mutrate_df["SAMPLE_ID"].map(old2new_sample_names)
    mutrate_df = mutrate_df[["CLEAN_SAMPLE_ID", "MUTRATE_MB"]].set_index("CLEAN_SAMPLE_ID")
    print("Mutation rate:")
    print(f'\tOverall median: {mutrate_df["MUTRATE_MB"].median().round(2)}')
    print(f'\tOverall mean: {mutrate_df["MUTRATE_MB"].mean().round(2)}')
    print(f'\tOverall range: {mutrate_df["MUTRATE_MB"].min().round(2)} - {mutrate_df["MUTRATE_MB"].max().round(2)}')
    print(f'\tOverall 95% interval: {mutrate_df["MUTRATE_MB"].quantile(0.05).round(2)} - {mutrate_df["MUTRATE_MB"].quantile(0.95).round(2)}')
    
    print(f'\tPA median: {mutrate_df_PA["MUTRATE_MB"].median().round(2)}')
    print(f'\tPA mean: {mutrate_df_PA["MUTRATE_MB"].mean().round(2)}')
    print(f'\tPA range: {mutrate_df_PA["MUTRATE_MB"].min().round(2)} - {mutrate_df_PA["MUTRATE_MB"].max().round(2)}')
    print(f'\tPA 95% interval: {mutrate_df_PA["MUTRATE_MB"].quantile(0.05).round(2)} - {mutrate_df_PA["MUTRATE_MB"].quantile(0.95).round(2)}')
    

    # load clinvars
    clinvars_df = pd.read_csv(clinvars_file, sep = "\t")
    clinvars_df["CLEAN_SAMPLE_ID"] = clinvars_df["SAMPLE_ID"].map(old2new_sample_names)
    clinvars_df = clinvars_df.set_index("CLEAN_SAMPLE_ID")    

    # merge and order by age
    cohort_df = clinvars_df.merge(nmuts_df, right_index = True, left_index = True, how = "outer").merge(
        mutrate_df, right_index = True, left_index = True, how = "outer").reset_index()
    cohort_df = cohort_df.sort_values(by = ["AGE", "CLEAN_SAMPLE_ID"], ascending = True)

    return cohort_df

def get_cancer_maf_cohort(path_vep_out, lst_genes = None):
    """
    Gets SNVs and indels from a VEP annotated MAF file for
    a specific set of genes. If only_protein_pos = True,
    filters out mutations that are not protein coding
    """
    
    # load VEP annotated MAF
    columns_to_load = ["SYMBOL", "Location", "Protein_position", "Feature", "Consequence"]
    vep_cohort = pd.read_table(path_vep_out, usecols = columns_to_load, low_memory = False)

    # keep selected genes (if applicable)
    if lst_genes is not None:
        vep_cohort = vep_cohort[vep_cohort["SYMBOL"].isin(lst_genes)]
    vep_cohort = vep_cohort[vep_cohort["SYMBOL"] != "-"].reset_index(drop=True)
        
    cols = ["SYMBOL", "Feature", "Location", "Protein_position", "Consequence"]
    vep_cohort = vep_cohort[cols].rename(columns={"SYMBOL" : "Gene", 
                                                "Feature" : "Ens_transcript_ID",
                                                "Protein_position" : "Pos"})
    
    return vep_cohort

def get_cancer_maf_all(cohort_df, path_all_vep_out, lst_genes = None):
    """
    Gets SNVs and indels from VEP files
    of a set of IntoGen cohorts. Keeps only
    mutations in specific genes. If only_protein_pos = True,
    filters out mutations that are not protein coding
    """
    
    # iterate through intogen cohorts and load VEP annotated MAF
    lst_vep_out = []
    for cohort, ttype in cohort_df[["COHORT", "CANCER_TYPE"]].values:
        path_vep = f"{path_all_vep_out}/{cohort}.tsv.gz"
        if os.path.exists(path_vep):

            vep_cohort = get_cancer_maf_cohort(path_vep, lst_genes)
            vep_cohort["Cohort"] = cohort
            vep_cohort["Cancer_type"] = ttype
            lst_vep_out.append(vep_cohort)
        else:
            print(f"Path {path_vep} doesn't exist: Skipping..")
    
    df = pd.concat(lst_vep_out).reset_index(drop=True)
    
    # rename frameshift variants to indel
    df["Consequence"] = df["Consequence"].replace("frameshift_variant", "inframe_insertion") 
    df["Consequence"] = get_broad_consequence(df["Consequence"])
    
    # parse DNA location
    df["Chr"] = df.Location.apply(lambda x: f'chr{x.split(":")[0]}')
    df["Location"] = df.Location.apply(lambda x: x.split(":")[1])
    df["DNA_pos"] = df.Location.apply(lambda x: int(x.split("-")[0]))
    
    return df.drop(columns=["Location"])

def get_normal_maf(path_maf, lst_genes, only_protein_pos = True):
    """
    Gets SNVs and indels from a MAF file for
    a specific set of genes. If only_protein_pos = True,
    filters out mutations that are not protein coding 
    """
    
    # read MAF
    maf_df = pd.read_table(path_maf)
    print(f"Total number of somatic mutations in the consensus panel: {len(maf_df)}")

    # filter genes of interest
    maf_df_f = maf_df.loc[(maf_df["canonical_SYMBOL"].isin(lst_genes))].reset_index(drop = True)
    print(f"\tMutations in selected genes: {len(maf_df_f)}")

    # keep SNVs and indels only
    maf_df_f = maf_df_f[(maf_df_f["TYPE"].isin(["SNV", "INSERTION", "DELETION"]))
                    ].reset_index(drop = True)
    print(f"\tSNVs and indels: {len(maf_df_f)}")
    
    # rename some consequence types
    maf_df_f.loc[(maf_df_f["TYPE"].isin(["INSERTION", "DELETION"])), "canonical_Consequence_broader"] = "indel"
    maf_df_f["canonical_Consequence_broader"] = maf_df_f["canonical_Consequence_broader"].replace("splice_region_variant", "splicing")
    maf_df_f["canonical_Consequence_broader"] = maf_df_f["canonical_Consequence_broader"].replace("essential_splice", "splicing")
    
    cols = ["canonical_SYMBOL", "canonical_Feature", "canonical_Protein_position",
            "canonical_Consequence_broader", "CHROM", "POS", "DEPTH", "ALT_DEPTH"]
    maf_df_f = maf_df_f[cols].rename(columns={"canonical_SYMBOL" : "Gene", 
                                            "canonical_Feature" : "Ens_transcript_ID",
                                            "canonical_Protein_position" : "Pos",
                                            "canonical_Consequence_broader" : "Consequence",
                                            "CHROM" : "CHR",
                                            "POS" : "DNA_POS"})
    
    return maf_df_f