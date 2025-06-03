# Module with utility functions to import so that we don't have have them repeated for all the analysis

import pandas as pd
import numpy as np
import requests
import time
import os


# Load data
# =========


# Load normal tissue mut
# ----------------------

def get_normal_maf(path_maf, gene_list, only_protein_pos=True, truncating=True):
    
    # TODO: check for correct filtering
    
    maf_df = pd.read_table(path_maf)
    maf_df["CLEAN_SAMPLE_ID"] = maf_df["SAMPLE_ID"].apply(lambda x: "_".join(x.split("_")[1:3]))
    maf_df_f = maf_df
    maf_df_f = maf_df.loc[(~maf_df["FILTER.not_covered"])].reset_index(drop = True)
    
    # Final filter
    maf_df_f = maf_df_f[
        (maf_df_f["TYPE"].isin(["SNV", "INSERTION", "DELETION"])
        ) & (maf_df_f["canonical_SYMBOL"].isin(gene_list))
        ].reset_index(drop = True)
    
    if only_protein_pos:
        maf_df_f = maf_df_f[maf_df_f["canonical_Protein_position"] != '-' ]
    
    maf_df_f.loc[(maf_df_f["TYPE"].isin(["INSERTION", "DELETION"])), "canonical_Consequence_broader"] = "indel"
    maf_df_f["canonical_Consequence_broader"] = maf_df_f["canonical_Consequence_broader"].replace("splice_region_variant", "splicing")
    maf_df_f["canonical_Consequence_broader"] = maf_df_f["canonical_Consequence_broader"].replace("essential_splice", "splicing")
    
    if truncating:
        maf_df_f["canonical_Consequence_broader"] = maf_df_f["canonical_Consequence_broader"].replace(
            {"nonsense": "truncating", "splicing": "truncating"}
            )
    
    maf_df_f['Alt_amino_acid'] = np.where(
        maf_df_f['canonical_Consequence_broader'] == 'missense', 
        maf_df_f['canonical_Amino_acids'].str.split('/').str[-1],
        np.nan 
    )
    
    # Parse
    cols = ["canonical_SYMBOL", 
            "canonical_Feature", 
            "canonical_Protein_position", 
            "canonical_Consequence_broader", 
            "Alt_amino_acid",
            "CHROM", 
            "POS", 
            "DEPTH", 
            "ALT_DEPTH"]
    maf_df_f = maf_df_f[cols].rename(columns={
        "canonical_SYMBOL" : "Gene", 
        "canonical_Feature" : "Ens_transcript_ID",
        "canonical_Protein_position" : "Pos",
        "canonical_Consequence_broader" : "Consequence",
        "CHROM" : "CHR",
        "POS" : "DNA_POS"}
        )
    
    # Extract first pos if multiple ones are defined
    maf_df_f['Pos'] = maf_df_f['Pos'].str.extract(r'^(\d+)')            
    maf_df_f['Pos'] = pd.to_numeric(maf_df_f['Pos'], errors='coerce')
    
    if only_protein_pos:
        maf_df_f = maf_df_f.dropna(subset="Pos").reset_index(drop=True)
        maf_df_f.Pos= maf_df_f.Pos.astype(int)
    
    return maf_df_f


# Load cancer mut
# ---------------

def get_cancer_maf_cohort(path_vep_out, lst_genes=None, only_protein_pos=True):
    
    columns_to_load = ["SYMBOL", "Location", "Protein_position", "Feature", "Consequence"]
    vep_cohort = pd.read_table(path_vep_out, usecols=columns_to_load, low_memory=False)
    if lst_genes is not None:
        vep_cohort = vep_cohort[vep_cohort["SYMBOL"].isin(lst_genes)]
    vep_cohort = vep_cohort[vep_cohort["SYMBOL"] != "-"].reset_index(drop=True)
    
    
    # Here we are removing the mutations that are in splice sites (as well as intronic and non coding exon regions)
    if only_protein_pos:
        vep_cohort = vep_cohort[vep_cohort["Protein_position"] != "-"].reset_index(drop=True)
        
    cols = ["SYMBOL", "Feature", "Location", "Protein_position", "Consequence"]
    vep_cohort = vep_cohort[cols].rename(columns={
        "SYMBOL" : "Gene", 
        "Feature" : "Ens_transcript_ID",
        "Protein_position" : "Pos"}
        )
    
    return vep_cohort


def get_cancer_maf_all(cohort_df, path_all_vep_out, lst_genes=None, only_protein_pos=True, truncating=True):
    
    lst_vep_out = []
    for cohort, ttype in cohort_df[["COHORT", "CANCER_TYPE"]].values:
        print(cohort)
        path_vep = f"{path_all_vep_out}/{cohort}.tsv.gz"
        if os.path.exists(path_vep):

            vep_cohort = get_cancer_maf_cohort(path_vep, lst_genes, only_protein_pos)
            vep_cohort["Cohort"] = cohort
            vep_cohort["Cancer_type"] = ttype
            lst_vep_out.append(vep_cohort)
        else:
            print(f"Path {path_vep} doesn't exist: Skipping..")
    
    df = pd.concat(lst_vep_out).reset_index(drop=True)
    
    # Extract first pos if multiple ones are defined
    df['Pos'] = df['Pos'].str.extract(r'^(\d+)')    
    df['Pos'] = pd.to_numeric(df['Pos'], errors='coerce')
    if only_protein_pos:
        df = df.dropna(subset="Pos").reset_index(drop=True)
        df.Pos= df.Pos.astype(int)
    
    # Ensure that all indels (including frameshift) are mapped to indels and not nonsense
    df["INDEL_INFRAME"] = False
    df.loc[(df["Consequence"] == 'inframe_insertion') | (df["Consequence"] == 'inframe_deletion'), "INDEL_INFRAME"] = True
    
    # TODO: CHECK ME
    df["Consequence"] = get_broad_consequence(df["Consequence"])
    # df["Consequence"] = df["Consequence"].replace("splice_region_variant", "splicing")
    # df["Consequence"] = df["Consequence"].replace("essential_splice", "splicing")
    
    if truncating:
        df["Consequence"] = df["Consequence"].replace(
            {"nonsense": "truncating", "splicing": "truncating"}
            )
    
    # Parse DNA location
    df["Chr"] = df.Location.apply(lambda x: f'chr{x.split(":")[0]}')
    df["Location"] = df.Location.apply(lambda x: x.split(":")[1])
    df["DNA_pos"] = df.Location.apply(lambda x: int(x.split("-")[0]))
    
    return df.drop(columns=["Location"])


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


# Load cancer clustering
# ----------------------

def get_o3d_gene_data(
    gene, 
    seq_df, 
    o3d_pos_df
    ):
    
    # Subset gene
    seq_df_gene = seq_df[seq_df["Gene"] == gene]
    gene_len = len(seq_df_gene.Seq.values[0])
    gene_pos = pd.DataFrame({"Pos" : range(1, gene_len+1)})
    uni_id, af_f = seq_df_gene[["Uniprot_ID", "F"]].values[0]
    score_gene_df = o3d_pos_df[o3d_pos_df["Gene"] == gene].reset_index(drop=True)
    
    ## O3D score vector
    score_gene_df = gene_pos.merge(score_gene_df[["Pos", "Score_obs_sim", "C", "C_ext"]], how="left", on="Pos")

    # Don't include Extended clusters
    score_gene_df["C"] = (score_gene_df["C"] == 1) & (score_gene_df["C_ext"] == 0)
    score_gene_df["C"] = score_gene_df["C"].astype(int)
    score_gene_df = score_gene_df.drop(columns=["C_ext"])

    score_gene_df.columns = "Pos", "O3D_score", "Cluster"
    score_gene_df["O3D_score"] = score_gene_df["O3D_score"].fillna(0)
    score_gene_df["Cluster"] = score_gene_df["Cluster"].fillna(0)
    
    return score_gene_df


def get_o3d_data_cancer(
    cohort_df, 
    o3d_cancer_path, 
    o3d_seq_df, 
    lst_genes
    ):
    
    lst_gene_df = []
    for gene in lst_genes:
        top_gene_qval = np.inf
        top_gene_score = 0
        top_cohort = np.nan
        top_gene_df = np.nan
        for cohort in cohort_df.COHORT.unique():

            o3d_pos_cancer_path_cohort = f"{o3d_cancer_path}/{cohort}/{cohort}.3d_clustering_pos.csv"   
            o3d_gene_cancer_path_cohort = f"{o3d_cancer_path}/{cohort}/{cohort}.3d_clustering_genes.csv"
            if os.path.exists(o3d_pos_cancer_path_cohort):
                columns_to_load = ["Gene", "Pos", "Score", "Score_obs_sim", "pval", "C", "C_ext"]
                o3d_pos_df = pd.read_csv(o3d_pos_cancer_path_cohort, usecols=columns_to_load, low_memory=False)
                o3d_pos_df = o3d_pos_df[o3d_pos_df["Gene"] == gene]
                o3d_gene_df = pd.read_csv(o3d_gene_cancer_path_cohort)
                o3d_gene_df = o3d_gene_df[o3d_gene_df["Gene"] == gene]
                if len(o3d_gene_df) > 0:
                
                    gene_qval = o3d_gene_df.qval.values[0]
                    gene_score = o3d_gene_df.Score_obs_sim_top_vol.values[0]
                    if gene_qval < top_gene_qval or (gene_qval == top_gene_qval and gene_score > top_gene_score):
                        top_gene_qval = gene_qval
                        top_gene_score = gene_score
                        top_cohort = cohort
                        top_gene_df = get_o3d_gene_data(gene, o3d_seq_df, o3d_pos_df)
                        top_gene_df["Gene"] = gene
                        top_gene_df["Cohort"] = cohort
        if not pd.isnull(top_cohort):               
            lst_gene_df.append(top_gene_df)
                        
    return pd.concat(lst_gene_df).reset_index(drop=True)


# Depth
# =====


# Get exons coord
# ---------------

def get_tr_lookup(transcript_id):
    
    server = "https://rest.ensembl.org"
    ext = f"/lookup/id/{transcript_id}?expand=1"

    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})

    while not r.ok:
        print("Retrying lookup..")
        time.sleep(5)
        r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})

    return r.json()


def get_cds_coord(transcript_id, len_cds_with_utr):

    server = "https://rest.ensembl.org"
    ext = f"/map/cds/{transcript_id}/1..{len_cds_with_utr}?"

    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
    
    while not r.ok:
        print("Retrying CDS map..")
        time.sleep(5)
        r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})

    return r.json()["mappings"]


def parse_cds_coord(exon):
    
    strand = exon["strand"]

    if strand == 1:
        start = exon["start"]
        end = exon["end"]
    else:
        start = exon["end"]
        end = exon["start"]
    
    if "id" in exon:
        exon_id = exon["id"]
        chrom = f'chr{exon["seq_region_name"]}'
        
        return exon_id, [chrom, start, end, strand]
    
    else:
        chrom = exon["seq_region_name"]
        
        return [chrom, start, end, strand]


# Get Exon coord to protein pos
# -----------------------------

def get_dna_exon_pos(exon_range, strand):

    if strand == -1:
         return np.arange(exon_range[1], exon_range[0] + 1)[::-1]
    else:
        return np.arange(exon_range[0], exon_range[1] + 1)
    

def get_exon_ix(i, exon_range, strand):
    
    len_exon = len(get_dna_exon_pos(exon_range, strand))
    
    return np.repeat(i, len_exon)


def get_dna_map_to_protein(coord_df):
    
    strand = coord_df.Strand.unique()[0]

    exons_range = coord_df[["Start", "End"]].values
    exons = np.concatenate([get_dna_exon_pos(exon, strand) for exon in exons_range])
    exons_ix = np.concatenate([get_exon_ix(i, exon, strand) for i, exon in enumerate(exons_range)])
    prot_pos = np.arange(len(exons)) // 3 + 1  
    
    df = pd.DataFrame({"GENE" : coord_df.Gene.unique()[0],
                       "CHR" : f'chr{coord_df.Chr.unique()[0]}', 
                       "DNA_POS" : exons, 
                       "PROT_POS" : prot_pos,
                       "REVERSE_STRAND" : strand,
                       "EXON_RANK" : exons_ix,
                       "TRANSCRIPT_ID" : coord_df.Ens_transcript_ID.unique()[0]})
    
    return df


def get_prot_coverage(dna_prot_df, gene, filter_masked_depth=True):
    
    gene_dna_prot_df = dna_prot_df[dna_prot_df["GENE"] == gene]
    gene_dna_prot_df = gene_dna_prot_df.dropna(subset=["PROT_POS"])[["PROT_POS", "COVERED", "DEPTH"]].reset_index(drop=True)
    gene_dna_prot_df = gene_dna_prot_df.groupby("PROT_POS").sum().reset_index()
    gene_dna_prot_df.COVERED = (gene_dna_prot_df.COVERED > 0).astype(int)
    
    return gene_dna_prot_df


def get_exon_coord_wrapper(maf):
    
    # Init df for coordinates
    coord_df = maf[["Gene", "Ens_transcript_ID"]].drop_duplicates().reset_index(drop=True)

    # Get coord
    coord_df_lst = []
    exons_coord_df_lst = []
    for gene, transcript in coord_df.values:
        print("Processing gene:", gene)
        coord_lst = []
        
        # Get the coord of exons with CDS and UTR as well as the lenght with UTR and exons ID
        exons_lookup = get_tr_lookup(transcript) # We will use this to get Exons ID
        for i, exon in enumerate(exons_lookup["Exon"]):
            exon_id, exons_coord = parse_cds_coord(exon)
            exons_coord_df_lst.append([f"{gene}--{i+1}_{transcript}_{exon_id}"] + exons_coord)
   
        # Get the coord of the exons without UTR to map to protein positions
        for i, exon in enumerate(get_cds_coord(transcript, exons_lookup["length"])):
            coord_lst.append((parse_cds_coord(exon) + [i]))

        gene_coord_df = pd.DataFrame(coord_lst, columns = ["Chr", "Start", "End", "Strand", "Exon_rank"])
        gene_coord_df["Gene"] = gene
        gene_coord_df["Ens_transcript_ID"] = transcript
        coord_df_lst.append(gene_coord_df)

    coord_df = pd.concat(coord_df_lst)
    exons_coord_df = pd.DataFrame(exons_coord_df_lst, columns = ["ID", "Chr", "Start", "End", "Strand"])
    
    return coord_df, exons_coord_df


# Get a DNA to protein mapping and coverage info & DNA to GENE annotation
# -----------------------------------------------------------------------

def dna2prot_depth(maf, coord_df, dna_sites, depth_df):
    
    # Map DNA to protein pos, get exons index to protein pos, etc
    dna_prot_df_lst = []
    for gene in maf.Gene.unique():
        gene_coord_df = coord_df[coord_df["Gene"] == gene]
        dna_prot_df_lst.append(get_dna_map_to_protein(gene_coord_df))
    dna_prot_df = pd.concat(dna_prot_df_lst)

    # Merge CDS position with availble sites (not masked) and depth info
    # and any other site that was included in the panel (splicing sites out of the CDS)
    dna_prot_df = dna_sites.merge(dna_prot_df, on=["GENE", "CHR", "DNA_POS"], how="outer")
    dna_prot_df["COVERED"] = dna_prot_df["CONTEXT"].notnull().astype(int)
    dna_prot_df = dna_prot_df.merge(depth_df.rename(columns={"CHROM" : "CHR", 
                                                             "POS" : "DNA_POS"}), 
                                    how="left", on=["CHR", "DNA_POS"])
    dna_prot_df.loc[dna_prot_df["COVERED"] == 0, "DEPTH"] = 0
    
    return dna_prot_df


def get_dna2prot_depth(maf, depth_df, consensus_df):
    
    consensus_df = consensus_df.merge(depth_df[["CHROM", "POS", "CONTEXT"]], on = ["CHROM", "POS"], how = 'left')
    consensus_df = consensus_df.rename(columns={"CHROM" : "CHR", "POS" : "DNA_POS"})
    
    if "DEPTH" not in depth_df:
        depth_df["DEPTH"] = depth_df.drop(columns=["CHROM", "POS", "CONTEXT"]).mean(1)
    depth_df = depth_df[["CHROM", "POS", "DEPTH"]].rename(columns = {"CHROM" : "CHR", "POS" : "DNA_POS"})
    
    coord_df, exons_coord_df = get_exon_coord_wrapper(maf)
    dna_prot_df = dna2prot_depth(maf, coord_df, consensus_df, depth_df)
    
    return dna_prot_df, exons_coord_df


# Utils function to retrieve exon ID from coordinate

def find_exon(x_coord, exon_coord_df):
    
    dna_pos, chrom, strand = x_coord["DNA_POS"], x_coord["CHR"], x_coord["REVERSE_STRAND"]
    
    if strand == -1:
        matches = exon_coord_df[(exon_coord_df['Chr'] == chrom) & (exon_coord_df['End'] <= dna_pos) & (dna_pos <= exon_coord_df['Start'])]
    
    else:
        matches = exon_coord_df[(exon_coord_df['Chr'] == chrom) & (exon_coord_df['End'] >= dna_pos) & (dna_pos >= exon_coord_df['Start'])]
    
    return matches['ID'].values[0] if not matches.empty else np.nan