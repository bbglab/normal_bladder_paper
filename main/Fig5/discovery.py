import warnings
warnings.filterwarnings("ignore")

import os
import tqdm
import sys
import click

import numpy as np
import pandas as pd


# Constants
sys.path.append('../../')
from consensus_variables import *


subsampling_rates = np.logspace(-2, np.log10(0.9), num=20)


def prob_min_uniform_sample_below_cut(N, n, cut):
    """
    Probability that the minimum of a sampling of n elements from [N] is lower or equal than cut
    """
    arr = np.array([np.log(N - cut - k) - np.log(N - k) for k in range(n)])
    return 1 - np.exp(np.sum(arr))


def load_mutations():

    somatic_mutations_file = f'{deepcsa_run_dir}/clean_somatic/all_samples.somatic.mutations.tsv'
    somatic_mutations = pd.read_csv(somatic_mutations_file, sep='\t', low_memory=False)
    mutations = somatic_mutations[
        ~(somatic_mutations['FILTER'].str.contains("not_in_exons"))
        & (somatic_mutations['canonical_Protein_affecting'] == 'protein_affecting')
        & (somatic_mutations['TYPE'] == 'SNV')
    ]
    mutations_lite = mutations[['CHROM', 'POS', 'REF', 'ALT', 'SAMPLE_ID', 'ALT_DEPTH', 'ALT_DEPTH_AM', 'SYMBOL', 'canonical_Consequence_broader']]
    mutations_lite = mutations_lite.groupby(by=['CHROM', 'POS', 'REF', 'ALT']).agg({'ALT_DEPTH_AM': 'sum', 'ALT_DEPTH': 'sum'}).reset_index()
    return mutations_lite


def get_aachange_format(r):

    if r['Protein_position'] == '-':
        return '-'
    elif len(r['Amino_acids']) == 1:
        return r['Amino_acids'] + r['Protein_position'] + r['Amino_acids']
    else:
        return r['Amino_acids'].split('/')[0] + r['Protein_position'] + r['Amino_acids'].split('/')[1]

        
def collect_vep():
    
    fn = os.path.join(deepcsa_run_dir, 'createpanels', 'panelannotation', 'captured_panel.tab.gz')
    vep_panel = pd.read_csv(fn, sep='\t', skiprows=44)
    dg = vep_panel[['#Uploaded_variation', 'Protein_position', 'Amino_acids', 'SYMBOL', 'CANONICAL']]
    dg = dg[dg['CANONICAL'] == 'YES']
    dg['CHROM'] = dg['#Uploaded_variation'].apply(lambda s: s.split('_')[0])
    dg['POS'] = dg['#Uploaded_variation'].apply(lambda s: s.split('_')[1])
    dg['POS'] = dg['POS'].astype(int)
    dg['REF'] = dg['#Uploaded_variation'].apply(lambda s: s.split('_')[2].split('/')[0])
    dg['ALT'] = dg['#Uploaded_variation'].apply(lambda s: s.split('_')[2].split('/')[1])
    dg['AACHANGE'] = dg.apply(get_aachange_format, axis=1)
    return dg


def load_panel(vep_annotations):

    df_panel = pd.read_csv(os.path.join(deepcsa_run_dir, 'createpanels', 'consensuspanels', 'consensus.exons_splice_sites.tsv'), sep='\t')
    df_depth = pd.read_csv(os.path.join(deepcsa_run_dir, 'annotatedepths', f'all_samples.depths.annotated.tsv.gz'), sep='\t')
    df_panel = pd.merge(df_panel, df_depth[['CHROM', 'POS', 'all_samples']], on=['CHROM', 'POS'], how='left')
    df_panel.rename(columns={'all_samples': 'DEPTH'}, inplace=True)
    df_panel = pd.merge(df_panel, vep_annotations[['CHROM', 'POS', 'REF', 'ALT', 'AACHANGE', 'SYMBOL']], 
                              left_on=['CHROM', 'POS', 'REF', 'ALT', 'GENE'], 
                              right_on=['CHROM', 'POS', 'REF', 'ALT', 'SYMBOL'],
                              how='left')
    df_panel = df_panel[df_panel['AACHANGE'] != '-']
    df_panel.dropna(inplace=True)
    df_panel['RESIDUE'] = df_panel['AACHANGE'].apply(lambda s: s[:-1])
    return df_panel


@click.command()
@click.option("--residue", is_flag=True, show_default=True, default=False, help="either genomic or residue based sites")
def cli(residue):
    
    mutations = load_mutations()
    vep = collect_vep()
    df_panel = load_panel(vep)
    mutations_lite = pd.merge(mutations, 
                          df_panel[['CHROM', 'POS', 'REF', 'ALT', 'DEPTH', 'GENE', 'AACHANGE']],
                          on=['CHROM', 'POS', 'REF', 'ALT'], 
                          how='left')
    mutations_lite.dropna(axis=0, inplace=True)  # drop nan sites
    mutations_lite = mutations_lite[mutations_lite['AACHANGE'] != '-']
    mutations_lite['RESIDUE'] = mutations_lite['AACHANGE'].apply(lambda s: s[:-1])

    if residue:
        mutations_lite = mutations_lite.groupby(['GENE', 'RESIDUE']).agg({'ALT_DEPTH': 'sum', 'DEPTH': 'mean'}).reset_index()    
    else:
        mutations_lite = mutations_lite.groupby(['GENE', 'POS']).agg({'ALT_DEPTH': 'sum', 'DEPTH': 'mean'}).reset_index()
    
    subsampling_rates = np.logspace(-2, np.log10(0.9), num=20)

    for i, p in tqdm.tqdm(enumerate(subsampling_rates)):
        mutations_lite[f'UNIQUE_RATE_{i}'] = mutations_lite.apply(lambda r: prob_min_uniform_sample_below_cut(r['DEPTH'], int(p * r['DEPTH']), r['ALT_DEPTH']), axis=1)
    
    if residue:
        mutations_lite.to_csv('./mutations_residue_rates.tsv', sep='\t', index=False)
    else:
        mutations_lite.to_csv('./mutations_genomic_rates.tsv', sep='\t', index=False)
        

if __name__ == '__main__':
    
    cli()
