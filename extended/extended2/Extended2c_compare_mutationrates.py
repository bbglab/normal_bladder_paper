import argparse
import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sys

sys.path.append('../../')
from consensus_variables import *


def produce_age(r):
    return int(r['microbiopsy_id'].split('_')[1].rstrip('M').rstrip('F'))

def produce_donor(r):
    return r['microbiopsy_id'].split('_')[0]

def process_data(data_file):
    data = pd.read_excel(data_file,sheet_name=0,skiprows=1)
    data = data[data['histological_feature']=='Urothelium']
    data = data[(data['ref'].isin(['A','C','G','T']))&
                (data['mut'].isin(['A','C','G','T']))]
    return data


def process_coverage(coverage_file):
    data = pd.read_excel(coverage_file,sheet_name=0,skiprows=1)
    data['donor'] = data.apply(produce_donor,axis=1)
    data = data[data['exome_median_coverage']!='-']
    data['exome_median_coverage'] = data['exome_median_coverage'].astype('int')
    data = data[data['exome_median_coverage']>=80]
    return list(data['donor'])

def read_mutrate_sample(mutrate_file):
    mutrate_sample = pd.read_csv(mutrate_file,sep='\t')
    mutrate_sample = mutrate_sample[(mutrate_sample['GENE']=='ALL_GENES')&
                                    (mutrate_sample['MUTTYPES']=='all_types')&
                                    (mutrate_sample['REGIONS']=='all')]
    mutrate_sample = mutrate_sample.rename(columns={'SAMPLE_ID':'sample'})
    return mutrate_sample

def process_mutrate(mutrate_dir,samples):
    mutrates = []
    for sample in samples:
        mutrate_file = os.path.join(mutrate_dir,sample + '.all.mutrates.tsv')
        mutrate_sample = read_mutrate_sample(mutrate_file)
        mutrates.append(mutrate_sample)
    return pd.concat(mutrates)


def process_metadata(metadata_file):
    metadata = pd.read_csv(metadata_file,sep='\t')
    metadata = metadata.rename(columns={'SAMPLE_ID':'sample',
                                        'AGE':'age'})
    return metadata

def mutation_count(data,genome_length,donors_coverage):
    counts = pd.DataFrame({'total_muts':data.groupby(['microbiopsy_id'])['pos'].count()}).reset_index()
    counts['age'] = counts.apply(produce_age,axis=1)
    counts['donor'] = counts.apply(produce_donor,axis=1)
    counts = counts[counts['donor'].isin(donors_coverage)]
    median = pd.DataFrame({'median_muts':counts.groupby(['donor','age'])['total_muts'].median()}).reset_index()
    median['mutrate'] = median['median_muts']/genome_length
    return median


def plot_regression(median, mutrates, figure_output_dir, m, b):
    figure_output_file = os.path.join(figure_output_dir, 'ExtendedFig2c_mutrate_comparison_lawson.pdf')

    fig, axs = plt.subplots(1, 1, figsize=(4,2.5))
    sns.scatterplot(data=median,
                    x='age',
                    y='mutrate',
                    color="red",
                    alpha = 0.5,
                    ax=axs,
                    s=20,
                    label='WES LCMs')
    sns.scatterplot(data=mutrates,
                    x='age',
                    y='MUTRATE_MB',
                    alpha=0.5,
                    s=20,
                    ax=axs,
                    label='Duplex sequencing')

    x_vals = np.array(axs.get_xlim())
    y_vals = b + m * x_vals
    axs.plot(x_vals, y_vals, '-', color='black', linewidth=1)

    axs.spines[['top', 'right']].set_visible(False)
    axs.set_xlabel('Age (years)')
    axs.set_ylabel('Mutations per megabase')
    axs.legend(frameon = False)
    plt.title('Normal bladder urothelium')

    plt.savefig(figure_output_file, bbox_inches = 'tight', dpi = 300)

def linear_regression(median):
    x = median['age'].values
    y = median['mutrate'].values
    m, b = np.polyfit(x, y, 1)
    return m,b


def main(deepCSA_directory):
    coverage_file = './data/lawson_et_al_data/aba8347_tables2.xlsx'

    genomes_file = './data/lawson_et_al_data/aba8347_tables5.xlsx'
    exomes_file = './data/lawson_et_al_data/aba8347_tables4.xlsx'

    deepCSA_dir = deepCSA_directory
    mutrate_dir = os.path.join(deepCSA_dir,'mutrate')

    figure_output_dir = './plots'

    #File with cohort metadata
    metadata_file = clinvars_regr_file

    ref_genome_length = 3200
    ref_exome_length = 35

    samples = all_sample_names_dirty

    #Genome data
    # data = process_data(genomes_file)
    # median = mutation_count(data,ref_genome_length)
    # plot_regression(median)

    #Exome data
    data = process_data(exomes_file)
    donors_coverage = process_coverage(coverage_file)
    median_lawson = mutation_count(data,ref_exome_length,donors_coverage)

    #mutrates
    mutrates = process_mutrate(mutrate_dir,samples)

    #metadata
    metadata = process_metadata(metadata_file)
    mutrates = mutrates.merge(metadata,on='sample')

    #plot alongside Lawson cohort
    m,b = linear_regression(median_lawson)
    plot_regression(median_lawson, mutrates, figure_output_dir, m, b)


if __name__ == '__main__':
    # parser = argparse.ArgumentParser()

    # parser.add_argument('-N',
    #                     '--normal_version',
    #                     required=True,
    #                     help='Name of version of deepCSA run')

    # args = parser.parse_args()

    # main(args.normal_version)
    main(deepcsa_run_dir)
