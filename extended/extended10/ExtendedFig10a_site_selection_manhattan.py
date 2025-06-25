import pandas as pd
import os, sys
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

sys.path.append('../../')
from consensus_variables import deepcsa_run_dir


def numeralize_chromosomes(r):
    if r['CHROM'] == 'chrX':
        return 23
    else:
        return int(r['CHROM'].lstrip('chr'))

def process_site_selection(normal_path, regions):
    site_selection = pd.read_csv(os.path.join(normal_path + '/sitecomparisongloballoc',
                                            'all_samples.global_loc.site.comparison.tsv.gz'),
                                    sep='\t')
    site_selection = site_selection[site_selection['GENE'].isin(regions)]
    site_selection['chr'] = site_selection.apply(numeralize_chromosomes,axis=1)
    site_selection = site_selection[['chr','POS','REF','ALT','GENE','OBS/EXP','p_value']].rename(columns={'POS':'pos',
                                                                                                            'REF':'ref',
                                                                                                            'ALT':'alt',
                                                                                                            'OBS/EXP':'site_selection'})
    return site_selection


def add_color(r,color_def):
    return color_def[r['GENE']]

def compute_logp(r):
    if r['p_value'] == 0:
        return 20
    else:
        return -np.log10(r['p_value'])

def plot_legend(color_def,axis):
    patches = []

    for gene in color_def.keys():
        patches.append(mpatches.Patch(color=color_def[gene], label=gene))
    axis.legend(handles=patches,bbox_to_anchor=(1.05, 1),
                                loc='upper left', borderaxespad=0.)

def axis_manhattan(s0,s1,s2,y_value,axis):

    axis.scatter(s0['rel_pos'].values, s0[y_value].values, alpha=0.3, s=3, facecolor = 'none', color=s0['color'])
    axis.scatter(s1['rel_pos'].values, s1[y_value].values, alpha=0.5, s=6, facecolor = 'none', color=s1['color'])
    axis.scatter(s2['rel_pos'].values, s2[y_value].values, alpha=1, s=10, color=s2['color'])

    axis.spines[['top', 'right', 'bottom']].set_visible(False)

    if y_value == 'logp':
        axis.set_ylabel('-log(p-value)', fontsize=8)
        axis.set_xlabel('')
    else:
        axis.set_ylabel('Site selection', fontsize=9)
        axis.set_xlabel('', fontsize=8)


def plot_manhattan(figure_output_dir,
                    site_selection_data,
                    color_def):

    figure_output_file = os.path.join(figure_output_dir, 'ExtendedFig10a_site_selection_manhattan.png')
    fig, ax = plt.subplots(figsize=(6,3))  # Adjust figure size as needed

    s0 = site_selection_data[site_selection_data['p_value']>=1e-5]
    s1 = site_selection_data[(site_selection_data['p_value']<1e-5)&(site_selection_data['p_value']>=1e-15)]
    s2 = site_selection_data[site_selection_data['p_value']<1e-15]

    # axis_manhattan(s0, s1, s2, 'p_value', ax_manhattan1)
    axis_manhattan(s0, s1, s2, 'site_selection', ax)
    # axis_volcano(s0,s1,s2,ax_volcano,color_def)

    x_tick_positions = site_selection_data.groupby(by = ["GENE"])['rel_pos'].median()
    # print(x_tick_positions.index.values)
    # print(x_tick_positions.values)
    ax.set_ylim(-20)
    ax.set_xticks(x_tick_positions.values)
    ax.set_xticklabels(x_tick_positions.index.values, rotation = 90, fontsize = 8.5)

    plt.tight_layout()
    plt.savefig(figure_output_file, dpi=300)

def main(normal_path):

    # Directory to write figure file
    figure_output_dir = './plots'

    #Regions to include in plot
    regions = ['ARID1A','CDKN1A','CREBBP','EP300','FOXQ1','KDM6A','KMT2C',
                'NOTCH2','KMT2D','PIK3CA','RB1','RBM10','STAG2','TERTpromoter','TP53']

    color_def = {'ARID1A':'#7f46c9',
                'NOTCH2':'#c1e2f3',
                'PIK3CA': '#de5a8f',
                #  'FGFR3': '#c1e5d1',
                'TERTp':'#7f46c9',
                'FOXQ1': '#6e6e6e',
                'CDKN1A':'#b4aab4',
                'KMT2C':'#2ab967',
                'KMT2D': '#fe6a90',
                'RB1': '#c9d2f1',
                'CREBBP':'#efd477',
                'TP53': '#e99d91',
                'EP300':'#f4ad73',
                #  'KDM6A':'#f2fab4',
                'KDM6A':'#c1e5d1',
                'RBM10': '#a2e296',
                'STAG2':'#ed54d4'}

    site_selection = process_site_selection(normal_path,regions)

    custom_gene_order = ['CDKN1A', 'ARID1A', 'PIK3CA', 'CREBBP','EP300','FOXQ1','KDM6A','KMT2C','TERTp',
                        'NOTCH2','KMT2D','RB1','RBM10','STAG2','TP53']

    site_selection = process_site_selection(normal_path, regions)
    site_selection["GENE"] = site_selection["GENE"].replace("TERTpromoter", "TERTp")

    # Convert 'GENE' to categorical with a custom order
    site_selection["GENE"] = pd.Categorical(site_selection["GENE"], categories=custom_gene_order, ordered=True)

    # Now sort based on the custom order, then by chr, pos, alt
    site_selection = site_selection.sort_values(['GENE', 'chr', 'pos', 'alt'])


    gap_size = 1700  # Adjust based on how much separation you want

    # Compute relative positions with separation between genes
    rel_pos = 0
    prev_gene = None
    rel_pos_list = []

    for _, row in site_selection.iterrows():
        if prev_gene is not None and row['GENE'] != prev_gene:
            rel_pos += gap_size  # Add gap when the gene changes
        rel_pos_list.append(rel_pos)
        rel_pos += 1
        prev_gene = row['GENE']

    site_selection.insert(0, 'rel_pos', rel_pos_list)

    site_selection['color'] = site_selection.apply(add_color,args=(color_def,),axis=1)
    site_selection['logp'] = site_selection.apply(compute_logp,axis=1)

    plot_manhattan(figure_output_dir, site_selection, color_def)



if __name__ == '__main__':

    main(deepcsa_run_dir)
