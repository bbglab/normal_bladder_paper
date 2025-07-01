import sys
import pandas as pd
import os

import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

sys.path.append('../../')

from consensus_variables import *


def convert_pyrimidines(ref, alt, pyrimidine_ref, pyrimidine_alt):
    converted_ref = ''
    converted_alt = ''

    if len(ref) == len(alt):
        for i in range(len(ref)):
            new_ref = pyrimidine_ref[ref[i]]
            if new_ref != ref[i]:
                converted_ref += new_ref
                converted_alt += pyrimidine_alt[alt[i]]
            else:
                converted_ref += ref[i]
                converted_alt += alt[i]

    else:
        converted_ref = ref
        converted_alt = alt

    return converted_ref, converted_alt


def reformat_site_selection_data(data_point, pyrimidine_ref, pyrimidine_alt):
    converted_ref, converted_alt = convert_pyrimidines(data_point['ref'], data_point['alt'], pyrimidine_ref, pyrimidine_alt)
    return [data_point['chr'],data_point['pos'],converted_ref,converted_alt,data_point['site_selection'],data_point['p_value']]


def categorize_normal_mutations_counts(r, threshold=1):
    if r['count'] > threshold:
        return '>' + str(threshold)
    elif r['count'] == threshold:
        return str(threshold)
    else:
        return str(0)

def significance_groups(r, threshold=1e-5):
    if r['site_selection'] > 0:
        if r['p_value'] > threshold:
            return 'Not significant'
        return 'Significant'
    return 'Not observed'


def process_normal_TERT_mutations(normal_path, TERT_info_38):
    normal_muts = pd.read_csv(os.path.join(normal_path + '/clean_somatic',
                                            'all_samples.somatic.mutations.tsv'),
                                sep='\t')
    normal_muts = normal_muts[['CHROM','POS','REF','ALT','ALT_DEPTH','VAF','SAMPLE_ID', "TYPE"]]
    normal_TERT_muts = normal_muts[(normal_muts['CHROM']==('chr' + TERT_info_38['chr'])) &
                                    (normal_muts['POS'] >= int(TERT_info_38['start'])) &
                                    (normal_muts['POS'] <= int(TERT_info_38['end'])) &
                                    (normal_muts['TYPE'] == 'SNV')
                                    ]
    normal_TERT_muts = normal_TERT_muts.drop('TYPE', axis = 1)
    return normal_TERT_muts.rename(columns={'CHROM':'chr',
                                            'POS':'pos',
                                            'REF':'ref',
                                            'ALT':'alt',
                                            'SAMPLE_ID':'sample'})


def process_site_selection(normal_path):
    site_selection = pd.read_csv(os.path.join(normal_path + '/sitecomparisongloballoc',
                                            'all_samples.global_loc.site.comparison.tsv.gz'),
                                    sep='\t')
    site_selection = site_selection[site_selection['GENE']=='TERTpromoter']
    site_selection = site_selection[['CHROM','POS','REF','ALT','OBS/EXP','p_value']].rename(columns={'CHROM':'chr',
                                                                                            'POS':'pos',
                                                                                            'REF':'ref',
                                                                                            'ALT':'alt',
                                                                                            'OBS/EXP':'site_selection'})
    return site_selection

def convert_relative_position(r, TERT_gene_TSS):
    return float(r['pos']) - TERT_gene_TSS

def fetch_saturation_data(saturation_file):
    data = pd.read_csv(saturation_file, sep='\t')
    return data[["chr", "pos", "ref", "alt", "value", "p", "system"]]


def draw_needleplot(axis,muts_data,ymax,factor,color_def,inverted=False):
    if inverted:
        axis.vlines(muts_data['rel_pos'], ymin=ymax-muts_data['count'], ymax=ymax+ymax/factor, lw=1, zorder=1, alpha=0.5,color='black')
        axis.scatter(muts_data['rel_pos'], ymax-muts_data['count'], color='white', zorder=3, lw=1, ec="none", s=20)
        axis.scatter(muts_data['rel_pos'].values, ymax-muts_data['count'].values, zorder=4,alpha=0.7, lw=0.1, ec="none", s=20, color=color_def['needle_point_count'])
    else:
        axis.vlines(muts_data['rel_pos'], ymin=0-ymax/factor, ymax=muts_data['count'], lw=1, zorder=1, alpha=0.5,color='black')
        axis.scatter(muts_data['rel_pos'], muts_data['count'], color='white', zorder=3, lw=1, ec="none", s=20)
        axis.scatter(muts_data['rel_pos'].values, muts_data['count'].values, zorder=4,alpha=0.7, lw=0.1, ec="none", s=20, color=color_def['needle_point_count'])


def draw_site_selection(axis,data,value,color_def):
    print(data)
    axis.plot(data['rel_pos'],
                data[value],
                zorder=1,
                color=color_def[value],
                lw=0.5)
    axis.fill_between(data['rel_pos'],
                        0,
                        data[value],
                        color=color_def[value],
                        alpha=0.4,
                        zorder=0,
                        lw=0.5)

    signif_data = data.copy()

    significant_rows = signif_data[signif_data['p_value'] < 1e-5]
    significant_indices = significant_rows.index
    left_indices = significant_indices - 1
    right_indices = significant_indices + 1
    all_indices = set(significant_indices) | set(left_indices) | set(right_indices)
    all_indices = [idx for idx in all_indices if 0 <= idx < len(signif_data)]  # Ensure valid indices

    # Keep only the selected rows
    signif_data = signif_data.loc[all_indices].sort_index()
    axis.fill_between(signif_data['rel_pos'],
                      0,
                      signif_data[value],
                      color='blue',#color_def[value],
                      alpha=1,
                      zorder=2,
                      lw=1)



def draw_continuous_line(axis,data,value,color_def):
    axis.plot(data['rel_pos'],
              data[value],
              zorder=2,
              color=color_def[value],
              lw=0.5)
    axis.fill_between(data['rel_pos'],
                      0,
                      data[value],
                      color=color_def[value],
                      alpha=0.4,
                      zorder=0,
                      lw=0.5)


def draw_stripplot(axis, data, color_def, value):
    data = data.sort_values('category_count',ascending=True)
    sns.stripplot(data,
                    x='category_count',
                    y=value,
                    ax=axis,
                    jitter=0.1,
                    hue='category_count',
                    palette=color_def['palette_' + value],
                    alpha = 0.7,
                    size=4)


# Reusable plotting function adapted for grouped data
def plot_half_violin_group(ax, data, x_col, y_col, group_order, palette, title, ylabel=None):
    sns.violinplot(
        data=data,
        ax=ax,
        x=x_col,
        y=y_col,
        hue=x_col,
        hue_order=group_order,
        order=group_order,
        palette=palette,
        inner=None,
        linewidth=0,
        bw=0.2,
        cut=0,
        density_norm='width',
        zorder = 1
    )

    ymin, ymax = ax.get_ylim()

    # Cover left half
    for i, group in enumerate(group_order):
        ax.add_patch(plt.Rectangle((i - 0.41, ymin), 0.39, ymax - ymin, color='white', zorder=2))

    x_offset = -0.25  # fixed offset if you want all dots shifted
    jitter_width = 0.1  # match seaborn's jitter
    dot_size = 10

    for i, group in enumerate(group_order):
        group_data = data[data[x_col] == group]
        print(group_data.shape)
        y_vals = group_data[y_col]
        # x_vals = [i + x_offset] * len(y_vals)# + np.random.uniform(-jitter_width, jitter_width, size=len(y_vals))
        x_vals = i + x_offset + np.random.uniform(-jitter_width, jitter_width, size=len(y_vals))
        ax.scatter(x_vals,
                    y_vals,
                color=palette[group],
                alpha=0.6,
                s=dot_size,
                edgecolor='none',
                linewidth=0,
                label=group, zorder = 3)

    ax.set_xlabel('')
    ax.set_title(title)
    ax.spines[['top', 'right']].set_visible(False)
    ax.legend([], [], frameon=False)


def draw_stripplot_half_violin(axis, data, color_def, value):

    group_order = ['Not observed', 'Not significant', 'Significant']
    custom_palette = dict(zip(group_order, color_def['palette_' + value]))
    plot_half_violin_group(axis, data, 'category_count', value, group_order, custom_palette, '', ylabel='Experimental\nfunctional impact')

    # data = data.sort_values('category_count',ascending=True)
    # sns.stripplot(data,
    #                 x='category_count',
    #                 y=value,
    #                 ax=axis,
    #                 jitter=0.1,
    #                 hue='category_count',
    #                 palette=color_def['palette_' + value],
    #                 alpha = 0.7,
    #                 size=4)


def draw_scatter(axis_left,site_selection_data,saturation_data,color_def):
    selection_saturation = pd.merge(site_selection_data, saturation_data, on=['pos','rel_pos','ref','alt'], how = 'inner')
    sns.scatterplot(selection_saturation,
                    y='site_selection',
                    x='value',
                    ax=axis_left,
                    color=color_def['value'],
                    s=10)

def tidy_axis_needleplot(axis, ylabel, ylimit, factor, xmin, xmax, inverted=False):
    if inverted:
        axis.spines[['right']].set_visible(False)
        axis.set_yticks(range(ylimit,-100,-200))
        axis.set_yticklabels(range(ylimit-100,-100,-200)[::-1])
        axis.set_ylim([0,ylimit+ylimit/factor])
        # axis.vlines(0, 0, ylimit + ylimit / factor, colors='gray', linestyles='dashed', linewidth=1)
    else:
        axis.spines[['top','right']].set_visible(False)
        axis.set_ylim([0-ylimit/factor,ylimit])
        # axis.vlines(0, 0 - ylimit / factor, ylimit, colors='gray', linestyles='dashed', linewidth=1)

    axis.set_xlim([xmin, xmax])
    axis.set_ylabel(ylabel, rotation=0, horizontalalignment='right',verticalalignment='center')
    # axis.tick_params(axis='both', labelsize=6)



def tidy_axis_continuous(axis, ylabel, ymin, ymax):
    axis.spines[['top','right']].set_visible(False)
    axis.set_ylabel(ylabel, rotation=0, horizontalalignment='right',verticalalignment='center')
    # axis.tick_params(axis='both', labelsize=6)
    axis.set_ylim([ymin,ymax])
    # axis.vlines(0,ymin,ymax,colors='gray', linestyles='dashed', linewidth=1)
    # axis.set_yscale('log')

def tidy_axis_stripplot(axis, xlabel, ylabel):
    axis.spines[['top','right']].set_visible(False)
    axis.set_ylabel(ylabel)
    axis.set_xlabel(xlabel)
    # axis.tick_params(axis='x', labelsize=6)
    # axis.tick_params(axis='y', labelsize=6)

def tidy_axis_scatterplot(axis, xlabel, ylabel):
    axis.spines[['top','right']].set_visible(False)
    axis.set_ylabel(ylabel)
    axis.set_xlabel(xlabel)
    # axis.tick_params(axis='both', labelsize=6)

def plot_TERT_panel(figure_output_dir,
                    unique_normal_TERT_muts,
                    site_selection_data_format,
                    saturation_data,
                    unique_tumor_TERT_muts,
                    color_def,
                    metadata):

    site_selection_data_format = site_selection_data_format.merge(saturation_data,on=['chr','pos','ref','alt','rel_pos'],how='outer')
    site_selection_data_format = site_selection_data_format[['chr','pos','ref','alt','site_selection','rel_pos', 'p_value']]
    site_selection_data_format['site_selection'] = site_selection_data_format['site_selection'].fillna(0)

    figure_output_file = os.path.join(figure_output_dir,'Fig5f_TERT_needleplots.pdf')
    fig = plt.figure(figsize=(7, 3.5))
    gs = GridSpec(20, 11, figure=fig, wspace=0, hspace=0)
    ax_tumor = fig.add_subplot(gs[13:19, 0:8])
    ax_normal = fig.add_subplot(gs[0:6, 0:8], sharex=ax_tumor)
    ax_selection = fig.add_subplot(gs[7:9, 0:8], sharex=ax_tumor)
    ax_saturation = fig.add_subplot(gs[10:12, 0:8], sharex=ax_tumor)
    ax_stripplot_saturation = fig.add_subplot(gs[0:8, 9:])
    ax_scatter = fig.add_subplot(gs[11:19, 9:])
    plt.setp(ax_normal.get_xticklabels(), visible=False)
    plt.setp(ax_saturation.get_xticklabels(), visible=False)
    plt.setp(ax_selection.get_xticklabels(), visible=False)

    min_TERT_site = unique_normal_TERT_muts['rel_pos'].min() - 5
    max_TERT_site = unique_normal_TERT_muts['rel_pos'].max() + 5

    #Draw panel with needleplot of mutations in normal
    draw_needleplot(ax_normal, unique_normal_TERT_muts, 25, 60, color_def, inverted=False)
    ylabel_normal = 'Mutations\nin normals\n(N=' + str(metadata['total_normal']) + ')'
    tidy_axis_needleplot(ax_normal, ylabel_normal, 25, 60, min_TERT_site, max_TERT_site)

    #Draw panel of site selection
    site_selection_position = pd.DataFrame({'site_selection':site_selection_data_format.groupby(by='rel_pos')['site_selection'].max()})
    p_value_position = pd.DataFrame({'p_value':site_selection_data_format.groupby(by='rel_pos')['p_value'].min()})
    site_selection_position = pd.concat((site_selection_position, p_value_position), axis = 'columns').reset_index()
    draw_site_selection(ax_selection, site_selection_position, 'site_selection', color_def)
    ylabel_selection = 'Site\nselection'
    tidy_axis_continuous(ax_selection, ylabel_selection, 0, 900)

    #Draw panel of saturation
    extended_saturation_data = pd.merge(saturation_data,unique_normal_TERT_muts,on='rel_pos',how='outer')
    extended_saturation_data[['count','value']] = extended_saturation_data[['count','value']].fillna(0)
    draw_continuous_line(ax_saturation, extended_saturation_data, 'value', color_def)
    ylabel_saturation = 'Experimental\nfunctional\nimpact'
    tidy_axis_continuous(ax_saturation, ylabel_saturation, extended_saturation_data['value'].min()-0.5,extended_saturation_data['value'].max()+0.5)

    #Draw panel with needleplot of mutations in tumors
    draw_needleplot(ax_tumor, unique_tumor_TERT_muts, 700, 25, color_def, inverted=True)
    ylabel_tumor = f'Mutations\nin tumors\n(N={metadata["total_tumors"]:,})'
    tidy_axis_needleplot(ax_tumor, ylabel_tumor, 700, 25, min_TERT_site, max_TERT_site, inverted=True)
    ax_tumor.set_xlabel('TERT_promoter sites')

    site_saturation = pd.merge(saturation_data,site_selection_data_format,on=['chr','pos','ref', 'alt','rel_pos'], how='left')

    site_saturation['p_value'] = site_saturation['p_value'].fillna(1)
    site_saturation['site_selection'] = site_saturation['site_selection'].fillna(0)
    site_saturation['category_count'] = site_saturation.apply(significance_groups,axis=1)

    draw_stripplot_half_violin(ax_stripplot_saturation, site_saturation, color_def,'value')

    xlabel_strip = ''
    ylabel_strip_saturation = 'Experimental\nfunctional impact'

    tidy_axis_stripplot(ax_stripplot_saturation, xlabel_strip, ylabel_strip_saturation)

    ax_stripplot_saturation.set_xticklabels([f"Not\nobserved\n(N={site_saturation[site_saturation['category_count'] == 'Not observed'].shape[0]:,})",
                                                f"Not\nsignificant\n(N={site_saturation[site_saturation['category_count'] == 'Not significant'].shape[0]:,})",
                                                f"Significant\n(N={site_saturation[site_saturation['category_count'] == 'Significant'].shape[0]:,})"])

    #Draw panel with scatterplot of site selection vs saturation
    draw_scatter(ax_scatter, site_selection_data_format, saturation_data, color_def)
    xlabel_scatter = 'Experimental\nfunctional impact'
    ylabel_scatter = 'Site selection'
    tidy_axis_scatterplot(ax_scatter, xlabel_scatter, ylabel_scatter)

    #Save plot
    print(figure_output_file)
    plt.savefig(figure_output_file,dpi=300)


def main(normal_path, figure_output_dir):


    #Processed data files
    saturation_file = f"{input_data_dir}/TERT_merged_data.tsv"

    #Sample sizes
    metadata = {'total_hartwig': 5582,
                'total_pcawg': 2554,
                'total_tumors': 8136,
                'total_normal': 79}

    #TERT genomic coordinates
    TERT_info_38 = {'chr':'5',
                    'start':1294942,
                    'end':1295289}

    TERT_gene_TSS = 1295068

    #Pyrimidine conversion data
    pyrimidine_ref = {'A':'T',
                      'C':'C',
                      'G':'C',
                      'T':'T'}

    pyrimidine_alt = {'A':'T',
                      'C':'G',
                      'G':'C',
                      'T':'A'}

    #Fetch normal data
    normal_TERT_muts = process_normal_TERT_mutations(normal_path, TERT_info_38)

    #Fetch tumor data
    unique_tumor_TERT_muts = pd.read_table(f"{input_data_dir}/TERT_merged_data.tsv")[["pos", "count_WGS_tumor"]]
    unique_tumor_TERT_muts.columns = ["pos", "count"]
    unique_tumor_TERT_muts = unique_tumor_TERT_muts[unique_tumor_TERT_muts["count"]>0].reset_index(drop = True)

    #Fetch saturation data
    saturation_data = fetch_saturation_data(saturation_file)
    saturation_data['rel_pos'] = saturation_data.apply(convert_relative_position,args=(TERT_gene_TSS,),axis=1)


    #Fetch site_selection data
    site_selection_data = process_site_selection(normal_path)
    site_selection_data_list = [reformat_site_selection_data(r, pyrimidine_ref, pyrimidine_alt) for k,r in site_selection_data.iterrows()]
    site_selection_data_format = pd.DataFrame(site_selection_data_list,columns=['chr','pos','ref','alt','site_selection', 'p_value'])
    site_selection_data_format['rel_pos'] = site_selection_data_format.apply(convert_relative_position,args=(TERT_gene_TSS,),axis=1)

    #Process normal data
    #Take mutation with max VAF at every position
    filtered_normal_TERT_muts = normal_TERT_muts.groupby(by=['chr','pos','sample'],
                                                            as_index = False).max()

    #Count unique mutations and sum reads per position
    unique_normal_TERT_muts = pd.DataFrame({'count':filtered_normal_TERT_muts.groupby(['pos'])['VAF'].count(),
                                            'sum_reads':filtered_normal_TERT_muts.groupby(['pos'])['ALT_DEPTH'].sum()}).reset_index()
    unique_normal_TERT_muts['rel_pos'] = unique_normal_TERT_muts.apply(convert_relative_position,args=(TERT_gene_TSS,),axis=1)

    #Processed tumor data
    unique_tumor_TERT_muts['rel_pos'] = unique_tumor_TERT_muts.apply(convert_relative_position,args=(TERT_gene_TSS,),axis=1)

    #####################Draw Figure#########################
    color_def = {'needle_point_count':'#a65628',
                 'needle_point_reads':'#a65628',
                 'site_selection':'#b3b3ff',
                 'value':'#df0086',
                 'palette_site_selection':['#d9d9ff',
                                           '#c2c2ff',
                                           '#b3b3ff'],
                 'palette_value':['gray',
                                 '#6baed6',
                                 'blue']
                 }

    plot_TERT_panel(figure_output_dir,
                    unique_normal_TERT_muts,
                    site_selection_data_format,
                    saturation_data,
                    unique_tumor_TERT_muts,
                    color_def,
                    metadata)


if __name__ == '__main__':

    input_data_dir = '../../data/tert_data'

    #Directory to write figure file
    figure_output_dir = './plots'

    main(deepcsa_run_dir, figure_output_dir)
