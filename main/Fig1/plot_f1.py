import pandas as pd
import os
import re
import numpy as np
import sys
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
import matplotlib.patches as mpatches
from matplotlib.font_manager import FontProperties
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.patches as patches


# sys.path.append('../..') 
sys.path.append('/data/bbg/projects/bladder_ts/notebooks/manuscript_figures_vMarch2025') 
from consensus_variables import * 

## -- Auxiliary -- ##
def modif_index(idx):
    match = re.match(r"[A|T|C|G]\[[A|T|C|G]>[A|T|C|G]\][A|T|C|G]", idx) #A[C>A]A
    if match:
        return f"{idx[0]}{idx[2]}{idx[-1]}>{idx[4]}" #ACA>A
    return idx 

def minor_tick_labels():
    major_labels = ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']
    flanks = ['AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT',
            'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT']
    minor_labels = []
    for subs in major_labels:
        for flank in flanks:
            minor_labels.append(flank[0] + subs[0] + flank[1])
    return minor_labels

contexts_ordered = ['ACA>A', 'ACC>A', 'ACG>A', 'ACT>A', 'CCA>A', 'CCC>A', 'CCG>A', 'CCT>A', 'GCA>A', 'GCC>A', 'GCG>A', 'GCT>A', 'TCA>A', 'TCC>A', 'TCG>A', 'TCT>A',
                        'ACA>G', 'ACC>G', 'ACG>G', 'ACT>G', 'CCA>G', 'CCC>G', 'CCG>G', 'CCT>G', 'GCA>G', 'GCC>G', 'GCG>G', 'GCT>G', 'TCA>G', 'TCC>G', 'TCG>G', 'TCT>G',
                        'ACA>T', 'ACC>T', 'ACG>T', 'ACT>T', 'CCA>T', 'CCC>T', 'CCG>T', 'CCT>T', 'GCA>T', 'GCC>T', 'GCG>T', 'GCT>T', 'TCA>T', 'TCC>T', 'TCG>T', 'TCT>T',
                        'ATA>A', 'ATC>A', 'ATG>A', 'ATT>A', 'CTA>A', 'CTC>A', 'CTG>A', 'CTT>A', 'GTA>A', 'GTC>A', 'GTG>A', 'GTT>A', 'TTA>A', 'TTC>A', 'TTG>A', 'TTT>A',
                        'ATA>C', 'ATC>C', 'ATG>C', 'ATT>C', 'CTA>C', 'CTC>C', 'CTG>C', 'CTT>C', 'GTA>C', 'GTC>C', 'GTG>C', 'GTT>C', 'TTA>C', 'TTC>C', 'TTG>C', 'TTT>C',
                        'ATA>G', 'ATC>G', 'ATG>G', 'ATT>G', 'CTA>G', 'CTC>G', 'CTG>G', 'CTT>G', 'GTA>G', 'GTC>G', 'GTG>G', 'GTT>G', 'TTA>G', 'TTC>G', 'TTG>G', 'TTT>G']

## -- Main -- ##

def cohort_descrp_plot(data_df, save_file, plot_config = {}):
    """
    """

    # plot config (check default for missing values)
    plot_config_dflt = {
        "nrows": 10,
        "sample_id": "CLEAN_SAMPLE_ID",
        "sample_highlight_color": "#EDEDED",
        "rows_ordered": ["N_MUTS", "MUTRATE_MB", 
                        "white_space", "is_dome",
                        "AGE", "BMI", "is_male", 
                        "smoking_status", "history_drinking",
                        "had_prior_chemotherapy"],
        "legend_order": ["is_dome","is_male", 
                        "smoking_status", "history_drinking",
                        "had_prior_chemotherapy",
                        "AGE", "BMI"],
        "barplot_vars": ["N_MUTS", "MUTRATE_MB"],
        "clinvar2color": clinvar2color,
        "clinvar_id2name": clinvar_id2name,
        "metric_id2name": metric_id2name_ylabel,
        "barplot_color": "#FFC5B5",
        "figsize": (18, 4),
        "height_ratios": [4, 4, 0.3, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6], #TODO: this can be better
        "legend_fontsize": 8,
        "legend_fontweight": "bold",
        "xylabels_fontsize": 10,
        "ylabels_pad": {"is_dome": 80, "AGE": 17,
                        "BMI": 17, "is_male": 17,
                        "smoking_status": 41,
                        "history_drinking": 40,
                        "had_prior_chemotherapy": 73},
        "legend_loc": {"is_dome": [1.14, 1],
                        "AGE": [1.017, -3, 1.3, 9],
                        "BMI": [1.017, -5.8, 1.3, 9], 
                        "is_male": [1.135, 1.7],
                        "smoking_status": [1.145, 19],
                        "history_drinking": [1.11, 7],
                        "had_prior_chemotherapy": [1.11, 5.5]},
        "legend_title_loc": {"AGE": -32,
                            "BMI": -33}
    }
    plot_config_dflt["samples_ordered"] = list(data_df[plot_config_dflt["sample_id"]].values)
    plot_config_dflt["dummy2clinvar"] = {v: k for k, v in clinvar2dummy.items()}

    for k in plot_config_dflt:
        if k not in plot_config:
            plot_config[k] = plot_config_dflt[k]

    # plot
    fig, axs = plt.subplots(plot_config["nrows"], 1,
                        figsize = plot_config["figsize"],
                        sharex = False,
                        gridspec_kw = {"height_ratios": plot_config["height_ratios"]},
                        constrained_layout = False)
    legends_lst = []

    for i, var in enumerate(plot_config["rows_ordered"]):

        # barplots
        if var in plot_config["barplot_vars"]:
            sns.barplot(data = data_df,
            x = plot_config["sample_id"], y = var,
            ax = axs[i], color = plot_config["barplot_color"],
            dodge = True, order = plot_config["samples_ordered"])
            
            axs[i].set_xticklabels([])
            axs[i].tick_params(axis = 'x', which = 'both', length = 0)
            axs[i].set_xlabel("")
            axs[i].set_ylabel(plot_config["metric_id2name"][var],
                            fontsize = plot_config["xylabels_fontsize"])
            axs[i].spines['top'].set_visible(False)
            axs[i].spines['right'].set_visible(False)
        
        elif var == "white_space":
            axs[i].spines['top'].set_visible(False)
            axs[i].spines['right'].set_visible(False)
            axs[i].spines['bottom'].set_visible(False)
            axs[i].spines['left'].set_visible(False)
            axs[i].set_xticklabels([])
            axs[i].set_yticklabels([])
            axs[i].tick_params(axis = 'both', which = 'both', length = 0)
            axs[i].set_xlabel("")
            axs[i].set_ylabel("")

        # heatmaps
        else:
            ## cmap depending on categorical or continuous var
            if var in plot_config["dummy2clinvar"]:
                var_original = plot_config["dummy2clinvar"][var]
                cmap = ListedColormap(plot_config["clinvar2color"][var_original].values())

            else:
                var_original = var
                cmap = LinearSegmentedColormap.from_list("", list(plot_config["clinvar2color"][var_original].values()))

            ## draw heatmap
            xticklabels = False
            if i == len(plot_config["rows_ordered"])-1:
                xticklabels = True
            sns.heatmap(data = data_df[[plot_config["sample_id"], var]].set_index(plot_config["sample_id"]).T,
            ax = axs[i], cbar = False, xticklabels = xticklabels, cmap = cmap)
            axs[i].set_xlabel("")
            axs[i].set_ylabel(plot_config["clinvar_id2name"][var_original], 
                            fontsize = plot_config["xylabels_fontsize"], rotation = 0, 
                            labelpad = plot_config["ylabels_pad"][var], loc = 'bottom')
            axs[i].set_yticklabels([])


    # sample names
    axs[i].set_xticklabels(axs[i].get_xticklabels(), fontfamily = "monospace", fontsize = plot_config["xylabels_fontsize"])
    ## highlight samples from the same donor
    color = plot_config["sample_highlight_color"]
    subject_id = ""
    for label in axs[i].get_xticklabels():
        if label.get_text().split("_")[0] != subject_id:
            if color == plot_config["sample_highlight_color"]:
                color = "white"
                label.set_bbox(dict(facecolor = color, edgecolor = 'none', boxstyle = 'square,pad=0.3'))
            elif color == "white":
                color = plot_config["sample_highlight_color"]
                label.set_bbox(dict(facecolor = color, edgecolor = 'none', boxstyle = 'square,pad=0.3'))
            subject_id = label.get_text().split("_")[0]
        else:
            label.set_bbox(dict(facecolor = color, edgecolor = 'none', boxstyle = 'square,pad=0.3'))
            subject_id = label.get_text().split("_")[0]

    # legend depending on categorical or continuous var (different order so it fits)
    for i, var in enumerate(plot_config["legend_order"]):
        if var in plot_config["dummy2clinvar"]:
            var_original = plot_config["dummy2clinvar"][var]
            cmap = ListedColormap(plot_config["clinvar2color"][var_original].values())

            patches = []
            for categ, color in plot_config["clinvar2color"][var_original].items():
                patches.append(mpatches.Patch(color = color,
                                            label = plot_config["clinvar_id2name"][categ]))
            legend = axs[i].legend(handles = patches, bbox_to_anchor = plot_config["legend_loc"][var], 
                    title = plot_config["clinvar_id2name"][var_original], ncols = 2,
                    frameon = False, prop = {'size': plot_config["legend_fontsize"]},
                    title_fontproperties = FontProperties(weight = plot_config["legend_fontweight"],
                                                        size = plot_config["legend_fontsize"]))
            legend._legend_box.align = "left"
            if i > 0:
                axs[i].add_artist(legend)
            legends_lst.append(legend)

        else:
            var_original = var
            cmap = LinearSegmentedColormap.from_list("", list(plot_config["clinvar2color"][var_original].values()))

            sm = plt.cm.ScalarMappable(cmap = cmap)
            sm.set_clim(vmin = data_df[var].min(), vmax = data_df[var].max())
            inset_ax = inset_axes(axs[i], width = "5%", height = "10%", loc = 'center left', 
            bbox_to_anchor = plot_config["legend_loc"][var], bbox_transform = axs[i].transAxes,
            borderpad = 0)
            cbar = plt.colorbar(sm, cax = inset_ax, orientation = 'horizontal')
            cbar.set_label(plot_config["clinvar_id2name"][var], 
                        labelpad = plot_config["legend_title_loc"][var],
                        fontproperties = FontProperties(weight = plot_config["legend_fontweight"],
                                                        size = plot_config["legend_fontsize"]),
                        loc = 'left')
            cbar.ax.tick_params(labelsize = plot_config["legend_fontsize"])
            cbar.outline.set_visible(False)


    # save 
    plt.savefig(save_file, dpi = 300, bbox_inches = 'tight', bbox_extra_artists = legends_lst) #TODO: resolution to general config

    return fig

def plot_mut_count_comparison(df_count,
                            figsize=(5, 4),
                            color="#E0E0E0",
                            save=False,
                            filename=None):

    pivoted = df_count.pivot(index="Gene", columns="Type", values="SNVs").fillna(0)
    pivoted = pivoted.sort_values(by="Cancer", ascending=True)

    fig, (ax_left, ax_right) = plt.subplots(
        1, 2, figsize=figsize, sharey=True, gridspec_kw={"width_ratios": [1, 1]}
    )

    # Barplots
    bars_left = ax_left.barh(pivoted.index, pivoted["Normal"], color="#E0E0E0")
    bars_right = ax_right.barh(pivoted.index, pivoted["Cancer"], color="#E0E0E0")
    ax_left.invert_xaxis() 
    ax_left.tick_params(axis="y", labelleft=False, labelright=True, pad=35)

    # Add values to the bars
    for bar in bars_left:
        ax_left.text(
            bar.get_width() + 80,
            bar.get_y() + bar.get_height() / 2, 
            f'{int(bar.get_width())}',
            va='center', 
            ha='right',  
            fontsize=8,
            color="black"
        )
    for bar in bars_right:
        ax_right.text(
            bar.get_width() + 80,  
            bar.get_y() + bar.get_height() / 2,  
            f'{int(bar.get_width())}',  
            va='center',  
            ha='left',  
            fontsize=8,
            color="black"
        )

    # Details
    for tick in ax_left.get_yticklabels():
        tick.set_horizontalalignment('center')
        
    ax_left.tick_params(left=False, right=True)
    ax_left.set_xlabel("Number of mutations")
    ax_right.set_xlabel("Number of mutations")
    ax_right.set_xlim([0, pivoted["Normal"].max()])

    ax_left.spines['left'].set_visible(False)
    ax_left.spines['top'].set_visible(False)
    ax_right.spines['top'].set_visible(False)
    ax_right.spines['right'].set_visible(False)

    # ax_left.set_title(f"$\\mathbf{{Normal\\ bladder}}$\n{pivoted['Normal'].sum():,}", fontsize=10)
    # ax_right.set_title(f"$\\mathbf{{Bladder\\ tumors}}$\n{pivoted['Cancer'].sum():,}", fontsize=10)
    ax_left.set_title(f"$\\mathbf{{Normal\\ bladder}}$\nDonors: 45", fontsize=10)
    ax_right.set_title(f"$\\mathbf{{Bladder\\ cancer}}$\nTumors: 892", fontsize=10)
    plt.tight_layout()
    
    if save and filename is not None:
        fig.savefig(filename, dpi=300, bbox_inches='tight')
    plt.show()

def plot_signature(profile_df, ax, ttype = "Frequency", add_contexts = False, text = None, textloc_x = 1,
                fontfactor = 0, ylabel = ""):
    """
    Args:
        profile_df: 96-array dataframe with contexts as index and 
        values to plot as column
        ax: matplotlib axis object to plot on
        ttype: type of values to be plotted (default: mutation frequencies)
        add_contexts: boolean, whether to add contexts as tick labels
        text: additional text annotation
        fontfactor: factor to scale all the fonts in the plot

    Returns:
        Modifies the provided ax with the signature bar plot.
    """

    # correct contexts' format and order
    profile_df.index = [modif_index(idx) for idx in profile_df.index]
    profile_df = profile_df.reindex(contexts_ordered)
    profile = profile_df.iloc[:, 0]

    # if ttype is "frequency", ensure profile sums to one
    if ttype == "Frequency":
        total = np.sum(profile)
        if abs(total - 1) > 0.01:
            profile = profile / total
            print("profile transformed to mutation probabilities (sum to one)")
    elif ttype == "Percentage":
        total = np.sum(profile)
        if abs(total - 1) > 0.01:
            profile = profile / total
            print("profile transformed to mutation probabilities (sum to one)")
        profile = profile * 100
        ttype = "% SBS"    
    
    if ylabel:
        ttype = ylabel

    # bar plot: one color per change
    barlist = ax.bar(range(96), profile, width= 0.5)
    mut2color = {'C>A': '#1ebff0',
                    'C>G': '#050708',
                    'C>T': '#e62725',
                    'T>A': '#cbcacb',
                    'T>C': '#a1cf64',
                    'T>G': '#edc8c5'}

    for i, mut in enumerate(mut2color.keys()):
        ## 16 trinucleotides per substitution ttype
        for j in range(16):  
            barlist[i * 16 + j].set_color(mut2color[mut])

    ax.set_xlim([-0.5, 96])
    ymax = np.max(profile) * 1.2
    ax.set_ylim(0, ymax)  # Ensure space for labels above bars
    ax.tick_params(axis = 'y', labelsize = 5+fontfactor)
    ax.set_ylabel(ttype, fontsize = 6+fontfactor)

    # add labels substitution ttype on top of the corresponding bars
    for i, mut in enumerate(mut2color.keys()):

        x_pos = i * 16 + 8  # text: center above the group
        y_pos = ymax * 1.05  # text: slightly above the highest bar
        rect = patches.Rectangle((x_pos - 7.5, y_pos - y_pos/20), 15.5, y_pos/15,
                                color = mut2color[mut], clip_on = False)
        ax.add_patch(rect)
        ax.text(x_pos, y_pos+y_pos/12, mut, color = 'black',
                ha = 'center', va = 'center', fontsize = 7+fontfactor, 
                fontweight = 'bold')

    # add contexts if needed
    if add_contexts:
        minor_ticks = np.arange(0.2, 96.2, 1)
        ax.set_xticks(minor_ticks)
        ax.set_xticklabels(minor_tick_labels(), rotation = 90, fontsize = 5+fontfactor)
    else:
        ax.set_xticklabels([])
        ax.tick_params(axis = 'x', length = 0)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # add extra info if needed
    if text:
        ax.text(textloc_x, ymax - ymax / 10, text, color = 'black',
                ha = 'left', va = 'top', fontsize = 6+fontfactor)

    return ax