import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
import sys
import numpy as np
import pandas as pd
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test, multivariate_logrank_test
from matplotlib.gridspec import GridSpec
import matplotlib.patches as patches


sys.path.append('../..') #TODO: pass this file to repo dir
from consensus_variables import * #TODO: pass this file to repo dir

## -- Auxiliary functions -- ##
def add_age_triangle(ax, width, add_text = False, text_fontsize = None):
    # Triangle growing from left to right, centered above the heatmap
    triangle = patches.Polygon(
        xy = [[0, 1.02], [width, 1.02], [width, 1.08]],  # flat base from 0 to width, height at right edge
        closed = True,
        transform = ax.get_xaxis_transform(),
        color = "#C5E8C3",
        alpha = 0.8,
        clip_on = False
    )
    ax.add_patch(triangle)
    if add_text:
        ax.text(1.04, 1.05, "Age",
            transform=ax.transAxes,
            ha='center', va='center',
            fontsize=text_fontsize, fontweight = "bold",
            color='black', clip_on=False)

def add_legend_patch(ax, fontsize):
    x_offset = 1.02  # slight right of the plot
    y_center = 0.95  # top of the legend area
    
    # Box height
    rect_height = 0.03

    # Draw the two stacked boxes
    ax.add_patch(patches.Rectangle((x_offset, y_center + rect_height/2), 0.05, rect_height, clip_on=False,
                                   transform=ax.transAxes, fill=True, facecolor="#C6D2EA", edgecolor='black', linewidth=0.3))
    ax.text(x_offset + rect_height - 0.005, y_center + rect_height, "N",
            transform=ax.transAxes,
            ha='center', va='center',
            fontsize=fontsize,
            color='black', clip_on=False)
    
    ax.add_patch(patches.Rectangle((x_offset, y_center - rect_height/2), 0.05, rect_height, clip_on=False,
                                   transform=ax.transAxes, fill=True, facecolor="#C6D2EA", edgecolor='black', linewidth=0.3))
    ax.text(x_offset + rect_height - 0.005, y_center, "N",
            transform=ax.transAxes,
            ha='center', va='center',
            fontsize=fontsize,
            color='black', clip_on=False)

    # Add text labels for the rectangles
    ax.text(x_offset + 0.06, y_center + rect_height/2, "Driver SNVs",
            transform=ax.transAxes, va='bottom', ha='left', fontsize=fontsize)
    ax.text(x_offset + 0.06, y_center + 0.01, "Protein-affecting\nindels",
            transform=ax.transAxes, va='top', ha='left', fontsize=fontsize)

    # Add black dot and label for "Donor with 2 samples"
    dot_y = y_center - 0.1  # place below the boxes
    ax.plot(x_offset + 0.015, dot_y, 'o', transform=ax.transAxes, color='black', markersize=4, clip_on=False)
    ax.text(x_offset + 0.03, dot_y, "Donor with\n2 samples",
            transform=ax.transAxes, va='center', ha='left', fontsize=fontsize, clip_on=False)

def add_dots_below_xticks(ax, dot_samples, offset = 28, size = 1.5, color = "black"):
    """
    Add dots below x-axis tick labels for specific samples.
    """

    xticks = ax.get_xticks()
    xticklabels = [label.get_text() for label in ax.get_xticklabels()]

    for tick, label in zip(xticks, xticklabels):
        if label in dot_samples:
            ax.plot(tick, offset, marker = 'o', color = color, markersize = size, 
                    clip_on = False, zorder = 10)
    
## -- Main functions -- ##

def plot_double_heatmap(data, heatmap_config, save_file):
    """
    """

    # initialize figure
    fig, axs = plt.subplots(1, 2, figsize = heatmap_config["figsize"], sharey = False,
                            gridspec_kw = {'width_ratios': heatmap_config["width_ratios"]})

    # normalize data for heatmap (gene-wise)
    data = data.reindex(heatmap_config["genes_order"])
    data_norm = (data.T / data.max(axis='columns')).T

    # plot h1
    sns.heatmap(
    data_norm[heatmap_config["h1_subset"]], 
    annot = data[heatmap_config["h1_subset"]],
    fmt = "", 
    cmap = heatmap_config["cmap"], 
    ax = axs[0],
    annot_kws = {"size": heatmap_config["annot_fontsize"]},
    cbar = False)
    add_age_triangle(axs[0], len(heatmap_config["h1_subset"]))
    duplicated_samples = [s for s in heatmap_config["duplicated_samples"] if s in heatmap_config["h1_subset"]]
    add_dots_below_xticks(axs[0], duplicated_samples)

    ## redefine gene names
    genes = [gene.split("_")[0] for gene in data.index[::2]]
    yticks_positions = [i + 1 for i in range(0, len(data), 2)]
    axs[0].set_yticks(yticks_positions)
    axs[0].set_yticklabels(genes, rotation = 0, va = "center",
                            fontsize = heatmap_config["xyticks_fontsize"])
    axs[0].tick_params(axis = 'x', labelsize = heatmap_config["xyticks_fontsize"])
    axs[0].set_xlabel(heatmap_config["h1_title"], fontsize = heatmap_config["title_fontsize"], fontweight = "bold")

    
    # plot h2
    sns.heatmap(
    data_norm[heatmap_config["h2_subset"]], 
    annot = data[heatmap_config["h2_subset"]],
    fmt = "", 
    cmap = heatmap_config["cmap"], 
    ax = axs[1],
    annot_kws = {"size": heatmap_config["annot_fontsize"]},
    cbar = False)
    add_age_triangle(axs[1], len(heatmap_config["h2_subset"]), 
                    add_text = True, text_fontsize = heatmap_config["annot_fontsize"])
    duplicated_samples = [s for s in heatmap_config["duplicated_samples"] if s in heatmap_config["h2_subset"]]
    add_dots_below_xticks(axs[1], duplicated_samples)

    axs[1].set_yticks([])
    axs[1].set_yticklabels([])
    axs[1].tick_params(axis = 'x', labelsize = heatmap_config["xyticks_fontsize"])
    axs[1].set_xlabel(heatmap_config["h2_title"], fontsize = heatmap_config["title_fontsize"], fontweight = "bold")

    # separate genes
    for row in heatmap_config["rows2separate"]:
        axs[0].hlines(y = row, xmin = 0, xmax = len(heatmap_config["h1_subset"]),
                    colors = "white", linewidth = 0.5)
        axs[1].hlines(y = row, xmin = 0, xmax = len(heatmap_config["h2_subset"]),
                    colors = "white", linewidth = 0.5)

    cbar_ax = fig.add_axes(heatmap_config["cmap_loc"])  # [left, bottom, width, height]
    norm = plt.Normalize(vmin = data_norm.min().min(), vmax = data_norm.max().max())
    sm = plt.cm.ScalarMappable(cmap = heatmap_config["cmap"], norm = norm)
    cbar = fig.colorbar(sm, cax = cbar_ax, orientation = heatmap_config["cmap_orient"])
    cbar.set_label(heatmap_config["cmap_label"], fontsize = heatmap_config["cmap_label_fontsize"])
    cbar.outline.set_visible(False)
    cbar.ax.tick_params(labelsize = heatmap_config["cmap_ticks_fontsize"])

    add_legend_patch(axs[1], fontsize = heatmap_config["annot_fontsize"])

    axs[0].set_ylabel("")
    axs[1].set_ylabel("")

    plt.tight_layout()
    # save
    plt.savefig(save_file, dpi = 300, #TODO: resolution to general config
                bbox_inches = 'tight')

    return None

def km_plot(data, time_var, res_var, categs_dict, save_file,
            comp_var = None, plot_config = {}, do_logrank = True):
    """
    """

    fig, ax = plt.subplots(1, 1, figsize = plot_config["figsize"])

    # calculate KM curves and plot
    kmf = KaplanMeierFitter()
    for categ in categs_dict:
        if categ == "no_categ":
            data_categ = data
        else:
            data_categ = data.loc[data[comp_var] == categ]

        kmf.fit(durations = data_categ[time_var], event_observed = data_categ[res_var],
                label = f'{categs_dict[categ][1]} (n = {len(data_categ)})')
        kmf.plot_cumulative_density(ax = ax, ci_show = False, at_risk_counts = False,
                                    color = categs_dict[categ][0])
    
    # plot config
    ax.set_xlabel(plot_config["xlabel"], fontsize = plot_config["xlabel_fontsize"])
    ax.set_ylabel(plot_config["ylabel"], fontsize = plot_config["ylabel_fontsize"])
    ax.tick_params(axis = "x", labelsize = plot_config["xticks_fontsize"])
    ax.tick_params(axis = "y", labelsize = plot_config["yticks_fontsize"])
    ax.legend(fontsize = plot_config["legend_fontsize"], frameon = False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    # if indicated, do log-rank test to compare the categories of the comp_var
    if do_logrank and comp_var is not None:

        if len(data[comp_var].unique()) == 2:
            pval = logrank_test(data.loc[data[comp_var] == list(categs_dict.keys())[0]][time_var],
                                data.loc[data[comp_var] == list(categs_dict.keys())[1]][time_var],
                                event_observed_A = data.loc[data[comp_var] == list(categs_dict.keys())[0]][res_var], 
                                event_observed_B = data.loc[data[comp_var] == list(categs_dict.keys())[1]][res_var]).p_value

        elif len(data[comp_var].unique()) > 2:
            pval = multivariate_logrank_test(data[time_var], data[comp_var], data[res_var]).pvalue

        ax.add_artist(AnchoredText(f'p-val < {pval:.2e}',
                    frameon = False, loc = 'lower right', 
                    prop = dict(size = plot_config["text_fontsize"])))

    # save
    plt.savefig(save_file, dpi = plot_config["figsave_resolution"], 
                bbox_inches = 'tight')

    return None

def regr_res_scatterplot(metric, clinvar, data_df, regrres_df, plot_config, plot_config2, 
                        plots_general_config, save_file, add_legend = True, cut_at_zero = False,
                        common_ylabel = False, text2add = None):
    """
    """
    if cut_at_zero:
        regrres_df = regrres_df.loc[regrres_df["predicted"] >= 0]

    nplots = 1
    if len(plot_config) == 2:
        nplots = 2
    elif len(plot_config) > 2:
        print("cannot make more than 2 plots")
        return None
        
    fig, axs = plt.subplots(nplots, 1, figsize = plot_config2["figsize"], sharex = True)

    if nplots == 1:
        n = list(plot_config.keys())[0]
        responses = plot_config[n].keys() 
        sns.scatterplot(
            data = data_df.loc[data_df["response"].isin(responses)],
            x = clinvar, y = metric, hue = "response",
            palette = plot_config[n], ax = axs,
            legend = False, s = plots_general_config["dot_size_scplot"],
            linewidth = 0.1  # Reduced edge width
        )
        sns.lineplot(data = regrres_df.loc[regrres_df["response"].isin(responses)],
                    x = clinvar, y = "predicted", hue = "response",
                    palette = plot_config[n], ax = axs, legend = False)

        axs.set_ylabel(plot_config2["ylabel"], fontsize = plots_general_config["xylabel_fontsize"]) #TODO: fontsize to general config
        axs.set_xlabel(plot_config2["xlabel"], fontsize = plots_general_config["xylabel_fontsize"]) #TODO: fontsize to general config
        if add_legend:
            legend_handles = [plt.Line2D([0], [0], marker = 'o', color = 'w', markerfacecolor = color, markersize = 7) for color in plot_config[n].values()]
            axs.legend(handles = legend_handles, labels = [" ".join(g.split("_")) for g in plot_config[n].keys()], 
                    title = "", fontsize = plots_general_config["legend_fontsize"], frameon = False, loc = "upper left") #TODO: fontsize to general config
        axs.spines['right'].set_visible(False)
        axs.spines['top'].set_visible(False)
        axs.tick_params(axis = 'x', labelsize = plots_general_config["xyticks_fontsize"])  
        axs.tick_params(axis = 'y', labelsize = plots_general_config["xyticks_fontsize"])  
        axs.set_title(plot_config2["title"], fontsize = plots_general_config["title_fontsize"]) #TODO: fontsize to general config

        if text2add != None: # only implemented for 1 plot mode!
            x_pos = regrres_df[clinvar].min()
            y_pos = max(regrres_df[metric].max(), data_df[metric].max())
            axs.text(x_pos, y_pos, text2add, color = 'black',
                ha = 'left', va = 'center', fontsize = plots_general_config["legend_fontsize"])

    elif nplots == 2:
        
        axs = axs.flatten()
        for i, n in enumerate(plot_config.keys()):
            responses = plot_config[n].keys() 
            sns.scatterplot(data = data_df.loc[data_df["response"].isin(responses)],
                            x = clinvar, y = metric, hue = "response",
                            palette = plot_config[n], ax = axs[i],
                            legend = False, s = plots_general_config["dot_size_scplot"])
            sns.lineplot(data = regrres_df.loc[regrres_df["response"].isin(responses)],
                        x = clinvar, y = "predicted", hue = "response",
                        palette = plot_config[n], ax = axs[i], legend = False)
    
            axs[i].set_ylabel(plot_config2["ylabel"], fontsize = plots_general_config["xylabel_fontsize"]) #TODO: fontsize to general config
            axs[i].set_xlabel(plot_config2["xlabel"], fontsize = plots_general_config["xylabel_fontsize"]) #TODO: fontsize to general config
            if add_legend:
                legend_handles = [plt.Line2D([0], [0], marker = 'o', color = 'w', markerfacecolor = color, markersize = 7) for color in plot_config[n].values()]
                axs[i].legend(handles = legend_handles, labels = [" ".join(g.split("_")) for g in plot_config[n].keys()], 
                        title = "", fontsize = plots_general_config["legend_fontsize"], frameon = False, loc = "upper left") #TODO: fontsize to general config
            axs[i].spines['right'].set_visible(False)
            axs[i].spines['top'].set_visible(False)
            axs[i].tick_params(axis = 'x', labelsize = plots_general_config["xyticks_fontsize"])  
            axs[i].tick_params(axis = 'y', labelsize = plots_general_config["xyticks_fontsize"]) 
            if i == 0:
                axs[i].set_title(plot_config2["title"], fontsize = plots_general_config["title_fontsize"]) #TODO: fontsize to general config

    if common_ylabel:
        axs[i-1].set_ylabel("")
        axs[i].set_ylabel("")
        fig.text(-0.05, 0.5, # TODO: pass coords to config
                plot_config2["common_ylabel"], va = 'center', ha = 'center', rotation = 'vertical', fontsize = plots_general_config["xylabel_fontsize"])

        
    # save
    plt.savefig(save_file, dpi = 300, #TODO: resolution to general config
                bbox_inches = 'tight')

    return None

def regr_res_coeffplot(regrres_df, plot_config, plots_general_config, responses, 
                    save_file, regrres2compare_df = None, figsize = (3, 5), remove_ylabels = False,
                    colored = False, xscale_log = False, write_coeff = False, write_dndstype = False,
                    add_arrow = True):
    """
    """

    fig, ax = plt.subplots(1, 1, figsize = figsize) 
    regrres_df = regrres_df.set_index("gene") #TODO: maybe is better to have/not have gene set as index differently in the functions of this nb
    if regrres2compare_df is not None:
        regrres2compare_df = regrres2compare_df.set_index("gene")
        main_dot_color = plots_general_config["dot_colorabove_coeffplot"]
    else:
        main_dot_color = plots_general_config["dot_color_coeffplot"]

    # plot one line per coefficient
    go_back = 0
    for i, resp in enumerate(responses):

        ## check whether the coeff is significant to highlight the point in black
        if regrres_df.loc[resp]["qval"] < plot_config["sign_threshold"]:
            edgecolors = "black"
        else:
            edgecolors = "face"

        ## plot dot + confidence interval
        if colored:
            main_dot_color = gene2color[resp]
        if resp == "CDKN1At":
            i -= 0.85
        ax.scatter(regrres_df.loc[resp]["coeff"], i-go_back, 
                color = main_dot_color,
                s = plots_general_config["dot_size_coeffplot"],
                edgecolors = edgecolors, linewidths = plots_general_config["dot_edgewidth_coeffplot"])
        if write_coeff:
            ax.text(regrres_df.loc[resp]["coeff"], i + plot_config["writecoeff_offset"],
                    # f"{regrres_df.loc[resp]['coeff']:.2f} ({regrres_df.loc[resp]['qval']:.3f})",
                    f"{regrres_df.loc[resp]['qval']:.3f}", 
                    ha = "center", va = "bottom", 
                    fontsize = plot_config["coeff_fontsize"], color = "black")
        ax.hlines(i-go_back, regrres_df.loc[resp]["lowci"], regrres_df.loc[resp]["highci"],
                color = main_dot_color)
        if write_dndstype:
            if resp == "CDKN1Am":
                dndstype = "Missense"
                extra_right = 3
            else:
                dndstype = "Truncating"
                extra_right = 2
            ax.text(regrres_df.loc[resp]["highci"]+extra_right, i-go_back,
                    dndstype, 
                    ha = "left", va = "center", 
                    fontsize = 8, color = "black")
        if resp == "CDKN1At":
            go_back = 1

        ## plot an additional coefficient line below to compare
        if regrres2compare_df is not None:
            
            if regrres2compare_df.loc[resp]["qval"] < plot_config["sign_threshold"]:
                edgecolors = "black"
            else:
                edgecolors = "face"
            
            offset = 0.3
            ax.scatter(regrres2compare_df.loc[resp]["coeff"], i - offset, 
                color = plots_general_config["dot_colorbelow_coeffplot"],
                s = plots_general_config["dot_size_coeffplot"],
                edgecolors = edgecolors, linewidths = plots_general_config["dot_edgewidth_coeffplot"])
            if write_coeff:
                ax.text(regrres2compare_df.loc[resp]["coeff"], i - offset + plot_config["writecoeff_offset"],
                        # f"{regrres2compare_df.loc[resp]['coeff']:.2f} ({regrres2compare_df.loc[resp]['qval']:.3f})", 
                        f"{regrres2compare_df.loc[resp]['qval']:.3f}", 
                        ha = "center", va = "bottom", 
                        fontsize = plot_config["coeff_fontsize"], color = "black")
            ax.hlines(i - offset, regrres2compare_df.loc[resp]["lowci"], regrres2compare_df.loc[resp]["highci"],
                    color = plots_general_config["dot_colorbelow_coeffplot"])

    # add labels to the plotted coefficients
    y_ticks = [i for i in range(len(responses))]
    if regrres2compare_df is not None:
        y_labels = [plot_config["ylabels"][0] for i in range(len(responses))]
        for i in range(len(responses)):
            y_ticks.append(i - offset)
            y_labels.append(plot_config["ylabels"][1])
        for i, g in enumerate(responses):
            y_coord = i - offset
            x_coord = plot_config["ylabels_gene_xcoord"]
            if g == "total":
                g = "All genes"
            ax.text(x_coord, y_coord, g, fontsize = plots_general_config["xyticks_fontsize"], ha = "right")


    else:
        y_labels = responses

    y_labels = ["All genes" if l == "total" else l for l in y_labels]
    
    if remove_ylabels:
        ax.set_yticks([], [])
    else:
        ax.set_yticks(y_ticks, y_labels, fontsize = plots_general_config["xyticks_fontsize"])

    # set zero effect in x axis
    ax.vlines(plot_config["null_effect"], -0.7, len(responses)-0.5, ls = '--', color = 'grey')

    # subplot config
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    # ax.set_xlabel('Effect size', fontsize = 8.5)
    if add_arrow:
        ax.arrow(plot_config["null_effect"],  plot_config["arrow_yloc"], plot_config["arrow_xlim"], 0, clip_on = False,
                head_width = plot_config["arrow_head_width"], head_length = plot_config["arrow_head_length"],
                color = "black", alpha = 0.65)
        ax.text(plot_config["effectsize_text_loc"], plot_config["effectsize_text_yloc"], plot_config["effectsize_text"],
                fontsize = plots_general_config["xylabel_fontsize"])
    ax.set_ylim(-0.7, len(responses)-0.5-go_back)
    ax.set_title(plot_config["title"], fontsize = plots_general_config["title_fontsize"])
    ax.tick_params(axis = 'x', labelsize = plots_general_config["xyticks_fontsize"])  
    # ax.tick_params(axis = 'y', labelsize = plots_general_config["xyticks_fontsize"]) 
    if xscale_log:
        ax.set_xscale('log')

    # save
    plt.savefig(save_file, dpi = 300, bbox_inches = 'tight') #TODO: resolution to general config

    return None

def regr_res_coeffplot_multi(regrres_df_dict, plot_config, plots_general_config, responses, 
                    save_file, regrres2compare_df = None, figsize = (3, 5), remove_ylabels = False,
                    colored = False, xscale_log = False, write_coeff = False):
    """
    """

    fig, axs = plt.subplots(1, len(regrres_df_dict), figsize = figsize, sharey = True) 

    for j, regrres_df in enumerate(regrres_df_dict.values()):
        regrres_df = regrres_df.set_index("gene") #TODO: maybe is better to have/not have gene set as index differently in the functions of this nb
        if regrres2compare_df is not None:
            regrres2compare_df = regrres2compare_df.set_index("gene")
            main_dot_color = plots_general_config["dot_colorabove_coeffplot"]
        else:
            main_dot_color = plots_general_config["dot_color_coeffplot"]

        # plot one line per coefficient
        for i, resp in enumerate(responses):

            ## check whether the coeff is significant to highlight the point in black
            if regrres_df.loc[resp]["qval"] < plot_config["sign_threshold"]:
                edgecolors = "black"
            else:
                edgecolors = "face"

            ## plot dot + confidence interval
            if colored:
                main_dot_color = gene2color[resp]
            axs[j].scatter(regrres_df.loc[resp]["coeff"], i, 
                    color = main_dot_color,
                    s = plots_general_config["dot_size_coeffplot"],
                    edgecolors = edgecolors, linewidths = plots_general_config["dot_edgewidth_coeffplot"])
            if write_coeff:
                axs[j].text(regrres_df.loc[resp]["coeff"], i + plot_config["writecoeff_offset"],
                        # f"{regrres_df.loc[resp]['coeff']:.2f} ({regrres_df.loc[resp]['qval']:.3f})", 
                        f"{regrres_df.loc[resp]['qval']:.3f}", 
                        ha = "center", va = "bottom", 
                        fontsize = plot_config["coeff_fontsize"], color = "black")
            axs[j].hlines(i, regrres_df.loc[resp]["lowci"], regrres_df.loc[resp]["highci"],
                    color = main_dot_color)

            ## plot an additional coefficient line below to compare
            if regrres2compare_df is not None:
                
                if regrres2compare_df.loc[resp]["qval"] < plot_config["sign_threshold"]:
                    edgecolors = "black"
                else:
                    edgecolors = "face"
                
                offset = 0.3
                axs[j].scatter(regrres2compare_df.loc[resp]["coeff"], i - offset, 
                    color = plots_general_config["dot_colorbelow_coeffplot"],
                    s = plots_general_config["dot_size_coeffplot"],
                    edgecolors = edgecolors, linewidths = plots_general_config["dot_edgewidth_coeffplot"])
                if write_coeff:
                    axs[j].text(regrres2compare_df.loc[resp]["coeff"], i - offset + plot_config["writecoeff_offset"],
                        # f"{regrres2compare_df.loc[resp]['coeff']:.2f} ({regrres2compare_df.loc[resp]['qval']:.3f})", 
                        f"{regrres2compare_df.loc[resp]['qval']:.3f}", 
                        ha = "center", va = "bottom", 
                        fontsize = plot_config["coeff_fontsize"], color = "black")
                axs[j].hlines(i - offset, regrres2compare_df.loc[resp]["lowci"], regrres2compare_df.loc[resp]["highci"],
                        color = plots_general_config["dot_colorbelow_coeffplot"])

        # add labels to the plotted coefficients
        y_ticks = [i for i in range(len(responses))]
        if regrres2compare_df is not None:
            y_labels = [plot_config["ylabels"][0] for i in range(len(responses))]
            for i in range(len(responses)):
                y_ticks.append(i - offset)
                y_labels.append(plot_config["ylabels"][1])
            for i, g in enumerate(responses):
                y_coord = i - offset/2
                x_coord = plot_config["ylabels_gene_xcoord"]
                if g == "total":
                    g = "All genes"
                axs[j].text(x_coord, y_coord, g, fontsize = plots_general_config["xyticks_fontsize"], ha = "right")


        else:
            y_labels = responses

        y_labels = ["All genes" if l == "total" else l for l in y_labels]
        
        if remove_ylabels:
            axs[j].set_yticks([], [])
        else:
            axs[j].set_yticks(y_ticks, y_labels, fontsize = plots_general_config["xyticks_fontsize"])

        # set zero effect in x axis
        axs[j].vlines(plot_config["null_effect"], -0.7, len(responses)-0.5, ls = '--', color = 'grey')

        # subplot config
        axs[j].spines['top'].set_visible(False)
        axs[j].spines['right'].set_visible(False)
        # ax.set_xlabel('Effect size', fontsize = 8.5)
        axs[j].arrow(plot_config["null_effect"],  plot_config["arrow_yloc"][j], plot_config["arrow_xlim"][j], 0, clip_on = False,
                head_width = plot_config["arrow_head_width"][j], head_length = plot_config["arrow_head_length"][j],
                color = "black", alpha = 0.65)
        axs[j].text(plot_config["effectsize_text_loc"][j], plot_config["effectsize_text_yloc"][j], plot_config["effectsize_text"][j],
                fontsize = plots_general_config["xylabel_fontsize"])
        axs[j].set_ylim(-0.7, len(responses)-0.5)
        axs[j].set_title(plot_config["title"][j], fontsize = plots_general_config["title_fontsize"])
        axs[j].tick_params(axis = 'x', labelsize = plots_general_config["xyticks_fontsize"])  
        # ax.tick_params(axis = 'y', labelsize = plots_general_config["xyticks_fontsize"]) 
        if xscale_log:
            axs[j].set_xscale('log')

    # save
    plt.savefig(save_file, dpi = 300, bbox_inches = 'tight') #TODO: resolution to general config

    return None

# below: Abel code to create sigmoid plots

def prepare_twin_plot_data(data, clinvar,genes, def_colors, metric):
    data = data[['gene','sample',metric,clinvar]].dropna()
    data['color'] = data.apply(do_colors,args=(def_colors,clinvar,),axis=1)
    return prepare_twin_sigmoid_plot(data, genes,clinvar,metric,def_colors)

def prepare_twin_sigmoid_plot(data,genes,clinvar,metric,def_colors):
    prepared_data = []
    clinvar_categs = def_colors.keys()

    for gene in genes:
        for i, categ in enumerate(clinvar_categs):
            if i == 0:
                extra = 0
            elif i == 1:
                extra = 1
            elif i == 2:
                extra = 2.5

            gene_data = data[(data['gene']==gene)&(data[clinvar]==categ)].sort_values(metric)
            gene_data['y'] = gene_data[metric]
            gene_data['x'] = np.linspace(0.1+extra,1+extra,len(gene_data))
            prepared_data.append(gene_data)

    return pd.concat(prepared_data)

def do_colors(r,def_colors,clinvar):
    return def_colors[r[clinvar]]

def draw_vertical_sigmoid(axis,data,clinvar,def_colors,dot_size,linewidth):
    clinvar_categs = def_colors.keys()
    axis.scatter(data['x'], data['y'], color='white', zorder=3, lw=1, ec="white", s=dot_size)
    axis.scatter(data['x'].values, data['y'].values, zorder=4,alpha=0.5, lw=0.1, ec=data['color'], s=dot_size*2, color=data['color'])

    for categ in clinvar_categs:
        axis.hlines(np.median(data[data[clinvar]==categ]['y'].values),
                np.min(data[data[clinvar]==categ]['x'].values),
                np.max(data[data[clinvar]==categ]['x'].values),
                linewidth=linewidth,
                color=def_colors[categ])

def tidy_sigmoid_axis(axis,i,genes,gene_data,clinvar,def_colors,plot_config,plots_general_config,mode):
    clinvar_categs = list(def_colors.keys())
    if i == 0:
        axis.spines[['top','right']].set_visible(False)
        axis.set_ylabel(plot_config["ylabel"], fontsize = plots_general_config["xylabel_fontsize"])
    else:
        axis.spines[['top', 'right']].set_visible(False)
        axis.set_ylabel('')

    N_1 = len(gene_data[gene_data[clinvar]==clinvar_categs[0]])
    N_2 = len(gene_data[gene_data[clinvar]==clinvar_categs[1]])

    if mode == 'twin':
        if len(clinvar_categs) == 2:
            axis.set_xticks([0.5,1.5],[f'{clinvar_categs[0]}\n(N={str(N_1)})', f'{clinvar_categs[1]}\n(N={str(N_2)})'], fontsize = plots_general_config["xyticks_fontsize"]-2)
            axis.set_xlim([-0.1,2.1])
        elif len(clinvar_categs) == 3:
            N_3 = len(gene_data[gene_data[clinvar]==clinvar_categs[2]])
            axis.set_xticks([0.5, 1.5, 3],
                            [f'{clinvar_categs[0]}', f'{clinvar_categs[1]}', f'{clinvar_categs[2]}'],
                            fontsize = plots_general_config["xyticks_fontsize"])
            axis.text(0.5, axis.get_ylim()[1] + 0.05, f'N = {str(N_1)}', ha='center', fontsize=plots_general_config["xyticks_fontsize"] - 2)
            axis.text(1.5, axis.get_ylim()[1] + 0.05, f'N = {str(N_2)}', ha='center', fontsize=plots_general_config["xyticks_fontsize"] - 2)
            axis.text(3, axis.get_ylim()[1] + 0.05, f'N = {str(N_3)}', ha='center', fontsize=plots_general_config["xyticks_fontsize"] - 2)
            axis.set_xlim([-0.1,3.6])
        else:
            print("this function only works with variables with two or three categories")
    axis.set_xlabel(plot_config["xlabel"], fontsize = plots_general_config["xylabel_fontsize"])
    if genes[i] in plot_config["titles"]:
        gene4title = plot_config["titles"][genes[i]]
    else:
        gene4title = genes[i]

    axis.set_title(gene4title, fontsize=plots_general_config["title_fontsize"])
    axis.tick_params(axis='y', labelsize=plots_general_config["xyticks_fontsize"])

def plot_sigmoids(figure_output_file,
                omega_data,
                clinvar,
                genes,
                def_colors,
                plot_config,
                plots_general_config,
                mode='together'):

    fig = plt.figure(figsize=plot_config["figsize"])
    gs = GridSpec(plot_config["grid_nrows"], plot_config["grid_ncols"], figure=fig, wspace=10, hspace=0)

    ax = {}

    for i in range(len(genes)):
        ax[genes[i]] = fig.add_subplot(gs[0:plot_config["grid_nrows"]-2, 9*i:9*i+8])
        draw_vertical_sigmoid(ax[genes[i]], omega_data[omega_data['gene']==genes[i]], clinvar, def_colors,
                            plot_config["dot_size"], plot_config["line_width"])
        tidy_sigmoid_axis(ax[genes[i]], i, genes, omega_data[omega_data['gene']==genes[i]], clinvar, def_colors, plot_config,
        plots_general_config, mode=mode)

    plt.savefig(figure_output_file, bbox_inches='tight', dpi = 300)
