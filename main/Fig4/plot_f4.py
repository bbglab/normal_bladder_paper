import matplotlib.pyplot as plt

def plot_count_track_tert(count_df,
                        axes, 
                        colors_dict,
                        ax=0, 
                        alpha=1,
                        dot_size = 60):
    
    
    for cnsq in ['missense', 'synonymous']:
    
        count_cnsq_df = count_df[count_df["Consequence"] == cnsq].reset_index(drop=True)

        axes[ax].vlines(count_cnsq_df["Pos"], ymin=0, ymax=count_cnsq_df["Count"], lw=1, zorder=1, alpha=0.5, color=colors_dict["hv_lines_needle"])
        axes[ax].scatter(count_cnsq_df["Pos"].values, count_cnsq_df["Count"].values, zorder=4,
                        alpha=alpha, lw=0, ec="black", s=dot_size,
                        label='Mutation observed in tumors (activating)' if cnsq == 'missense' else 'Mutation not observed in tumors',
                        color=colors_dict["omega_synon_tert" if cnsq == 'synonymous' else "omega_miss_tert"]) 

def needle_plot(data, subset, subset_name, plot_config, plots_general_config, save_file,
                add_legend = True, add_xticks = True, add_xlabel = True):
    """
    """
    
    # create count df
    counts_per_position = data[data["SAMPLE_ID"].isin(subset)].groupby(
        by = ['MUT_ID', 'Consequence', 'POS'])['ALT_DEPTH'].size().to_frame('Count').reset_index()
    counts_per_position.columns = ['Mutation', 'Consequence', 'Pos', 'Count']

    # plot needle
    fig, ax = plt.subplots(1,1, figsize = plot_config["figsize"])
    plot_count_track_tert(counts_per_position, axes=[ax], ax=0,
                        colors_dict= plot_config["colors"], 
                        alpha = 1, dot_size = plot_config["dot_size"],
                        )
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_ylim(0, 22)
    if add_xlabel:
        ax.set_xlabel("Position in chromosome 5", 
                    fontsize = plots_general_config["xylabel_fontsize"])
    else:
        ax.set_xlabel("")
    ax.set_ylabel("Mutation count", fontsize = plots_general_config["xylabel_fontsize"])
    ax.text(1294951, 17, f"{subset_name} (N={len(subset)})", fontsize = plots_general_config["xylabel_fontsize"],
            fontweight = "bold")

    ytickss = [x for x in ax.get_yticks()[0:-1:1] if not (x % 1) ] 
    ax.set_yticks(ytickss)
    ax.set_yticklabels([f"{int(x):,}" for x in ytickss], 
                    fontsize = plots_general_config["xyticks_fontsize"])

    ax.set_xlim(1294942, 1295292)

    if add_xticks:
        xtickss = ax.get_xticks()[1:-1:1]
        ax.set_xticks(xtickss)
        ax.set_xticklabels([f"{int(x):,}" if not x%100 else "" for x in xtickss], 
                        fontsize = plots_general_config["xyticks_fontsize"])

    else:
        ax.tick_params(axis='x', labelbottom=False)

    if add_legend:
        fig.legend(frameon = False, bbox_to_anchor = (0.9, 0.9),
                fontsize = plots_general_config["legend_fontsize"])
    fig.savefig(save_file, bbox_inches = 'tight', dpi = 300)
    plt.show()