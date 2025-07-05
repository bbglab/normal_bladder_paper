def plot_count_track_tert(count_df,
                        axes, 
                        colors_dict,
                        ax=0, 
                        negative=False, 
                        label_pos_track=None,
                        label_neg_track=None,
                        ymargin=None,
                        alpha=1,
                        indel=False,
                        dot_size = 60):
    
    
    for cnsq in ['missense', 'synonymous']:
    
        count_cnsq_df = count_df[count_df["Consequence"] == cnsq].reset_index(drop=True)

        axes[ax].vlines(count_cnsq_df["Pos"], ymin=0, ymax=count_cnsq_df["Count"], lw=1, zorder=1, alpha=0.5, color=colors_dict["hv_lines_needle"])
        axes[ax].scatter(count_cnsq_df["Pos"].values, count_cnsq_df["Count"].values, zorder=4,
                        alpha=alpha, lw=0, ec="black", s=dot_size,
                        label='Mutation observed in tumors (activating)' if cnsq == 'missense' else 'Mutation not observed in tumors',
                        color=colors_dict["omega_synon_tert" if cnsq == 'synonymous' else "omega_miss_tert"]) 