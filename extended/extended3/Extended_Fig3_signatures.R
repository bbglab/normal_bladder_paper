library(reshape2)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(grid)
library(cowplot)

deepCSA_run <- "/data/bbg/nobackup/bladder_ts/results/2025-05-14_deepCSA_45_donors"
outplot <- file.path("./plots/Extended_Fig3_signatures.png")

plot_sig_profiles <- function(deepCSA_run, extraction_method){

    if (extraction_method == "SigProfiler"){
        sig_path <- paste("../../supplementary/supp_note_5/output/SigProfilerExtractor_out.med.all/SBS96/Suggested_Solution/SBS96_De-Novo_Solution/Signatures/SBS96_De-Novo_Signatures.txt", sep="")
        df <-  read.table(sig_path, sep="\t", header=TRUE, row.names=1)
        sigs_names <- c("SBS96A","SBS96B","SBS96C")
    }else{
        sig_path <- paste(deepCSA_run, "/signatures_hdp/samples_matrix.all.compared_output_dir/signatures.txt", sep="")
        df <- read.table(sig_path, sep="\t", header=TRUE, row.names=1)
        sigs_names <- c("N1","N3","N2", "N4","N5")

    }
    sigs <- t(df)
    # plot components
    # SBS #
    sigs_df <- reshape2::melt(sigs)
    colnames(sigs_df) <- c('component', 'channel', 'value')
    sigs_df$group <- substr(as.character(sigs_df$channel), start = 3, stop = 5)
    sigs_df <- sigs_df %>%
        mutate(channel = factor(channel, levels =  colnames(sigs)),
        group = factor(group, levels = unique(sigs_df$group)),
        value = value * 100)
  
    #set colours
    colours <- setNames(c('#1ebff0', '#050708', '#e62725', '#cbcacb', '#a1cf64', '#edc8c5'), unique(sigs_df$group))
    strip_name_colours <- c('black','white','white','black','black','black')
    ####

    plots <- list()
    print(unique(sigs_df$component))
    for (x in sigs_names) {
        plot_data <- sigs_df[sigs_df$component == x,]
        plot_data$xlabels <- paste0(substr(as.character(plot_data$channel), 1, 1),
            substr(as.character(plot_data$channel), 3, 3),
            substr(as.character(plot_data$channel), 7, 7))
  
    p <- ggplot(plot_data, aes(x = xlabels, y = value, fill = group)) +
        geom_bar(stat = 'identity') + 
        facet_grid(. ~ group, space = 'free_x', scales = 'free_x') +
        scale_fill_manual(name = '', values = colours, guide = 'none') +
        xlab('') +
        ylab('% SBS') +
        scale_y_continuous(expand = c(0,0,0.05,0)) +
        ggtitle(paste0(x)) +
        theme_classic() + 
        theme(plot.title = element_text(hjust = 1, size = 15,margin=margin(t=30,b=-40)), 
            #axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
            strip.text = element_text(face = 'bold'),
            axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank())

    # Change colors and labels of strips
    g <- ggplot_gtable(ggplot_build(p))
    striprt <- which(grepl('strip-r', g$layout$name) | grepl('strip-t', g$layout$name))
    k <- 1
    for (i in striprt) {
        j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
        g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- colours[k]
        g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$col  <- colours[k]
    
        t <- which(grepl('text', g$grobs[[i]]$grobs[[1]]$childrenOrder))
        g$grobs[[i]]$grobs[[1]]$children[[t]]$children[[1]]$gp$col <- strip_name_colours[k]
    
        k <- k + 1
    }

    # Add the plot to the list
    plots[[length(plots) + 1]] <- g
    }

    # Determine the number of rows and columns for the grid
    n_col <- 1  # Number of columns in the grid
    n_row <- ceiling(length(plots) / n_col)

    profile_plot <- grid.arrange(grobs = plots, ncol = n_col, nrow = n_row, top = extraction_method)
    return(profile_plot)

}

plot_exposures_HDP <- function(deepCSA_run){
    colors <-c("#005f73","#94d2bd","#0a9396","#e9d8a6", "#ee9b00","#f37674", "#9b2226")
    all_mut_df <- read.table(paste("../..//supplementary/supp_note_5/data/Mutation_stats_per_genome_per_sample.med.all.csv", sep=""), sep=",", header =TRUE)
    samples_order = all_mut_df[order(all_mut_df$Total_mut_genome, decreasing=TRUE),]$SAMPLE_ID
    exposures_file = read.table(paste(deepCSA_run, "/signatures_hdp/samples_matrix.all.compared_output_dir/signatureExposures_counts.txt", sep=""), sep="\t", header=T)
    names_list <- c("N1", "N2", "N3", "N4", "N5", "N6", "N7", "N8", "N9", "N10", "N11", "N12", "N13", "N14", "N15", "N16")
    colnames(exposures_file) <- names_list[1:ncol(exposures_file)]
    n_sigs <- ncol(exposures_file)
    exposures_file$Samples <- rownames(exposures_file)
    exposures_file_melted = melt(exposures_file)
    colnames(exposures_file_melted)<- c("sample", "signature", "N_mutations")
    exposures_file$Total_mut = rowSums(exposures_file[,c(1:ncol(exposures_file)-1)])
    print(head(exposures_file))
    plot_exposures_hdp <- (ggplot(exposures_file_melted, aes(x=factor(sample, levels=samples_order), y=N_mutations, fill=signature))
        + geom_bar(stat="identity")
        + theme_bw()
        + theme(axis.text.x = element_blank())
        + xlab("")
        + scale_fill_manual(values = colors[1:n_sigs])
        + ggtitle("HDP signatures exposures")
    )
    return(plot_exposures_hdp)
}

plot_exposures_sigprofiler <- function(deepCSA_run){
    colors <-c("#005f73","#0a9396", "#94d2bd","#e9d8a6", "#ee9b00","#f15855", "#9b2226")
    all_mut_df <- read.table(paste("../..//supplementary/supp_note_5/data/Mutation_stats_per_genome_per_sample.med.all.csv", sep=""), sep=",", header =TRUE)
    samples_order = all_mut_df[order(all_mut_df$Total_mut_genome, decreasing=TRUE),]$SAMPLE_ID
    exposures_file = read.table(paste0("../../supplementary/supp_note_5/output/SigProfilerExtractor_out.med.all/SBS96/Suggested_Solution/SBS96_De-Novo_Solution/Activities/SBS96_De-Novo_Activities_refit.txt"), sep="\t", header=T)
    n_sigs = ncol(exposures_file)-1
    exposures_file_melted = melt(exposures_file)
    colnames(exposures_file_melted)<- c("sample", "signature", "N_mutations")
    exposures_file$Total_mut = rowSums(exposures_file[,c(2:ncol(exposures_file))])
    plot_exposures_sigprofiler <- (ggplot(exposures_file_melted, aes(x=factor(sample, levels=samples_order), y=N_mutations, fill=signature))
        + geom_bar(stat="identity")
        + theme_bw()
        + theme(axis.text.x = element_blank())
        + xlab("")
        + scale_fill_manual(values = colors[1:n_sigs])
        + ggtitle("SigProfiler signatures exposures")
    )
    return(plot_exposures_sigprofiler)
}


profile_hdp <- arrangeGrob(plot_sig_profiles(deepCSA_run, "HDP"), top = textGrob("C", x = unit(0, "npc"), y = unit(1, "npc"), just=c("left","top"), hjust=-1.5, vjust=2,gp = gpar(fontface = "bold")))
exposures_hdp <- arrangeGrob(plot_exposures_HDP(deepCSA_run), top = textGrob("D", x = unit(0, "npc"), hjust=-1.5, vjust=2,gp = gpar(fontface = "bold")))
profile_sigprofiler <- arrangeGrob(plot_sig_profiles(deepCSA_run, "SigProfiler"), top = textGrob("A", x = unit(0, "npc"), hjust=-1.5, vjust=2,gp = gpar(fontface = "bold")))
exposures_sigprofiler <- arrangeGrob(plot_exposures_sigprofiler(deepCSA_run), top = textGrob("B", x = unit(0, "npc"), hjust=-1.5, vjust=2,gp = gpar(fontface = "bold")))


hlay <- rbind(c(1,1,2,2),
                c(1,1,2,2),
                c(1,1,2,2),
                c(1,1,2,2),
                c(NA,NA,2,2),
                c(NA,NA,2,2),
                c(3,3,4,4),
                c(3,3,4,4),
                c(3,3,4,4)
                )

png(filename = outplot, width = 35, height = 25, units = "cm", res = 300)
grid.arrange(profile_sigprofiler,profile_hdp,exposures_sigprofiler,exposures_hdp,layout_matrix=hlay)
dev.off()



