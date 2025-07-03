library(ggplot2)
library(dplyr)
library(lsa)
library(gridExtra)

# Load mutation_counts from LCMs Lawson et al.
samples_lawson = read.csv("./data/lawson_et_al_data/Lawson_bladder_processed_muts.tsv", header=T)
lawson_data = unique(samples_lawson[,c("sampleID", "type")])
genomes_lawson <- lawson_data[lawson_data$type == "wgs", ]$sampleID
mutation_counts_lawson <- read.table("./data/lawson_et_al_data/Lawson_bladder_mutation_matrix.tsv", header = T, sep="\t")
# select only wgs samples
mutation_counts_lawson <- mutation_counts_lawson[,genomes_lawson] 
mutation_counts_lawson <- mutation_counts_lawson[,c(2:ncol(mutation_counts_lawson))]
print(mutation_counts_lawson[1:5,1:5])

# Load duplex data
deepCSA_run_dir <- "/data/bbg/nobackup/bladder_ts/results/2025-05-14_deepCSA_45_donors"
output_dir <- paste0("plots/")
df <- as.data.frame(t(read.table(paste(deepCSA_run_dir, "/signatures_hdp/input/samples_matrix.all.hdp.csv", sep=""), sep=" ", header=TRUE, row.names=1)))
rownames(df) <- paste0(substr(rownames(df),1,1),"[", substr(rownames(df),3,3), ">", substr(rownames(df),5,5), "]",substr(rownames(df),7,7))
print(head(df[,1:5]))

#merge duplex data and lawson data
df_full <- merge(mutation_counts_lawson,df,by=0)
rownames(df_full) = df_full$Row.names
df_full <-  df_full[,!(names(df_full) %in% "Row.names")]
df_full <- as.data.frame(t(df_full))
df_upd <- df_full %>% mutate(dataset=case_when(grepl("^P19", rownames(df_full)) ~ "duplex panel all mutations", !grepl("^P19", rownames(df_full)) ~ "WGS LCMs"))

sigs <- df_upd %>%
    group_by(dataset) %>%
    summarise_all(sum)

sigs <- as.data.frame(sigs)
print(head(sigs))
sigs_freq <- sigs %>% ungroup() %>%  mutate(across(where(is.numeric))/rowSums(across(where(is.numeric))))
print(sigs_freq[1:3,1:5])

print("cosine lcm vs all panel")
cosine_lcm_all <- (cosine(as.numeric((sigs_freq[1,2:ncol(sigs_freq)])), as.numeric((sigs_freq[2,2:ncol(sigs_freq)]))))
print(cosine_lcm_all)


# plot components
# SBS #
sigs_df <- reshape2::melt(sigs_freq)
colnames(sigs_df) <- c('component', 'channel', 'value')
print(head(sigs_df))
sigs_df$group <- paste0(substr(as.character(sigs_df$channel), start = 3, stop = 3), ">", substr(as.character(sigs_df$channel), start = 5, stop = 5))
sigs_df <- sigs_df %>%
    mutate(channel = factor(channel, levels =  colnames(sigs)),
    group = factor(group, levels = unique(sigs_df$group)),
    value = value * 100)
  
#set colours
colours <- setNames(c('#1ebff0', '#050708', '#e62725', '#cbcacb', '#a1cf64', '#edc8c5'), unique(sigs_df$group))
strip_name_colours <- c('black','white','white','black','black','black')
####

lcm_profile <- sigs_df[sigs_df$component == "WGS LCMs",]$value

plots <- list()

for (x in unique(sigs_df$component)) {
  plot_data <- sigs_df[sigs_df$component == x,]
  plot_data$xlabels <- paste0(substr(as.character(plot_data$channel), 1, 1),
                              substr(as.character(plot_data$channel), 3, 3),
                              substr(as.character(plot_data$channel), 7, 7))
  sample_profile = plot_data$value
  cosine_val=cosine(lcm_profile, sample_profile)
  if (x=="WGS LCMs"){
    title_text = x
  }else{
    title_text = paste(x, "\nCosine similarity: ", round(cosine_val,2), sep="")
  }
  
  p <- ggplot(plot_data, aes(x = xlabels, y = value, fill = group)) +
    geom_bar(stat = 'identity') + 
    facet_grid(. ~ group, space = 'free_x', scales = 'free_x') +
    scale_fill_manual(name = '', values = colours, guide = 'none') +
    xlab('') +
    ylab('% SBS') +
    scale_y_continuous(expand = c(0,0,0.05,0)) +
    ggtitle(title_text) +
    theme_classic() + 
    theme(plot.title = element_text(hjust = 1, vjust=-10, size = 12,margin=margin(t=25,b=-20)), 
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

# Save the combined plot as a PNG file
png(filename = paste0(output_dir, 'ExtendedFig2b_LCMs_lawson_vs_duplex_all_profiles.png'), width = 25, height = 10, units = "cm", res = 300)
grid.arrange(grobs = plots, ncol = n_col, nrow = n_row)
dev.off()

