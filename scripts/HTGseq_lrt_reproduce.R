library(GSVA)
library(dplyr)
library(SummarizedExperiment)
library(DESeq2)
library(DEGreport)
library(tidyverse)
library(msigdbr)
library(fgsea)
library(data.table)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(ggplot2)
library(DOSE)
library(org.Hs.eg.db)
library(GGally)

# directories
# path_to_input <- "/input/data/folder/"
# path_to_output <- "/ouptut/results/folder/"
path_to_input <- "/home/soterinosogo/2026.01.Reproduce_GeoMx/data/"
path_to_output <- "/home/soterinosogo/2026.01.Reproduce_GeoMx/results/"

# build color coding dataframe
patient_id <- c("S2210959", "S2210105", "S2210320",
               "S2210388", "S2203938", "S2202105", "S2216712", "S2228101", "S2202680")
patient_pub_id <- c("Patient 1", "Patient 2", "Patient 3", 
                   "Patient 4", "Patient 5", "Patient 2", "Patient 6", "Patient 7", "Patient 8")
# #000000, #E69F00, #56B4E9, #009E73, #F0E442, #0072B2, #D55E00, #CC79A7
color <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                    "#F0E442", "#E69F00","#0072B2", "#D55E00", "#CC79A7")
                   
patient_df <- data.frame("patient_id" = patient_id,
                        "patient_pub_id" = patient_pub_id,
                        "color" = color)

# load data
data <- readxl::read_excel(file.path(path_to_input, "21019221207_HTGseq_rawdata.xlsx"), skip = 9, sheet = 2) %>%
  as.data.frame()
# remove total counts
data <- data[2:dim(data)[1],]
rownames(data) <- data[,1]
data[,1] <- NULL

# create metadata fromm data names
# create vector of column names from the count matrix
colnames_data <- colnames(data)

# function to classify tissues based on suffix
classify_tissue <- function(suffix) {
  if (suffix == "SAIN") {
    return("Healthy")
  } else if (suffix == "TUM") {
    return("Tumor")
  } else if (suffix == "PRE") {
    return("AK/In situ")
  } else {
    return(NA)
  }
}

# extract patient ID and tissue classification
meta <- data.frame(
  Patient = sapply(strsplit(colnames_data, "-"), function(x) x[2]),
  tissue = sapply(strsplit(colnames_data, "-"), function(x) classify_tissue(x[3])),
  row.names = colnames_data
)
meta$tissue <- factor(meta$tissue, levels = c("Healthy", "AK/In situ", "Tumor"))

# introduce patient publication name for metadata
meta <- meta %>%
  left_join(patient_df, by = c("Patient" = "patient_id"))
rownames(meta) <- colnames_data

# create dds object
dds <- DESeqDataSetFromMatrix(countData = data,
                             colData = meta,
                             design = ~ tissue)

# transform counts for data visualization
dds_rlog <- rlog(dds, blind=FALSE)
dds_rlog <- assay(dds_rlog)


# perform DE using LRT
dds <- DESeq(dds, test="LRT", reduced = ~ 1)
# extract results
dds_results <- results(dds)
# subset the LRT results to return genes with padj < 0.01
dds_results <- dds_results %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(padj < 0.01)


# obtain rlog values for those significant genes and perform clustering for DE
dds_clustering <- dds_rlog[dds_results$gene, ]
# use the `degPatterns` function from the 'DEGreport' package to show gene clusters across sample groups
dds_clustering <- degPatterns(dds_clustering, metadata = meta, time = "tissue", col = NULL)

# reorder the plot to match group order with layout
# extract the data from the plot object
plot_data <- dds_clustering$plot$data
# extract group numbers from the 'title' column using regex
plot_data$group_number <- as.numeric(gsub("Group:\\s*(\\d+)\\s*-.*", "\\1", plot_data$title))
# reorder the 'title' column based on the numeric group number
plot_data$title <- factor(plot_data$title, levels = plot_data$title[order(plot_data$group_number)] %>% unique)
# remove the temporary 'group_number' column
plot_data$group_number <- NULL
# update the plot data in the object
dds_clustering$plot$data <- plot_data

# save a more visual representation of the plot
modified_plot <- dds_clustering$plot +
  aes(col = tissue) +
  scale_color_manual(values=c("#FDE725", "#21908C", "#3B528B"))

# remove previous line (only way to make it work)
modified_plot$layers[[3]] <- ggplot2::geom_smooth(
  mapping = ggplot2::aes(x = tissue, y = value, group = 1),
  method = "loess", 
  color = NA, 
  se = FALSE, 
  size = 1.5
)
# include black line
modified_plot + ggplot2::geom_smooth(
  mapping = ggplot2::aes(x = tissue, y = value, group = 1),
  method = "loess", 
  color = "black", 
  se = FALSE, 
  size = 1.5)
ggsave(file.path(path_to_output, "Publication_LRT_DEGclustering_HTGEdgeSeq.png"),
       device = "png", units = "px", width = 750 *6, height = 1000 *6, dpi = 600)


# save gene lists to csv file
# extract unique titles (groups) and associated genes
gene_group_data <- dds_clustering$plot$data

# create a list of genes for each group
grouped_genes <- split(gene_group_data$genes %>% unique(), gene_group_data$title)

# convert the list into a data frame with each group as a column
# make the list elements equal in length by padding with NA
max_length <- max(sapply(grouped_genes, length))
gene_df <- as.data.frame(do.call(cbind, lapply(grouped_genes, function(x) {
  c(x, rep("", max_length - length(x))) # Pad with NA
})))

# set column names as group names
colnames(gene_df) <- names(grouped_genes)

# write the data frame to a CSV file
file_name <- "Publication_HTG_Edgeseq_gene_groups.csv"
write.csv(gene_df, file.path(path_to_output, file_name), row.names = FALSE)


######################################################
########## DIFFERENTIAL EXPRESSION MARKERS ###########
######################################################

# store results in different list items
list_sig_results = list()

for (stage in unique(meta$tissue)){
  # create a temporary variable to loop through
  meta$tissue_alternative = ifelse(meta$tissue == stage, stage, "rest")
  meta$tissue_alternative = factor(meta$tissue_alternative, levels = c(stage, "rest"))
  # create DESeq object
  dds_stage = DESeqDataSetFromMatrix(countData = data %>% as.matrix(),
                                          colData = meta,
                                          design = ~ tissue_alternative)
  # run DESeq2 DE testing
  dds_stage = DESeq(dds_stage)
  # store results in list
  list_sig_results[[stage]] = results(dds_stage, contrast = c("tissue_alternative", stage, "rest"))
  # turn gene name into column
  list_sig_results[[stage]]$gene = rownames(list_sig_results[[stage]])
  # turn into dataframe
  list_sig_results[[stage]] = list_sig_results[[stage]] %>% as.data.frame()
  # subset padj < 0.05
  list_sig_results[[stage]] = list_sig_results[[stage]][list_sig_results[[stage]]$padj < 0.05,]
  # sort in decreasing order
  list_sig_results[[stage]] = list_sig_results[[stage]][order(list_sig_results[[stage]]$log2FoldChange, decreasing = T),]
  print(dim(list_sig_results[[stage]]))
}

# create dir for results
dir.create(file.path(path_to_output, "DE_genes_by_stage"), showWarnings = FALSE)

# combine results
combined_results <- dplyr::bind_rows(
  lapply(names(list_sig_results), function(stage) {
    df <- list_sig_results[[stage]]
    df$stage <- stage
    df
  })
)

# write results into a single Excel document
write.csv(combined_results, file.path(path_to_output, "DE_genes_by_stage", "HTGseq_DE_stage.csv"))

##############################################
############## PATHWAY ANALYSES ##############
##############################################

# pathway analyses was performed via gprofiler and ORSUM summarization
# see /Pathway_analyses folder


#################################################
########### PLOT INDIVIDUAL GENES ###############
#################################################

# plot the genes for Macrophages and PanCK
genes <- c("MYBL2", "PLOD2", "CD36", "SERPINA12", "SERPINB12")

# build named vector: names = patient_pub_id, values = hex colors
color_map <- setNames(patient_df$color, patient_df$patient_pub_id)
color_map_stages <- setNames(c("#FDE725", "#21908C", "#440154"), levels(meta$tissue))

# create directory for plots
dir.create(file.path(path_to_output, "Single_gene_plots"), showWarnings = FALSE)

for (gene in genes){
  # generate plots
  counts <- plotCounts(dds, gene = gene, intgroup=c("tissue", "patient_pub_id"), normalized = TRUE, returnData=TRUE)
  ggplot(counts, aes(x=tissue, y=count, color =patient_pub_id)) +
    geom_jitter(size = 2, position = position_jitter(width = 0.3)) +
    scale_color_manual(values = color_map) +
    theme_bw() +
    theme(text = element_text(size = 20),
          axis.text.x = element_text(angle=45, hjust=1),
          legend.title = element_blank()) +
    ggtitle(gene) +
    xlab("Tissue of origin") +
    ylab("Normalized counts")
  
  # save plot
  file_name <- paste0("HTGseq_", gene, ".png")
  ggsave(file.path(path_to_output, "Single_gene_plots", file_name),
         device = "png", units = "px", width = 750 *4, height = 450 *4, dpi = 400)
}

######################################
####### CREATE HEAMAPS HTGseq  #######
######################################

# load file
heatmap_lists <- readxl::read_excel(file.path(path_to_input, "Heatmaps_gene_list.xlsx"), sheet = 3)
colnames(heatmap_lists) <- gsub(" - ", "-", colnames(heatmap_lists))

# extract counts
htgseq_mat <- counts(dds, normalized = T)

# prepare prefixes for groups of genes
group_prefix <- c("Group 1-", "Group 2-")

# reorder matrix and metadata
meta <- meta[order(meta$tissue),]
htgseq_mat <- htgseq_mat[, match(rownames(meta), colnames(htgseq_mat))]
# since ids in meta and expression are in the same order, annotations can be made straight from meta_macro
identical(rownames(meta), colnames(htgseq_mat))

# create heatmap annotations
ha <- HeatmapAnnotation(
  Patient = meta$patient_pub_id,
  Tissue = meta$tissue,
  col = list(
    Patient = color_map,
    Tissue = color_map_stages),
  na_col = "white"
)

# create directory for plots
dir.create(file.path(path_to_output, "Heatmap_plots"), showWarnings = FALSE)

# plot heatmaps
for (i in 1:length(group_prefix)){
  col_subset <- startsWith(colnames(heatmap_lists), group_prefix[i])
  tmp_data <- heatmap_lists[,col_subset]
  colnames(tmp_data) <- gsub(group_prefix[i], "", colnames(tmp_data))
  # melt the gene groups data frame into long format
  gene_long <- tmp_data %>%
    pivot_longer(everything(), names_to = "Pathway", values_to = "Gene") %>%
    filter(!is.na(Gene))
  
  # subset exp mat to genes present in pathways
  gene_long <- gene_long[order(gene_long$Pathway),]
  genes_in_order <- gene_long$Gene
  expr_subset <- htgseq_mat[genes_in_order, , drop = FALSE]
  
  # factor row split levels for proper order
  row_split <- factor(gene_long$Pathway, levels = unique(gene_long$Pathway))
  
  # plot to extract gene order
  ht_plot <- Heatmap(
    t(scale(t(expr_subset))),
    name = "Expression",
    top_annotation = ha,
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    split = row_split,
    row_title_gp = gpar(fontsize = 12),
    row_title_rot = 0,
    show_row_names = FALSE,
    show_column_names = FALSE
  )
  
  # get parameters to adjust heatmaps height based on gene length
  scale_factor <- length(genes_in_order)
  print(3.6 * scale_factor)
  
  # save heatmap
  png(file.path(path_to_output, "Heatmap_plots", paste0("HTGseq_heatmap_", group_prefix[i], ".png")), width = 4* 1100, height = (4* (500 + 3.6 * scale_factor)), units = "px", res = 400)
  draw(ht_plot)
  dev.off()
  
  # extract row order
  ht <- draw(ht_plot)
  list_genes <- row_order(ht)
  
  # generate gene annotation for the text box
  text <- split(gene_long, gene_long$Pathway)
  class(text)
  text_ann <- list()
  for (j in 1:length(list_genes)){
    text_ann[[j]] <- data.frame("text" = rownames(expr_subset)[list_genes[[j]]],
                               "col" = "black",
                               "fontsize" = 12)
    names(text_ann)[j] <- names(list_genes)[j]
  }
  
  # plot with gene order next to the pathways
  ht_plot <- Heatmap(
    t(scale(t(expr_subset))),
    name = "Expression",
    top_annotation = ha,
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    split = row_split,
    row_title_gp = gpar(fontsize = 12),
    right_annotation = rowAnnotation(textbox = anno_textbox(row_split, text_ann)),
    row_title_rot = 0,
    show_row_names = FALSE,
    show_column_names = FALSE
  )
  
  # save heatmap
  png(file.path(path_to_output, "Heatmap_plots", paste0("HTGseq_heatmap_", group_prefix[i], "_w_genes.png")), width = 4* 1500, height = (4* (500 + 3.6 * scale_factor)), units = "px", res = 400)
  draw(ht_plot)
  dev.off()
}
