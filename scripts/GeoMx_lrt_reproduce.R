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
library(writexl)

# directories
# path_to_input <- "/input/data/folder/"
# path_to_output <- "/ouptut/results/folder/"
path_to_input <- "/home/soterinosogo/2026.01.Reproduce_GeoMx/data/"
path_to_output <- "/home/soterinosogo/2026.01.Reproduce_GeoMx/results/"

# load data and metadata
data <- read.csv(file.path(path_to_input, "21019220628_rawdata.csv"), header = T, row.names = 1, sep = ";")
meta <- read.csv(file.path(path_to_input, "21019220628_metadata.csv"), header = T, row.names = 1)

# NOTE PanCK.Vimentin is only in Mixed tisse
# this tumor area was initially (histologically) included as an intermediate stage between tumor core and invasive front
# finally, it was removed from the presented analysis since we found it to be (transcriptomically) very similar to the tumor core
meta$tissue <- ifelse(meta$segment == "PanCK.Vimentin", "Invasive front", meta$tissue)
meta$tissue <- ifelse(meta$tissue == "Mixed", "Invasive", meta$tissue)
meta$tissue <- ifelse(meta$tissue == "in situ", "In situ", meta$tissue)
meta$segment_manipulated <- ifelse(meta$segment == "PanCK.Vimentin", "PanCK", meta$segment)

# format columns
data <- data[rownames(data) != "NegProbe-WTX",]
colnames(data) <- gsub("\\.", "-", colnames(data))
colnames(data) <- gsub("-dcc", "\\.dcc", colnames(data))
meta$segment_manipulated <- factor(meta$segment_manipulated, levels = unique(meta$segment_manipulated))
meta$Sample_ID <- rownames(meta)
meta$segment <- factor(meta$segment, levels = unique(meta$segment))
meta$DetectionThreshold <- factor(meta$DetectionThreshold, levels = c("1-5%", "5-10%", "10-15%", ">15%"))
meta$Patient <- factor(meta$Patient, levels = unique(meta$Patient))

# patient S2014485 was removed from the analysis, since it was processed in a different batch, which skewed the results from our analyses
meta <- meta[meta$Patient != "S2014485",]
# PanCK Invasive state was removed, but Macrophages Invasive state was kept, in the absence of an "Invasive front"
meta <- meta[!(meta$tissue == "Invasive" & meta$segment == "PanCK"), ]
# In situ stage is actually AK/In situ, since some cases posed difficulties to distinguish between the two
meta$tissue <- ifelse(meta$tissue == "In situ", "AK/In situ", meta$tissue)
# rename macrophage invasive
meta$tissue <- ifelse(meta$tissue == "Invasive", "Invasive front", meta$tissue)
# convert into factor
meta$tissue <- factor(meta$tissue, levels = c("Healthy", "Peritumoral", "AK/In situ", "Tumor", "Invasive front"))
meta$segment_manipulated <- factor(meta$segment_manipulated, levels = unique(meta$segment_manipulated))
data <- data[, colnames(data) %in% rownames(meta)]
identical(rownames(meta), colnames(data))

# create separate objects for macrophages vs tumor cells
list_data <- list()
segment_manipulated_levels <- levels(meta$segment_manipulated) %>% as.character()

for (i in segment_manipulated_levels){
  list_data[[i]] <- DESeqDataSetFromMatrix(countData = data[,colnames(data) %in% meta[meta$segment_manipulated == i,]$Sample_ID] %>% as.matrix(),
                                          colData = meta[meta$segment_manipulated == i,],
                                          design = ~ tissue)
}

# transform counts for data visualization
list_rlog <- list()

for (i in segment_manipulated_levels){
  # perform rlog transformation
  list_rlog[[i]] <- rlog(list_data[[i]], blind=FALSE)
  # Extract the rlog matrix from the object
  list_rlog[[i]] <- assay(list_rlog[[i]])
}

# perform DE using LRT
list_sig_results <- list()

for (i in segment_manipulated_levels){
  list_data[[i]] <- DESeq(list_data[[i]], test="LRT", reduced = ~ 1)
  # extract results
  list_sig_results[[i]] <- results(list_data[[i]])
  # subset the LRT results to return genes with padj < 0.01
  list_sig_results[[i]] <- list_sig_results[[i]] %>%
    data.frame() %>%
    rownames_to_column(var="gene") %>% 
    as_tibble() %>% 
    filter(padj < 0.01)
}


# bbtain rlog values for those significant genes and perform clustering for DE
list_clustering <- list()

for (i in segment_manipulated_levels){
  # extract the clustering data for Macrophages and PanCK
  list_clustering[[i]] <- list_rlog[[i]][list_sig_results[[i]]$gene, ]
  # use the `degPatterns` function from the 'DEGreport' package to show gene clusters across sample groups
  list_clustering[[i]] <- degPatterns(list_clustering[[i]], metadata = meta[meta$segment_manipulated == i,], time = "tissue", col = NULL)
  
  # reorder the plot to number groups correctly and adapt aesthetics
  plot_data <- list_clustering[[i]]$plot$data
  
  if (i == "PanCK"){
    # reorder the group numbers and skip the missing ones (adapted to PanCK case)
    plot_data <- plot_data %>%
      mutate(
        group_num = as.numeric(str_extract(title, "\\d+")),
        new_group = case_when(
          group_num >= 12 & group_num <= 18 ~ group_num - 1,
          group_num == 20                       ~ group_num - 2,
          TRUE                                  ~ group_num
        ),
        titles_corrected = str_replace(title, "Group:\\s*\\d+", paste0("Ker-G", new_group))
      )
    plot_data$title <- plot_data$titles_corrected
  } else if (i == "Macrophages") {
    # reorder the group numbers and skip the missing ones (adapted to Macrophages case)
    plot_data <- plot_data %>%
      mutate(
        group_num = as.numeric(str_extract(title, "\\d+")),
        new_group = case_when(
          group_num >= 2 & group_num <= 5 ~ group_num - 1,
          group_num == 6                        ~ group_num - 4,
          TRUE                                  ~ group_num
        ),
        titles_corrected = str_replace(title, "Group:\\s*\\d+", paste0("Mac-G", new_group))
      )
    plot_data$title = plot_data$titles_corrected
  }
  
  # extract group numbers from the 'title' column using regex
  plot_data$group_number <- as.numeric(gsub(".*:\\s*(\\d+)\\s*-.*", "\\1", plot_data$title))
  # reorder the 'title' column based on the numeric group number
  plot_data$title <- factor(plot_data$title, levels = plot_data$title[order(plot_data$group_number)] %>% unique)
  # remove the temporary 'group_number' column
  plot_data$group_number <- NULL
  
  # update the plot data in the object
  list_clustering[[i]]$plot$data <- plot_data
  
  # adapt color coding
  modified_plot <- list_clustering[[i]]$plot +
    aes(col = tissue) +
    scale_color_manual(values=c("#FDE725", "#5DC863", "#21908C", "#3B528B", "#440154"))
  
  # change font size
  if (i == "PanCK"){
    modified_plot <- modified_plot +
    theme(
      strip.text = element_text(size = 12)
    )
  } else if (i == "Macrophages") {
        modified_plot <- modified_plot +
    theme(
      strip.text = element_text(size = 24)
    )
  }

  # include trend line
  modified_plot$layers[[3]] <- ggplot2::geom_smooth(
    mapping = ggplot2::aes(x = tissue, y = value, group = 1),
    method = "loess", 
    color = NA, 
    se = FALSE, 
    size = 1.5
  )
  modified_plot + ggplot2::geom_smooth(
    mapping = ggplot2::aes(x = tissue, y = value, group = 1),
    method = "loess", 
    color = "black", 
    se = FALSE, 
    size = 1.5)

  # save plot
  plot_name <- paste0("Publication_LRT_DEGclustering_", i,".png")
  ggsave(file.path(path_to_output, plot_name),
         device = "png", units = "px", width = 1300 *6, height = 1000 *6, dpi = 600)
}

# save gene lists to csv file
for (i in segment_manipulated_levels){
  # extract unique titles (groups) and associated genes
  gene_group_data <- list_clustering[[i]]$plot$data
  
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
  file_name <- paste0(file = "Publication_", i, "_gene_groups.csv")
  write.csv(gene_df, file.path(path_to_output, file_name), row.names = FALSE)
  
}

##############################################
############## PATHWAY ANALYSES ##############
##############################################

# pathway analyses was performed via gprofiler and ORSUM summarization
# see /Pathway_analyses folder

#################################################
########### PLOT INDIVIDUAL GENES ###############
#################################################

# color consistency between GeoMx and HTG-seq
# shared between GeoMx and HTG-seq
# "S2210105" --> typo in patient metadata, correct label is S2202105
# "S2210320"
# "S2210959"
# 
# only HTG-seq
# "S2203938"
# "S2210388"
# "S2202680"
# "S2228101"
# "S2216712"

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

# introduce patient publication name for metadata
meta <- meta %>%
  left_join(patient_df, by = c("Patient" = "patient_id"))

# create separate objects for macrophages vs tumor cells
list_data <- list()
segment_manipulated_levels <- levels(meta$segment_manipulated) %>% as.character()

for (i in segment_manipulated_levels){
  list_data[[i]] <- DESeqDataSetFromMatrix(countData = data[,colnames(data) %in% meta[meta$segment_manipulated == i,]$Sample_ID] %>% as.matrix(),
                                          colData = meta[meta$segment_manipulated == i,],
                                          design = ~ tissue)
}

# create directory for plots
dir.create(file.path(path_to_output, "Single_gene_plots"), showWarnings = FALSE)

# plot the genes for Macrophages and PanCK
genes <- c("AQP3", "IFI6", "UBE2L6", "TYMP", "PLEK2", "LAMB3", "VIM")

# build named vector: names = patient_pub_id, values = hex colors
color_map <- setNames(patient_df$color, patient_df$patient_pub_id)

groups <- names(list_data)
for (group in groups){
  for (gene in genes){
    # generate plots
    counts <- plotCounts(list_data[[group]], gene = gene, intgroup=c("tissue", "patient_pub_id"), normalized = TRUE, returnData=TRUE)
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
    file_name = paste0("GeoMx_", gene, "_", group, ".png")
    ggsave(file.path(path_to_output, "Single_gene_plots", file_name),
           device = "png", units = "px", width = 750 *4, height = 450 *4, dpi = 400)
  }
}

#################################################
########## CREATE HEAMAPS GeoMx PANCK ###########
#################################################

# split matrix for macrophages vs panck
panck <- DESeq(list_data[["PanCK"]])
panck_mat <- counts(panck, normalized = T)
macrophages <- DESeq(list_data[["Macrophages"]])
macrophages_mat <- counts(macrophages, normalized = T)

# load file
heatmap_lists <- readxl::read_excel(file.path(path_to_input, "Heatmaps_gene_list.xlsx"))
colnames(heatmap_lists) <- gsub(" - ", "-", colnames(heatmap_lists))

# annotation for columns
meta_panck <- meta[meta$segment_manipulated == "PanCK",]

# reorder matrix and metadata
meta_panck <- meta_panck[order(meta_panck$tissue),]
panck_mat <- panck_mat[, match(meta_panck$Sample_ID, colnames(panck_mat))]
# since ids in meta and expression are in the same order, annotations can be made straight from meta_panck
identical(meta_panck$Sample_ID, colnames(panck_mat))

# build named vector: names = patient_pub_id, values = hex colors
color_map <- setNames(patient_df$color, patient_df$patient_pub_id)
# build named vector: names = Tissue, values = hex colors
color_map_stages <- setNames(c("#FDE725", "#5DC863", "#21908C", "#3B528B", "#440154"), levels(meta$tissue))

# create heatmap annotations
ha <- HeatmapAnnotation(
  Patient = meta_panck$patient_pub_id,
  Tissue = meta_panck$tissue,
  col = list(
    Patient = color_map,
    Tissue = color_map_stages),
  na_col = "white"
)

# create directory for plots
dir.create(file.path(path_to_output, "Heatmap_plots"), showWarnings = FALSE)

# prepare prefixes for groups of genes
group_prefix <- c("Group 1-", "Group 3-", "Group 2-", "Group 4-", "Group 6-", "Group 10-",  "Group B-")

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
  expr_subset <- panck_mat[genes_in_order, , drop = FALSE]
  
  # factor row split levels for proper order
  row_split <- factor(gene_long$Pathway, levels = unique(gene_long$Pathway))
  
  # heatmap plot to extract gene order
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
  png(file.path(path_to_output, "Heatmap_plots", paste0("GeoMx_PanCK_heatmap_", group_prefix[i], ".png")), width = 4* 900, height = (4* (500 + 3.6 * scale_factor)), units = "px", res = 400)
  draw(ht_plot)
  dev.off()
  
  # extract row order
  ht <- draw(ht_plot)
  list_genes <- row_order(ht)
  # slice heatmap in pathways
  split_levels <- levels(row_split)

  # exception for single slice groups
  if (!is.list(list_genes)) {
    # if single slice, convert vector to named list
    list_genes <- setNames(list(list_genes), split_levels[1])
  }
  if (is.null(names(list_genes))) {
    names(list_genes) <- split_levels[seq_along(list_genes)]
  }

  # build textbox annotation
  text_ann <- vector("list", length(split_levels))
  names(text_ann) <- split_levels

  for (lvl in names(list_genes)) {
    idx <- list_genes[[lvl]]

    text_ann[[lvl]] <- data.frame(
      text = rownames(expr_subset)[idx],
      col = "black",
      fontsize = 12
    )
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
  png(file.path(path_to_output, "Heatmap_plots", paste0("GeoMx_PanCK_heatmap_", group_prefix[i], "_w_genes.png")), width = 4* 1500, height = (4* (500 + 3.6 * scale_factor)), units = "px", res = 400)
  draw(ht_plot)
  dev.off()
}


#######################################################
########## CREATE HEAMAPS GeoMx Macrophages ###########
#######################################################

# load file
heatmap_lists <- readxl::read_excel(file.path(path_to_input, "Heatmaps_gene_list.xlsx"), sheet = 2)
colnames(heatmap_lists) <- gsub(" - ", "-", colnames(heatmap_lists))

# annotation for columns
meta_macro <- meta[meta$segment_manipulated == "Macrophages",]

# reorder matrix and metadata
meta_macro <- meta_macro[order(meta_macro$tissue),]
macrophages_mat <- macrophages_mat[, match(meta_macro$Sample_ID, colnames(macrophages_mat))]
# since ids in meta and expression are in the same order, annotations can be made straight from meta_macro
identical(meta_macro$Sample_ID, colnames(macrophages_mat))

# build named vector: names = patient_pub_id, values = hex colors
color_map <- setNames(patient_df$color, patient_df$patient_pub_id)
# build named vector: names = Tissue, values = hex colors
color_map_stages <- setNames(c("#FDE725", "#5DC863", "#21908C", "#3B528B", "#440154"), levels(meta$tissue))

# create heatmap annotations
ha <- HeatmapAnnotation(
  Patient = meta_macro$patient_pub_id,
  Tissue = meta_macro$tissue,
  col = list(
    Patient = color_map,
    Tissue = color_map_stages),
  na_col = "white"
)

# prepare prefixes for groups of genes
group_prefix <- c("Group 1-")

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
  expr_subset <- macrophages_mat[genes_in_order, , drop = FALSE]
  
  # factor row split levels for proper order
  row_split <- factor(gene_long$Pathway, levels = unique(gene_long$Pathway))
  
  # heatmap plot to extract gene order
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
  png(file.path(path_to_output, "Heatmap_plots", paste0("GeoMx_Macrophages_heatmap_", group_prefix[i], ".png")), width = 4* 900, height = (4* (500 + 3.6 * scale_factor)), units = "px", res = 400)
  draw(ht_plot)
  dev.off()
  
  # extract row order
  ht <- draw(ht_plot)
  list_genes <- row_order(ht)
  # slice heatmap in pathways
  split_levels <- levels(row_split)

  # exception for single slice groups
  if (!is.list(list_genes)) {
    # if single slice, convert vector to named list
    list_genes <- setNames(list(list_genes), split_levels[1])
  }
  if (is.null(names(list_genes))) {
    names(list_genes) <- split_levels[seq_along(list_genes)]
  }

  # build textbox annotation
  text_ann <- vector("list", length(split_levels))
  names(text_ann) <- split_levels

  for (lvl in names(list_genes)) {
    idx <- list_genes[[lvl]]

    text_ann[[lvl]] <- data.frame(
      text = rownames(expr_subset)[idx],
      col = "black",
      fontsize = 12
    )
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
  png(file.path(path_to_output, "Heatmap_plots", paste0("GeoMx_Macrophages_heatmap_", group_prefix[i], "_w_genes.png")), width = 4* 1500, height = (4* (500 + 3.6 * scale_factor)), units = "px", res = 400)
  draw(ht_plot)
  dev.off()
}


######################################################
########## DIFFERENTIAL EXPRESSION MARKERS ###########
######################################################

# perform DE for each of the groups
list_sig_results = list("Macrophages" = list(),
            "PanCK" = list())
# store PanCK and Macrophages in different list items
list_data = list()

for (i in segment_manipulated_levels){
 # define design based on the specific level of the tissue variable
 # extract levels of tissue variable to create an alternative variable
 groups = meta[meta$segment_manipulated == i,]$tissue %>% unique()
 for (group in groups){
  # create a temporary variable to loop through
  meta$tissue_alternative = ifelse(meta$tissue == group, group, "rest")
  meta$tissue_alternative = factor(meta$tissue_alternative, levels = c(group, "rest"))
  # create DESeq object
  list_data[[i]] = DESeqDataSetFromMatrix(countData = data[,colnames(data) %in% meta[meta$segment_manipulated == i,]$Sample_ID] %>% as.matrix(),
                      colData = meta[meta$segment_manipulated == i,],
                      design = ~ tissue_alternative)
  # run DESeq2 DE testing
  list_data[[i]] = DESeq(list_data[[i]])
  # store results in list
  list_sig_results[[i]][[group]] = results(list_data[[i]], contrast = c("tissue_alternative", group, "rest"))
  # turn gene name into column
  list_sig_results[[i]][[group]]$gene = rownames(list_sig_results[[i]][[group]])
  # turn into dataframe
  list_sig_results[[i]][[group]] = list_sig_results[[i]][[group]] %>% as.data.frame()
  # sort in decreasing order
  list_sig_results[[i]][[group]] = list_sig_results[[i]][[group]][order(list_sig_results[[i]][[group]]$log2FoldChange, decreasing = T),]
 }
}

# create dir for results
dir.create(file.path(path_to_output, "DE_genes_by_stage"), showWarnings = FALSE)
# remove "/" from "AK/In situ" to be able to save the xlsx files
names(list_sig_results[["Macrophages"]]) <- gsub("/", "_", groups)
names(list_sig_results[["PanCK"]]) <- gsub("/", "_", groups)

# write results into a single Excel document
write_xlsx(list_sig_results[["Macrophages"]],file.path(path_to_output, "DE_genes_by_stage", "Macrophages_DE_stage.xlsx"))
write_xlsx(list_sig_results[["PanCK"]], file.path(path_to_output, "DE_genes_by_stage", "PanCK_DE_stage.xlsx"))