# IMPORTANT NOTES: 
# see lines 73-77 to see GMT and organism dependent parameters
# the script will loop through input file headers and generate an output for each of them
# output dirs will be named after input file header names
# see gost() documentation for supported gene input formats https://www.rdocumentation.org/packages/gprofiler2/versions/0.2.3/topics/gost

#############################
####### LOAD ARGUMENTS ######
#############################

# parse arguments from bash console to R language
parseArgs <- function(defaults = list()) {
  args <- commandArgs(trailingOnly = TRUE)
  parsed_args <- list()
  
  for (arg in args) {
    if (grepl("^--", arg)) {
      key_value <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
      parsed_args[[key_value[1]]] <- key_value[2]
    }
  }
  
  # Combine defaults with parsed arguments, giving priority to parsed ones
  args_combined <- modifyList(defaults, parsed_args)
  return(args_combined)
}

# default values for arguments
defaults <- list(
  orsum_dir = NULL,
  gmt_dir = NULL,
  output_dir = NULL,
  input_filepath = NULL,
  ordered_query = FALSE,
  orsum_term_size_limit = 'none',
  organism = 'hsapiens'
)

# load parsed arguments
args <- parseArgs(defaults)

# load individual arguments into the environment
orsum_dir <- args$orsum_dir
gmt_dir <- args$gmt_dir
output_dir <- args$output_dir
input_filepath <- args$input_filepath
ordered_query <- as.logical(gsub(" ", "", args$ordered_query))
orsum_term_size_limit <- args$orsum_term_size_limit
organism <- args$organism

# ensure parameters have been introduced correctly
if (is.null(gmt_dir) || is.null(output_dir) || is.null(input_filepath)) {
  stop("Error: --gmt_dir, --output_dir and --input_filepath are required arguments.")
}

# ensure that ordered_query is logical
if (class(ordered_query) != "logical") {
  stop("Error: --ordered_query is not logical.")
}

# ensure that orsum_term_size_limit is a number if inputed
if (orsum_term_size_limit != "none") {
  orsum_term_size_limit <- as.numeric(orsum_term_size_limit)
  if (class(orsum_term_size_limit) != "numeric") {
      stop("Error: --orsum_term_size_limit is not 'none' nor numeric.")
  }
}

# define repositories to run gprofiler on
sources <- c("GO:BP", "GO:MF", "GO:CC", "REAC", "KEGG", "TF", "MIRNA", "CORUM", "HP", "WP")
# Parameters for ORSUM
patterns <- c("GO_BP", "GO_MF", "GO_CC", "MIRNA", "WP", "REAC", "CORUM", "HP")
ORSUM_output_dir <- file.path(output_dir, "ORSUM_results")
# create output directory for ORSUM script
dir.create(ORSUM_output_dir)

###############################
####### CUSTOM FUNCTIONS ######
###############################

# load libraries
library(dplyr)
library(ggplot2)
library(openxlsx)
library(gprofiler2)
library(glue)
library(gtools)

# Function 1/3
# saving results into Excel and tsv, and save .txt column used as input for ORSUM
# Excel files are formatted specifically to show a barplot with pathway enrichment
export_gprofiler_excel <- function(gost.res, sample_name, output_dir, build_report = FALSE) {
  # Ensure output directories exist
  enrichment_dir <- file.path(output_dir, sample_name)
  orsum_dir <- file.path(output_dir, "ORSUM")
  dir.create(enrichment_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(orsum_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Iterate over unique sources in the data
  for (source in unique(gost.res$source)) {
    # Filter data for the current source
    source_data <- gost.res[gost.res$source == source, ]
    filtered_data <- source_data[, c(
      "source", "term_name", "term_id", "p_value", 
      "negative_log10_of_adjusted_p_value", "term_size", 
      "query_size", "intersection_size", "intersection"
    )]
    
    # IMPORTANT NOTE:
    # gprofiler reports directly padj values in the output when specifying correction method
    # however, column name still read "p_value", as if it was uncorreted
    # thus, column name is changed to avoid result misinterpretation
    colnames(filtered_data) <- gsub("p_value", "padj", colnames(filtered_data))
    
    # Create and save the tsv file
    # Save results as TSV
    tsv_file <- file.path(enrichment_dir, paste0(
      "gprofiler2_", gsub(":", "_", source), ".tsv"
    ))
    write.table(filtered_data, file = tsv_file, row.names = FALSE, sep = "\t", quote = FALSE)

    # Create and save the Excel file
    excel_file <- file.path(enrichment_dir, paste0(
      "gprofiler2_", gsub(":", "_", source), ".xlsx"
    ))
    wb <- createWorkbook()
    addWorksheet(wb, "enrichment_results")
    writeData(wb, "enrichment_results", filtered_data, 
              headerStyle = createStyle(textDecoration = "Bold", border = "Bottom"))
    conditionalFormatting(wb, "enrichment_results", 
                          cols = grep("negative_log10_of_adjusted_padj", colnames(filtered_data)), 
                          rows = 1:(nrow(filtered_data) + 1), 
                          type = "databar", 
                          rule = c(0, 50), 
                          style = c("#0BA31B", "#3FE075"), 
                          gradient = FALSE)
    setColWidths(wb, "enrichment_results", cols = 1:ncol(filtered_data), widths = "auto")
    saveWorkbook(wb, excel_file, overwrite = TRUE)
    
    # Extract and save the "term_id" column as a .txt file for ORSUM
    # Remove terms bigger than "orsum_term_size_limit" if specified
    if (orsum_term_size_limit == "none") {
    term_id_data <- source_data[, "term_id", drop = FALSE]
    txt_file <- file.path(orsum_dir, paste0("gprofilerres__", sample_name, "__", gsub(":", "_", source), ".txt"))           # Double underscore helps dealing with "_" in group naming
    write.table(term_id_data, file = txt_file, sep = "\t", 
                row.names = FALSE, col.names = FALSE, quote = FALSE)
    } else {
    source_data_subset = source_data[source_data$term_size < orsum_term_size_limit,]
    term_id_data <- source_data_subset[, "term_id", drop = FALSE]
    txt_file <- file.path(orsum_dir, paste0("gprofilerres__", sample_name, "__", gsub(":", "_", source), ".txt"))           # Double underscore helps dealing with "_" in group naming
    write.table(term_id_data, file = txt_file, sep = "\t", 
                row.names = FALSE, col.names = FALSE, quote = FALSE)  
    }
  }
}

# Function 2/3
# Plotting gprofiler results
plot_gprofiler2 <- function(df, output_dir, col_name, nterms = 15, title = 'gprofiler enrichment') {
  # Ensure output directory exists
  plot_dir <- file.path(output_dir, col_name, "Figures")
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Filter sources of interest
  df <- df[df$source %in% sources, ]
  
  # Iterate through sources
  for (source in unique(df$source)) {
    # Subset data for this source
    df_source <- df[df$source == source, ]
    
    # Sort and select the top terms
    df_source <- df_source[order(df_source$negative_log10_of_adjusted_p_value, decreasing = TRUE), ]
    df_source <- head(df_source, nterms)
    
    # Skip if there are no terms to plot
    if (nrow(df_source) == 0) next
    
    # Prepare the plot
    plot <- ggplot(df_source, aes(x = reorder(term_name, negative_log10_of_adjusted_p_value),
                                  y = negative_log10_of_adjusted_p_value)) +
      geom_point(aes(size = intersection_size, color = source)) +
      labs(title = glue("{title} - {col_name} ({source})"),
           subtitle = "Functional enrichment analysis",
           x = "Term Name",
           y = "-log10(Adjusted P-Value)",
           size = "Gene Count") +
      theme_minimal() +
      theme(axis.text.y = element_text(size = 10, hjust = 1),
            axis.title.y = element_blank()) +
      coord_flip()
    
    # Save the plot
    source_safe <- gsub(":", "_", source)  # Replace ":" to avoid file name issues
    plot_file <- file.path(plot_dir, glue("{col_name}_{source_safe}.png"))
    ggsave(plot, filename = plot_file, width = 10, height = 7)
    
    print(glue("Plot saved: {plot_file}"))
  }
}

# Function 3/3
# generate ORSUM command lines and save them into a .txt
# the .txt will next be parsed to ORSUM via the bash script
generate_orsum_commands <- function(ORSUM_output_dir = ORSUM_output_dir, ORSUM_input_dir = file.path(output_dir, "ORSUM"), gmt_dir, patterns) {
  # Output directory for ORSUM commands
  orsum_command_file <- file.path(ORSUM_output_dir, "Orsum_command_lines.txt")
  if (file.exists(orsum_command_file)) {
    file.remove(orsum_command_file)  # Remove old file if it exists
  }
  
  # Helper function to generate the command line for ORSUM
  generate_command_line <- function(pattern, files, ORSUM_input_dir) {
    command_line <- paste0(
      "orsum.py --gmt '", gmt_dir, "/", organism, ".", pattern, ".name.gmt' --files ",
      paste(paste0("'", ORSUM_input_dir, "/", files, "'"), collapse = " "),
      " --fileAliases ",
      paste(paste0(sub("gprofilerres__(.*?)__.*", "\\1", files), collapse = " ")),       # Names for heatmap
      " --outputFolder '", ORSUM_output_dir, "/Output_", pattern, "'\n"
    )
    return(command_line)
  }
  
  # Loop over patterns
  for (pattern in patterns) {
    print(glue("Processing pattern: {pattern}"))
    
    # Search for files matching the pattern
    files <- list.files(path = ORSUM_input_dir, pattern = pattern)
    # Reorder files correctly in case they have a number (i.e.: Group_1, Group_2, ...)
    files <- gtools::mixedsort(files)
    print(glue("Found files: {paste(files, collapse = ', ')}"))
    
    # If sufficient files are found, generate the command line
    if (length(files) >= 1) {
      command_line <- generate_command_line(pattern, files, ORSUM_input_dir)
      
      # Append the command line to the text file
      write(command_line, file = orsum_command_file, append = TRUE)
    }
  }
  print(glue("ORSUM command lines saved to {orsum_command_file}"))
}

############################
####### RUN GPROFILER ######
############################

# Load input data
resSig_df <- read.table(input_filepath, header = TRUE, sep = "\t", quote = "")

# Loop through each column (except the first column if it's an ID or metadata column)
for (col_name in colnames(resSig_df)) {
  # Extract the gene list from the current column
  gene_list <- resSig_df[[col_name]]
  
  # Skip if the column is not relevant (e.g., non-numeric or metadata)
  if (!is.character(gene_list) && !is.factor(gene_list)) next
  
  # Remove NA or invalid entries
  gene_list <- na.omit(as.character(gene_list))
  # Remove empty strings
  gene_list2 <- gene_list[!(gene_list == "")]


  if (length(gene_list) == 0) next
  
  # Perform enrichment analysis using gprofiler
  gostres <- gost(query = gene_list,
                  organism = organism,
                  user_threshold = 0.05,
                  significant = TRUE,
                  ordered_query = ordered_query,
                  domain_scope = "annotated",
                  correction_method = "fdr",
                  sources = sources,
                  evcodes = TRUE,
                  exclude_iea = TRUE)
  
  # Skip column if genes are not represented in pathway databases
  if (is.null(gostres)) {
    print(glue("No pathways were found enriched for: {col_name}"))
  } else {
  # Extract and format results
  df <- as.data.frame(gostres$result)
  df$negative_log10_of_adjusted_p_value <- -log10(df$p_value)
  
  # Generate and save plots
  plot_gprofiler2(df = df, col_name = col_name, output_dir = output_dir, nterms = 15)
  
  # Export results to Excel (and tsv), and save ORSUM file
  export_gprofiler_excel(gost.res = df, sample_name = col_name, output_dir = output_dir)

  # Write ORSUM command lines
  generate_orsum_commands(ORSUM_output_dir = ORSUM_output_dir, gmt_dir = gmt_dir, patterns = patterns)
  }
}
