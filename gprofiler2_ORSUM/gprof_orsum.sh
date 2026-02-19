#!/bin/bash
#SBATCH --job-name=gprof_orsum
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --error=gprof_orsum_%j.err
#SBATCH --output=gprof_orsum_%j.out

# EXAMPLE USE CASE
# sbatch gprof_orsum.sh -i Publication_Macrophages_gene_groups.tsv -o /results -gmt /gprofiler/gmt/hsapiens -os '500' -org 'hsapiens'

# IMPORTANT NOTES
# gmt files for your species of interest can be downloaded from "Data sources" here https://biit.cs.ut.ee/gprofiler/gost
# within the R script, see lines 73-77 to see GMT and organism dependent parameters
# gmt file names need to be renamed by substituing the ":" for "_" (i.e.: GO:BP to GO_BP)
#
# input file needs to be a TSV, where column names will be used for directory organization and plot labeling
# it is recommended to use underscores ("_") to separate words int eh tsv headers. Spaces will be showed as "." in the output folder and ORSUM heatmaps
#
# it is recommended to introduce the absolute paths in the bash console (-i, -o and -gmt), since they will be saved in a file at the results directory
#

# definition of usage instructions in case a parameter is not defined correctly
usage ()
{
  echo ' Usage : gprof_orsum_merged_modules.sh -i <input_file> -o <out_dir> -gmt <gmt_dir> -or <orsum_dir> -oq <ordered_query> -os <orsum_term_size_limit> -org <organism>'
  echo " Options:"
  echo "            -i|--input_file                     Absolute path to input file. Input file should be a tab-separated text file. Note that column names will be used as comparison names, and will be used for results organization."
  echo "            -o|--out_dir                        Absolute path to the output directory used to store gprofiler and ORSUM results."
  echo "            -gmt|--gmt_dir                      Absolute path to the directory with the GMT files necessary for ORSUM."
  echo "            -oq|--ordered_query                 Boolean (TRUE, FALSE) gprofiler parameter to obtain GSEA-like adjusted p-values if your gene lists are ranked. Note that p-value should NOT be used as a ranking metric in gprofiler. Default is FALSE."
  echo "            -os|--orsum_term_size_limit         Maximal gprofiler term size to be used by ORSUM to summarize results. Default is 'none' (unfiltered gprofiler results), but can be set to any number (i.e.: '500')."
  echo "            -org|--organism                     Supported organism by gost(). If unespecified, default is 'hsapiens'."
  echo
  exit
}

POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

# read console input
case $key in
    -i|--input_file)
    input_file="$2"
    shift # past argument
    shift # past value
    ;;
    -o|--out_dir)
    out_dir="$2"
    shift # past argument
    shift # past value
    ;;
    -gmt|--gmt_dir)
    gmt_dir="$2"
    shift # past argument
    shift # past value
    ;;
    -oq|--ordered_query)
    ordered_query="$2"
    shift # past argument
    shift # past value
    ;;
    -os|--orsum_term_size_limit)
    orsum_term_size_limit="$2"
    shift # past argument
    shift # past value
    ;;
    -org|--organism)
    organism="$2"
    shift # past argument
    shift # past value
    ;;
esac
done

# restore positional parameters
set -- "${POSITIONAL[@]}" 

if [[ -n $1 ]]; then
    echo "Last line of file specified as non-opt/last argument:"
    tail -1 "$1"
fi

if [ "$input_file" = "" ]
then
    usage
fi

if [ "$out_dir" = "" ]
then
    usage
fi

if [ "$gmt_dir" = "" ]
then
    usage
fi

if [ "$ordered_query" = "" ]
then
    ordered_query='FALSE'
fi

if [ "$orsum_term_size_limit" = "" ]
then
    orsum_term_size_limit='none'
fi

if [ "$organism" = "" ]
then
    organism='hsapiens'
fi

# load necessary environment for the R script to run, and for Python3 to run ORSUM
source ~/.bashrc
conda activate gprof_orsum_env

# create ORSUM results output directory
orsum_dir="${out_dir}/ORSUM_results/"

# ensure output and orsum dir exists
if [ ! -d "$orsum_dir" ]; then
    echo "Output directory does not exist. Creating it now..."
    mkdir -p "$orsum_dir"
else
    echo "Output directory already exists."
fi

# find bash script directory to run R script (needs to be located in the same directory)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# initiallize count for runtime
START_TIME=$(date +"%Y-%m-%d %H:%M:%S")
SECONDS=0

# save run parameters into a text (run summary) file at the output directory
PARAMS_FILE="${out_dir}/run_parameters.txt"
echo "Scripts directory: ${SCRIPT_DIR}" > "$PARAMS_FILE"
echo "Input File: ${input_file}" >> "$PARAMS_FILE"
echo "Output Directory: ${out_dir}" >> "$PARAMS_FILE"
echo "GMT Directory: ${gmt_dir}" >> "$PARAMS_FILE"
echo "Ordered Query: ${ordered_query}" >> "$PARAMS_FILE"
echo "ORSUM Term Size Limit: ${orsum_term_size_limit}" >> "$PARAMS_FILE"
echo "Organism: ${organism}" >> "$PARAMS_FILE"

# call R script
Rscript "${SCRIPT_DIR}/gprofiler_standalone.R" \
         "--orsum_dir=${orsum_dir}" \
         "--gmt_dir=${gmt_dir}" \
         "--output_dir=${out_dir}" \
         "--input_filepath=${input_file}" \
         "--ordered_query=${ordered_query}" \
         "--orsum_term_size_limit=${orsum_term_size_limit}" \
         "--organism=${organism}"

echo "${orsum_dir}"
# parse commands to ORSUM from the text file generated in the R script
while IFS= read -r line; do
    eval "$line"
done < "${orsum_dir}/Orsum_command_lines.txt"

# save runtime in the run summary file
END_TIME=$(date +"%Y-%m-%d %H:%M:%S")
ELAPSED_TIME=$SECONDS

# convert elapsed time to human-readable format (hh:mm:ss)
HOURS=$((ELAPSED_TIME / 3600))
MINUTES=$(( (ELAPSED_TIME % 3600) / 60 ))
SECONDS=$((ELAPSED_TIME % 60))
RUNTIME=$(printf "%02d:%02d:%02d" $HOURS $MINUTES $SECONDS)

# save run time in the run summary file
echo "Run Start Time: ${START_TIME}" >> "$PARAMS_FILE"
echo "Run End Time: ${END_TIME}" >> "$PARAMS_FILE"
echo "Total Runtime: ${RUNTIME}" >> "$PARAMS_FILE"
