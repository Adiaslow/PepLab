#!/bin/bash
#SBATCH --job-name=gen_pep_lib
#SBATCH --output=peptide_library_%j.out
#SBATCH --error=peptide_library_%j.err
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1

# Check if correct number of arguments provided
if [ "$#" -ne 2 ]; then
    echo "Usage: sbatch $0 input_file.csv config.json"
    exit 1
fi

# Get input arguments
INPUT_FILE=$1
CONFIG_FILE=$2

# Validate input files exist
if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: Input file $INPUT_FILE not found"
    exit 1
fi

if [ ! -f "$CONFIG_FILE" ]; then
    echo "Error: Config file $CONFIG_FILE not found"
    exit 1
fi

# Load required modules (adjust based on your cluster setup)
module load rdkit

# Activate virtual environment if needed
# source /path/to/your/venv/bin/activate

# Run the peptide library generator
python generate_peptide_library.py \
    --input "${INPUT_FILE}" \
    --config "${CONFIG_FILE}" \
    --verbose

# Check exit status
if [ $? -eq 0 ]; then
    echo "Peptide library generation completed successfully"
    echo "Library saved to library.csv"
else
    echo "Peptide library generation failed"
    exit 1
fi
