#!/bin/bash
#SBATCH --partition="cascade"
#SBATCH --nodes=1
#SBATCH --account="punim2251"
#SBATCH --ntasks=1
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-user=mschauhan@student.unimelb.edu.au
#SBATCH --cpus-per-task=4
#SBATCH --mem=350gb
#time in days-hours:mins:sec
#SBATCH --time=1-00:00:00

source /home/mschauhan/anaconda3/bin/activate /home/mschauhan/.conda/envs/bbmap_env

echo "Conda environment initialized"

# Load modules
module load GCC/11.3.0  
module load OpenMPI/4.1.4 
module load R/4.4.0

fastq="/data/gpfs/projects/punim2251/LongBench_data/flames_out/ont_sn/matched_reads.fastq"
outputDir="/data/gpfs/projects/punim2251/Aim1_LongBench/ReadRarefaction_wFixedCells/data/LongBench_All/ont_sn"
targetMedians=(1291 517) # List of target medians to iterate over


for targetMedian in "${targetMedians[@]}"; do
    echo "Processing target median: ${targetMedian}"
    
    # Run the R script with the current target median
    Rscript subsampleToExactMedian.R ${fastq} ${targetMedian} ${outputDir}
    
    # Rename the log file to avoid errors
    mv ${outputDir}/fine_MedianBarcodeReads_log.csv ${outputDir}/${targetMedian}_intermediate_log.csv
done