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
#SBATCH --mem=450gb 
#time in days-hours:mins:sec
#SBATCH --time=0-10:30:00  

# Load modules 
module load GCC/11.3.0   
module load OpenMPI/4.1.4  
module load R/4.4.0 
module load Anaconda3/2024.02-1  

# Activate conda env with python and bbmap 
eval "$(conda shell.bash hook)" 
conda activate bbmap_env 
echo "Conda environment initialized"  

fastq="/data/gpfs/projects/punim2251/Neurodevelopmental_Models/iPSC_Neuron_Diff/matched_reads.fastq" 
outputDir="/data/gpfs/projects/punim2251/Aim1_LongBench/ReadRarefaction_wFixedCells/data/iPSC_D55_Rarefaction"  
./1_coarseSubsamplingScript.sh ${fastq} ${outputDir}