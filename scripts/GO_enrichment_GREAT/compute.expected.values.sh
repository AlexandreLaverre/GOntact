#!/bin/bash

export sp=$1
export method=$2 ## classical_upstream5kb_downstream1kb_extend1Mb
export space=$3 ## biological_process, molecular_function, cellular component

####################################################################################

export path=/beegfs/data/necsulea/GOntact
export pathGO=${path}/data/GeneOntology
export pathGenomes=${path}/data/genome_sequences
export pathResults=${path}/results/GO_enrichment_GREAT
export pathScripts=${path}/scripts/GO_enrichment_GREAT

export ensrelease=94

####################################################################################

echo "#SBATCH --job-name=exp_${sp}_${method}_${space}" >>  ${pathScripts}/bsub_script_expected
echo "#SBATCH --partition=normal" >>  ${pathScripts}/bsub_script_expected
echo "#SBATCH --output=${pathScripts}/std_out_exp_${sp}_${method}_${space}" >>  ${pathScripts}/bsub_script_expected
echo "#SBATCH --error=${pathScripts}/std_err_exp_${sp}_${method}_${space}" >>  ${pathScripts}/bsub_script_expected
echo "#SBATCH --cpus-per-task=1" >>  ${pathScripts}/bsub_script_expected ## 1 CPUs
echo "#SBATCH --time=3:00:00" >>  ${pathScripts}/bsub_script_expected ## 3 hours
echo "#SBATCH --mem=20G" >>  ${pathScripts}/bsub_script_expected ## 20g per CPU

echo "perl ${pathScripts}/compute.expected.values.pl --pathGOCategories=${pathGO}/GOCategories.txt --pathGOAnnotations=${pathGO}/${sp}.gene.annotation.txt --pathRegulatoryRegions=${pathResults}/${sp}/${method}/regulatory_regions_Ensembl${ensrelease}_${space}.txt --pathGenomeSequence=${pathGenomes}/${sp}/genome_Ensembl${ensrelease}.fa --pathOutput=${pathResults}/${sp}/${method}/expected_values_Ensembl${ensrelease}_${space}.txt " >>  ${pathScripts}/bsub_script_expected 

sbatch ${pathScripts}/bsub_script_expected

####################################################################################
