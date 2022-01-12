#!/bin/bash

export sp=$1
export dataset=$2
export background=$3
export method=$4
export space=$5

####################################################################################

export path=/beegfs/data/necsulea/GOntact
export pathGO=${path}/data/GeneOntology
export pathEnhancers=${path}/data/enhancers
export pathResults=${path}/results/GO_enrichment_GREAT
export pathScripts=${path}/scripts/GO_enrichment_GREAT

export ensrelease=94

####################################################################################

if [ -e ${pathResults}/${sp}/${method}/${dataset} ]; then
    echo "path output exists"
else
    mkdir -p ${pathResults}/${sp}/${method}/${dataset} 
fi

####################################################################################

echo "#!/bin/bash" > ${pathScripts}/bsub_script_observed
echo "#SBATCH --job-name=exp_${sp}_${method}_${space}" >>  ${pathScripts}/bsub_script_observed
echo "#SBATCH --partition=normal" >>  ${pathScripts}/bsub_script_observed
echo "#SBATCH --output=${pathScripts}/log/std_out_obs_${sp}_${method}_${space}" >>  ${pathScripts}/bsub_script_observed
echo "#SBATCH --error=${pathScripts}/log/std_err_obs_${sp}_${method}_${space}" >>  ${pathScripts}/bsub_script_observed
echo "#SBATCH --cpus-per-task=1" >>  ${pathScripts}/bsub_script_observed ## 1 CPUs
echo "#SBATCH --time=3:00:00" >>  ${pathScripts}/bsub_script_observed ## 3 hours
echo "#SBATCH --mem=10G" >>  ${pathScripts}/bsub_script_observed ## 10g per CPU

echo "perl ${pathScripts}/compute.observed.values.pl --pathInputElements=${pathEnhancers}/${sp}/${dataset}.bed --pathBackgroundElements=${pathEnhancers}/${sp}/${background}.bed --pathGOCategories=${pathGO}/GOCategories.txt --pathGOAnnotations=${pathGO}/${sp}.simplified.gene.annotation.${space}.txt --pathRegulatoryRegions=${pathResults}/${sp}/${method}/regulatory_regions_Ensembl${ensrelease}_${space}.txt --pathOutput=${pathResults}/${sp}/${method}/${dataset}/observed_values_Ensembl${ensrelease}_background${background}_${space}.txt  --pathOutputAssociation=${pathResults}/${sp}/${method}/${dataset}/enhancer_GO_association_Ensembl${ensrelease}_${space}.txt " >>  ${pathScripts}/bsub_script_observed

####################################################################################

## do the test

echo "Rscript ${pathScripts}/test.enrichment.R ${sp} ${dataset} ${background} ${method} ${space} > ${pathScripts}/log/test.enrichment.${sp}.${dataset}.${method}.${space}.Rout " >>  ${pathScripts}/bsub_script_observed

####################################################################################

sbatch ${pathScripts}/bsub_script_observed

####################################################################################
