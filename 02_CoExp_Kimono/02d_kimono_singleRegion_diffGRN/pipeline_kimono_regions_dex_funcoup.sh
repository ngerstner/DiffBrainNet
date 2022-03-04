#!/bin/bash
#
#SBATCH --job-name=kimono_region
#SBATCH --output=err_out/kimono_region_%A_%a.out
#SBATCH --error=err_out/kimono_region_%A_%a.err
#SBATCH --array=0-25
#SBATCH --mem=5000
#SBATCH --cpus-per-task=12
#SBATCH --partition=hp
#SBATCH --exclude=hp11

region="dCA1"
dex=("0" "1")
startnodes=("1" "1001" "2001" "3001" "4001" "5001" "6001" "7001" "8001" "9001" "10001" "11001" "12001") 

nstartnodes=${#startnodes[@]}
ndex=${#dex[@]}

#get region and dex index for each job id
istartnode=$((SLURM_ARRAY_TASK_ID / ndex)) #divide task id by number of dex status
istartnode=$iregion|cut -f1 -d"." #take floor of the index
idex=$(($SLURM_ARRAY_TASK_ID%$ndex))

echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "${startnodes[$istartnode]}"
echo "${dex[$idex]}"
echo $istartnode
echo $idex

Rscript --vanilla 04_runKimono_funcoup.R $region ${dex[$idex]} ${startnodes[$istartnode]}
