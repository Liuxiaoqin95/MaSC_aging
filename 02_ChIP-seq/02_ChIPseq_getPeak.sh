###get binding peaks from mapped bam file

#!/bin/bash

#slurm options
#SBATCH -p intel-debug                  
#SBATCH -q debug                             
#SBATCH -J getBED                      
#SBATCH -c 10                              
#SBATCH -o getBED.log                         

macs2=$path_to_macs2
treated=$path_to_treatment
control=$path_to_IgG
id=$result_file_name

$macs2 callpeak -B -t $treated -c $control -n $id -f BAMPE -g mm --outdir $resultdir --keep-dup all --nomodel 

