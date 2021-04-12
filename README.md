HRS

srun --account=bakulski1 --partition=largemem --nodes=1 --cpus-per-task=10 --time=300 --mem-per-cpu=30G --pty /bin/bash
module load Rtidyverse/3.6.1
Rscript main.R
