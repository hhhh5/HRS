Analyses for pre-print "Computational deconvolution of fifteen leukocyte subtypes from DNA methylation microarrays trained on flow cytometry data in the Health and Retirement Study"

Code was run on the UofM Great Lakes cluster

```
srun --account=bakulski1 --partition=largemem --nodes=1 --cpus-per-task=10 --time=300 --mem-per-cpu=30G --pty /bin/bash
module load Rtidyverse/4.2.0
Rscript main.R
```
