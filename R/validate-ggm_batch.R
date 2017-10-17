# submits the ggm calculation for each availablie sentinel to the server
# parses the sentinels from the available data files
load("results/data.processed.RData")

for(sentinel in names(data)){
  cmd <- paste0("qsub -cwd -V -q long_fed25 -pe smp 1 -hard -l job_mem=4G -b y ")
  cmd <- paste0(cmd, "-N ", sentinel," -o ", sentinel, ".validation..out -e ", sentinel, ".validation.out ")
  cmd <- paste0(cmd, "Rscript 2-validate-ggm.R ", sentinel)
  system(cmd)
}
