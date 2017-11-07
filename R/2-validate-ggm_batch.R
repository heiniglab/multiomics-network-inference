# submits the ggm calculation for each availablie sentinel to the server
# parses the sentinels from the available data files
load("results/current/data.processed.RData")

for(sentinel in c("rs3809627","rs17850402","rs12412214","rs11793438")) { #names(data)){
  cmd <- paste0("qsub -cwd -V -q long_fed25 -pe smp 1 -hard -l job_mem=4G -b y ")
  cmd <- paste0(cmd, "-N ", sentinel, " ",
                "-o results/current/", sentinel, ".validation.out ", 
                "-e results/current/", sentinel, ".validation.out ")
  cmd <- paste0(cmd, "Rscript R/2-validate-ggm.R ", sentinel)
  system(cmd)
}
