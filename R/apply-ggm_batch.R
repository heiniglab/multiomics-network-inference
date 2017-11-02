# submits the ggm calculation for each availablie sentinel to the server
# parses the sentinels from the available data files
load("results/current/data.processed.RData")

cores <- 10

for(sentinel in names(data)){
  cmd <- paste0("qsub -cwd -V -q long_fed25 -pe smp ", cores, " -hard -l job_mem=2G -b y ")
  cmd <- paste0(cmd, "-N ", sentinel," -o results/current/fits/", sentinel, 
                ".out -e results/current/fits/", sentinel, ".out ")
  cmd <- paste0(cmd, "Rscript R/apply-ggm_server.R ", sentinel, " ", cores)
  system(cmd)
}
