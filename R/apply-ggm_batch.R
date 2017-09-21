# submits the ggm calculation for each availablie sentinel to the server
# parses the sentinels from the available data files
files <- list.files("data/kora/", ".*.data.RData", full.names = F)
for(f in files){
  # get the snp
  snp <- strsplit(f,"\\.")[[1]][1]
  cmd <- paste0("qsub -cwd -V -q long_fed25 -pe smp 8 -hard -l job_mem=2G -b y ")
  cmd <- paste0(cmd, "-o ", snp, ".out -e ", snp, ".out ")
  cmd <- paste0(cmd, "Rscript R/apply-ggm_server.R ", snp)
  system(cmd)
}
