# submits the ggm calculation for each availablie sentinel to the server
# parses the sentinels from the available data files
sentinels <- NULL

# check whether arguments provide a list of sentinels
args <- commandArgs(trailingOnly =T )
if(!is.null(args[1]) & !is.na(args[1])){
  if(!file.exists(args[1])){
    stop("File with sentinels does not exist.")  
  }
  sentinels <- read.table(args[1], header=F)[,1]
} else {
  # instead of argument, load sentinel names from the preprocessed data
  load("results/current/data.processed.RData")
  sentinels <- names(data)
}
cat("Loaded", length(sentinels), "sentinels.\n")

for(sentinel in sentinels) {
  # check whether output files exist
  rdir <- "results/current/fits/"
  f1 <- paste0(rdir, sentinel, ".lolipop.RData")
  f2 <- paste0(rdir, sentinel, ".kora.RData")
  if(file.exists(f1) & file.exists(f2)) {
    cmd <- paste0("qsub -cwd -V -q long_fed25 -pe smp 1 -hard -l job_mem=4G -b y ")
    cmd <- paste0(cmd, "-N ", sentinel, " ",
                  "-o results/current/validation/", sentinel, ".validation.out ", 
                  "-e results/current/validation/", sentinel, ".validation.out ")
    cmd <- paste0(cmd, "Rscript R/2-validate-ggm.R ", sentinel)
    system(cmd)
  } else {
    cat(sentinel, file="results/current/validation.call.failed.txt\n", append=T)
  }
}
