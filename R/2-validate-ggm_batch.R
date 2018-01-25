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

rdir <- "results/current/fits_nopriors/"
vdir <- "results/current/validation_nopriors/"
dir.create(vdir)

for(sentinel in sentinels) {
  # check whether output files exist
  f1 <- paste0(rdir, sentinel, ".lolipop.RData")
  f2 <- paste0(rdir, sentinel, ".kora.RData")
  f3 <- paste0(vdir, sentinel, ".validation.txt")
  if(file.exists(f1) & file.exists(f2) & !file.exists(f3)) {
    cmd <- paste0("qsub -cwd -V -q long_fed25 -pe smp 1 -hard -l job_mem=4G -b y ")
    cmd <- paste0(cmd, "-N ", sentinel, " ",
                  "-o ", vdir, sentinel, ".validation.out ", 
                  "-e ", vdir, sentinel, ".validation.out ")
    cmd <- paste0(cmd, "Rscript R/2-validate-ggm.R ", sentinel)
    system(cmd)
  } else {
    if(exists(f3)) {
      cat("Results for sentinel", sentinel, "already exist.\n", 
          file="results/current/validation.call.failed.txt", append=T)
    } else {
      cat(sentinel, "\n", file="results/current/validation.call.failed.txt", append=T)
    }
  }
}
