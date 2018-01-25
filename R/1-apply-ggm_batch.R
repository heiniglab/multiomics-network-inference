library(optparse)
# submits the ggm calculation for each availablie sentinel to the server
# parses the sentinels from the available data files

sentinels <- NULL

# check whether arguments provide a list of sentinels
opt = parse_args(OptionParser(option_list=list(
  make_option(c("-f", "--file"), type="character", default=NULL),
  make_option("--nopriors", type="logical", action="store_true", default=FALSE),
  make_option(c("--cores", "-c"), type="integer", default=10))))

ifile <- opt$file
nopriors <- opt$nopriors
cores <- opt$cores

if(!is.null(ifile)){
  if(!file.exists(ifile)){
   stop("File with sentinels does not exist.")  
  }
  sentinels <- read.table(ifile, header=F)[,1]
} else {
  # instead of argument, load sentinel names from the preprocessed data
  load("results/current/data.processed.RData")
  sentinels <- sort(names(data), T)
}

cat("Loaded", length(sentinels), "sentinels.\n")

# create output dirs
if(nopriors) {
  fodir <- "results/current/fits_nopriors/"
  podir <- "results/current/plots_nopriors/"
} else {
  fodir <- "results/current/fits/"
  podir <- "results/current/plots/"

}
dir.create(fodir)
dir.create(podir)

for(sentinel in sentinels) {
  if(!file.exists(paste0(fodir,sentinel,".RData"))) {
    cmd <- paste0("qsub -cwd -V -q long_fed25 -pe smp ", cores, " -hard -l job_mem=1G -b y ")
    cmd <- paste0(cmd, "-N ", sentinel," -o ", fodir, sentinel, 
                  ".out -e ", fodir, sentinel, ".out ")
    cmd <- paste0(cmd, "Rscript R/1-apply-ggm.R -s ", sentinel, " -c ", cores)
    # check whether to add nopriors flag
    if(nopriors) {
      cmd <- paste0(cmd, " --nopriors")
    }
    system(cmd)
  }
}
