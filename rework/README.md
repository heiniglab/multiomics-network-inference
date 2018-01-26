###
# data prep
###
We first need to gather some data from other projects:

cohort data:
TODO

get the randwom_walk data:
```{bash}
d="data/current/networks/"
mkdir $d
cd $d
for i in /storage/groups/groups_epigenereg/analyses/meQTLs/results/20170517/networks/random_walk/rw_string_v9_ld_wb_plots/*.RData ; do 
  ln -s $i ; 
done
```

###
# pipeline
###

The pipeline is described as a snakemake pipeline.
However, the already existing scripts merely have been altered such that they allow execution
using snakemake and likely more improvements could/should be done.

Modifications too the pipeline can be done in

./Snakfile

A new run of the pipeline can be performed (with cluster execution) using the following command:

```{bash}
nohup nice snakemake -u cluster.config --jobs=100 --local-cores=10 \
           --cluster "qsub -pe smp {threads} -hard -l job_mem={resources.mem_mb}M \
           -q {cluster.q} -cwd -V -o {cluster.o} -e {cluster.e} -N {cluster.N}" all & 
```
