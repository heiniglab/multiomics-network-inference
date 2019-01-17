## Get the randwom_walk data:
Currently we are relying on some precalculated information which we atm
cannot reproduce using our workflow (TODO). Therefore in a first step
we just link the needed data to our working directory.

```{bash}
d="data/current/networks/"
mkdir $d
cd $d
for i in /storage/groups/groups_epigenereg/analyses/meQTLs/results/20170517/networks/random_walk/rw_string_v9_ld_wb_plots/*.RData ; do 
  ln -s $i ; 
done
cd $WDIR
```

## Snakemake pipeline
The analysis is implemented as a snakemake pipeline.
Both the simulation and the cohort-data based analyses are integrated in a single workflow,
although split into separate snakemake files, and both rely on the subworkflow located under
'workflows/1_extract_hotspts.sm'. This workflow needs to be run before everything else, since
it extracts the hotspot loci (as dummy sentinel files) on which the rest of the pipeline
is based. Although this is automatically done when calling any downstream rules, the hotspot
extraction can also be done individually:

```{bash}
snakemake -s workflows/1_extract_hotspots.sm all
```

There are several meta targets to obtain intermediate results.

### Generate ranges overview
The call below creates all hotspots networks as 'ranges' objects/files and generates
a summary plot.

```{bash}
nohup nice snakemake -j 10 -k all_ranges &
```

### Generate data overview
This call collects and normalized all cohort data for the created ranges objects and
generates an overview plot.

```{bash}
nohup nice snakemake -j 10 -k all_data &
```

### Generate overview on GGM fits
This generates all ranges and data collections and fits different GGMs to the data

```{bash}
nohup nice snakemake -j 10 --res mem_mb=80000 -k all_ggm &
```

### Running the complete cohort data pipeline and cluster execution
The below code can be used to run the full network inference pipeline on a SGE cluster,
althrough the cluster configuration (cluster.config) should be adjusted before 
executing. Cohort and simulation study can be run separately by using either the 
*all_cohort* or *all_simulation* rule, respectively. To call the full pipeline 
including both studies the *all* rule is used:

```{bash}
nohup nice snakemake -u cluster.config --jobs=100 --local-cores=10 \
           --cluster "qsub -pe smp {threads} -hard -l job_mem={resources.mem_mb}M \
           -q {cluster.q} -cwd -V -o {log} -e {log} -N {cluster.N}" all & 
```

