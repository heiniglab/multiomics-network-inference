## Table of Contents

[[_TOC_]]

## Online collaboration resources:

- [Paper on overleaf](https://www.overleaf.com/2873145782fcvqtppkhgdb)
- [Reviewer comments on google sheets](https://docs.google.com/spreadsheets/d/1rtgCK2rgXGBGJyT-Xg43TGysO_0IBbDzELT9FNgs7do/edit?usp=sharing)
- [Response to reviewer comments on google docs](https://docs.google.com/document/d/1x9KLl1AKVdHe64g4aDI2arbfn_9Jjp2zWgFmfdooiMs/edit?usp=sharing)

-----

- [seal python package](https://github.com/muhanzhang/SEAL/tree/master/Python)
- [scanpy in R](https://lazappi.id.au/project/scanpy-in-r/)
- [reticulate R package](https://cran.r-project.org/web/packages/reticulate/index.html)

## Getting the randwom_walk data:

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

## Snakemake pipeline: cohort study, simulation and benchmark

The analysis is implemented as a snakemake pipeline.
Both the simulation and the cohort-data based analyses are integrated in a single workflow,
although split into separate snakemake files, and both rely on the subworkflow located under
'workflows/1_extract_hotspts.sm'. This workflow needs to be run before everything else, since
it extracts the hotspot loci (as dummy sentinel files) on which the rest of the pipeline
is based. Although this is automatically done when calling any downstream rules, the hotspot
extraction can also be done individually:

```{bash}
snakemake --profile default -s workflows/1_extract_hotspots.sm all
```

> NOTE: For some of the hotspots, none of the respective SNP genes have any expression
> probes in our data. We have to remove them manually. The list of SNPs is:
> rs57743634,rs17420384,rs2295981,rs7924137,rs1570038,rs57743634,rs2685252

### Snakemake profiles

We use profiles to handle local and cluster execution (saved under [./profiles/](./profiles)).

- `default` profile: available under `./profiles/default/`, use `--profile profiles/default` for snakemake
- `slurm` profile: available under `./profiles/slurm/`, use `--profile profiles/slurm` for snakemake

Both profiles specify some default snakemake parameters. The SLURM profile also specifies necessary SLURM parameters specific to the ICB queue. The SLURM profile was installed using cookiecutter (compare [this link](https://github.com/Snakemake-Profiles/slurm)) then modified to fit our queue:

```
# install cookiecutter if not yet available
# conda install -c conda-forge cookiecutter

# install SLURM profile (use appropriate path under profiles/)
cookiecutter https://github.com/Snakemake-Profiles/slurm.git

# move to ./profiles/ 
# mkdir profiles
# mv slurm/ profiles/

# modify profiles/slurm/* as needed...
```

### Charliecloud image and conda environment usage

At the moment of creating the workflow, the only option of using containers
on the cluster was via charliecloud, which unfortunately is not supported by snakemake.
To make it work, we manually edited the `jobscript.sh` from snakemake to extract the
image if needed. In addition, we edited the `script.py` from snakemake to wrap any Rscript
calls in a `ch-run` call with our specific charliecloud container. This is not ideal and not
very portable, but allows as to have a simple software container in place which we can use on 
other systems, too (such as MARCC from JHU).

In addition, a general conda environment is defined in *envs/bioR.yaml". This environment
has been linked to all rules which are based on R scripts. In principle, this conda env could be
used instead of the charliecloud image, but it is **deprecated** now.
To execute the pipelin using conda, you have to specify the *--use-conda*
parameter in the snakemake call (done by default in the `./profiles/default` snakemake profile).

> NOTE: Some packages were not able to be installed via Conda. It is absolutely
> necessary to make an initial call to snakemake which installs remaining 
> packages like so: **snakemake --use-conda config_r** .
> NOTE 2: Be aware that if you set any repository paths on startup of R 
> you might want to adjust that for usage with the conda environment

### Meta target rules

There are some 'meta' targets which can be used to run only parts of the pipeline for all hotspots.

#### Generate ranges overview

The call below creates all hotspots networks as 'ranges' objects/files and generates
a summary plot.

```{bash}
snakemake --profile=./profiles/default all_ranges
```

#### Generate data overview

This call collects and normalized all cohort data for the created ranges objects and
generates an overview plot.

```{bash}
snakemake --profile=./profiles/default all_data
```

#### Generate overview on GGM fits

This generates all ranges and data collections and fits different GGMs to the data. Note that we use the SLURM profile, as running this locally is not entirely sensible.

```{bash}
snakemake --profile=./profiles/slurm/ all_ggm &
```

### Cohort data pipeline

The below code can be used to run the full network inference pipeline on our SLURM cluster on the cohort data.

> NOTE: Cohort and simulation study can be run separately by using either the 
*all_cohort* or *all_simulation* (see below) rule/commands, respectively.

To execute the cohort study using the SLURM cluster (using the `./profiles/slurm` profile), call:

```{bash}
./submit_slurm_cohort_study.sh
```

#### Replication analysis with prior noise

Based on cohort replication, we also investigate the effect of noisy priors.
To run all needed model fits including systematic creation of noise priors, simply call:

```{bash}
./submit_slurm_cohort_replication_prior_noise.sh
```

### Method benchmarking

We implemented a benchmarking procedure to formally test runtimes for all applied models.
The complete benchmark can be run using the following call:

```
./submit_slurm_benchmark.sh
```

The above script implicetely calls snakemake with the `gather_benchmark_results` rule and submits all SLURM jobs to the same compute node (currently: `icb-neu-001`) to achieve comparable performance estimates.
Results are summarized under `results/current/benchmark/summary.pdf`.

### Simulation study

We implemented a simulation study, where we generate ground truth graphs as well
as noisy priors in order to compare inference methods with respect to 1) their 
general performance to recover the ground truth network 2) how much they are 
influenced by the exisiting prior information and noise therein (in addition to their computational complexity, see above). 
The target rule for this simulation study is `all_simulation` and is implicetely used when calling:

```{bash}
./submit_slurm_simulation.sh
```

The above command will execute the simulation study using the SLURM cluster.

## DEPRECATED: LRZ specific submit (simulations)

At some point we had to run our simulation studies on the LRZ cluster. As this is not necessary anymore, below code is merely here for documentation purposes.

```{bash}
module load python/3.6_intel
nohup nice snakemake --use-conda -u configs/slurm.json --jobs=1000 --local-cores=1 --cluster \
  "sbatch --time=24:00:00 --ntasks 1 --clusters=mpp2 -c {cluster.cpu} --mem-per-cpu={cluster.mem} \
      -o {cluster.log} -e {cluster.log}" all_simulation &> all_simulation.out &
```

