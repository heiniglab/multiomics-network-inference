# Hawe et al. 2022 - Network reconstruction for trans acting genetic loci using multi-omics data and prior information

## Snakemake pipeline: cohort study, simulation and benchmark

The [snakemake workflow](Snakefile) first defines a set of trans-QTL hotspots from available QTL results (i.e. genetic variants with 5+ trans associations).
For each hotspot, a 'locus set' of entities is defined (see Methods in the manuscript) and cohort/simulated data collected.
On each of these datasets, different network inference methods are applied and evaluated.

### Main scripts

The most important script for the analyses are listed here (located under `scripts/`):

- [reg_net.R](scripts/reg_net.R): Main script containing methods for inferring regulatory networks given the collected data and locus set information. Relies on helper methods defined in [reg_net_utils.R](scripts/reg_net_utils.R)
- [collect_ranges.R](scripts/collect_ranges.R): The main script to collect all entities in a locus set for a given hotspot.
- [collect_ranges_methods.R](scripts/collect_ranges_methods.R): Defines helper methods for collecting locus sets.
- [collect_priors.R](scripts/collect_priors.R): Main script to collect all prior information for a given hotspot locus. Utilizes helper methods defined in [scripts/priors.R](scripts/priors.R)
- [benchmark.R](scripts/benchmark.R): This script and its corresponding [methods script](scripts/benchmark_methods.R) gathers all code for generating the runtime benchmarks.
- [create_priors.R](scripts/create_priors.R): Main script to create the global eQTL and gene-gene prior informartion
- [apply_ggm.R](scripts/apply_ggm.R): Main script to run GGM inference for a specific hotspot locus set.
- [simulation/simulate_data.R](scripts/simulation/simulate_data.R): Main script for generating simulated data for individual hotspots, including erroneous graphs/priors.
- [simulation/run_ggm.R](scripts/simulation/run_ggm.R): Main script to run the GGM inference on simulated data.

> NOTE: general helper methods are defined in separate `lib.R` scripts (main [lib.R](scripts/lib.R), simulation [lib.R](scripts/simulation/lib.R))

Generally, to implement your own analysis you'd want to

1) collect a set of trans QTLs
2) curate prior information as needed (e.g. utilizing the [priors.R](scripts/priors.R) script)
3) use the 'collect_ranges' methods to define locus sets
4) on the defined locus sets, follow the 'apply_ggm' script to infer GGM networks (this heavily utilizes the functions in [reg_net.R](scripts/reg_net.R) and [reg_net_utils.R](scripts/reg_net_utils.R)) 
  
### Configuration

There are few configuration options for the workflow. These can be adjusted 
in the (configs/workflow.json)[configs/workflow.json] file. 
Here is the list of options:

* `hots_thres` - Threshold for the number of trans entities to define a hotspot
* `ppi_db` - Either 'string', 'biogrid' or 'biogrid_stringent'
* `suffix_tfa_expr` - Whether to use expression or TF activities for TFs
* `eqtl_prior_type` - eQTL prior to be used, either 'eqtlgen' or 'gtex'

### Preprocessing 

Both the simulation and the cohort-data based analyses are integrated in a single workflow,
although split into separate snakemake files, and both rely on the subworkflow located under
'workflows/1_extract_hotspts.sm'. 
This workflow needs to be run before everything else as it extracts the hotspot loci (as dummy sentinel files) on which the rest of the pipeline
is based. 
Although this is automatically done when calling any downstream rules, the hotspot
extraction can also be done individually:

```{bash}
snakemake --profile default -s workflows/1_extract_hotspots.sm all
```

> NOTE: For some of the hotspots, none of the respective SNP genes have any expression
> probes in our data. We have to remove them manually. The list of SNPs is:
> rs57743634,rs17420384,rs2295981,rs7924137,rs1570038,rs57743634,rs2685252

### Cohort data pipeline

The below code can be used to run the full network inference pipeline on a SLURM cluster on the cohort data.
The snakemake rules for this part of the pipeline are defined in (`workflow/2_1_cohort_data.sm`)(workflow/2_1_cohort_data.sm).

> NOTE: Cohort and simulation study can be run separately by using either the 
*all_cohort* or *all_simulation* (see below) rule/commands, respectively.

To execute the cohort study using the SLURM cluster (using the `./profiles/slurm` profile, see below), call:

```{bash}
./submit_slurm_cohort_study.sh
```

See script contents for more details.

#### Replication analysis with prior noise

Based on cohort replication, we also investigate the effect of erroneous priors.
To run all needed model fits including systematic creation of erroneous priors, simply call:

```{bash}
./submit_slurm_cohort_replication_prior_noise.sh
```

To investigate the results, you can have a look at the Rmarkdown document under [validation_reviewer_comments.Rmd](validation_reviewer_comments.Rmd).

> The markdown is not yet implemented in the snakemake workflow

### Method benchmarking

We implemented a benchmarking procedure to formally test runtimes for all applied models.
The complete benchmark can be run using the following call:

```
./submit_slurm_benchmark.sh
```

The above script implicetely calls snakemake with the `all_benchmark` rule and submits all SLURM jobs to the same compute node (needs to be adjusted for different clusters) to achieve comparable performance estimates.
Results are summarized under `results/current/benchmark/summary.pdf`.

### Simulation study

We implemented a simulation study, were we generate ground truth graphs as well
as erroneous priors in order to compare inference methods with respect to 
  1) their general performance to recover the ground truth network 
  2) how much they are influenced by the exisiting prior information and noise therein and 
  3) the effect of sample size on inference performance. 
As we simulate data such that they match the properties of the original data (genotype frequencies,
sample size and number of nodes per locus), this is coupled to the original data 
workflow. 
Simulation specific rules can be found under [`workflows/2_2_simulation.sm`](workflows/2_2_simulation.sm).
The target rule for this simulation study is `all_simulation`.

```{bash}
./submit_slurm_simulation.sh
```

The above command will execute the simulation study using the SLURM cluster.

### Available target rules

There are several meta targets in the workflow to obtain intermediate results.
Call `snakemake <target-name>` to execute the specific part of the workflow.

* all_ranges - The call below creates a locus set for each hotspot as GenomicRanges objects/files and generates
a summary plot.
* all_data - This call collects and normalized all cohort data for the created ranges objects
* all_ggm - This generates all ranges and data collections and fits different network models to the data.
* plot_relevant_graphs - Generates *dot* based plots for some of the networks with desireable properties.
* all_cohort - Run the full cohort data analysis.
* all_simulation - Run full simulation study (see also below)

**Example: Generate ranges overview**

Use snakemake to create all 'ranges' (hotspot based locus sets) collections:

```{bash}
snakemake --profile=./profiles/default all_ranges
```

## Additional information

### Snakemake profiles

We use profiles to handle local and cluster execution (saved under [./profiles/](./profiles)).

- `default` profile: available under `./profiles/default/`, use `--profile profiles/default` for snakemake
- `slurm` profile: available under `./profiles/slurm/`, use `--profile profiles/slurm` for snakemake

Both profiles specify some default snakemake parameters. The SLURM profile also specifies necessary SLURM parameters specific to the ICB queue. The SLURM profile was installed using cookiecutter (compare [this link](https://github.com/Snakemake-Profiles/slurm)), then modified to fit our queue:

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

We were only able to use charliecloud on our HPC which is unfortunately not supported by snakemake.
To make it work, we manually edited the `jobscript.sh` from `snakemake` to extract the
image if needed. In addition, we edited the `script.py` from snakemake to wrap any Rscript
calls in a `ch-run` call with our specific charliecloud container. This is not ideal and not
very portable, but allows as to have a simple software container in place which we can use on 
other systems, too.

The [Dockerfile](Dockerfile) for this container is provided in this repository.
Generally, if e.g. Docker or Singularity is available on the system the workflow is executed, a standard 
installation of *Snakemake* can be used to exeucted the workflow in the respective environment.
Otherwise, one needs to modify the *script.py* script of *Snakemake* to execute all scripts within the
provided container.

**Deprecated** In addition, a general conda environment is defined in *envs/bioR.yaml*. This environment
has been linked to all rules which are based on R scripts. In principle, this conda env could be
used instead of the charliecloud image.
To execute the pipelin using conda, you have to specify the *--use-conda*
parameter in the snakemake call (done by default in the `./profiles/default` snakemake profile).

> NOTE: Some packages were not able to be installed via Conda. It is absolutely
> necessary to make an initial call to snakemake which installs remaining 
> packages like so: **snakemake --use-conda config_r** .
> NOTE 2: Be aware that if you set any repository paths on startup of R 
> you might want to adjust that for usage with the conda environment

### Getting the random walk data:

We use some previously established results (results from the random walk analysis) from the [Hawe et al. 2022](https://www.nature.com/articles/s41588-021-00969-x) publication.
The data are provided in this file](rw_string_v9_ld_wb_prioritize_full_with_empirical_p_lte_0.05.txt) of this repository.
