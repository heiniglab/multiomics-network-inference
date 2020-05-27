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

> NOTE: For some of the hotspots, none of the respective SNP genes have any expression
> probes in our data. We have to remove them manually. The list of SNPs is:
> rs57743634,rs17420384,rs2295981,rs7924137,rs1570038,rs57743634,rs2685252

There are few configuration options for the workflow. These can be adjusted 
in the (configs/workflow.json)[configs/workflow.json] file. 
Here is the list of options:

* hots_thres - Threshold for the number of trans entities to define a hotspot
* ppi_db - Either 'string', 'biogrid' or 'biogrid_stringent'
* suffix_tfa_expr - Whether to use expression or TF activities for TFs
* eqtl_prior_type - eQTL prior to be used, either 'eqtlgen' or 'gtex'

There are several meta targets in the workflow to obtain intermediate results.
Call `snakemake <target-name>` to execute the specific part of the workflow.

* all_ranges - The call below creates all hotspots networks as 'ranges' objects/files and generates
a summary plot.
* all_data - This call collects and normalized all cohort data for the created ranges objects and
generates an overview plot.
* all_ggm - This generates all ranges and data collections and fits different network models to the data.
* plot_relevant_graphs - Generates *dot* based plots for some of the networks with desireable properties.
* all_cohort - Run the full cohort data analysis.
* all_simulation - Run full simulation study (see also below)

An example for a call of the Snakemake pipeline using SLURM cluster execution is given below (some cluster
settings can be set in (configs/slurm.json)[configs/slurm.json]).

```{bash}
nohup nice snakemake -u configs/slurm.json --jobs=200 -k --local-cores=1 \
    --latency-wait 40 --cluster "sbatch --nice=10000 -t {cluster.time} \
    -c {cluster.cpu} --mem-per-cpu {cluster.mem} -p {cluster.partition}" \
    ${target} &> ${log} &
```

### Executing the simulation study
We implemented a simulation study, were we generate ground truth graphs as well
as noisy priors in order to compare inference methods with respect to 1) their 
general performance to recover the ground truth network 2) how much they are 
influenced by the exisiting prior information and noise therein and 3) the effect of
sample size on inference performance. 
As we simulate data such that they match the properties of the original data (genotype frequencies,
sample size and number of nodes per locus), this is coupled to the original data 
workflow. Simulation specific rules can be found under [workflows/2_2_simulation.sm](workflows/2_2_simulation.sm) 
in the repository. The target rule for this simulation study is `all_simulation`.

## Charliecloud container
We use a charliecloud container instead of Docker/Singularity for our computations.
To this end, we extended *Snakemake* to execute all RScript calls within a specific container.
The [Dockerfile](Dockerfile) for this container is provided in this repository.
Generally, if e.g. Docker or Singularity is available on the system the workflow is executed, a standard 
installation of *Snakemake* can be used to exeucted the workflow in the respective environment.
Otherwise, one needs to modify the *script.py* script of *Snakemake* to execute all scripts within the
provided container.

## Get the randwom_walk data:
Currently we are using some precalculated information which we at the moment
cannot reproduce in full using this specific workflow (as the data originate from the original meQTL
project from which we obtained the trans-meQTL hotspots). While this is not strictly
necessary, we use these data to be 100% in line with the previous study.
the data is contained [this file](rw_string_v9_ld_wb_prioritize_full_with_empirical_p_lte_0.05.txt)
