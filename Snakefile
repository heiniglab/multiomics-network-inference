#------------------------------------------------------------------------------
# Master Snakefile for the bayesian GGM project.
# So far we have essentially to sub-workflows: application of ggm on cohort data
# and a simulation study. Rules specific to these workflows are in separte
# snakemake files (under snakemake_rules/).
#
# @author: Johann Hawe <johann.hawe@helmholtz-muenchen.de>
# ------------------------------------------------------------------------------

configfile: "configs/workflow.json"
include: "workflows/common.sm"

# set global wildcard constraints
wildcard_constraints:
    sentinel="rs\d+",
    seed="eqtlgen|meqtl"

rule test_tfa:
	conda:
		"envs/bioR.yaml"
	threads: 12
	log:
		"logs/test_tfa.out"
	resources:
		mem_mb=3000
	script:
		"scripts/test_tfa_inference.R"

# -----------------------------------------------------------------------------
# include the hotspot extraction workflow which extracts all the sentinels
# for eqtlgen and meqtl as well as performs additional data prep
# -----------------------------------------------------------------------------
subworkflow preprocess:
    workdir:
        "./"
    snakefile:
        "workflows/1_preprocess.sm"
    configfile:
        "./configs/workflow.json"

# get the loci available in eQTLgen
# the preprocess() function ensures the dependency on the subworkflow, 
EQTLGEN = glob_wildcards(preprocess(DHOTSPOTS + "eqtlgen_thres" + 
                         config["hots_thres"] + "/loci/{sentinel}.dmy"))

# same for meQTL loci
MEQTL = glob_wildcards(preprocess(DHOTSPOTS + "meqtl_thres" + 
                         config["hots_thres"] + "/loci/{sentinel}.dmy"))

# -----------------------------------------------------------------------------
# Define rules which should only be executed locally
# -----------------------------------------------------------------------------
localrules:
        all, preprocess_kora_individuals, summarize_validation_meqtl,
        all_simulation, all_cohort,
        summarize_simulation, render_validation, summarize_validation_eqtlgen,
        create_stringdb, create_cosmo_splits, convert_cpg_context,
	all_priors, all_data, create_biogrid

# ------------------------------------------------------------------------------
# Include the rule-sets for the two individual analyses (cohort, simulation)
# ------------------------------------------------------------------------------
include: "workflows/2_1_cohort_data.sm"
include: "workflows/2_2_simulation.sm"

#-------------------------------------------------------------------------------
# Overall target rule (both simulation and cohort study are executed).
# Note: this runs a virtually forever and should definitely be executed on a cluster
#-------------------------------------------------------------------------------
rule all:
        input:
                DCOHORT_VAL + "stat_overview_meqtl.pdf",
                DRANGES + "meqtl_summary.pdf",
                DCOHORT_VAL + "stat_overview_eqtlgen.pdf",
                DRANGES + "eqtlgen_summary.pdf",
                "results/current/simulation/simulation.html"

#-------------------------------------------------------------------------------
# General rules used in both simulation and cohort studies
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Convert the precalculated cpg context to RDS
# TODO: calculate the actual context on our own
#-------------------------------------------------------------------------------
rule convert_cpg_context:
	input:
		"data/current/cpgs_with_chipseq_context_100.RData"
	output:
		"results/current/cpg_context.rds"
	conda:
		"envs/bioR.yaml"
	script:
		"scripts/convert_cpg_context.R"

#-------------------------------------------------------------------------------
# Create the GTeX based prior information
# The 45gb mem requirement is not a joke, unfortunately.
# we should try and reduce this somehow...
#-------------------------------------------------------------------------------
rule create_priors:
	input:	
		eqtl="data/current/gtex/Whole_Blood_Analysis.v6p.all_snpgene_pairs.txt.gz",
		snpinfo="data/current/gtex/GTEx_Analysis_v6_OMNI_genot_1KG_imputed_var_chr1to22_info4_maf01_CR95_CHR_POSb37_ID_REF_ALT.txt",
		expr="data/current/gtex/GTEx_Analysis_v6_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct.gz",
		sampleinfo="data/current/gtex/GTEx_Data_V6_Annotations_SampleAttributesDS.txt",
		pheno="data/current/gtex/GTEx_Data_V6_Annotations_SubjectPhenotypesDS.txt",
		ppi=PPI_DB
	output:	
		gene_priors=protected("results/current/" + PPI_NAME + "/gtex.gg.cors.rds"),
		eqtl_priors=protected("results/current/" + PPI_NAME + "/gtex.eqtl.priors.rds")
	threads: 1
	params:
		plot_dir = "results/current/" + PPI_NAME + "/plots/",
		time = "16:00:00"
	resources:
		mem_mb=55000
	conda:
		"envs/bioR.yaml"
	log:
		"logs/create_priors.log"
	benchmark:
		"benchmarks/create_priors.bmk"
	script:
		"scripts/create_priors.R"

#------------------------------------------------------------------------------
# eQTLgen based SNP-gene priors. Currently used for meQTL analysis instead of
# the GTEx based priors
#------------------------------------------------------------------------------
rule create_priors_eqtlgen:
	input:	
		eqtl="data/current/eqtl_gen/cis-eQTLs_full_20180905.txt.gz",
	output:	
		eqtl_priors=protected("results/current/" + PPI_NAME + "/eqtlgen_eqtl_priors.rds")
	threads: 1
	params:
		plot_dir = "results/current/" + PPI_NAME + "/plots/",
		time = "03:00:00"
	resources:
		mem_mb=55000
	conda:
		"envs/bioR.yaml"
	log:
		"logs/create_priors_eqtlgen.log"
	benchmark:
		"benchmarks/create_priors_eqtlgen.bmk"
	script:
		"scripts/create_priors_eqtlgen.R"

#------------------------------------------------------------------------------
# Preprocess PPI network
#------------------------------------------------------------------------------
rule create_ppi_db:
	input:
		string="data/current/string/human_gene_hgnc_symbol.links.detailed.v9.0.txt",
		biogrid="data/current/biogrid/3.5.166/by_organism/BIOGRID-ORGANISM-Homo_sapiens-3.5.166.tab2.txt",
		gtex="data/current/gtex/GTEx_Analysis_v6_RNA-seq_RNA-SeQCv1.1.8_gene_median_rpkm.gct",
		gene_annot = GENE_ANNOT
	output:
		PPI_DB
	params:
		ppi_name = config["ppi_db"]
	conda:
		"envs/bioR.yaml"
	script:
		"scripts/create_ppi_db.R"

#------------------------------------------------------------------------------
# Split the cosmo object in cis, longrange and trans
#------------------------------------------------------------------------------
rule create_cosmo_splits:
	input: 
		cosmo="data/current/meQTLs/cosmopairs_combined_151216.RData"
	output:
		trans="results/current/trans-cosmopairs_combined_151216.rds",
		cis="results/current/cis-cosmopairs_combined_151216.rds",
		longrange="results/current/longrange-cosmopairs_combined_151216.rds"
	conda:
		"envs/bioR.yaml"
	script:
		"scripts/create-cosmo-splits.R"

#------------------------------------------------------------------------------
# Collect range information per sentinel snp (snp genes, cpg genes, etc)
#------------------------------------------------------------------------------
rule collect_ranges:
	input: 
		cpgcontext="results/current/cpg_context.rds",
		ppi_db=PPI_DB,
		meqtl="data/current/meQTLs/transpairs_r02_110117_converted_1MB.txt",
		tcosmo="results/current/trans-cosmopairs_combined_151216.rds",
		priorization="data/current/rw_string_v9_ld_wb_prioritize_full_with_empirical_p_lte_0.05_bck.txt",
		gene_annot = GENE_ANNOT
	output: 
		DRANGES + "{sentinel}_meqtl.rds"
	threads: 1
	resources:
		mem_mb=3000
	params:
		time="02:00:00"
	conda:
		"envs/bioR.yaml"
	log: 
		"logs/collect_ranges/{sentinel}_meqtl.log"
	benchmark: 
		"benchmarks/collect_ranges/{sentinel}_meqtl.bmk"
	script:
		"scripts/collect_ranges.R"

#------------------------------------------------------------------------------
# Collect range information per sentinel snp (snp genes, connecting genes, 
# eqtl genes etc)
#------------------------------------------------------------------------------
rule collect_ranges_eqtlgen:
        input:
                eqtl="data/current/eqtl_gen/trans-eQTL_significant_20181017.txt",
                ppi=PPI_DB,
                tfbs_annot="results/current/tfbs_tss_annot.rds",
                gene_annot = GENE_ANNOT
        output:
                ranges = DRANGES + "{sentinel}_eqtlgen.rds",
                plot = DRANGES + "{sentinel}_eqtlgen.pdf"
        threads: 1
        resources:
                mem_mb=3000
        params:
                time="02:00:00"
        conda:
                "envs/bioR.yaml"
        log:
                "logs/collect_ranges/{sentinel}_eqtlgen.log"
        benchmark:
                "benchmarks/collect_ranges/{sentinel}_eqtlgen.bmk"
        script:
                "scripts/collect_ranges_eqtl.R"

# -----------------------------------------------------------------------------
# Target rule to generate all hotspot ranges collections for eqtlgen and to 
# create a summary plot.
# -----------------------------------------------------------------------------
rule all_ranges:
	input:
		expand(DRANGES + "{sentinel}_eqtlgen.rds", sentinel=EQTLGEN.sentinel),
		expand(DRANGES + "{sentinel}_meqtl.rds", sentinel=MEQTL.sentinel)
	output:
		DRANGES + "summary.pdf"
	threads: 1
	resources:
		mem_mb=1000
	params:
		time="00:10:00"
	script:
		"scripts/create_locus_summary.R"

# -----------------------------------------------------------------------------
# Annotate TSS with TFBS information
# -----------------------------------------------------------------------------
rule annotate_tss_with_tf:
        input:
                tfbs_remap="data/current/tfbs/filPeaks_public.bed",
                tfbs_encode="data/current/tfbs/wgEncodeRegTfbsClusteredWithCellsV3.bed",
                gene_annot = GENE_ANNOT
        output:
                tfbs_annot="results/current/tfbs_tss_annot.rds"
        conda:
                "envs/bioR.yaml"
        script:
                "scripts/annotate_tss_with_tf.R"

# -----------------------------------------------------------------------------
# Calculate the TF activities. Results are data matrices containing residual
# corrected expression/TFA values for all available genes
# -----------------------------------------------------------------------------
rule calculate_tfa:
	input:
		cohort_data="results/current/ggmdata_{cohort}.RData",
		tfbs_annot="results/current/tfbs_tss_annot.rds"
	output:
		plot="results/current/tfa/plot_{cohort}.pdf",
		tfa="results/current/tfa/activities_{cohort}.rds",
		expr="results/current/tfa/expression_{cohort}.rds",
	resources:
		mem_mb=10000
	params:
		time="02:00:00"
	log:
		"logs/calculate_tfa_{cohort}.log"
	benchmark:
		"benchmarks/calculate_tfa_{cohort}.bmk"
	script:
		"scripts/calculate_tfa.R"

#------------------------------------------------------------------------------
# Collect cohort data for a single sentinel locus
#------------------------------------------------------------------------------
rule collect_data:
	input: 
		ranges=DRANGES + "{sentinel}_{seed}.rds",
		kora=preprocess("results/current/ggmdata_kora.RData"),
		lolipop=preprocess("results/current/ggmdata_lolipop.RData"),
		kora_activities="results/current/tfa/activities_kora.rds",
		lolipop_activities="results/current/tfa/activities_lolipop.rds",
		ceqtl="data/current/kora/eqtl/kora-cis-eqtls.csv",
		ccosmo="results/current/cis-cosmopairs_combined_151216.rds"
	output:
		DCOHORT_DATA + "{cohort}/{sentinel}_{seed}.rds",
		DCOHORT_DATA + "{cohort}/{sentinel}_raw_{seed}.rds"
	threads: 1
	resources:
		mem_mb=25000
	conda:
		"envs/bioR.yaml"
	params:
		time="01:00:00",
		tfa_or_expr=config["suffix_tfa_expr"]
	log:
		"logs/collect_data/{cohort}/{sentinel}_{seed}.log"
	benchmark:
		"benchmarks/collect_data/{cohort}/{sentinel}_{seed}.bmk"
	script:
		"scripts/collect_data_tfa.R"

#------------------------------------------------------------------------------
# Meta rule to collect data for all sentinels and plot some information. 
# For eQTLgen we currently only have the KORA data available.
#------------------------------------------------------------------------------
rule all_data:
        input:
                expand(DCOHORT_DATA + "{cohort}/{sentinel}_meqtl.rds", 
                             sentinel=MEQTL.sentinel, cohort=COHORTS),
                expand(DCOHORT_DATA + "{cohort}/{sentinel}_eqtlgen.rds", 
                               sentinel=EQTLGEN.sentinel, cohort=COHORTS)
        output:
                DCOHORT_DATA + "summary.pdf"
        script:
                "scripts/create_data_summary.R"

#------------------------------------------------------------------------------
# Collect prior information for a sentinel locus
#------------------------------------------------------------------------------
rule collect_priors:
	input:
		gg_priors="results/current/" + PPI_NAME + "/gtex.gg.cors.rds", 
		gtex_eqtl_priors="results/current/" + PPI_NAME + "/gtex.eqtl.priors.rds",
		eqtlgen_eqtl_priors="results/current/" + PPI_NAME + "/eqtlgen_eqtl_priors.rds",
		ranges=DRANGES + "{sentinel}_{seed}.rds",
		ppi=PPI_DB,
		cpg_context="data/current/cpgs_with_chipseq_context_100.rds",
		cpg_annot="data/current/epigenetic_state_annotation_weighted_all_sentinels.txt"
	output: 
		DPRIORS + "{sentinel}_{seed}.rds",
		DPRIORS + "{sentinel}_{seed}.pdf"
	threads: 1
	resources:
		mem_mb=20000
	conda:
		"envs/bioR.yaml"
	params:
		time="02:00:00",
		eqtl_prior_type = config["eqtl_prior_type"]
	log:
		"logs/collect_priors/{sentinel}_{seed}.log"
	benchmark:
		"benchmarks/collect_priors/{sentinel}_{seed}.bmk"
	script:
		"scripts/collect_priors.R"

#------------------------------------------------------------------------------
# Meta rule to get priors for all sentinels
#------------------------------------------------------------------------------
rule all_priors:
	input: 
		expand(DPRIORS + "{sentinel}_meqtl.pdf", sentinel=MEQTL.sentinel),
		expand(DPRIORS + "{sentinel}_eqtlgen.pdf", sentinel=EQTLGEN.sentinel)
