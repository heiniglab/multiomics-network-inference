	
#------------------------------------------------------------------------------
# Master Snakefile for the bayesian GGM project.
# So far we have essentially to sub-workflows: application of ggm on cohort data
# and a simulation study. Rules specific to these workflows are in separte
# snakemake files (under snakemake_rules/).
#
# @author: Johann Hawe <johann.hawe@helmholtz-muenchen.de>
# ------------------------------------------------------------------------------

configfile: "config.json"

# -----------------------------------------------------------------------------
# incude the hotspot extraction workflow which extracts all the sentinels
# for eqtlgen and meqtl
# -----------------------------------------------------------------------------
subworkflow hotspots:
    workdir:
        "./"
    snakefile:
        "workflows/1_extract_hotspots.sm"
    configfile:
        "./config.json"

# -----------------------------------------------------------------------------
# Insert global vars
# -----------------------------------------------------------------------------
include: "workflows/common.sm"

# set global wildcard constraints
wildcard_constraints:
    sentinel="rs\d+",
    seed="eqtlgen|meqtl"

# the hotspots() funciton ensures the dependency on the subworkflow, 

# get the loci available in eQTLgen
EQTLGEN = glob_wildcards(hotspots(DHOTSPOTS + "eqtlgen_thres" + 
                         config["hots_thres"] + "/loci/{sentinel}.dmy"))

# get the meQTL loci
MEQTL = glob_wildcards(hotspots(DHOTSPOTS + "meqtl_thres" + 
                         config["hots_thres"] + "/loci/{sentinel}.dmy"))

# -----------------------------------------------------------------------------
# Define rules which should only be executed locally
# -----------------------------------------------------------------------------
localrules:
        all, preprocess_kora_individuals, summarize_validation_meqtl,
        all_simulation, all_cohort, validate_ggm_simulation, create_priors,
        summarize_simulation, render_validation, summarize_validation_eqtlgen,
        create_stringdb, create_cosmo_splits, convert_cpg_context, all_ranges,
	all_priors, all_data

# ------------------------------------------------------------------------------
# Include the rule-sets for the two individual analyses (cohort, simulation)
# ------------------------------------------------------------------------------
include: "snakemake_rules/cohort_data.sm"
include: "snakemake_rules/simulation.sm"

#-------------------------------------------------------------------------------
# Overall target rule (both simulation and cohort study are executed).
# Note: this runs a long time and should definitely be executed on a cluster
#-------------------------------------------------------------------------------
rule all:
        input:
                DCOHORT_VAL + "stat_overview_meqtl.pdf",
                DRANGES + "overview_meqtl.pdf",
                DCOHORT_VAL + "stat_overview_eqtlgen.pdf",
                DRANGES + "overview_eqtlgen.pdf",
                "results/current/simulation/simulation.html"

################################################################################
# General rules used in both simulation and cohort studies
################################################################################

#-------------------------------------------------------------------------------
# Convert the precalculated cpg context to RDS
# TODO: calculate it on our own
#-------------------------------------------------------------------------------
rule convert_cpg_context:
	input:
		"data/current/cpgs_with_chipseq_context_100.RData"
	output:
		"results/current/cpg_context.rds"
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
		gene_priors=protected("results/current/gtex.gg.cors.rds"),
		eqtl_priors=protected("results/current/gtex.eqtl.priors.rds")
	threads: 1
	params:
		plot_dir = "results/current/plots/"
	log:
		"logs/create_priors.log"
	benchmark:
		"benchmarks/create_priors.bmk"
	resources:
		mem_mb=45000
	script:
		"scripts/create_priors.R"

#------------------------------------------------------------------------------
# Preprocess stringdb PPI network
#------------------------------------------------------------------------------
rule create_stringdb:
	input:
		string="data/current/string/human_gene_hgnc_symbol.links.detailed.v9.0.txt",
		gtex="data/current/gtex/GTEx_Analysis_v6_RNA-seq_RNA-SeQCv1.1.8_gene_median_rpkm.gct"
	output:
		PPI_DB_STRING
	log:
		"logs/create_stringdb.log"
	benchmark:
		"benchmarks/create_stringdb.bmk"
	script:
		"scripts/create-stringdb.R"


#------------------------------------------------------------------------------
# Preprocess biogrid PPI network
#------------------------------------------------------------------------------
rule create_biogrid:
	input:
		biogrid="data/current/biogrid/3.5.166/by_organism/BIOGRID-ORGANISM-Homo_sapiens-3.5.166.tab2.txt",
		gtex="data/current/gtex/GTEx_Analysis_v6_RNA-seq_RNA-SeQCv1.1.8_gene_median_rpkm.gct"
	output:
		PPI_DB_BIOGRID
	script:
		"scripts/create_biogrid.R"

#------------------------------------------------------------------------------
# Preprocess stringdb PPI network
#------------------------------------------------------------------------------
rule create_cosmo_splits:
	input: 
		cosmo="data/current/meQTLs/cosmopairs_combined_151216.RData"
	output:
		trans="results/current/trans-cosmopairs_combined_151216.rds",
		cis="results/current/cis-cosmopairs_combined_151216.rds",
		longrange="results/current/longrange-cosmopairs_combined_151216.rds"
	script:
		"scripts/create-cosmo-splits.R"

#------------------------------------------------------------------------------
# Preprocess the sample-mapping sheet for kora
#------------------------------------------------------------------------------
rule preprocess_kora_individuals:
	input:
		"data/current/kora/individuals.csv"
	output:
		"results/current/kora_individuals.csv"
	shell:
		"""
		awk 'BEGIN {{ FS = ";" }} ; {{ if ($1 != "" && $5 != "" && $6 != "" ) print }}' {input} | sed s/zz_nr_//g > {output}
		"""

#------------------------------------------------------------------------------
# Gather kora data in single file for more convenient postprocessing
#------------------------------------------------------------------------------
rule prepare_kora_data:
	input:
		genotypes="data/current/kora/snp_dosages/MAF001/full_sorted.bgz",
		expression="data/current/kora/expression/kora_f4_normalized.Rdata",
		expression_cov="data/current/kora/expression/technical_covariables_kora_f4.Rdata",
		methylation="data/current/kora/methylation/KF4_beta_qn_bmiq.RData",
		methylation_cov="data/current/kora/methylation/control_probe_pcs_n1727.RData",
		individuals="results/current/kora_individuals.csv",
		impute_indiv="data/current/kora/imputation_individuals",
		trans_meqtl="data/current/meQTLs/transpairs_r02_110117_converted_1MB.txt",
		houseman="data/current/kora/methylation/Houseman/KF4_QN_estimated_cell_distribution_meanimpute473_lessThanOneTRUE.csv",
		kora_ceqtl="data/current/kora/eqtl/kora-cis-eqtls.csv",
		ccosmo="results/current/cis-cosmopairs_combined_151216.rds",
		eqtl_gen="data/current/eqtl_gen/trans-eQTL_significant_20181017.txt.gz"
	output:
		"results/current/ggmdata_kora.RData"
	resources:
		mem_mb=23000
	log:
		"logs/prepare-kora-data.log"
	threads: 6
	benchmark:
		"benchmarks/prepare-kora-data.bmk"
	script:
		"scripts/prepare-kora-data.R"

#------------------------------------------------------------------------------
# Prepare lolipop data for more convenient postprocessing
#------------------------------------------------------------------------------
rule prepare_lolipop_data:
	input: 
		lolipop="data/current/meQTLs/ggmdata_201017.RData"
	output: "results/current/ggmdata_lolipop.RData"
	resources:
		mem_mb=2000
	log:
		"logs/prepare-lolipop-data.log"
	benchmark:
		"benchmarks/prepare-lolipop-data.bmk"
	script:
		"scripts/prepare-lolipop-data.R"
		
#------------------------------------------------------------------------------
# Collect range information per sentinel snp (snp genes, cpg genes, etc)
#------------------------------------------------------------------------------
rule collect_ranges:
	input: 
		cpgcontext="results/current/cpg_context.rds",
		ppi_db=PPI_DB,
		meqtl="data/current/meQTLs/transpairs_r02_110117_converted_1MB.txt",
		tcosmo="results/current/trans-cosmopairs_combined_151216.rds",
		priorization="data/current/rw_string_v9_ld_wb_prioritize_full_with_empirical_p_lte_0.05.txt"
	output: 
		DRANGES + "{sentinel}_meqtl.rds"
	log: 
		"logs/collect_ranges/{sentinel}_meqtl.log"
	benchmark: 
		"benchmarks/collect_ranges/{sentinel}_meqtl.bmk"
	threads: 1
	resources:
		mem_mb=2300
	script:
		"scripts/collect_ranges.R"

#------------------------------------------------------------------------------
# Meta rule to collect ranges for all sentinels and plot some information
#------------------------------------------------------------------------------
rule create_meqtl_locus_summary:
	input:
		expand(DRANGES + "{sentinel}_meqtl.rds", zip, sentinel=MEQTL.sentinel)
	output:
		DRANGES + "meqtl_summary.pdf"
	script:
		"scripts/create_locus_summary.R"


#------------------------------------------------------------------------------
# Collect range information per sentinel snp (snp genes, connecting genes, 
# eqtl genes etc)
#------------------------------------------------------------------------------
rule collect_ranges_eqtlgen:
        input:
                eqtl="data/current/eqtl_gen/trans-eQTL_significant_20181017.txt.gz",
                ppi=PPI_DB,
                tfbs_annot="results/current/tfbs_tss_annot.rds"
        output:
                DRANGES + "{sentinel}_eqtlgen.rds"
        log:
                "logs/collect_ranges_egen/{sentinel}.log"
        benchmark:
                "benchmarks/collect_ranges_egen/{sentinel}.bmk"
        threads: 1
        resources:
                mem_mb=2300
        script:
                "scripts/collect_ranges_eqtl.R"

# -----------------------------------------------------------------------------
# Target rule to generate all hotspot ranges collections for eqtl gen and to 
# create a summary plot.
# -----------------------------------------------------------------------------
rule create_eqtlgen_locus_summary:
	input:
		expand(DRANGES + "{sentinel}_eqtlgen.rds", sentinel=EQTLGEN.sentinel)
	output:
		DRANGES + "eqtlgen_summary.pdf"
	script:
		"scripts/create_locus_summary.R"

# -----------------------------------------------------------------------------
# Meta rule to create all ranges and their summaries
# -----------------------------------------------------------------------------
rule all_ranges:
	input:
		DRANGES + "eqtlgen_summary.pdf",
		DRANGES + "meqtl_summary.pdf"
	
# -----------------------------------------------------------------------------
# Annotate TSS with TFBS information
# -----------------------------------------------------------------------------
rule annotate_tss_with_tf:
        input:
                tfbs_remap="data/current/tfbs/filPeaks_public.bed",
                tfbs_encode="data/current/tfbs/wgEncodeRegTfbsClusteredWithCellsV3.bed"
        output:
                tfbs_annot="results/current/tfbs_tss_annot.rds"
        script:
                "scripts/annotate_tss_with_tf.R"

#------------------------------------------------------------------------------
# Collect cohort data for a single sentinel locus
#------------------------------------------------------------------------------
rule collect_data:
	input: 
		ranges=DRANGES + "{sentinel}_{seed}.rds",
		kora="results/current/ggmdata_kora.RData",
		lolipop="results/current/ggmdata_lolipop.RData",
		ceqtl="data/current/kora/eqtl/kora-cis-eqtls.csv",
		ccosmo="results/current/cis-cosmopairs_combined_151216.rds"
	output:
		DCOHORT_DATA + "{cohort}/{sentinel}_{seed}.rds",
		DCOHORT_DATA + "{cohort}/{sentinel}_raw_{seed}.rds"
	threads: 2
	params:
		cohort="{cohort}",
		sentinel="{sentinel}"
	log:
		"logs/collect_data/{cohort}/{sentinel}_{seed}.log"
	benchmark:
		"benchmarks/collect_data/{cohort}/{sentinel}_{seed}.bmk"
	resources:
		mem_mb=20000
	script:
		"scripts/collect_data.R"

#------------------------------------------------------------------------------
# Meta rule to collect data for all sentinels and plot some information. 
# For eQTLgen we currently only have the KORA data available.
#------------------------------------------------------------------------------
rule all_data:
        input:
                meqtl=expand(DCOHORT_DATA + "{cohort}/{sentinel}_meqtl.rds", 
                             sentinel=MEQTL.sentinel, cohort=COHORTS),
                eqtlgen=expand(DCOHORT_DATA + "kora/{sentinel}_eqtlgen.rds", 
                               sentinel=EQTLGEN.sentinel)
        output:
                DCOHORT_DATA + "summary.pdf"
        script:
                "scripts/create_data_summary.R"

#------------------------------------------------------------------------------
# Collect prior information for a sentinel locus
#------------------------------------------------------------------------------
rule collect_priors:
	input:
		gg_priors="results/current/gtex.gg.cors.rds", 
		eqtl_priors="results/current/gtex.eqtl.priors.rds",
		ranges=DRANGES + "{sentinel}_{seed}.rds",
		ppi=PPI_DB,
		cpg_context="data/current/cpgs_with_chipseq_context_100.RData",
		cpg_annot="data/current/epigenetic_state_annotation_weighted_all_sentinels.txt"
	output: 
		DPRIORS + "{sentinel}_{seed}.rds",
		DPRIORS + "{sentinel}_{seed}.pdf"
	log:
		"logs/collect_priors/{sentinel}_{seed}.log"
	threads: 1
	benchmark:
		"benchmarks/collect_priors/{sentinel}_{seed}.bmk"
	resources:
		mem_mb=9000
	script:
		"scripts/collect_priors.R"

#------------------------------------------------------------------------------
# Meta rule to get priors for all sentinels
#------------------------------------------------------------------------------
rule all_priors:
	input: 
		expand(DPRIORS + "{sentinel}.pdf", sentinel=MEQTL.sentinel)
