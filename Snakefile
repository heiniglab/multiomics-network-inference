LISTS = glob_wildcards("data/current/sentinels/{sentinel}.dummy")
COHORTS = ["lolipop", "kora"]


localrules: 
	all, preprocess, summarize_validation, 
	all_sim, validate_ggm_simulation, create_priors,
	summarize_simulation, render_validation,
	create_stringdb, create_cosmo_pairs, all_ranges, all_data

###############################################################################
####### General rules used in both simulation and cohort studies ##############
###############################################################################

#------------------------------------------------------------------------------
# Create the GTeX based prior information
# The 45gb mem requirement is not a joke, unfortunately.
# we should try and reduce this somehow...
#------------------------------------------------------------------------------
rule create_priors:
	input:	
		eqtl="data/current/gtex/Whole_Blood_Analysis.v6p.all_snpgene_pairs.txt.gz",
		snpinfo="data/current/gtex/GTEx_Analysis_v6_OMNI_genot_1KG_imputed_var_chr1to22_info4_maf01_CR95_CHR_POSb37_ID_REF_ALT.txt",
		expr="data/current/gtex/GTEx_Analysis_v6_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct.gz",
		sampleinfo="data/current/gtex/GTEx_Data_V6_Annotations_SampleAttributesDS.txt",
		pheno="data/current/gtex/GTEx_Data_V6_Annotations_SubjectPhenotypesDS.txt"
	output:	
		gene=protected("results/current/gtex.gg.cors.rds"),
		eqtl=protected("results/current/gtex.eqtl.priors.rds")
	threads: 1
	params:
		plot_dir = "results/current/plots/"
	log:
		"logs/create-priors.log"
	benchmark:
		"benchmarks/create-priors.bmk"
	resources:
		mem_mb=45000
	script:
		"R/create-priors.R"

#------------------------------------------------------------------------------
# Preprocess stringdb PPI network
#------------------------------------------------------------------------------
rule create_stringdb:
	input:
		string="data/current/string/human_gene_hgnc_symbol.links.detailed.v9.0.txt",
		gtex="data/current/gtex/GTEx_Analysis_v6_RNA-seq_RNA-SeQCv1.1.8_gene_median_rpkm.gct"
	output:
		"results/current/string.v9.expr.rds"
	script:
		"scripts/create-stringdb.R"

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
# Creates a list of cis-eqtl snp locations for the gtex and kora eqtls
#------------------------------------------------------------------------------
rule get_snp_locations:
	input: 
		"data/current/kora/KORAF4-cis-eqtls_schramm_journal.pone.0093844.s005.csv",
		"data/current/gtex/GTEx_Analysis_v6_OMNI_genot_1KG_imputed_var_chr1to22_info4_maf01_CR95_CHR_POSb37_ID_REF_ALT.txt"
	output:
		"results/current/cis-eqtl-snp-locations.txt"
	shell:
		"""
		cut -f 1 -d ";" data/current/kora/KORAF4-cis-eqtls_schramm_journal.pone.0093844.s005.csv  | sort | uniq > results/current/cis-eqtl-snps.txt
	        cut -f 1,2,7 data/current/gtex/GTEx_Analysis_v6_OMNI_genot_1KG_imputed_var_chr1to22_info4_maf01_CR95_CHR_POSb37_ID_REF_ALT.txt > results/current/gtex-snp-locations.txt
	        fgrep -w -f results/current/cis-eqtl-snps.txt results/current/gtex-snp-locations.txt > results/current/cis-eqtl-snp-locations.txt
		rm results/current/cis-eqtl-snps.txt results/current/gtex-snp-locations.txt
		"""

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
		ceqtl_locations="results/current/cis-eqtl-snp-locations.txt",
		ccosmo="results/current/cis-cosmopairs_combined_151216.rds"
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
		
###############################################################################
###################### Rules for cohort GGM study #############################
###############################################################################

#------------------------------------------------------------------------------
# Collect range information per sentinel snp (snp genes, cpg genes, etc)
#------------------------------------------------------------------------------
rule collect_ranges:
	input: 
		cpgcontext="data/current/cpgs_with_chipseq_context_100.RData",
		string="results/current/string.v9.expr.rds",
		meqtl="data/current/meQTLs/transpairs_r02_110117_converted_1MB.txt",
		tcosmo="results/current/trans-cosmopairs_combined_151216.rds"
	output: "results/current/ranges/{sentinel}.rds"
	log: 
		"logs/collect-ranges/{sentinel}.log"
	benchmark: 
		"benchmarks/collect-ranges/{sentinel}.bmk"
	threads: 1
	resources:
		mem_mb=2300
	params:
		sentinel="{sentinel}",
		window=1e6
	script:
		"scripts/collect-ranges.R"

#------------------------------------------------------------------------------
# Pseudo rule to collect ranges for all sentinels
#------------------------------------------------------------------------------
rule all_ranges:
	input:
		expand("results/current/ranges/{sentinel}.rds", zip, sentinel=LISTS.sentinel)

#------------------------------------------------------------------------------
# Collect cohort data for a single sentinel locus
#------------------------------------------------------------------------------
rule collect_data:
	input: 
		ranges="results/current/ranges/{sentinel}.rds",
		kora="results/current/ggmdata_kora.RData",
		lolipop="results/current/ggmdata_lolipop.RData",
		ceqtl="data/current/kora/KORAF4-cis-eqtls_schramm_journal.pone.0093844.s005.csv",
		ccosmo="results/current/cis-cosmopairs_combined_151216.rds"
	output:
		"results/current/cohort-data/{cohort}/{sentinel}.rds"
	threads: 1
	params:
		cohort="{cohort}",
		sentinel="{sentinel}"
	log:
		"logs/collect-data/{cohort}/{sentinel}.log"
	benchmark:
		"benchmarks/collect-data/{cohort}/{sentinel}.log"
	resources:
		mem_mb=20000
	script:
		"scripts/collect-data.R"

#------------------------------------------------------------------------------
# Collect prior information for a sentinel locus
#------------------------------------------------------------------------------
rule collect_priors:
	input:
		gg_priors="results/current/gtex.gg.cors.rds", 
		eqtl_priors="results/current/gtex.eqtl.priors.rds",
		ranges="results/current/ranges/{sentinel}.rds",
		string="results/current/string.v9.expr.rds"
	output: 
		"results/current/priors/{sentinel}.rds"
	log:
		"logs/collect-priors/{sentinel}.log"
	threads: 1
	params:
		sentinel="{sentinel}",
		plot_file="results/current/priors/{sentinel}.pdf"
	benchmark:
		"benchmarks/collect-priors/{sentinel}.bmk"
	resources:
		mem_mb=9000
	script:
		"scripts/collect-priors.R"

#------------------------------------------------------------------------------
# Apply ggm on collected data and priors for a sentinel
#------------------------------------------------------------------------------
rule apply_ggm:
	input: data="results/current/cohort-data/{cohort}/{sentinel}.rds", 
		priors="results/current/priors/{sentinel}.rds",
		ranges="results/current/ranges/{sentinel}.rds"
	output: 
		"results/current/fits/{cohort}/{sentinel}.rds"
	threads: 10
	params:
		nriter=20000,
		burnin=5000,
		plot_file="results/current/fits/{cohort}/{sentinel}.pdf"
	resources:
		mem_mb=1000
	benchmark:
		"benchmarks/apply-ggm/{sentinel}-{cohort}.bmk"
	log:
		"logs/apply-ggm/{sentinel}-{cohort}.log"
	script:
		"scripts/apply-ggm.R"

#------------------------------------------------------------------------------
# Validate calculated ggms
#------------------------------------------------------------------------------
rule validate_ggm:
	input: 
		expand("results/current/fits/{{sentinel}}-{cohort}.rds", cohort=COHORTS)
	output: 
		"results/current/validation/{sentinel}.txt"
	log: 
		"logs/validate-ggm/{sentinel}.log"
	threads: 1
	benchmark:
		"benchmarks/validate-ggm/{sentinel}.bmk"
	resources:
		mem_mb=3000
	script:
		"R/validate-ggm.R"

#------------------------------------------------------------------------------
# Collect validation information in single file
#------------------------------------------------------------------------------
rule summarize_validation:
	input: expand("results/current/validation/{sentinel}.txt", zip, sentinel=LISTS.sentinel)
	output: "results/current/validation/validation.all.txt"
	shell:
		"cat {input} | sort -r | uniq > {output}"

#------------------------------------------------------------------------------
# Create summary report
#------------------------------------------------------------------------------
rule render_validation:
	input: "results/current/validation/validation.all.txt"
	output: "results/current/validation.html"
	log:	
		"logs/validation.log"
	script:
		"R/render-validation.Rmd"

#------------------------------------------------------------------------------
# Target rule for complete cohort study
#------------------------------------------------------------------------------
rule all:
	input: "results/current/validation.html"


###############################################################################
###################### Rules for simulation study #############################
###############################################################################

#------------------------------------------------------------------------------
# Simulate ground truth and data for simulation study
#------------------------------------------------------------------------------
RUNS = 100
rule simulate_data:
	input: 
		data="results/current/cohort-data/lolipop/{sentinel}.rds",
		ranges="results/current/ranges/{sentinel}.rds",
		priors="results/current/priors/{sentinel}.rds"
	output: 
		"results/current/simulation/data/{sentinel}.RData"
	threads: 10
	resources:
		mem_mb=1200
	params:
		sentinel="{sentinel}",
		runs=RUNS
	log:	
		"logs/simulation/simulate-data/{sentinel}.log"
	benchmark: 
		"benchmarks/simulation/simulate-data/{sentinel}.bmk"
	script:
		"scripts/simulate-data.R"

#------------------------------------------------------------------------------
# Apply ggm on simulated data
#------------------------------------------------------------------------------
ITERATIONS = range(1,RUNS+1)

rule apply_ggm_simulation:
	input: 
		"results/current/simulation/data/{sentinel}.RData",
	output: 
		"results/current/simulation/fits/{sentinel}-iter{iteration}.RData"
	params:
		nriter=20000,
		burnin=5000,
		iteration="{iteration}"
	threads: 16
        resources:
                mem_mb=3000
	log: 
		"logs/simulation/apply-ggm/{sentinel}-iter{iteration}.log"
	benchmark: 
		"benchmarks/simulation/apply-ggm/{sentinel}-iter{iteration}.bmk"
	script:
		"scripts/apply-ggm-simulation.R"

#------------------------------------------------------------------------------
# Validate a simulation run
#------------------------------------------------------------------------------
rule validate_ggm_simulation:
	input: 
		"results/current/simulation/fits/{sentinel}-iter{iteration}.RData"
	output: 
		"results/current/simulation/validation/{sentinel}-iter{iteration}.txt"
	log: 
		"logs/simulation/validate-ggm/{sentinel}-iter{iteration}.log"
	script:
		"scripts/validate-ggm-simulation.R"


#------------------------------------------------------------------------------
# Target rule to validate all simulation runs
#------------------------------------------------------------------------------
rule validate_all:
	input: expand("results/current/simulation/validation/{sentinel}-iter{iter}.txt", sentinel=LISTS.sentinel, iter=ITERATIONS)

#------------------------------------------------------------------------------
# Create simulation report
#------------------------------------------------------------------------------
rule summarize_simulation:
	input: 
		expand("results/current/simulation/validation/{sentinel}-iter{iter}.txt", zip, sentinel=LISTS.sentinel, iter=ITERATIONS)
	params:
		indir="../results/current/simulation/validation/"
	output: 
		"results/current/simulation/simulation.html"
	log: 
		"logs/simulation/summarize-simulation.log"
	script:
		"scripts/summarize-simulation.Rmd"

#------------------------------------------------------------------------------
# Target rule for complete simulation study
#------------------------------------------------------------------------------
rule all_sim:
	input: "results/current/simulation/simulation.html"
