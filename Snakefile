SENT_COHO = glob_wildcards("data/current/cohorts/{sentinel}.{cohort}.adjusted.data.RData")
LISTS = glob_wildcards("data/current/sentinels/{sentinel}.dummy")


localrules: 
	all, preprocess, summarize_validation, 
	all_sim, validate_ggm_simulation, create_priors,
	summarize_simulation, render_validation, simulate_data,
	create_stringdb, create_cosmo_pairs, all_ranges


#####
# Rule used in both 'subworkflows'
# The 45gb mem requirement is not a joke, unfortunately.
# we should try and reduce this somehow...
#####
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
	log:
		"logs/create-priors.log"
	benchmark:
		"benchmarks/create-priors.bmk"
	resources:
		mem_mb=45000
	script:
		"R/create-priors.R"
		
rule create_stringdb:
	input:
		string="data/current/string/human_gene_hgnc_symbol.links.detailed.v9.0.txt",
		gtex="data/current/gtex/GTEx_Analysis_v6_RNA-seq_RNA-SeQCv1.1.8_gene_median_rpkm.gct"
	output:
		"results/current/string.v9.expr.rds"
	script:
		"scripts/create-stringdb.R"

rule create_cosmo_splits:
	input: 
		cosmo="data/current/meQTLs/cosmopairs_combined_151216.RData"
	output:
		trans="results/current/trans-cosmopairs_combined_151216.rds",
		cis="results/current/cis-cosmopairs_combined_151216.rds",
		longrange="results/current/longrange-cosmopairs_combined_151216.rds"
	script:
		"scripts/create-cosmo-splits.R"

#####
# Rules for application of GGM on real data
#####

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
		"scripts/collect-ranges.R &> {log}"


rule all_ranges:
	input:
		expand("results/current/ranges/{sentinel}.rds", zip, sentinel=LISTS.sentinel)

rule collect_data:
	input: 
		"results/current/ranges/{sentinel}.rds"
	output:
		"results/current/cohort-data/{cohort}/{sentinel}.rds"
	log:
		""
	benchmark:
		""
	resources:
		mem_mb=2000
	script:
		"scripts/collect-data.R &> {log}"

rule collect_priors:
	input:
		gg_priors="results/current/gtex.gg.cors.rds", 
		eqtl_priors="results/current/gtex.eqtl.priors.rds",
		ranges="results/current/ranges/{sentinel}.rds"
	output: 
		"results/current/priors/{sentinel}.rds"
	log:
		""
	benchmark:
		""
	resources:
		mem_mb=5000
	script:
		"scripts/collect-priors.R"

rule apply_ggm:
	input: "results/current/cohort-data/{cohort}/{sentinel}.rds", 
		"results/current/priors/{sentinel}.rds",
		"results/current/ranges/{sentinel}.rds"
	output: 
		"results/current/fits/{sentinel}-{cohort}.rds"
	threads: 10
	params:
		plotdir="results/current/plots/{sentinel}-{cohort}/",
		nriter=20000,
		burnin=5000
	resources:
		mem_mb=1000
	benchmark:
		"benchmarks/apply-ggm/{sentinel}-{cohort}.bmk"
	log:
		"logs/apply-ggm/{sentinel}.log"
	script:
		"R/apply-ggm.R"

COHORTS = [{"lolipop", "kora"}]

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

rule summarize_validation:
	input: expand("results/current/validation/{sentinel}.txt", zip, sentinel=SENT_COHO.sentinel)
	output: "results/current/validation/validation.all.txt"
	shell:
		"cat {input} | sort -r | uniq > {output}"

rule render_validation:
	input: "results/current/validation/validation.all.txt"
	output: "results/current/validation.html"
	log:	
		"logs/validation.log"
	script:
		"R/render-validation.Rmd"

rule all:
	input: "results/current/validation.html"


####################################################################
# Rules for simulation study
####################################################################

SENT = glob_wildcards("data/current/sentinels/{sentinel}.dummy")

rule simulate_data:
	input: 
		"results/current/data.processed.RData"
	output: 
		expand("results/current/simulation/data/{sentinel}.RData", zip, sentinel=SENT.sentinel),
		odir="results/current/simulation/data/"
	threads: 10
	resources:
		mem_mb=10000
	params:
		random_walk_results="data/current/networks/",
		plotdir="results/current/simulation/plots/"
	log:	
		"logs/simulation/simulate_data.log"
	benchmark: 
		"benchmarks/simulation/simulate_data.bmk"
	script:
		"R/simulate-data.R"

#rule check_simulated_graphs:
#	input: 

rule apply_ggm_simulation:
	input: "results/current/simulation/data/{sentinel}.RData",
		"results/current/gtex.gg.cors.rds",
                "results/current/gtex.eqtl.priors.rds"
	output: 
		"results/current/simulation/fits/{sentinel}.RData"
	params:
		plotdir="results/current/simulation/plots/",
		nriter=20000,
		burnin=5000
	threads: 10
        resources:
                mem_mb=1300
	log: 
		"logs/simulation/apply-ggm.{sentinel}.log"
	benchmark: "benchmarks/simulation/apply-ggm.{sentinel}.bmk"
	script:
		"R/apply-ggm-simulation.R"

rule validate_ggm_simulation:
	input: 
		"results/current/simulation/fits/{sentinel}.RData"
	output: 
		"results/current/simulation/validation/{sentinel}.txt"
	log: 
		"logs/simulation/validate-ggm.{sentinel}.log"
	script:
		"R/validate-ggm-simulation.R"

rule summarize_simulation:
	input: 
		expand("results/current/simulation/validation/{sentinel}.txt", zip, sentinel=SENT.sentinel)
	params:
		indir="../results/current/simulation/validation/"
	output: 
		"results/current/simulation/simulation.html"
	log: 
		"logs/simulation/summarize-simulation.log"
	script:
		"R/summarize-simulation.Rmd"

rule all_sim:
	input: "results/current/simulation/simulation.html"
