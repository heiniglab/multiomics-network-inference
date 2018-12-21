# general params
COHORTS = ["lolipop", "kora"] # available cohorts

# define the available PPI networks and set the active one
PPI_DB_BIOGRID = "results/current/biogrid.rds"
PPI_DB_STRING = "results/current/string.v9.expr.rds"

# -----------------------------------------------------------------------------
# most files depend on the type of PPI db used.
# so we define a nice name in accordance to the used
# DB to be added to our directory definitions below
PPI_NAME = config["ppi_db"]
if PPI_NAME == "string":
  PPI_DB = PPI_DB_STRING
else:
  PPI_DB = PPI_DB_BIOGRID
# -----------------------------------------------------------------------------

# output directories
DCOHORT_VAL = "results/current/" + PPI_NAME + "/validation/"
DPRIORS = "results/current/" + PPI_NAME + "/priors/"
DRANGES = "results/current/" + PPI_NAME + "/ranges/"
DCOHORT_DATA = "results/current/" + PPI_NAME + "/cohort_data/"
DCOHORT_FITS = "results/current/" + PPI_NAME + "/fits/"
DMEDIATION = "results/current/" + PPI_NAME + "/mediation/"

# simulation study specific output directories
DSIM_DATA = "results/current/" + PPI_NAME + "/simulation/data/"
DSIM_VALIDATION = "results/current/" + PPI_NAME + "/simulation/validation/"
DSIM_FITS = "results/current/" + PPI_NAME + "/simulation/fits/"

# Helper method to get list of eQTL gen hotspot snps
def get_eqtlgen_hotspots(file):
        out = {}
        slist = []
        with open(file, "r") as f:
                for line in f:
                        slist.append(line.strip())
        out["sentinel"] = slist
        return(out)

# get the loci available in eQTLgen
# this is how we got ALL loci. it turned out that for some we currently
# dont have all trans.genes in our annotation, which resulted in fewer
# ranges. we then defined a new set based only on these ranges
EQTLGEN_ALL = get_eqtlgen_hotspots("data/current/eqtl_gen/trans_hotspots.txt")
EQTLGEN = glob_wildcards(DRANGES + "{sentinel}_eqtlgen.rds")

# get the meQTL loci
MEQTL = glob_wildcards("data/current/sentinels/{sentinel}.dummy")

