# general params
COHORTS = ["lolipop", "kora"] # available cohorts

# define the available PPI networks and set the active one
PPI_DB_BIOGRID = "results/current/biogrid.rds"
PPI_DB_STRING = "results/current/string.v9.expr.rds"


# -----------------------------------------------------------------------------
# most files depend on the type of PPI db used.
# so we define a nice name in accordance to the used
# DB to be added to our directory definitions below
PPI_DB = PPI_DB_STRING # the file containing the prepared PPI network
PPI_NAME = "string" # alternative: biogrid
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
        # remove some loci for which trans genes are not
        # available in our data
        slist.remove("rs3130573")
        slist.remove("rs114105355")
        slist.remove("rs114905291")
        out["sentinel"] = slist
        return(out)

# get the loci available in eQTLgen
EQTLGEN = get_eqtlgen_hotspots("data/current/eqtl_gen/trans_hotspots.txt")

# get the meQTL loci
MEQTL = glob_wildcards("data/current/sentinels/{sentinel}.dummy")

