
library(graphframes)
library(sparklyr)
library(dplyr)

# Set environ vars
Sys.setenv(SPARK_HOME = '/home/chw925/spark/spark-2.1.0-bin-hadoop2.7/')

# Configure cluster
conf <- spark_config()
conf$spark.executor.cores <- 4
conf$spark.executor.memory <- "20G"
conf$spark.yarn.am.cores  <- 2
conf$spark.yarn.am.memory <- "20G"

# Start spark session
sc <- spark_connect(master = "yarn-client", version = "2.1.0", config = conf)

infile_293t <- "/n/shared_db/ccb/bioplex/data/protein_interactions/BioPlex_293T_Network_10K_Dec_2019.tsv"
df_293t <- spark_read_csv(sc, infile_293t, delimiter="\t")

# "from" vertices
from_tbl <- df_293t %>%
    distinct(SymbolA) %>%
    transmute(id = SymbolA)

# "to" vertices
to_tbl <- df_293t %>%
    distinct(SymbolB) %>%
    transmute(id = SymbolB)

# Union of to and from veritices
vertices_tbl <- from_tbl %>%
    sdf_bind_rows(to_tbl)

# Edges
edges_tbl <- df_293t %>%
    transmute(src = SymbolA, dst = SymbolB)
#sdf_num_partitions(edges_tbl)

# Create graphframe
gf <- gf_graphframe(vertices_tbl, edges_tbl)

# Calculate page rank
pgrank <- gf %>%
    gf_pagerank(reset_prob = 0.15, max_iter = 10L, source_id = "ACADVL")
head(pgrank)

# Calculate vertex degrees
#degrees <- gf %>%
#    gf_degrees()
#head(degrees)

