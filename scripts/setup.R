library(duckdb)
library(DBI)
library(here)

db_path <- here("IC50_expr_cnv.duckdb")

# Delete if exists
if (file.exists(db_path)) {
  unlink(db_path)
}

# Load IC50_expr_cnv.csv into a table named `IC50_expr_cnv`
conn <- dbConnect(duckdb(), dbdir = db_path)
duckdb_read_csv(conn, "IC50_expr_cnv", here("IC50_expr_cnv.csv"))
dbDisconnect(conn)
