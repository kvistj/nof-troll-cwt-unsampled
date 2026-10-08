# 00_export_public_data.R
#
# AUTHOR-ONLY. Run once, in a session with access to the PSC CAMP cache
# (Data/camp_tables.rds), to build the public stock lookup in data/derived/.
# Users reproducing the analysis do not need to run this script.
#
# Only the three columns the analysis uses (tag_code, stock, brood_year) are
# exported, and only for tag codes actually recovered in the North of Falcon
# troll recovery file. Nothing else from CAMP is written out.

library(tidyverse)

pad_tag <- function(x) str_pad(as.character(x), width = 6, pad = "0")

camp_tables <- readRDS("Data/camp_tables.rds")   # private cache, not in the repo

recovered_tags <- read.table("data/raw/NOF Troll Recoveries 2005-2025.txt",
                             sep = ",", header = TRUE,
                             colClasses = c(tag_code = "character")) |>
  mutate(tag_code = pad_tag(tag_code)) |>
  distinct(tag_code)

lookup <- camp_tables$camp_cwt_release |>
  mutate(tag_code = pad_tag(tag_code)) |>
  select(tag_code, stock, brood_year) |>
  semi_join(recovered_tags, by = "tag_code") |>
  distinct() |>
  arrange(tag_code)

# One row per tag code, or the join in the analysis will duplicate recoveries
stopifnot(!anyDuplicated(lookup$tag_code))

dir.create("data/derived", showWarnings = FALSE, recursive = TRUE)
write.csv(lookup, "data/derived/indicator_stock_lookup.csv", row.names = FALSE)

message("Wrote ", nrow(lookup), " tag codes to data/derived/indicator_stock_lookup.csv")
