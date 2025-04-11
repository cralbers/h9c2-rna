library(gprofiler2)
library(dplyr)
library(plotly)
library(htmlwidgets)
library(docstring)

setwd("/Users/coraalbers/Documents/smith_lab/h9c2_RNA_seq/RNAlysis/cab_only_cluster")



load("gprofiler_all.RData")

tbl <- as.data.frame(gprofiler_results_oe$result)


# tbl <- apply(tbl,2,as.character)
# write.table(tbl, file = "gprofiler_all_results.tsv", sep = "\t", quote = FALSE)


find_term <- function(term, df) {
  
  if (term %in% df$term_name) {
    
    term_df <- df[df$term_name == term,]
    print(term_df)
  }
}

find_term("muscle tissue development", tbl)

find_term('cytoplasm', tbl)

find_term("actin cytoskeleton organization", tbl)
