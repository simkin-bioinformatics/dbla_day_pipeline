#---
#title: "tutorialDBLa_combine"
#author: "Kathryn Tiedje"
#date: "July, 7, 2021"
#output: html_document
#---

### Load Libraries
#```{r}
library(data.table)
library(tidyverse)
#```

### Load data
#```{r}
binary_example <-fread(snakemake@input[["binary_file"]], data.table = FALSE)
ups_example <-fread(snakemake@input[["ups_file"]], data.table = FALSE)
#```

### Binary Example: Rename OTU to DBLa type
#```{r}
binary_example<-binary_example%>% rename_at("#OTU ID", ~"DBLa_type")
#```

### Ups Example: Categorize as DBLa domains as upsA, non-upsA, and other
#```{r}
ups_example $read <- gsub(";.*", "", ups_example$read)
ups_example  <- ups_example  %>% rename_at("read", ~"DBLa_type")
ups_example$domain <- as.factor(ups_example$domain)

#Create an ups categorical variable to stratify as upsA and non-upsA.
ups_example $Ups <- ifelse(grepl("^DBLa1", ups_example$domain, ignore.case = T), "upsA",
                                  ifelse(grepl("^DBLa0", ups_example$domain, ignore.case = T), "non-upsA", 
                                         ifelse(grepl("^DBLa2", ups_example$domain, ignore.case = T), "non-upsA", "Other")))
ups_example$Ups <- as.factor(ups_example$Ups)
#```

## FINAL Combined: Combine ups grouping to the DBLa types in the isolate_example. This file used for the combined DBLa type analyses.
#```{r}
isolate_ups_example<-binary_example %>% inner_join(dplyr::select(ups_example, DBLa_type, Ups), by ="DBLa_type")
isolate_ups_example <- isolate_ups_example %>% relocate(Ups, .after = DBLa_type)

write.table(isolate_ups_example, file=snakemake@output[["combined_dbla_file"]], sep=",", row.names=F)
