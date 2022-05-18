library(dplyr)
library(tidyverse)
library

#Access numbers were collected using pysradb containing the following filters: query="epilepsy", 
#platform	= "Illumina", organism = "Homo sapiens",layout = "PAIRED"
ena.search <- read.csv("~/pysra_search/ena.csv")
geo.search <- read.csv("~/pysra_search/geo.csv")
sra.search <- read.csv("~/pysra_search/sra.csv")

#Filter results of ENA search to only keep transcriptomic experiments
ena.filtered <- filter(ena.search, library_source == "TRANSCRIPTOMIC", read_count > 30000000)
#Filter results of SRA search to only keep transcriptomic experiments
#at least 30-40M reads
sra.filtered <- filter(sra.search, experiment_library_source == "TRANSCRIPTOMIC", run_1_size > 30000000)
geo.filtered <- filter(geo.search, run_1_size > 30000000)

#Curate data by unique study accession, if they attend to the further requirements and if the data is publicly available

unique(ena.filtered$study_accession)
#"PRJNA563467" "PRJNA589589" "PRJNA280563" "PRJNA389517"

unique(geo.filtered$study_accession) 
#All GEO studies are contained in the SRA database, no need to curate further.
unique(sra.filtered$study_accession)
#"SRP360371" "SRP329510" "SRP318206" "SRP288732" "SRP273213" "SRP269629" "SRP246291" "SRP229967" "SRP150178" "ERP117136"
#"SRP220383" "SRP216012" "SRP132816" "SRP187565" "SRP108851" "SRP062617" "SRP061286" "SRP056956"

#After curation, most suitable datasets

PRJNA563467 - SRP220383
PRJNA280563 - SRP056956
SRP150178 - PRJNA369732
ERP117136 - PRJEB34260
SRP216012
SRP187565 - only MTLE samples
SRP108851 - only MTLE samples
SRP061286 - MTLE with and without HS