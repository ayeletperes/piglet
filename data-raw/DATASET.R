### prepare the datasets for the package.

## germline
HVGERM <- tigger::readIgFasta("data-raw/HVGERM.fasta")

usethis::use_data(HVGERM, overwrite = TRUE)

## functionality
hv_functionality <- data.table::fread("data-raw/ighv_functionality.tsv")

usethis::use_data(hv_functionality, overwrite = TRUE)

