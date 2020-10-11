library(testthat)
pbo <- pbapply::pboptions(type = "none") #turn progress bars off
test_check("linking-data")
pbapply::pboptions(pbo) #reset progress bar options
