library(COVID19)

for (country in c("AUS", "BRA", "CHL", "COL", "CZE" ,"DEU", "JPN", "LTU",  "ZAF") ) {
    x <- covid19(country, level = 1, verbose = FALSE, start = "2021-03-31", end = "2021-11-22")
    write.csv(x, file=paste0(country, ".csv"), row.names=FALSE)
}