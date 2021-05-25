# system("wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel")
populations <- read.table("integrated_call_samples_v3.20130502.ALL.panel", header = T)
EUR <- as.character(populations[populations$super_pop == "EUR","sample"])
write.table(cbind(EUR, EUR), "EUR.txt", col.names = F, row.names = F, quote = F)
