thai_summary_stats <- data.table::fread('../../data/thailand_all.summary_stats.txt') %>% dplyr::filter(Chrom == 21 & Pos > 43289997-50000 & Pos < 43289997+50000)

query_SNP <- 'rs12151990'

library(LDlinkR)

LD_Mat <- LDmatrix(snps = c(query_SNP,thai_summary_stats$RS), 
         pop = "CHB", 
         r2d = "r2", 
         token = 'b9f8399d9bab',
         genome_build = "grch38"
)
