############################
####       MBI 17Z      ####
####  Artur M. Brodzki  ####
####  Adam Małkowski    ####
############################

# Download gdsfmt and SNPRelate packages
source("http://bioconductor.org/biocLite.R")
biocLite("gdsfmt")
biocLite("SNPRelate")

library(gdsfmt)
library(SNPRelate)
library(shiny)

# load VCF file

vcf = "~/Projekty/MBI/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf"

vcf.fn <- "~/Projekty/MBI/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf"
snpgdsVCF2GDS(vcf.fn, "chr22.gds", method="biallelic.only")

# load GDFS file

chr22 <- snpgdsOpen("Projekty/MBI/chr22.gds")

# load pruns chr22

snpsetcsv <- read.csv("~/Projekty/MBI/chr22-pruned.csv")
snpsetcsv.id = unlist(snpsetcsv$chr22)

pca <- snpgdsPCA(chr22, snp.id = snpsetcsv.id, num.thread=4)

write.csv(pca, "~/Projekty/MBI/pca.csv")

pop_code <- read.table("~/Projekty/MBI/integrated_call_samples_v3.20130502.ALL.panel", sep = "", header = T, nrows = 2504)

# AFR = 1
# AMR = 2
# EAS = 3
# EUR = 4
# SAS = 5
tab <- data.frame(sample.id = pca$sample.id,
                  pop = pop_code$super_pop,
                  EV1 = pca$eigenvect[,3],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)

learn_pop = tab[1:1504, ]
test_pop = tab[1505:2504, ]

test_list = as.vector(test_pop)

geneticDistance <- function (p1, p2) {
  (p1$EV1 - p2$EV1)^2 + (p1$EV2 - p2$EV2)^2
}

print(learn_pop[1, ])

closestRelative <- function(pcatab, x) {
  closeRelative = NULL
  closeDistance = 100000
  for(i in 1:nrow(pcatab)) {
    if(geneticDistance(pcatab[i, ], x) < closeDistance) {
      closeDistance = geneticDistance(pcatab[i, ],x)
      closeRelative = pcatab[i, ]
    }
  }
  closeRelative
}

closestRelative1 = function(x) {
  closestRelative(learn_pop, x)$pop
}

test_pop2 = unlist(test_pop)

r = lapply(test_list, closestRelative1)


closestRelative1(test_pop[200, ])
                

pop_code$super_pop_int = as.integer(tab$pop)





View(pop_code)

plot(tab$EV2, tab$EV1, xlab="eigenvector 2", ylab="eigenvector 1", col=as.integer(tab$pop))
?legend
palette()


runApp('mbi-app')



