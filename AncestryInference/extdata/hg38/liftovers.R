library(data.table)
library(rtracklayer)

chain19to38 <- import.chain("~/programs/liftOver/hg19ToHg38.over.chain")
chain38to19 <- import.chain("~/programs/liftOver/hg38ToHg19.over.chain")

# aims
aims19 <- import.bed("../aims.hg19.bed")
aims38 <- liftOver(aims19, chain19to38)
unlifted <- which(sapply(aims38, isEmpty))
aims38 <- unlist(aims38) # GRangesList to GR
aims38_noChr <- aims38
seqlevelsStyle(aims38_noChr) <- "Ensembl"

# Freqs
freq19 <- fread("../freqs.hg19.txt")
freq38 <- freq19[!unlifted] #Two unlifted positions
dt19 <- as.data.table(aims19)
all(dt19$seqnames == freq19$CHR)
all(dt19$end == freq19$POS)
freq38[, POS := end(aims38)]
freq38_noChr <- copy(freq38)
freq38_noChr[, CHR := gsub("chr", "", CHR)]

# Exports
export.bed(aims38, "aims.hg38.bed")
export.bed(aims38_noChr, "aims.hg38.noChr.bed")
fwrite(freq38, "freqs.hg38.txt", sep = " ")
fwrite(freq38_noChr, "freqs.hg38.noChr.txt", sep = " ")

