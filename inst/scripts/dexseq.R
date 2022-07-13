############################################################
# 
# author: Ludwig Geistlinger
# date: 2022-07-13 15:24:52
# 
# descr: construction of HEK293T_HCT116_exon_counts 
#        DEXSeqDataSet from GEO 
# 
############################################################

# RNA-seq data for HEK293T cells was obtained from GEO accession GSE122633
# (runs labeled UnTfx1, UnTfx2 and UnTfx3 were used as those assayed untransformed
# wild type HEK293T cells)

# RNA-seq data for HCT116 cells was obtained from GEO accession GSE52429 
# (runs labeled repPF and repPJ were used as those assayed untransformed wild type
# HCT116 cells)

# The reads were aligned to the GRCh38.104 genome by STAR 2.7.9a 
# and counted by DEXseq 1.42.0.

# The resulting DEXSeqDataSet was filtered to exclude genes for which all exons
# have zero read counts using the following code:
# ("dxd" is the name of the DEXSeqDataSet)

library(DEXSeq)
rs <- rowSums(assay(dxd)[,1:5])
all.zero <- unname(rs) == 0
spl <- split(all.zero, rowData(dxd)$groupID)
all.zero <- vapply(spl, all, logical(1))
all.zero <- names(spl)[all.zero]

ind <- rowData(dxd)$groupID %in% all.zero
dxd <- dxd[!ind,]
