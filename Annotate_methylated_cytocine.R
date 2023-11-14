library(methylKit)
library(genomation)

db.filename <- "./split_cg_context_2/methylDiff_CG.txt.bgz"

bed.filename <- "annotate.bed"

qvalue <- 1e-5


methDiff <- readMethylDB(db.filename)

gene.obj <- readTranscriptFeatures(bed.filename,
                                   remove.unusual=FALSE)

# Hyper

Diff25p.hyper <- getMethylDiff(methDiff,
                               difference=25,
                               qvalue=qvalue,
                               type="hyper")

diffAnn <- annotateWithGeneParts(as(Diff25p.hyper,
                                    "GRanges"),
                                 gene.obj)

print(diffAnn@precedence)
print(diffAnn@num.precedence)

promoter.idx <- which((diffAnn@members[, 1]) > 0)
exon.idx <- which((diffAnn@members[, 2] - diffAnn@members[, 1]) > 0)
intron.idx <- which((diffAnn@members[, 3] - diffAnn@members[, 2] - diffAnn@members[, 1]) > 0)

tss <- getAssociationWithTSS(diffAnn)

write.table(tss[(tss[, 1] %in% promoter.idx), ],
            file="dmr_CG_hyper_promoter.csv",
            row.names=F,
            col.names=T,
            sep=",",
            quote=F)

write.table(tss[(tss[, 1] %in% exon.idx), ],
            file="dmr_CG_hyper_exon.csv",
            row.names=F,
            col.names=T,
            sep=",",
            quote=F)

write.table(tss[(tss[, 1] %in% intron.idx), ],
            file="dmr_CG_hyper_intron.csv",
            row.names=F,
            col.names=T,
            sep=",",
            quote=F)

# Hypo

Diff25p.hypo <- getMethylDiff(methDiff,
                              difference=25,
                              qvalue=qvalue,
                              type="hypo")

diffAnn <- annotateWithGeneParts(as(Diff25p.hypo,
                                    "GRanges"),
                                 gene.obj)

print(diffAnn@precedence)
print(diffAnn@num.precedence)

promoter.idx <- which((diffAnn@members[, 1]) > 0)
exon.idx <- which((diffAnn@members[, 2] - diffAnn@members[, 1]) > 0)
intron.idx <- which((diffAnn@members[, 3] - diffAnn@members[, 2] - diffAnn@members[, 1]) > 0)

tss <- getAssociationWithTSS(diffAnn)

write.table(tss[(tss[, 1] %in% promoter.idx), ],
            file="dmr_CG_hypo_promoter.csv",
            row.names=F,
            col.names=T,
            sep=",",
            quote=F)

write.table(tss[(tss[, 1] %in% exon.idx), ],
            file="dmr_CG_hypo_exon.csv",
            row.names=F,
            col.names=T,
            sep=",",
            quote=F)

write.table(tss[(tss[, 1] %in% intron.idx), ],
            file="dmr_CG_hypo_intron.csv",
            row.names=F,
            col.names=T,
            sep=",",
            quote=F)