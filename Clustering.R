library(pvclust)
library(parallel)
number_cluster <- 8
cl <- makeCluster(number_cluster)

context <- "CG"

in_dir <- "./split_cg_context/"
prefix <- "methylBase_"
suffix <- ".txt"
out_dir <- "./cor_clus_figs/"

in_filename <- paste(in_dir, prefix, context, "_methMath", suffix, sep="")

meth <- read.table(in_filename,
                   sep="\t",
                   header=T,
                   row.names=1)

clust_dist <- "correlation"
clust_method <- "complete"
png(paste(out_dir, context, "_clustering_", clust_method, "_pval.png",
          sep=""),
    res=300,
    width=2450,
    height=1400)
result <- pvclust(meth,
                  method.dist=clust_dist,
                  method.hclust=clust_method,
                  use.cor="all.obs",
                  parallel=cl,
                  nboot=1000)

plot(result)
dev.off()
print(paste(out_dir, context, "_clustering_", clust_method, "_pval.png", sep=""))