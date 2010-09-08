wtccc2.read.samplefile <- function(coh, platform) {
    file <- file.path("data/wtccc2", coh, platform, "calls", sprintf("%s_%s.sample", coh, platform))
    d <- read.table(file, skip=2)[, 1:3, drop=FALSE]
    colnames(d) <- c("id1", "id2", "missing")
    d
}

wtccc2.read.chiamo <- function(platform, chrom, cols) {
    p <- pipe(sprintf("zcat data/wtccc2/POBI/%s/calls/POBI_%s_%s.gen.gz | cut -d' ' -f%s",
                      platform, chrom, platform, paste(cols, collapse=",")))
    on.exit(close(p))
    read.table(p, as.is=TRUE)
}
