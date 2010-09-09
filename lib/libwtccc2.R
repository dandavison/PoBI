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

is.same.strand <- function(alleles1, alleles2) {
    ## At this stage the allele frequencies correspond; although the
    ## two data sets may be referring to different strands. Return
    ## logical vector of same strand indicators

    ## alleles1 is 2xL character matrix

    flip <- function(pair) c(A="T", C="G", G="C", T="A")[pair]
    alleles1.f <- apply(alleles1, 2, flip)

    id <- alleles1 == alleles2
    id.f <- alleles1.f == alleles2
    noflip <- colSums(id) == 2
    flip <- colSums(id.f) == 2
    stopifnot(xor(flip, noflip))
    noflip
}
