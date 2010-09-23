wtccc2.read.samplefile <- function(coh, platform) {
    file <- file.path("data/wtccc2", coh, platform, "calls", sprintf("%s_%s.sample", coh, platform))
    d <- read.table(file, skip=2)[, 1:3, drop=FALSE]
    colnames(d) <- c("ID_1", "ID_2", "missing")
    rownames(d) <- d$ID_1
    d
}

wtccc2.write.samplefile <- function(d, file, datatypes) {
    stopifnot(ncol(d) >= 3, colnames(d)[1:3] == c("ID_1","ID_2","missing"))
    if(missing(datatypes))
        datatypes <- c(rep("0",3), rep("C", ncol(d)-3))
    else
        stopifnot(datatypes[1:3] == "0")
    write.table(d[FALSE,], file, col.names=TRUE, row.names=FALSE, quote=FALSE)
    cat(datatypes, "\n", file=file, append=TRUE)
    write.table(d, file, col.names=FALSE, row.names=FALSE, quote=FALSE, append=TRUE)
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
