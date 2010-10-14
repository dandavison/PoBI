cp.read.chunk.table <- function(chrom, dir, subdir, id, what, exclude)
    cp.read.chunk.table(cp.list.files(chrom, file.path(dir, subdir), id=id, what=what), exclude)

cp.read.donor.chunk.counts <- function(chrom, dir, exclude)
    cp.read.chunk.table(chrom, dir, "LeaveOneOutNeResults", "LeaveOneOutNeEst", "chunkcounts.out", exclude)

cp.read.recipient.chunk.counts <- function(chrom, dir, exclude)
    cp.read.chunk.table(chrom, dir, "PaintingNeResults", "LeaveOneOutNeEst", "chunkcounts.out", exclude)

cp.read.recipient.chunk.lengths <- function(chrom, dir, exclude)
    cp.read.chunk.table(chrom, dir, "PaintingNeResults", "LeaveOneOutNeEst", "chunklengths.out", exclude)

cp.read.chunk.table <- function(files, exclude) {
    names(files) <- substr(basename(files), 1, 3)
    if(!missing(exclude)) files <- files[! names(files) %in% exclude]

    x <- sapply(files, read.table, row.names=1, header=TRUE, simplify=FALSE)
    z <- sapply(x, colnames)
    stopifnot(z == z[,1])
    lapply(x, as.matrix)
}

cp.list.files <- function(chrom, path, id, what)
    list.files(path,
               pattern=sprintf("[A-Z][A-Z][A-Z]Chrom%s%s\\.%s", chrom, id, what),
               full.names=TRUE)
