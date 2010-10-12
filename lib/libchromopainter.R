chromopainter.read.chunk.table <- function(files) {
    names(files) <- substr(basename(files), 1, 3)
    x <- sapply(files, read.table, row.names=1, header=TRUE, simplify=FALSE)
    z <- sapply(x, colnames)
    stopifnot(z == z[,1])
    lapply(x, as.matrix)
}

chromopainter.read.donor.chunk.counts <- function(dir, chrom)
    chromopainter.read.chunk.table(list.files(path=file.path(dir, "LeaveOneOutNeResults"),
                                              pattern=sprintf("...Chrom%sLeaveOneOutNeEst\\.chunkcounts\\.out", chrom),
                                              full.names=TRUE))

chromopainter.read.recipient.chunk.counts <- function(dir, chrom)
    chromopainter.read.chunk.table(list.files(path=file.path(dir, "PaintingSamples"),
                                              pattern=sprintf("...Chrom%sPaintingNeEstSamples\\.chunkcounts\\.out", chrom),
                                              full.names=TRUE))

chromopainter.read.recipient.chunk.lengths <- function(dir, chrom)
    chromopainter.read.chunk.table(list.files(path=file.path(dir, "PaintingSamples"),
                                              pattern=sprintf("...Chrom%sPaintingNeEstSamples\\.chunklengths\\.out", chrom),
                                              full.names=TRUE))
