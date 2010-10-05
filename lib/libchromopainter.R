chromopainter.read.donor.copy.counts <- function(dir) {
    files <- list.files(path=file.path(dir, "LeaveOneOutNeResults"),
                        pattern=sprintf("...Chrom%sLeaveOneOutNeEst.chunkcounts.out", chrom),
                        full.names=TRUE)
    names(files) <- substr(basename(files), 1, 3)
    x <- sapply(files, read.table, row.names=1, header=TRUE, simplify=FALSE)
    lapply(x, as.matrix)
}

chromopainter.read.recipient.copy.counts <- function(dir) {
    files <- list.files(path=file.path(dir, "PaintingSamples"),
                        pattern=sprintf("...Chrom%sPaintingNeEstSamples.chunklengths.out", chrom),
                        full.names=TRUE)
    names(files) <- substr(basename(files), 1, 3)
    x <- sapply(files, read.table, row.names=1, header=TRUE, simplify=FALSE)
    lapply(x, as.matrix)
}
