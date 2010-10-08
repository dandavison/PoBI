chromopainter.read.counts <- function(files) {
    names(files) <- substr(basename(files), 1, 3)
    x <- sapply(files, read.table, row.names=1, header=TRUE, simplify=FALSE)
    lapply(x, as.matrix)
}

chromopainter.read.donor.copy.counts <- function(dir)
    chromopainter.read.counts(list.files(path=file.path(dir, "LeaveOneOutNeResults"),
                                         pattern=sprintf("...Chrom%sLeaveOneOutNeEst.chunkcounts.out", chrom),
                                         full.names=TRUE))

chromopainter.read.recipient.copy.counts <- function(dir)
    chromopainter.read.counts(list.files(path=file.path(dir, "PaintingSamples"),
                                         pattern=sprintf("...Chrom%sPaintingNeEstSamples.chunkcounts.out", chrom),
                                         full.names=TRUE))
