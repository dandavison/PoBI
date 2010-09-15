pobi.read.samplefile <- function(platform="illumina") {
    d <- wtccc2.read.samplefile("POBI", platform)
    manifest <- pobi.read.manifest()
    reg <- manifest[d$id2,"GEOGRAPHICAL.REGION"]
    reg[is.na(reg)] <- "NA"
    d$reg <- factor(reg)
    d
}

pobi.read.manifest <- function() {
    d <- read.delim("data/POBI/sanger-sample-manifest.tsv", na.strings="", as.is=TRUE)
    d <- d[c("SANGER.SAMPLE.ID", "GEOGRAPHICAL.REGION")]
    rownames(d) <- d[,1]

    ## Edit region names
    z <- d[,"GEOGRAPHICAL.REGION"]
    z[z == "Nottinghampshire"] <- "Nottinghamshire"
    z[z == "Buckinghamshire"] <- "Oxfordshire"
    d[,"GEOGRAPHICAL.REGION"] <- z
    
    d
}

ms.read.samplefile <- function(platform="illumina") {
    d <- wtccc2.read.samplefile("MS", platform)
    MS.info.file <- "data/MS/MS_illumina.sample.geoinfoIII"
    MS.info <- read.table(MS.info.file, header=TRUE, row.names=2)[-1]
    map <- match(rownames(d), rownames(MS.info))
    if(any(noinfo <- is.na(map)))
        warning(sum(noinfo), " individuals are absent from ", MS.info.file)
    d$reg <- factor(MS.info$geo.info[map])
    d
}

pobi.ms.read.samplefile <- function()
    rbind(pobi.read.samplefile(),
          ms.read.samplefile())

extract.haplotypes <- function(pop, haps, indivs.d) {
    ## Extract all haplotypes from region POP.
    ids <- colnames(haps)
    stopifnot(ids %in% rownames(indivs.d))
    select <- indivs.d[ids, "reg"] == pop
    haps[, select, drop=FALSE]
}

read.haplotypes <- function(hapfile, ids) {
    ## Return L x 2n binary matrix
    p <- pipe(paste("cut -d ' ' -f 6- <", hapfile))
    matrix(scan(p, what=integer()), ncol=2*length(ids),
           dimnames=list(NULL, rep(ids, each=2)))
}

read.haplotypes.legend <- function(hapfile, ids) {
    ## Return L x 2n binary matrix
    p <- pipe(paste("cut -d ' ' -f 6- <", hapfile))
    matrix(scan(p, what=integer()), ncol=2*length(ids),
           dimnames=list(NULL, rep(ids, each=2)))
}

read.chiamo.legend <- function(chiamofile)
    read.table(pipe(paste("cut -d' ' -f1-5 <", chiamofile)),
               col.names=c("id1","id2","pos", "allele1", "allelel2"),
               as.is=TRUE)

make.intervals <- function(points, width, overlap) {
    ## Return intervals of width WIDTH overlapping by OVERLAP
    i <- 1
    starts <- min(points)
    ends <- starts[i] + width
    while(ends[i] < max(points)) {
        i <- i+1
        starts[i] <- ends[i-1] - overlap
        ends[i] <- starts[i] + width
    }
    ends[i] <- max(points)
    cbind(start=starts, end=ends)
}

qsub.script <- function(cmd, name, outfile, errfile)
    paste("#$ -N ", name, "\n",
          "#$ -o ", outfile, "\n",
          "#$ -e ", errfile, "\n",
          "#$ -cwd", "\n",
          "#$ -V", "\n",
          "#$ -pe level2.pe 1", "\n",
          "#$ -S /bin/bash", "\n",
          cmd, "\n",
          sep="")

stitch.haplotypes <- function(files, ids, outfile, thresh=.9) {
    write.haplotypes <- function(leg, h, file, append)
        write.table(cbind(leg, h), file=file, append=append,
                    quote=FALSE, row.names=FALSE, col.names=FALSE)

    match.haplotypes <- function(haps1, leg1, haps2, leg2) {
        ## Return haps2 haplotypes matched to haps1 haplotypes
        ## Haplotypes are in columns

        n1 <- nrow(leg1)

        ## Identify overlapping SNPs
        stopifnot(leg2$id2[1] %in% leg1$id2)
        olap <- which(leg2$id2[1] == leg1$id2)
        o1.idx <- olap:n1
        o2.idx <- 1:(n1-olap+1)
        stopifnot(leg1$id2[o1.idx] == leg2$id2[o2.idx])
        
        ## Form column index vector matching haplotypes
        o1 <- haps1[o1.idx,,drop=FALSE]
        o2 <- haps2[o2.idx,,drop=FALSE]

        stopifnot(dim(o1) == dim(o2), ncol(o1) %% 2 == 0)
        n <- ncol(o1)/2
        idx <- rep(NA, 2*n)

        for(i in 1:n) {
            ii <- 2*(i-1)
            p11 <- mean(o1[,ii + 1] == o2[,ii + 1])
            p12 <- mean(o1[,ii + 1] == o2[,ii + 2])
            p21 <- mean(o1[,ii + 2] == o2[,ii + 1])
            p22 <- mean(o1[,ii + 2] == o2[,ii + 2])

            if(p11 > thresh) {
                stopifnot(p22 > thresh)
                idx[ii + (1:2)] <- ii + (1:2)
            }
            else if(p12 > thresh) {
                stopifnot(p21 > thresh)
                idx[ii + (1:2)] <- ii + (2:1)
            }
            else {
                msg <- paste("Failed to find match for individual", i)
                cat(msg, "\n")
                warning(msg)
            }
        }
        haps2[-o2.idx,idx,drop=FALSE]
    }

    
    hprev <- read.haplotypes(files[1], ids)
    lprev <- read.chiamo.legend(files[1])
    write.haplotypes(lprev, hprev, outfile, append=FALSE)
    for(f in files[-1]) {
        h <- read.haplotypes(f, ids)
        l <- read.chiamo.legend(f)
        h <- match.haplotypes(hprev, lprev, h, l)
        write.haplotypes(l, h, outfile, append=TRUE)
        hprev <- h
        lprev <- l
    }
}
