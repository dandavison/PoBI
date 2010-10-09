pobi.read.samplefile <- function(platform="illumina") {
    d <- wtccc2.read.samplefile("POBI", platform)
    manifest <- pobi.read.manifest()
    reg <- manifest[d$ID_2,"GEOGRAPHICAL.REGION"]
    reg[is.na(reg)] <- "NA"
    d$reg <- factor(reg)
    d
}

pobi.read.manifest <- function(file="data/POBI/sanger-sample-manifest-gender-updated.csv", sep=",") {
    d <- read.table(file, header=TRUE, sep=sep, na.strings="", as.is=TRUE)
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
    rn <- scan(pipe(paste("cut -d ' ' -f 2 <", hapfile)), what="")
    p <- pipe(paste("cut -d ' ' -f 6- <", hapfile))
    matrix(scan(p, what=integer()), ncol=2*length(ids), byrow=TRUE,
           dimnames=list(rn, rep(ids, each=2)))
}

read.chiamo.legend <- function(chiamofile)
    read.table(pipe(paste("cut -d' ' -f1-5 <", chiamofile)),
               col.names=c("ID_1","ID_2","pos", "allele1", "allele2"),
               as.is=TRUE)

read.chiamo <- function(genfile, ids=1:n, gz=FALSE, thresh=.9) {
    leg <- read.chiamo.legend(genfile)
    L <- nrow(leg)
    p <- scan(pipe(paste("cut -d' ' -f 6- <", genfile)), what=double())
    stopifnot(length(p) %% (3*L) == 0)
    n <- length(p) / (3*L)
    stopifnot(length(ids) == n)
    dim(p) <- c(3,n,L)
    w <- p > thresh
    ## g <- apply(w, c(2,3), function(s) ifelse(any(s), which(s), NA)) - 1
    g <- array(as.integer(NA), dim=c(n, L))
    g[w[1,,]] <- 0
    g[w[2,,]] <- 1
    g[w[3,,]] <- 2
    dimnames(g) <- list(ids, leg$ID_2)
    t(g)
}

compare.haplotypes <- function(h1, h2, imgfile, twon=ncol(h1)) {

    sum.adjacent.columns <- function(x, colnames=NULL) {
        stopifnot((twon <- ncol(x)) %% 2 == 0)
        n <- twon / 2
        xo <- x[,seq.int(1, twon-1, 2)]
        xe <- x[,seq.int(2, twon,   2)]
        matrix(xo + xe, nrow(x), n, dimnames=list(rownames(x), colnames))
    }
    
    ## L x 2n
    stopifnot(dim(h1) == dim(h1), twon %% 2 == 0)
    stopifnot(identical(dimnames(h1), dimnames(h2)))
    h1 <- h1[,1:twon] ; h2 <- h2[,1:twon]
    n <- twon / 2
    
    ## Check genotypes are the same
    odds <- seq(1, twon-1, by=2)
    ids <- colnames(h1)[odds]
    g1 <- sum.adjacent.columns(h1, colnames=ids)
    g2 <- sum.adjacent.columns(h2, colnames=ids)
    if(any(bad <- g1 != g2)) {
        bad.indiv1 <- which(colSums(bad) > 0)[1]
        bad.snp1 <- which(bad[,bad.indiv1])[1]
        cat(sum(bad), "genotypes disagree out of", length(bad), "\n")
        cat("These are in", sum(colSums(bad) > 0), "out of", ncol(bad), "individuals\n")
        cat("and involve", sum(rowSums(bad) > 0), "out of", nrow(bad), "SNPs\n")
        cat("The first bad genotype is individual", bad.indiv1, "SNP", bad.snp1, ":\n")
        print(cbind(h1[bad.snp1,2*(bad.indiv1-1)+(1:2),drop=FALSE],
                    h2[bad.snp1,2*(bad.indiv1-1)+(1:2),drop=FALSE]))
        pdf(file="genotype-agreements.pdf")
        image(z=(g1 == g2), x=c(0,seq(nrow(g1))), y=c(0,1:n))
        dev.off()
    }

    ## Compare haplotypes
    cat("Comparing haplotypes\n")
    png(file=imgfile, width=1000, height=500)
    homo <- g1 %in% c(0,2) | g2 %in% c(0,2)
    dim(homo) <- dim(g1)
    hap.homo <- array(dim=c(nrow(homo),2,ncol(homo)))
    hap.homo[,1,] <- hap.homo[,2,] <- homo
    dim(hap.homo) <- c(nrow(homo), 2*ncol(homo))
    h1[hap.homo] <- h2[hap.homo] <- NA
    
    image(z=(h1[,odds] == h2[,odds]), x=c(0,seq(nrow(h1))), y=c(0,1:n), col=c("red","blue"))
    dev.off()
}

make.intervals <- function(n, width, overlap) {
    ## There are n points on a line. Return start and end indices of
    ## intervals of width WIDTH, overlapping by OVERLAP points.
    i <- 1
    starts <- 1
    ends <- integer()
    while(TRUE) {
        ends[i] <- starts[i] + width - 1
        if(ends[i] > n) { ends[i] <- n ; break }
        i <- i+1
        starts[i] <- ends[i-1] - overlap + 1
    }
   cbind(start=starts, end=ends)
}

qsub.script <- function(cmd, name, outfile, errfile, level=2)
    paste("#$ -N ", name, "\n",
          "#$ -o ", outfile, "\n",
          "#$ -e ", errfile, "\n",
          "#$ -cwd", "\n",
          "#$ -V", "\n",
          "#$ -pe level", level, ".pe 1", "\n",
          "#$ -S /bin/bash", "\n",
          cmd, "\n",
          sep="")

stitch.haplotypes <- function(files, ids, outfile, nolap, thresh=.9, olap.nuse=10) {
    write.haplotypes <- function(leg, h, file, append)
        write.table(cbind(leg, h), file=file, append=append,
                    quote=FALSE, row.names=FALSE, col.names=FALSE)

    match.haplotypes <- function(haps1, leg1, haps2, leg2) {
        ## Return haps2 haplotypes matched to haps1 haplotypes
        ## Haplotypes are in columns

        n1 <- nrow(leg1)

        ## Identify overlapping SNPs
        stopifnot(leg2$ID_2[1] %in% leg1$ID_2)
        ostart1 <- which(leg2$ID_2[1] == leg1$ID_2)
        o1.idx <- ostart1:n1
        print(length(o1.idx))

        ##stopifnot(length(o1.idx) == nolap)
        if(length(o1.idx) != nolap) {
            msg <- paste(length(o1.idx), "in overlapping region yet nolap=", nolap, "\n")
            cat(msg)
            warning(msg)
        }

        nolap <- length(o1.idx)
        o2.idx <- 1:nolap
        stopifnot(leg1$ID_2[n1] == leg2$ID_2[nolap])
        stopifnot(leg1$ID_2[o1.idx] == leg2$ID_2[o2.idx])
        stopifnot(olap.nuse <= nolap)
        
        ## Form column index vector matching haplotypes
        o1 <- haps1[o1.idx,,drop=FALSE]
        o2 <- haps2[o2.idx,,drop=FALSE]

        stopifnot(dim(o1) == dim(o2), ncol(o1) %% 2 == 0)
        n <- ncol(o1)/2
        idx <- rep(NA, 2*n)

        for(i in 1:n) {
            ii <- 2*(i-1)

            ## Check that genotypes are identical
            if(any(bad <- rowSums(o1[,ii+(1:2)]) != rowSums(o2[,ii+(1:2)]))) {
                msg <- paste(sum(bad), "genotypes disagree")
                cat(msg, "\n")
                warning(msg)
            }
            
            p <- array(dim=c(2,2))
            w <- olap.nuse
            while(TRUE) {
                p[1,1] <- mean(o1[1:w,ii + 1] == o2[1:w,ii + 1])
                p[1,2] <- mean(o1[1:w,ii + 1] == o2[1:w,ii + 2])
                p[2,1] <- mean(o1[1:w,ii + 2] == o2[1:w,ii + 1])
                p[2,2] <- mean(o1[1:w,ii + 2] == o2[1:w,ii + 2])
                if(any(p < 1) || w == nolap) break
                w <- min(w + olap.nuse, nolap)
                cat("Extending overlap region to", w, "\n")
            }

            ## print(round(p, 2))
            p11 <- p[1,1] + p[2,2]
            p12 <- p[1,2] + p[2,1]
            if(p11 == p12)
                pair <- sample(1:2)
            else if(p11 > p12)
                pair <- 1:2
            else
                pair <- 2:1

            idx[ii + (1:2)] <- ii + pair
            
            if(!(p11 >= 2*thresh || p12 >= 2*thresh)) {
                msg <- paste("No clear resolution for individual", i)
                cat(msg, "\n")
                warning(msg)
            }
        }
        list(h=haps2[-o2.idx,idx,drop=FALSE], l=leg2[-o2.idx,,drop=FALSE])
    }

    
    hprev <- read.haplotypes(files[1], ids)
    lprev <- read.chiamo.legend(files[1])
    write.haplotypes(lprev, hprev, outfile, append=FALSE)
    for(f in files[-1]) {
        h <- read.haplotypes(f, ids)
        l <- read.chiamo.legend(f)
        hl <- match.haplotypes(hprev, lprev, h, l)
##        stopifnot(nrow(hl$h)+nolap == nrow(h), nrow(hl$l)+nolap == nrow(l))
        if(nrow(hl$h)+nolap == nrow(h) || nrow(hl$l)+nolap == nrow(l)) {
            msg1 <- paste(f, "\n",nrow(hl$h), "+", nolap, "vs.", nrow(h), "\n")
            msg2 <- paste(nrow(hl$l),"+",nolap,"vs.", nrow(l))
            msg <- paste(msg1, msg2)
            cat(msg)
            warning(msg)
        }
        write.haplotypes(hl$l, hl$h, outfile, append=TRUE)
        hprev <- hl$h
        lprev <- hl$l
    }
}
