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
    ids <- row.names(haps)
    stopifnot(ids %in% rownames(indivs.d))
    haps[indivs.d[ids, pop] == pop,,drop=FALSE]
}

read.haplotypes <- function(hapfile, ids) {
    p <- pipe(paste("cut -d ' ' -f 6- <", hapfile))
    matrix(scan(p, what=integer()), ncol=2*length(ids),
           dimnames=list(NULL, rep(ids, each=2)))
}
