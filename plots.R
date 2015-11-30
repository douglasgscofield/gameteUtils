reload = function(doit=FALSE) if (doit) source("plots.R")
library(naturalsort)

read.data = function(file="Sample_17_all.singleParentTest.txt") {
    dat = read.delim(file, header=FALSE, na.strings=".", stringsAsFactors=FALSE)
    names(dat) = c("chrom", "pos", "totaltest", "totaltest_p", "pooltest", "twopool_p", "onepool_p")
    all.chromosomes = naturalsort(unique(sort(dat$chr)))
    attr(dat, "all.chromosomes") = all.chromosomes
    attr(dat, "chromosomes") = grep("chr", all.chromosomes, value = TRUE)
    ####
    dat
}

#dat = read.data()
#cat("chromosomes =", paste(collapse=" ", attr(dat, "chromosomes")), "\n\n")

do.plot = function(chromosomes = attr(dat, "chromosomes"), do.png=FALSE) {
    for (chr in chromosomes) {
        sdat = subset(dat, chrom == chr)
        xl = range(c(1, sdat$pos))
        if (do.png) png(paste0(chr, ".png"), width=1000, height=150)
        par(mar=c(3, 7, 0.5, 0), tcl=-0.3, mgp=c(1.5, 0.5, 0), las=1)
        cat("generating figure for", chr, "\n")
        opdat = subset(sdat, onepool_p <= 0.0001)
        plot(opdat$pos, rep(1, nrow(opdat)), pch="|", cex=1.2, xlim=xl, ylim=c(1,3), col="#0000FF80", xlab=paste("Chromosome", chr), ylab="", axes=FALSE)
        tpdat = subset(sdat, twopool_p <= 0.001)
        points(tpdat$pos, rep(2, nrow(tpdat)), pch="|", cex=1.2, col="#00b00040")
        tpdat2 = subset(tpdat, twopool_p <= 0.00005)
        points(tpdat2$pos, rep(3, nrow(tpdat2)), pch="|", cex=1.2, col="#FF0000FF")
        axis(1)
        axis(2, at=c(1,2,3), labels=c("one-p 0.0001", "two-p 0.001", "two-p 0.00005"))
        n.opp = nrow(opdat)
        n.tpp = nrow(tpdat)
        n.tpp2 = nrow(tpdat2)
        legend(x=mean(xl), y=2.5, col=c("#FF0000", "#00b000", "#0000FF"), pch=15, pt.cex = 1.5, horiz=TRUE, xjust = 0.5, yjust = 0.5, adj=c(0.05, 0.5),
               legend=c(paste0("two-p 0.00005 N = ", n.tpp2), paste0("two-p 0.001 N = ", n.tpp), paste0("one-p 0.0001 N = ", n.opp)),
               bty="n")
        if (do.png) dev.off()
    }
}
