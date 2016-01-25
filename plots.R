options(prompt='gamete> ', scipen=15, stringsAsFactors=FALSE)
reload = function(doit=FALSE) if (doit) source("plots.R")
library(naturalsort)
library(stringr)

site1.col = "violetred"
site2.col = "gray10"

plot.hist = function(x, bks=seq(-9, -2, length.out=200)) {
    pdf("hist.pdf", width=2.5, height=1)
    opa = par(mar=c(1.5, 1.2, 0.0, 0.0), mgp=c(0.65, 0.2, 0), tcl=-0.1, las=1, ps=10, cex.axis=0.6, cex.lab=1.2)
    h = hist(x, breaks=bks, plot=FALSE)
    cuts = as.character(cut(h$breaks, c(-Inf, log10(0.000001), log10(0.00005), Inf), labels=c("violetred","gray10","gray")))
    is.na(cuts)[h$counts == 0] = TRUE
    plot(h, col=cuts, border=cuts, main="", xlab="", ylab="", axes=FALSE)
    #plot.new()
    #plot.window(xlim=lim, ylim=lim, asp=1, bty="o")
    #qq = qqplot(x, y, plot.it = FALSE)
    #qq = as.list(as.data.frame(qq)[length(qq$x):1, ]) # reverse x and y so lower values overplot higher values
    ##segments(-6.5, -6, -5.5, -6, col="darkblue", lwd=1)
    #segments(-10, -10, -2, -2, col="grey", lwd=0.5)
    #point.col = ifelse(qq$y <= log10(0.000001), site1.col, ifelse(qq$y <= log10(0.00005), site2.col, "grey60"))
    #print(table(point.col))
    #points(qq$x, qq$y, pch=20, cex=0.8, col=point.col)
    #box(bty="o", lwd=0.5)
    #par(mgp=c(1.2, 0, 0))
    axis(1, lwd=0.5, mgp=c(0, -0.3, 0))
    title(xlab=expression("Observed two-pool log-"*italic(L)))
    #par(mgp=c(1.0, 0.3, 0))
    l = a = seq(0, 1500, 500)
    is.na(l[length(l)]) = TRUE
    axis(2, lwd=0.5, at=a, labels=l)
    #title(ylab="Counts")
    #text(-10.3, -2, pos=4, substitute(italic(N)~obs==NN, list(NN=length(y))))
    #text(-10.3, -3, pos=4, substitute(italic(N)~obs~log[10]~italic(P)*"< -6" == NN, list(NN=sum(y < (-6)))))
    text(-9.9, 1600, "A", xpd=NA, font=2)
    par(opa)
    dev.off()
}

plot.qq = function(x, y, lim=c(-9.2, -2)) {
    pdf("qq.pdf", width=2.5, height=2.5)
    opa = par(mar=c(2.2, 2.2, 0.1, 0.1), mgp=c(1.0, 0.3, 0), tcl=-0.2, las=1, ps=10, cex.axis=0.7, cex.lab=1.3)
    plot.new()
    plot.window(xlim=lim, ylim=lim, asp=1, bty="o")
    qq = qqplot(x, y, plot.it = FALSE)
    qq = as.list(as.data.frame(qq)[length(qq$x):1, ]) # reverse x and y so lower values overplot higher values
    #segments(-6.5, -6, -5.5, -6, col="darkblue", lwd=1)
    segments(-10, -10, -2, -2, col="grey", lwd=0.5)
    point.col = ifelse(qq$y <= log10(0.000001), site1.col, ifelse(qq$y <= log10(0.00005), site2.col, "grey60"))
    print(table(point.col))
    points(qq$x, qq$y, pch=20, cex=0.8, col=point.col)
    box(bty="o", lwd=0.5)
    par(mgp=c(1.2, 0, 0))
    axis(1, lwd=0.5)
    title(xlab=expression("Null two-pool log-"*italic(L)))
    par(mgp=c(1.0, 0.3, 0))
    axis(2, lwd=0.5)
    title(ylab=expression("Observed two-pool log-"*italic(L)))
    #text(-10.3, -2, pos=4, substitute(italic(N)~obs==NN, list(NN=length(y))))
    #text(-10.3, -3, pos=4, substitute(italic(N)~obs~log[10]~italic(P)*"< -6" == NN, list(NN=sum(y < (-6)))))
    par(opa)
    dev.off()
}

f = function(cov=50, n=1000000, p.thresh=0.1) {
    rbinom.1 = function(cov, n, p=0.5) { 
        ncov = rbinom(n, cov, p=0.98)
        x = rbinom(n, ncov, p)
        #ifelse(x > cov / 2.0, cov - x, x)
    }
    h = rbinom.1(cov, 50*n, 0.5) 
    p = pbinom(h, cov, 0.5)
    lp = length(p)
    p = p[p <= p.thresh]
    cat("p started with",lp,"elements and now has",length(p),"elements, using",2*n,"elements\n")
    stopifnot(length(p) >= 2*n)
    log10(p[1:n] * p[(n+1):(2*n)])
}

read.data = function(file="Sample_17_all.singleParentTest.binom.new.txt") {
    dat = read.delim(file, header=FALSE, na.strings=".", stringsAsFactors=FALSE, nrows = 2000000)
    names(dat) = c("chrom", "pos", "totaltest", "totaltest_p", "pooltest", "twopool_p", "onepool_p")
    if (any(dat$twopool_p < 0)) {
        cat("read.data: twopool_p appears to be logged, changing name to twopool_log10_p and creating twopool_p ...\n")
        names(dat)[6] = "twopool_log10_p"
        dat$twopool_p = 10 ** dat$twopool_log10_p
    }
    all.chromosomes = naturalsort(unique(sort(dat$chrom)))
    attr(dat, "all.chromosomes") = all.chromosomes
    attr(dat, "chromosomes") = grep("chr", all.chromosomes, value = TRUE)
    ####
    dat
}

get.depths = function(x) {
    d = str_split_fixed(x$totaltest, fixed(':'), 2)[, 2]
    d = str_split_fixed(d, fixed('/'), 2)
    storage.mode(d) = 'integer'
    d
}

load.null = function(file = 'null-uniform-100000.singleParentNull.txt') {
    dat.null = read.data(file)
    dat.null <<- dat.null
}

generate.null = function(d, n.iter = 3, method = c('binom', 'uniform'),
                         file = paste0('null-', match.arg(method), '-', n.iter, '.txt')) {
    method = match.arg(method)
    if (is.data.frame(d))
        d = get.depths(d)
    stopifnot(is.matrix(d))
    ans = matrix(0, n.iter, 5)
    randomrows = sample(nrow(d), n.iter, replace = TRUE)
    #randomrows = 1:n.iter
    redo = 0
    for (i in 1:n.iter) {
        if (i %% 50000 == 0) cat("iteration",i,"...\n")
        row = randomrows[i]
        tot = sum(d[row, ])  # total number of alleles
        fr = d[row, 1] / tot # frequency of allele 1
        ok = 0
        while (! ok) {
            # The number of alleles in each pool is always binomial
            tot.1 = rbinom(1, tot, 0.5) # number of alleles in pool 1
            if (method == 'binom') {
                # binomial allele counts
                p1.1 = rbinom(1, tot.1, fr) # number of allele 1 in pool 1
                p1.2 = tot.1 - p1.1  # number of allele 2 in pool 1
                p1 = c(p1.1, p1.2)
                p2 = d[row, ] - p1
            } else if (method == 'uniform') {
                # uniform allele counts, this is likely wrong
                p1 = c(sample.int(d[row, 1], 1), sample.int(d[row, 2], 1))
                p2 = d[row, ] - p1
            } else {
                stop('unknown method:', method)
            }
            if (all(c(p1, p2) >= 0) && any(p1 != 0) && any(p2 != 0))
                ok = 1
            else redo = redo + 1
        }
        ans[i, ] = c(i, p1, p2)
        #cat(paste(collapse=":", c(ans[i, ], p1 + p2)), "\n")
    }
    ans = cbind("chrnull", as.data.frame(ans))
    write.table(ans, file = file, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
    cat("wrote", n.iter, "rows to file", file, ", method =", method, ", redo =", redo, "\n")
    invisible(ans)
}


#dat = read.data()
#cat("chromosomes =", paste(collapse=" ", attr(dat, "chromosomes")), "\n\n")

do.chromosome.plot.simplified = function(dat, chromosomes = attr(dat, "chromosomes")) {
    #dat = subset(dat, ! is.na(twopool_log10_p))
    g = with(read.delim('danRer7.fa.genome', header=TRUE), setNames(length, name))
    pdf("simplified.pdf", width=2.5, height=2.5)
    opa=par(mar=c(2.0, 1.1, 0, 0.3), tcl=-0.2, mgp=c(0.9, 0, 0), las=1, ps=10, cex.axis=0.6, cex.lab=1.2)
    plot.new()
    plot.window(xlim=c(1, max(g)), ylim=c(length(g)-1, 0))
    for (ichr in 1:length(chromosomes)) {
        chr = chromosomes[ichr]
        #cat("generating line for", chr, "\n")
        ychr = ichr - 1#0.5
        sdat = subset(dat, chrom == chr)
        segments(1, ychr, g[chr], ychr, col="grey", lwd=0.5)
        p.neg5 = subset(sdat, twopool_log10_p <= log10(0.00005))$pos
        t.off = 0.20
        if (length(p.neg5)) segments(p.neg5, ychr - t.off, p.neg5, ychr + t.off, col=site2.col, lwd=0.8)
        #p.neg5a = subset(sdat, twopool_log10_p <= log10(0.00001))$pos
        #if (length(p.neg5a)) segments(p.neg5a, ychr - 0.12, p.neg5a, ychr + 0.12, col=site2.col)
        p.neg6 = subset(sdat, twopool_log10_p <= log10(0.000001))$pos
        #if (length(p.neg6)) segments(p.neg6, ychr - 0.25, p.neg6, ychr + 0.25, col=site1.col, lwd=2)
        #if (length(p.neg6)) points(p.neg6, rep(ychr, length(p.neg6)), pch=20, col=site1.col)
        if (length(p.neg6)) points(p.neg6, rep(ychr, length(p.neg6)), pch=21, col="black", bg=site1.col, cex=0.8, lwd=0.5)

        #opdat = subset(sdat, onepool_p <= 0.0001)
        #plot(opdat$pos, rep(1, nrow(opdat)), pch="|", cex=1.2, xlim=xl, ylim=c(1,3), col="#0000FF80", xlab=paste("Chromosome", chr), ylab="", axes=FALSE)
        #tpdat = subset(sdat, twopool_p <= 0.001)
        #points(tpdat$pos, rep(2, nrow(tpdat)), pch="|", cex=1.2, col="#00b00040")
        #tpdat2 = subset(tpdat, twopool_p <= 0.00005)
        #points(tpdat2$pos, rep(3, nrow(tpdat2)), pch="|", cex=1.2, col="#FF0000FF")
        title(ylab="", xlab="Chromosome position (Mbp)")
        xpos = seq(0, 80, 10)
        axis(1, at=xpos * 1000000, labels=xpos, mgp=c(0, -0.1, 0))
        axis(2, lwd=0, at=(1:length(chromosomes)) - 1, labels=substring(chromosomes, 4), mgp=c(0, -0.1, 0))
    }
    text(-9 * 1000000, 0, "B", xpd=NA, font=2)
    dev.off()
    par(opa)
}

do.chromosome.plot = function(chromosomes = attr(dat, "chromosomes"), do.png=FALSE) {
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
