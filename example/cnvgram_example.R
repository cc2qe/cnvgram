# THE SOURCE NEEDS TO POINT TO THE GENOME DRAW SOURCE FILE (genomeDrawSrc.R)
source("/Users/colby/Documents/chgr/code/cnv_draw/src/genomeDrawSrc.R")
setwd("/Users/colby/Documents/chgr/code/cnv_draw/example")

######################################

# YOU MAY HAVE TO ADJUST THE HEIGHT TO GET THE GRAPH TO LOOK NICE
dev.new(width=9, height=4)
par(mar=c(0,0,0,0), oma=c(0,0,0,0))
plotWidth <- 10000
plotHeight <- 100

# YOU MAY HAVE TO ADJUST THE LOWER/RAISE YLIM TO GET THE GRAPH TO LOOK NICE
plot(0, type='n', xlab='', ylab="", axes=F, xlim=c(0, plotWidth), ylim=c(30, plotHeight))

# INPUT IDEOGRAM SPECS HERE
ideo <- NULL
ideo$name <- "13q12.2 - 13q12.3"
ideo$y <- 97	# PROBABLY DON'T NEED TO CHANGE THIS
ideo$height <- 3	# OR THIS
ideo$chr <- "chr13"
ideo$start <- 26500000
ideo$end <- 31000000
ideo$plotWidth <- plotWidth	# DON'T CHANGE THIS

# ADJUST THE VALUES HERE TO CHANGE TICK START/END AND FREQUENCY
drawIdeogram(ideo, 26500000, 31000000, 500000)

# INPUT THE CNV FILE YOU WANT. CHANGE AUTOCOLLAPSE TO "T" TO MAKE IT GROUP CNVS THAT SPAN THE REGION
CnvArray <- importCnvs("example_cnvs.txt", autoCollapse=F)
for (p in CnvArray) {
	drawCnv(p, ideo)
}

# INPUT THE GENE FILE YOU WANT. ADJUST THE NUMBER (35) TO CHANGE THE Y POSITION OF THE FIRST ROW
geneArray <- importGenes("example_genes.txt", 50, sep=4)
for (q in geneArray) {
	drawGene(q)
}

# I'VE COMMENTED OUT THE CODE BELOW, BUT YOU CAN ADD TRANSLOCATIONS TOO, ALBEIT SOMEWHAT MANUALLY
#
#trx <- 72
#breakCoord <- (199976847+199985345)/2
#scaledBrk <- scaleCoord(breakCoord, ideo)
#
#arrows(scaledBrk, trx+1.5, scaledBrk, trx-1.5, length=0, col='black')
#text(scaledBrk+50, trx, "*transloc, t(2;6)", adj=c(0,.5), cex=0.8)
#
#text(scaledBrk-100, trx-5+1.5, '...ATTCAACCATCTATT', adj=c(1, 0.5), cex=0.7, col='royalblue4')
#text(scaledBrk+100, trx-5-1.5, 'TATCAAACCAGTAGG...', adj=c(0, 0.5), cex=0.7, col='red3')
#
#text(scaledBrk-100, trx-5-1.5, '...AAAGGGTTGGCCTC', adj=c(1, 0.5), cex=0.7, col='red3')
#text(scaledBrk+100, trx-5+1.5, 'GACAAATGATGAAGC...', adj=c(0, 0.5), cex=0.7, col='royalblue4')
#
#arrows(scaledBrk-100, trx-5+1.5, scaledBrk+100, trx-5-1.5, length=0)
#arrows(scaledBrk-100, trx-5-1.5, scaledBrk+100, trx-5+1.5, length=0)