options(scipen=3)

scaleCoord <- function(x, myIdeo) {
	ideoLength <- myIdeo$end-myIdeo$start
	ratio <- (x-myIdeo$start)/ideoLength	
	scaled <- ratio*plotWidth
	
	return(scaled)
}

# inverse of scaleCoord NOT CORRECT
unscaleCoord <- function(scaled, myChr) {
	chrLength <- myChr$end - myChr$start
	ratio <- scaled/drawLength
	unscaled <- (myChr$start + ratio*chrLength)
	
	return(unscaled)
}

drawIdeogram <- function(myIdeo, tickStart, tickEnd, tickSep) {
	rect(0, myIdeo$y-myIdeo$height/2, plotWidth, myIdeo$y+myIdeo$height/2, col='gray', border='gray')
	text((0+plotWidth)/2, myIdeo$y, myIdeo$name, cex=0.8)
	
	# add tick marks
	tick = tickStart
	while (tick <= tickEnd && tick <= myIdeo$end) {
		scaledTick <- scaleCoord(tick, myIdeo)
		tickText <- paste(format(tick/1000000, nsmall=2), " Mb", sep="")
		text(scaledTick, myIdeo$y+4, tickText, cex=0.7)
		arrows(scaledTick, myIdeo$y+myIdeo$height/2, scaledTick, myIdeo$y+myIdeo$height/2+.8, length=0)
		
		tick <- tick+tickSep
	}
}

importCnvs <- function(cnvFile, autoCollapse=FALSE, height=3, sep=4, fontsize=0.8, dupCol='mediumaquamarine', delCol='indianred1') {
	textTable <- matrix(scan(cnvFile, what="raw", sep="\t", skip=1), ncol=16, byrow=T)
	
	CnvArray <- vector(mode='list')
	for (i in 1:nrow(textTable)) {
		v <- NULL
		v$name <- textTable[i,1]
		v$align <- textTable[i,3]
		v$row <- as.numeric(textTable[i,13])
		v$y <- 97-(ideo$height-height)/2-sep*v$row
		v$height <- height
		v$chr <- textTable[i,7]
		v$start <- as.numeric(textTable[i,8])
		v$end <- as.numeric(textTable[i,9])
		v$whiskStart <- as.numeric(textTable[i,11])
		v$whiskEnd <- as.numeric(textTable[i,12])
		v$event <- textTable[i,4]
		v$geneSpecific <- as.numeric(textTable[i,16])
		v$fontsize <- fontsize
		
		if (v$event=="Copy Loss") {
			v$color=delCol
		}
		else if (v$event=="Copy Gain") {
			v$color=dupCol
		}
		else if (v$event=="Special") {
			v$color=c(dupCol, delCol)
			v$gain=as.numeric(textTable[i,14])
			v$loss=as.numeric(textTable[i,15])
		}
		else {
			v$color='steelblue3'
		}
		
		CnvArray[[i]] <- v
	}
	
	# if autoCollapse on, compress samples that extend beyond ideogram into Special split bar
	if (autoCollapse) {
		collapsedArray <- vector(mode='list')
		numGains <- 0
		numLosses <- 0

		for (i in 1:length(CnvArray)) {
			v <- CnvArray[[i]]			
			
			if (v$start >= ideo$start || v$end <= ideo$end) {
				collapsedArray[[length(collapsedArray)+1]] <- v
			}
			else {
				if (v$event == "Copy Loss") {
					numLosses <- numLosses+1
					#shift all the rows below this item
					for (j in 1:length(CnvArray)) {
						if (CnvArray[[j]]$row >= v$row) {
							CnvArray[[j]]$y <- CnvArray[[j]]$y+sep
						}
					}
					if (length(collapsedArray)>0) {
						for (j in 1:length(collapsedArray)) {
							if (collapsedArray[[j]]$row > v$row) {
								collapsedArray[[j]]$y <- collapsedArray[[j]]$y+sep
							}
						}
					}
				}
				else if (v$event == "Copy Gain") {
					numGains <- numGains+1
					#shift all the rows below this item
					for (j in 1:length(CnvArray)) {
						if (CnvArray[[j]]$row >= v$row) {
							CnvArray[[j]]$y <- CnvArray[[j]]$y+sep
						}
					}
					if (length(collapsedArray)>0) {
						for (j in 1:length(collapsedArray)) {
							if (collapsedArray[[j]]$row > v$row) {
								collapsedArray[[j]]$y <- collapsedArray[[j]]$y+sep
							}
						}
					}
				}
				else { collapsedArray[[length(collapsedArray)+1]] <- v }
			}
		}

		if (numGains==1) { englishGains <- "1 gain" }
		else { englishGains <- paste(numGains, "gains", sep=" ") }
		if (numLosses==1) { englishLosses <- "1 loss" }
		else { englishLosses <- paste(numLosses, "losses", sep=" ") }

		splitCnv <- NULL
		if (numLosses==0) { splitCnv$name <- englishGains }
		else if (numGains==0) { splitCnv$name <- englishLosses }
		else { splitCnv$name <- paste(englishGains, englishLosses, sep=", ") }
		splitCnv$align <- "center"
		splitCnv$y <- 97-(ideo$height-height)/2-sep*1
		splitCnv$height <- height
		splitCnv$chr <- ideo$chr
		splitCnv$start <- ideo$start
		splitCnv$end <- ideo$end
		splitCnv$whiskStart <- -1
		splitCnv$whiskEnd <- -1
		splitCnv$event <- "Special"
		splitCnv$geneSpecific <- 0
		splitCnv$color=c(dupCol, delCol)
		splitCnv$gain=numGains
		splitCnv$loss=numLosses
		
		
		if (numGains+numLosses>0) {
			#shift all the bars down before adding the splitCnv bar
			if (length(collapsedArray)>0) {
				for (j in 1:length(collapsedArray)) {
					collapsedArray[[j]]$y <- collapsedArray[[j]]$y-sep
				}
			}
			collapsedArray[[length(collapsedArray)+1]] <- splitCnv
			
			# rename it
			CnvArray <- collapsedArray
		}
	}
	
	return(CnvArray)
}

drawCnv <- function(myCnv, myIdeo, labels=TRUE, arrLength=0.08) {
	if (myCnv$chr == myIdeo$chr && myCnv$start <= myIdeo$end && myCnv$end >= myIdeo$start) {

		# if special, draw the split rectangle
		if (myCnv$event=="Special") {
			scaledStart <- scaleCoord(myIdeo$start, myIdeo)
			scaledEnd <- scaleCoord(myIdeo$end, myIdeo)
			gainHeight <- myCnv$gain/(myCnv$gain+myCnv$loss)*myCnv$height
			
			#gains
			if (myCnv$gain!=0) {
				rect(scaledStart, myCnv$y+myCnv$height/2, scaledEnd, myCnv$y+myCnv$height/2-gainHeight, col=myCnv$col[1], border=myCnv$col[1])
			}
			
			#losses
			if (myCnv$loss!=0) {
				rect(scaledStart, myCnv$y+myCnv$height/2-gainHeight, scaledEnd, myCnv$y-myCnv$height/2, col=myCnv$col[2], border=myCnv$col[2])
			}
			if (gainHeight > myCnv$height/2) {
				arrowCol=myCnv$col[1]
			}
			else { arrowCol=myCnv$col[2] }
			arrows(scaledStart, myCnv$y, scaledStart-200, myCnv$y, col=arrowCol, length=arrLength)
			arrows(scaledEnd, myCnv$y, scaledEnd+200, myCnv$y, col=arrowCol, length=arrLength)
		}

		else {
			# draw rectangle
			# if beyond bounds, crop
			if (myCnv$start < myIdeo$start) {
				scaledStart <- scaleCoord(myIdeo$start, myIdeo)
				arrows(scaledStart, myCnv$y, scaledStart-200, myCnv$y, col=myCnv$col, length=arrLength)
			}
			else { scaledStart <- scaleCoord(myCnv$start, myIdeo) }
			
			if (myCnv$end > myIdeo$end) {
				scaledEnd <- scaleCoord(myIdeo$end, myIdeo)
				arrows(scaledEnd, myCnv$y, scaledEnd+200, myCnv$y, col=myCnv$col, length=arrLength)
			}
			else { scaledEnd <- scaleCoord(myCnv$end, myIdeo) }
			
			rect(scaledStart, myCnv$y-myCnv$height/2, scaledEnd, myCnv$y+myCnv$height/2, col=myCnv$col, border=myCnv$col)
			if (myCnv$geneSpecific==1) {
				rect(scaledStart, myCnv$y-myCnv$height/2, scaledEnd, myCnv$y+myCnv$height/2, col='gray60', border=F, density=15)
			}
		
			# draw whiskers
			if (myCnv$whiskStart > -1 && myCnv$whiskStart > myIdeo$start) {
				scaledWhiskStart <- scaleCoord(myCnv$whiskStart, myIdeo)
				arrows(scaledStart, myCnv$y, scaledWhiskStart, myCnv$y, col=myCnv$col, angle=90, length=0)
				arrows(scaledWhiskStart, myCnv$y-myCnv$height/2*0.8, scaledWhiskStart, myCnv$y+myCnv$height/2*0.8, col=myCnv$col, angle=90, length=0)
			}
			# if exceeds bounds
			if (myCnv$whiskStart > -1 && myCnv$whiskStart < myIdeo$start) {
				scaledWhiskStart <- scaleCoord(myIdeo$start, myIdeo)
				arrows(scaledStart, myCnv$y, scaledWhiskStart, myCnv$y, col=myCnv$col, angle=90, length=0)
			}
			
			if (myCnv$whiskEnd > -1 && myCnv$whiskEnd < myIdeo$end) {
				scaledWhiskEnd <- scaleCoord(myCnv$whiskEnd, myIdeo)
				arrows(scaledEnd, myCnv$y, scaledWhiskEnd, myCnv$y, col=myCnv$col, angle=90, length=0)
				arrows(scaledWhiskEnd, myCnv$y-myCnv$height/2*0.8, scaledWhiskEnd, myCnv$y+myCnv$height/2*0.8, col=myCnv$col, angle=90, length=0)
			}
			# if exceeds bounds
			if (myCnv$whiskEnd > -1 && myCnv$whiskEnd > myIdeo$end) {
				scaledWhiskEnd <- scaleCoord(myIdeo$end, myIdeo)
				arrows(scaledEnd, myCnv$y, scaledWhiskEnd, myCnv$y, col=myCnv$col, angle=90, length=0)
			}
		}
		
		if (labels) {
			# label it
			if (myCnv$align=="center") {
				text((scaledStart+scaledEnd)/2, myCnv$y, myCnv$name, cex=myCnv$fontsize)
			}
			else if (myCnv$align=="right") {
				if (myCnv$whiskEnd == -1) { scaledWhiskEnd=scaledEnd }
				text(scaledWhiskEnd+40, myCnv$y, myCnv$name, cex=myCnv$fontsize, adj=c(0,0.5))
			}
			else if (myCnv$align=="left") {
				if (myCnv$whiskStart == -1) { scaledWhiskStart=scaledStart }
				text(scaledWhiskStart-40, myCnv$y, myCnv$name, cex=myCnv$fontsize, adj=c(1,0.5))
			}
		}
	}
}

importGenes <- function(geneFile, top, sep=4, height=3, fontsize=0.8) {
	textTable <- matrix(scan(geneFile, what="raw", sep="\t", skip=1), ncol=18, byrow=T)
	
	geneArray <- vector(mode='list')

	for (i in 1:nrow(textTable)) {	
		GENE <- NULL
		GENE$name <- textTable[i,13]
		GENE$y <- top-sep*as.numeric(textTable[i,17])
		GENE$height <- height
		GENE$start <- as.numeric(textTable[i,5])
		GENE$end <- as.numeric(textTable[i,6])
		GENE$strand <- textTable[i,4]
		GENE$exonStarts <- unlist(strsplit(textTable[i,10], split=","))
		GENE$exonEnds <- unlist(strsplit(textTable[i,11], split=","))
		GENE$align <- textTable[i,18]
		GENE$fontsize <- fontsize
			
		geneArray[[i]] <- GENE
	}
	return(geneArray)
}

drawGene <- function(myGene) {
	
	# draw line
	scaledStart <- scaleCoord(myGene$start, ideo)
	scaledEnd <- scaleCoord(myGene$end, ideo)
	arrows(scaledStart, myGene$y, scaledEnd, myGene$y, col="royalblue4", length=0)
	
	# draw arrows
	width <- 23.33*myGene$height
	height <- 0.4*myGene$height
	spacing <- 2*width
	
	# if neg strand, the flip the arrows
	if (myGene$strand == "-") {
		width = -width
	}
	
	x <- scaledStart+abs(width)
	while (x+abs(width/2) < scaledEnd) {
		polygon(c(x-width/2, x+width/2, x-width/2), c(myGene$y-height/2, myGene$y, myGene$y+height/2), col="royalblue4", border="royalblue4")
		x <- x+spacing
	}
	
	# draw exons
	for (i in 1:length(myGene$exonStarts)) {
#		exon <- myGene$exons[i,]
		exonStart <- scaleCoord(as.numeric(myGene$exonStarts[i]), ideo)
		exonEnd <- scaleCoord(as.numeric(myGene$exonEnds[i]), ideo)

		rect(exonStart, myGene$y-myGene$height/2, exonEnd, myGene$y+myGene$height/2, col='royalblue4', border='royalblue4')
	}
	
	# label it
	if (myGene$align=="center") {
		text((scaledStart+scaledEnd)/2, myGene$y-3, myGene$name, cex=myGene$fontsize)
	}
	else if (myGene$align=="right") {
		text(scaledEnd+50, myGene$y, myGene$name, cex=myGene$fontsize, adj=c(0,0.5))
	}
	else if (myGene$align=="left") {
		text(scaledStart-50, myGene$y, myGene$name, cex=myGene$fontsize, adj=c(1,0.5))
	}
}





