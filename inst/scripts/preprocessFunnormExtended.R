######################################################
## Extended Functional normalization of the 450k array
## Jean-Philippe Fortin 
## Apr 7 2014
#####################################################


library(matrixStats)

preprocessFunnormExtended <- function(rgSet, nPCs = 2, sex = NULL, verbose = TRUE, design = NULL) {
	minfi:::.isRG(rgSet)
	rgSet <- updateObject(rgSet) ## FIXM: might not KDH: technically, this should not be needed, but might be nice
	if (verbose) 
		cat("[preprocessFunnorm] Mapping to genome\n")
	gmSet <- mapToGenome(rgSet)
	subverbose <- max(as.integer(verbose) - 1L, 0)
	if (verbose) 
		cat("[preprocessFunnorm] Quantile extraction\n")
	extractedData <- .extractFromRGSet450k(rgSet)
	if (is.null(sex)) {
		gmSet <- addSex(gmSet, minfi::getSex(gmSet, cutoff = -3))
		sex <- rep(1L, length(gmSet$predictedSex))
		sex[gmSet$predictedSex == "F"] <- 2L
	}
	rm(rgSet)
	if (verbose) 
		cat("[preprocessFunnorm] Normalization\n")
	CN <- getCN(gmSet)
	gmSet <- .normalizeFunnorm450kExtended(object = gmSet, extractedData = extractedData, sex = sex, nPCs = nPCs, verbose = subverbose, design = design)
	grSet <- ratioConvert(gmSet, type = "Illumina")
	assay(grSet, "CN") <- CN
	grSet@preprocessMethod <- c(preprocessMethod(gmSet), mu.norm = sprintf("Funnorm, nPCs=%s", nPCs))
	return(grSet)
}


.normalizeFunnorm450kExtended <- function(object, extractedData, nPCs, sex, verbose = TRUE, design = NULL) {
	normalizeQuantilesExtended <- function(matrix, indices, sex = NULL, design = NULL) {
		matrix <- matrix[indices, , drop = FALSE]
		## uses probs, model.matrix, nPCS, through scoping)
		oldQuantiles <- t(colQuantiles(matrix, probs = probs))
		if (is.null(sex)) {
			newQuantiles <- .returnFitExtended(controlMatrix = model.matrix, quantiles = oldQuantiles, nPCs = nPCs, design = design)
		} else {
			newQuantiles <- .returnFitBySexExtended(controlMatrix = model.matrix, quantiles = oldQuantiles, nPCs = nPCs, sex = sex, design = design)
		}
		.normalizeMatrix(matrix, newQuantiles)
	}


	indicesList <- .getFunnormIndices(object)
	model.matrix <- .buildControlMatrix450k(extractedData)
	probs <- seq(from = 0, to = 1, length.out = 500)
	Meth <- getMeth(object)
	Unmeth <- getUnmeth(object)
	for (type in c("IGrn", "IRed", "II")) {
		indices <- indicesList[[type]]
		if (length(indices) > 0) {
			if (verbose) 
				cat(sprintf("[normalizeFunnorm450k] Normalization of the %s probes\n", type))
			Unmeth[indices, ] <- normalizeQuantilesExtended(Unmeth, indices = indices, sex = NULL, design = design)
			Meth[indices, ] <- normalizeQuantilesExtended(Meth, indices = indices, sex = NULL, design = design)
		}
	}

	indices <- indicesList[["X"]]
	if (length(indices) > 0) {
		if (verbose) 
			cat("[normalizeFunnorm450k] Normalization of the X-chromosome")
		Unmeth[indices, ] <- normalizeQuantilesExtended(Unmeth, indices = indices, sex = sex, design = design)
		Meth[indices, ] <- normalizeQuantilesExtended(Meth, indices = indices, sex = sex, design = design)
	}

	indices <- indicesList[["Y"]]
	if (length(indices) > 0) {
		if (verbose) 
			cat("[normalizeFunnorm450k] Normalization of the Y-chromosome")
		sex <- as.character(sex)
		levels <- unique(sex)
		nSexes <- length(levels)
		if (nSexes == 2) {
			level1 <- levels[1]
			level2 <- levels[2]
		}
		if (nSexes == 2) {
			if (sum(sex == level1) > 1) {
				Meth[indices, sex == level1] <- preprocessCore::normalize.quantiles(Meth[indices, sex == level1, drop = FALSE])
				Unmeth[indices, sex == level1] <- preprocessCore::normalize.quantiles(Unmeth[indices, sex == level1, drop = FALSE])
			}
			if (sum(sex == level2) > 1) {
				Meth[indices, sex == level2] <- preprocessCore::normalize.quantiles(Meth[indices, sex == level2, drop = FALSE])
				Unmeth[indices, sex == level2] <- preprocessCore::normalize.quantiles(Unmeth[indices, sex == level2, drop = FALSE])
			}
		} else {
			Meth[indices, ] <- preprocessCore::normalize.quantiles(Meth[indices, ])
			Unmeth[indices, ] <- preprocessCore::normalize.quantiles(Unmeth[indices, ])
		}
	}
	assay(object, "Meth") <- Meth
	assay(object, "Unmeth") <- Unmeth
	return(object)
}



.extractFromRGSet450k <- minfi:::.extractFromRGSet450k
.buildControlMatrix450k <- minfi:::.buildControlMatrix450k
.normalizeMatrix <- minfi:::.normalizeMatrix
.getFunnormIndices <- minfi:::.getFunnormIndices



### Return the normalized quantile functions
.returnFitExtended <- function(controlMatrix, quantiles, nPCs, design = NULL) {
	stopifnot(is.matrix(quantiles))
	stopifnot(is.matrix(controlMatrix))
	stopifnot(ncol(quantiles) == nrow(controlMatrix))
	## Fixing potential problems with extreme quantiles
	quantiles[1, ] <- 0
	quantiles[nrow(quantiles), ] <- quantiles[nrow(quantiles) - 1, ] + 1000
	meanFunction <- rowMeans(quantiles)
	res <- quantiles - meanFunction
	controlPCs <- prcomp(controlMatrix)$x[, 1:nPCs, drop = FALSE]
	design1 <- model.matrix(~controlPCs - 1)
	design1 <- cbind(design1, design)
	fits <- lm.fit(x = design1, y = t(res))
	newQuantiles <- meanFunction + t(fits$residuals)
	return(newQuantiles)
}

.returnFitBySexExtended <- function(controlMatrix, quantiles, nPCs, sex, design = NULL) {
	stopifnot(is.matrix(quantiles))
	stopifnot(is.matrix(controlMatrix))
	stopifnot(ncol(quantiles) == nrow(controlMatrix))
	sex <- as.character(sex)
	levels <- unique(sex)
	nSexes <- length(levels)
	if (nSexes == 2) {
		sex1 <- sum(sex == levels[1])
		sex2 <- sum(sex == levels[2])
	} else {
		sex1 <- sum(sex == levels[1])
		sex2 <- 0
	}
	## When normalization should not be performed by sex separately:
	if ((sex1 <= 10) | (sex2 <= 10)) {
		newQuantiles <- .returnFitExtended(controlMatrix = controlMatrix, quantiles = quantiles, nPCs = nPCs, design = design)
	} else {
		quantiles1 <- quantiles[, sex == levels[1]]
		controlMatrix1 <- controlMatrix[sex == levels[1], ]
		if (!is.null(design)) {
			designMatrix1 <- design[sex == levels[1], ]
		} else {
			designMatrix1 <- NULL
		}

		newQuantiles1 <- .returnFitExtended(controlMatrix = controlMatrix1, quantiles = quantiles1, nPCs = nPCs, design = designMatrix1)
		quantiles2 <- quantiles[, sex == levels[2]]
		controlMatrix2 <- controlMatrix[sex == levels[2], ]
		if (!is.null(design)) {
			designMatrix2 <- design[sex == levels[2], ]
		} else {
			designMatrix2 <- NULL
		}
		newQuantiles2 <- .returnFitExtended(controlMatrix = controlMatrix2, quantiles = quantiles2, nPCs = nPCs, design = designMatrix2)
		newQuantiles <- quantiles
		newQuantiles[, sex == levels[1]] <- newQuantiles1
		newQuantiles[, sex == levels[2]] <- newQuantiles2
	}
	return(newQuantiles)
}
