trans <- function(rawdata, method.transform){
	if (method.transform == "squareroot")
		return(sqrt(rawdata))
	else if (method.transform == "log")
		return(log(rawdata))
	else if (method.transform == "PA")
		return((rawdata > 0) + 0) # +0 converts T/F to 1/0
	else if (typeof(method.transform) == "double")
		return(rawdata^method.transform)
	}

computeAverage <- function(expectedprofile.simprof, num.expected){
	return(.Call("computeAverage", expectedprofile.simprof, num.expected))
	}
	
tsComparison <- function(simulatedprofile.simprof, expectedprofile.average, num.simulated, teststatistic){
	pi_matrix <- matrix(data = NA, nrow = 1, ncol = num.simulated)
	piAsExtremeAsTS <- 0
	for (i in 1:num.simulated){
		pi_matrix[i] <- computeTestStatistic(simulatedprofile.simprof[[i]], expectedprofile.average)
		if (pi_matrix[i] >= teststatistic)
			piAsExtremeAsTS <- piAsExtremeAsTS + 1
		}
	pval <- piAsExtremeAsTS/num.simulated
	return(pval)
	}

diveDeep <- function(rawdata, num.expected, num.simulated, method.cluster, method.distance, originaldata, 
						alpha, clust.order, startrow, pMatrix, side, const, silent, increment){
							
	side.startrow <- findNextStartRow(oldstartrow=startrow, mergematrix=clust.order, side=side)
	if(side.startrow > 0){
		if (toupper(side)=="LEFT")
			side.samples <- findleftgroup(startrow=startrow, mergematrix=clust.order)
		else if (toupper(side)=="RIGHT")
			side.samples <- findrightgroup(startrow=startrow, mergematrix=clust.order)
		side.rawdata <- matrix(data = NA, nrow = length(side.samples), ncol = ncol(originaldata))
		for (i in 1:ncol(side.rawdata))
			side.rawdata[,i] <- findsamples(rawdata=originaldata, samples=side.samples, column=i)		
		simprof.results <- simprof.body(rawdata=side.rawdata, num.expected=num.expected, num.simulated=num.simulated, 
			method.cluster=method.cluster, method.distance=method.distance,
			originaldata=originaldata, alpha=alpha, clust.order=clust.order, 
			startrow=side.startrow, pMatrix=pMatrix, currentsamples=side.samples, const=const, silent=silent, increment=increment)
		
		return(simprof.results)
		}
	else{
		simprof.results <- list()
		simprof.results[["samples"]] <- c(-1*side.startrow)
		simprof.results[["pval"]] <- pMatrix
		return(simprof.results)
		}
	}

findNextStartRow <- function(oldstartrow, mergematrix, side){
	if (toupper(side) == "LEFT")
		return (mergematrix[oldstartrow,1])
	else if (toupper (side) == "RIGHT")
		return (mergematrix[oldstartrow,2])
	}

findsamples <- function(rawdata, samples, column){
	temp <- matrix(data = NA, nrow=length(samples), ncol=1)
	for (i in 1:length(samples))
		temp[i,1] <- rawdata[samples[i],column]
	return(temp)
	}
	
computeTestStatistic <- function(rawdata.simprof, expectedprofile.average){
	return(.Call("computeTestStatistic", rawdata.simprof, expectedprofile.average, length(rawdata.simprof[1,])))
	}
	
findComplementaryIndices <- function(n, rawdata.thisgroup, rawdata.otherindices){
	### Generate a list of all possible sample IDs, NA those that are in the other group, remove the NAs, sample from that
	availableindices<-c(1:n)
	for(i in 1:length(rawdata.otherindices))
		availableindices[rawdata.otherindices[i]]<-NA
	availableindices<-na.omit(availableindices)
	return(sample(availableindices, length(rawdata.thisgroup), replace=FALSE))
	}
	
genProfile <- function(rawdata, originaldata, num.expected, method.distance, const, silent, increment, type){
	### "expectedprofile" is a misnomer for these variables:
	### both expected and simulated profiles are generated with genProfile()
	expectedprofile.indices <- list()
	expectedprofile.data 	<- list()
	expectedprofile.simprof <- list()

	for (i in 1:num.expected){
		if (!(silent)){
			if (i %% increment == 0){
				print(paste(type,"iteration",i))
				}
			}
		
		expectedprofile.data[[i]] <- columnPermuter(mat=rawdata)
		expectedprofile.simprof[[i]] <- genSimilarityProfile(expectedprofile.data[[i]], method.distance, const)
		if (!dim(expectedprofile.simprof[[i]])[2]==1)
			expectedprofile.simprof[[i]] <- expectedprofile.simprof[[i]][,order(expectedprofile.simprof[[i]][2,])] ## order the matrix by ranks and keep the distances attached
		}
	return(expectedprofile.simprof)
	}
	
genSimilarityProfile <- function(rawdata.samples, method.distance, const) {
	if (is.function(method.distance))
		rawdata.dist <- method.distance(rawdata.samples)
	else if (method.distance == "braycurtis")
		rawdata.dist <- braycurtis(rawdata.samples, const)
	else
		rawdata.dist <- dist(rawdata.samples, method.distance)
	rawdata.distvec <- as.vector(rawdata.dist)
	rawdata.ranks 	<- rank(rawdata.distvec)
	rawdata.simprof <- rbind(rawdata.distvec, rawdata.ranks)
	return(rawdata.simprof)
	}
	
columnPermuter <- function(mat){
	return(.Call("columnPermuter", mat, nrow(mat), ncol(mat)))
	}
	
### findleftgroup, findrightgroup, and treedive are used for making sense of hclust's $merge results	
findleftgroup <- function(startrow, mergematrix){
	if (mergematrix[startrow,1] < 0) # negative means we've encountered a singleton
		samples <- mergematrix[startrow,1]
	else if (mergematrix[startrow,1] > 0) # positive means we've encountered a subgroup
		samples <- treedive(mergematrix[startrow,1], mergematrix)
	return(-1*samples)
	}
			
findrightgroup <- function(startrow, mergematrix){
	if (mergematrix[startrow,2] < 0) # negative means we've encountered a singleton
		samples <- mergematrix[startrow,1]
	else if (mergematrix[startrow,2] > 0) # positive means we've encountered a subgroup
		samples <- treedive(mergematrix[startrow,2], mergematrix)
	return(-1*samples)
	}			
			
treedive <- function(startrow, mergematrix){
	
	## left
	if (mergematrix[startrow,1] < 0)
		samples.left <- mergematrix[startrow,1]
	else 
		samples.left <- treedive(mergematrix[startrow,1], mergematrix)
		
	## right
	if (mergematrix[startrow,2] < 0)
		samples.right <- mergematrix[startrow,2]
	else
		samples.right <- treedive(mergematrix[startrow,2], mergematrix)
		
	return(append(samples.left, samples.right))
	}

braycurtis <- function(data, const){
	return(as.dist(.Call("braycurtis", data, nrow(data), ncol(data), const)))
	}
	
### this function will put the p-value in a column on the same row as the split it is testing.
pTracker <- function(PMmat, startrow, pval){
	PMmat[startrow,3] <- pval
	return(PMmat)
	}
	
pTrackerMerge <- function(PMmatL, PMmatR){
	for (i in 1:nrow(PMmatL))
		if (is.na(PMmatL[i,3]) && !is.na(PMmatR[i,3]))
			PMmatL[i,3] <- PMmatR[i,3]
	return(PMmatL)
	}
