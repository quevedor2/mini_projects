#!/usr/bin/env Rscript

## Collect arguments
args <- commandArgs(TRUE)
 
## Default setting when no arguments passed
if(length(args) < 4) {
  args <- c("--help")
}
 
## Help section
if("--help" %in% args) {
  cat("
     Calculate area under the curve threshold for CUT&RUN peaks 
 
      Arguments:
			--exp=someValue   - Input AUC values from experiment CUT&RUN
			--ctrl=someValue   - Input AUC values from control CUT&RUN
			--norm=[yes|no]     - Whether to normalize control and experimental files
			--output=someValue   - Output prefix
")
 
  q(save="no")
}

## Parse arguments (we expect the form --arg=value)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
argsL <- as.list(as.character(argsDF$V2))
names(argsL) <- argsDF$V1
invis <- gc(verbose=FALSE) 

## Arg1 default
#if(is.null(args[1])){
if(is.null(argsL$exp) | is.null(argsL$ctrl) | is.null(argsL$output) | is.null(argsL$norm)) {
  stop("Argument is missing!
     Calculate area under the curve threshold for CUT&RUN peaks 
 
      Arguments:
			--exp=someValue   - Input AUC values from experiment CUT&RUN
			--ctrl=someValue   - Input AUC values from control CUT&RUN
			--norm=[yes|no]     - Whether to normalize control and experimental files
			--output=someValue   - Output prefix
")
 
  q(save="no")
}

## Read in the AUC,Num data from the TARGET/EXP file
exp<-read.table(argsL$exp)
expvec<-exp$V1  # AUC
expmax<-exp$V2  # Num
rm(exp)

## If second file is bedgraph, returns NA; else if a threshold, it will be the AUC threshold
suppressWarnings(numtest<-as.numeric(argsL$ctrl))
invis <- gc(verbose=FALSE)

## If 2nd field is a bedgraph, calculate empirical threshold
if(is.na(numtest)){ 
#	print("Ctrl is a file")
        # Read in the CONTROL AUC,Num data
	ctrl<-read.table(argsL$ctrl)
	ctrlvec<-ctrl$V1  # AUC
	ctrlmax<-ctrl$V2  # Num
	rm(ctrl)
	invis <- gc(verbose=FALSE)
	if(argsL$norm=="yes"){  
                ## Calculate peaks of density plots to generate normalization factor
		
                ## 2D distance (?) m = [b-c, a-b]; det(m) / sqrt(sum((b-c) * (b-c)))
		# a = c(nullAUC, obsAUC)
		# b = 0
		# c = 1
		#
		# if(a = c(0.9, 0.9)):  d = 0
		# if(a = c(0.9, 0.95)): d = -0.05
		# if(a = c(0.9, 0.5)):  d = 0.4
		#
                dist2d<-function(a,b,c){
                   v1<- b - c;
                   v2<- a - b;
                   m<-cbind(v1,v2); 
                   d<-det(m)/sqrt(sum(v1*v1))
                }
 		
		# dataframe of sorted AUC (decreasing) [experimental]
		expframe<-data.frame(count=seq(1,0,length=length(expvec)), 
                                     quant=sort(expvec,decreasing=TRUE)/max(expvec), 
				     value=sort(expvec,decreasing=TRUE))
		expframe$diff<-abs(expframe$count-expframe$quant)            		# Difference between null AUC to observed AUC
		expframe<-expframe[expframe$diff > 0.9*max(expframe$diff),]   		# subset where null-obs difference is > 90th quantile 
		expframe$dist<-apply(expframe,1,function(x) dist2d(c(x[1],x[2]),0,1))	# dist2d() determinant matrix to determine 2d distance between nullAUC and obsAUC [a=(nullAUC,obsAUC), b=0, c=1]
		
		# dataframe of sorted AUC (Decreasing) [control]
		ctrlframe<-data.frame(count=seq(1,0,length=length(ctrlvec)), 
				      quant=sort(ctrlvec,decreasing=TRUE)/max(ctrlvec), 
				      value=sort(ctrlvec,decreasing=TRUE))
		ctrlframe$diff<-abs(ctrlframe$count-ctrlframe$quant)
		ctrlframe<-ctrlframe[ctrlframe$diff > 0.9*max(ctrlframe$diff),]
		ctrlframe$dist<-apply(ctrlframe,1,function(x) dist2d(c(x[1],x[2]),0,1))
		
		# If the maxAUC for the maxDistance (nullAUC-obsAUC) is greater than the 90th percentile AUC 
		if(ctrlframe$value[ctrlframe$dist==max(ctrlframe$dist)][1] > sort(ctrlvec)[as.integer(0.9*length(ctrlvec))]){
		  # Set CtrlValue to the maxAUC
		  ctrlvalue<-ctrlframe$value[ctrlframe$dist==max(ctrlframe$dist)][1]
		} else {
		  # Set CtrlValue to the 90th percentile AUC
		  ctrlvalue<-sort(ctrlvec)[as.integer(0.9*length(ctrlvec))] ## Added 7/15/19 to improve memory performance
		}
		
		# Same as above
		if(expframe$value[expframe$dist==max(expframe$dist)][1] > sort(expvec)[as.integer(0.9*length(expvec))]){
		  expvalue<-expframe$value[expframe$dist==max(expframe$dist)][1]
		} else {
		  expvalue<-sort(expvec)[as.integer(0.9*length(expvec))] ## Added 7/15/19 to improve memory performance
		}
		
		# Create a density curve for EXP & CONTROL AUC values less than the ctrlvalue or expvalue
		# constant = expAUC_peak / ctrlAUC_peak
		ctrltest<-density(ctrlvec[ctrlvec <= ctrlvalue]) ## New for SEACR_1.1
		exptest<-density(expvec[expvec <= expvalue]) ## New for SEACR_1.1
		constant<-(exptest$x[exptest$y==max(exptest$y)])/(ctrltest$x[ctrltest$y==max(ctrltest$y)])
		ctrlvec<-ctrlvec*constant
	}
 
        ## Calculate total signal and max signal thresholds
	both<-c(expvec,ctrlvec)
	pctremain<-function(x) {
	  # calculate number of peaks above a given quantile in the expAUC and the combined [expAUC,ctrlAUC] vectors 
	  # and return the ratio between Exp/(Exp,Ctrl)
	  i <- length(expvec)-(ecdf(expvec)(x) * length(expvec))  
          j <- length(both)- (ecdf(both)(x) * length(both))
	  i/j
	}
	
	# Find the quantile for every AUC value
	# Return the AUC value that maximizes the difference between number of peaks for EXP compared to EXP+CTRL 
	x<-sort(unique(both)) ## New for SEACR_1.1
	x0<-x[which(na.omit(pctremain(x[pctremain(x) < 1])) == max(na.omit(pctremain(x[pctremain(x) < 1]))))]  ## New for SEACR_1.1

	# Return the AUC value corresponding to the minimum difference between number of peaks for EXP compared to EXP+CTRL	
	z<-x[x <= x0[1]]	## New for SEACR_1.1
	z2<-z[abs(((pctremain(x0)+min(pctremain(z)))/2)-pctremain(z))==min(abs(((pctremain(x0)+min(pctremain(z)))/2)-pctremain(z)))] ## New for SEACR_1.1

	# Set z0 for minDiff
	if(x0[1]!=z2[1]){  ## Added 7/15/19 to avoid omitting z when x0==z2
		# if minDiff does not equal maxDiff		
		# AUC[<= maxDiff && >= minDiff]
		z<-z[z > z2[1]] ## New for SEACR_1.1
		z0<-z[abs(z-(max(z)-((1/2)*(max(z)-min(z)))))==min(abs(z-(max(z)-((1/2)*(max(z)-min(z))))))] ## New for SEACR_1.1
	}else{  ## Added 7/15/19 to avoid omitting z when x0==z2
		z0<-x0  ## Added 7/15/19 to avoid omitting z when x0==z2
	}  ## Added 7/15/19 to avoid omitting z when x0==z2
	
	## The following code segment was added to avoid spurious high thresholding when the peak of a lower threshold is within 95% of the peak of the maximum threshold
	# Data frame of difference of peaks for each AUC value
	frame<-data.frame(thresh=x[1:(length(x)-1)], 
			  pct=pctremain(x[1:(length(x)-1)]), 
			  diff=abs(diff(pctremain(x))))
	frame<-na.omit(frame)
	
	# Get quantile at which there's a difference between differences
	i<-2
	output<-0
	while(output==0){
		test3<-as.numeric(paste(c(0,".",rep(9,i)),sep="",collapse=""))  # e.g. 0.911 when i=11
		output<-as.numeric(quantile(frame$diff, test3))			# Get quantile for difference between differences at test3 
#		print(output)
		i<-i+1
	}
	
	# Differences where the diff_of_diff is !=0 and less than the positive inflection value
	a<-frame$thresh[frame$diff != 0 & frame$diff < quantile(frame$diff, test3)]
	a0<-a[which(na.omit(pctremain(a[pctremain(a) < 1])) == max(na.omit(pctremain(a[pctremain(a) <  1]))))]
	
	# Repeat the same step as above, except for all values of a less than max(a)
	b<-a[a <= a0[1]]
	b2<-b[abs(((pctremain(a0)+min(pctremain(b)))/2)-pctremain(b))==min(abs(((pctremain(a0)+min(pctremain(b)))/2)-pctremain(b)))]
	
	#
	if(a0[1]!=b2[1]){  ## Added 7/15/19 to avoid omitting b when a0==b2
		b<-b[b > b2[1]]
		b0<-b[abs(b-(max(b)-((1/2)*(max(b)-min(b)))))==min(abs(b-(max(b)-((1/2)*(max(b)-min(b))))))]
	}else{  ## Added 7/15/19 to avoid omitting b when a0==b2
		b0<-a0  ## Added 7/15/19 to avoid omitting b when a0==b2
	}  ## Added 7/15/19 to avoid omitting b when a0==b2
	if(max(na.omit(pctremain(a[pctremain(a) < 1])))/max(na.omit(pctremain(x[pctremain(x) < 1]))) > 0.95){
		x0<-a0
		z0<-b0
	}
	both2<-c(expmax,ctrlmax)
	d<-sort(unique(both2))
	pctremain2<-function(x) 1-(ecdf(expmax)(x)-ecdf(ctrlmax)(x))
	if(length(d[pctremain2(d) > 1]) > 0){
		d0<-min(d[pctremain2(d) > 1])
	}else{
		d0<-1
	}
	invis <- gc(verbose=FALSE)
	fdr<-c(1-pctremain(x0[1]), 1-pctremain(z0[1])) ## New for SEACR_1.1
}else{ ## If 2nd field is numeric, calculate percentile threshold
#	print("Ctrl is numeric")
	test<-ecdf(expvec)(expvec)
	frame<-data.frame(values=expvec, percentile=1-test)
	test2<-ecdf(expmax)(expmax)
	frame2<-data.frame(values=expmax, percentile=1-test2)
	ctrl<-as.vector(as.numeric(paste(0,argsL$ctrl,sep="")))
	x0<-min(frame$values[frame$percentile <= ctrl[1]])
	z0<-min(frame2$values[frame2$percentile <= ctrl[1]])
	d0<-0
	fdr<-ctrl[1] ## New for SEACR_1.1
}
invis <- gc(verbose=FALSE)
write.table(c(x0[1],z0[1],d0[1]), file=paste(argsL$output, ".threshold.txt", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
if(argsL$norm=="yes"){
	write.table(constant, file=paste(argsL$output, ".norm.txt", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE) #Added 7/19/18 to ensure norm value is multiplied by ctrl
}
invis <- gc(verbose=FALSE)
write.table(fdr, file=paste(argsL$output, ".fdr.txt", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE) #Added 5/15/19 to report empirical FDR for threshold detection
