library(data.table)

args<-commandArgs(TRUE)
options(scipen=999)

library(dplyr)

getScriptPath <- function(){
    cmd.args <- commandArgs()
    m <- regexpr("(?<=^--file=).+", cmd.args, perl=TRUE)
    script.dir <- dirname(regmatches(cmd.args, m))
    if(length(script.dir) == 0) stop("can't determine script dir: please call the script with Rscript")
    if(length(script.dir) > 1) stop("can't determine script dir: more than one '--file' argument detected")
    return(normalizePath(script.dir))
}


scriptDir <- getScriptPath()

#print(scriptDir)

infilename <- args[1]

listname <- args[2]

mythreshold <- args[3]


admixtureResults <- read.table(infilename, header=FALSE, sep=" ", stringsAsFactors = FALSE)
inputlist <- read.table(listname, header=FALSE, sep="\t", stringsAsFactors = FALSE)


adminputonly <- merge(admixtureResults,inputlist,by=c("V1"))

tmplist <- inputlist

foundlist <- data.frame()
paste0("Original number of subjects: ", nrow(inputlist))

for(i in 2:6){

	#print(i)
	#mysubset <- subset(adminputonly, i>=0.95)	
	mysubset <- adminputonly[ which(adminputonly[,i]>=mythreshold),]

	if(nrow(mysubset)>0){


		tmplist <- anti_join(as.data.frame(tmplist), as.data.frame(mysubset), by="V1")
		outname<-paste0(infilename,".",i-1,".ids")
		#print(mysubset)
		#print(length(mysubset$V1))
		#print(nrow(tmplist))
		#print(nrow(inputlist))
		keep <- inputlist[inputlist$V1 %in% mysubset$V1, ]
		print(paste0("Ancestry group ", i-1, ", number of subjects: ", length(mysubset$V1)))
		write.table(file=outname, keep[,c(2,3)], row.names=F, col.names=F, sep="\t", quote=F)#FAM FILES
		#write.table(file=outname, keep, row.names=F, col.names=F, sep="\t", quote=F)
	}


}


if(nrow(tmplist)>0){
	print(paste0("Admixed number of subjects remaining: ", nrow(tmplist)))
	outname<-paste0(infilename,".mixed.ids")

	out <- file(outname, "w", encoding="ASCII")
		write.table(file=outname, tmplist[,c(2,3)], row.names=F, col.names=F, sep="\t", quote=F)
	close(out)
}

#head(inputlist)

#tail(adminputonly)

