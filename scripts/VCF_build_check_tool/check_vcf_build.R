library(data.table)

args<-commandArgs(TRUE)


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

#Right now the comparison uses only SNPs that belong to the positive strand of hgTable files, SNPs from negative strand of hgTables are not compared to avoid complications
#Todo: Development of function for detecting flips 
#When flip detection is done, use the function bellow to fix flip by generating reverse complement
convertToComplement<-function(x){
	bases=c("A","C","G","T")
	xx<-unlist(strsplit(toupper(x),NULL))
	paste(unlist(lapply(xx,function(bbb){
		if(bbb=="A") compString<-"T"
		if(bbb=="C") compString<-"G"
		if(bbb=="G") compString<-"C"
		if(bbb=="T") compString<-"A"
		if(!bbb %in% bases) compString<-"N"
		return(compString)
	})),collapse="")
}


filename <- args[1]

command <- paste0("grep -m 100000 -w \"^1\\|^chr1\" ",filename, " | cut -f 1,2,3,4,5")
#print(command) 

print(paste0("Reading first 100K markers from chr1 in input file: ", filename, "..."))
input <- fread(cmd=command)
input <- as.data.frame(input)
names(input) <- c("chr", "pos", "originalID", "ref", "alt")

#head(input)

hgtables <- list()
table_list <- c("chr1_GRCh37-hg19.hgTable.bi", "chr1_GRCh38-hg38.hgTable.bi", "chr1_NCBI34-hg16.hgTable.bi", "chr1_NCBI35-hg17.hgTable.bi", "chr1_NCBI36-hg18.hgTable.bi")
ref_names <- c("GRCh37/hg19", "GRCh38/hg38", "NCBI34/hg16", "NCBI35/hg17", "NCBI36/hg18")
chain_names1 <- c("none, already GRCh37", "GRCh38_to_GRCh37.chain", "NCBI34_to_GRCh37.chain", "NCBI35_to_GRCh37.chain", "NCBI36_to_GRCh37.chain")
chain_names2 <- c("hg19_to_GRCh37.chain", "hg38_to_GRCh37.chain", "hg16_to_hg19.chain -> hg19_to_GRCh37.chain", "hg17_to_hg19.chain -> hg19_to_GRCh37.chain", "hg18_to_hg19.chain -> hg19_to_GRCh37.chain")

#hg -> chr1
#GRCh -> 1

i<-1

Nmatches <- NULL

myheader <- c("chr", "start","pos", "rs", "strand", "ref", "alt")

for(hgtablename in table_list){

	print(paste0("Reading ", hgtablename, " from ", scriptDir))

	hgtablename <- paste0(scriptDir,"/",hgtablename)

	hgtables[[i]] <- read.table(hgtablename, header=FALSE, sep="\t")
	names(hgtables[[i]]) <- myheader

	hgtables[[i]] <- subset(hgtables[[i]], strand == "+")

	print(paste0("Identifying matches between ", hgtablename, " and ", filename))
	matches <- merge(input, hgtables[[i]], by=c("pos", "ref", "alt"))

	Nmatches[i] <- nrow(matches)

	print(paste0("Exact pos/ref/alt matches found: ", filename, " versus ", hgtablename, " = ", Nmatches[i]))

	i<-i+1

}

maxindex <- which(Nmatches==max(Nmatches))

others <- sum(Nmatches[-maxindex])

print(paste0("Dataset is based on build: ", ref_names[maxindex], ". Exact pos/ref/alt matches: ", Nmatches[maxindex], ". Total accidental matches from other builds, summed together: ", others) )


#hg -> chr1
#GRCh -> 1


if(input$chr[1]=="1"){

print(paste0("Use chain file(s): ", chain_names1[maxindex]))


}

if(input$chr[1]=="chr1"){

print(paste0("Use chain file(s): ", chain_names2[maxindex]))

}

