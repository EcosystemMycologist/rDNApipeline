
## Open two terminals: one for computer, one for abacus
# ssh iad31@abacus

##One-off database prep
## downloaded from https://maarjam.botany.ut.ee/ganja/?action=bFiles
## note that there is an excel file with the taxonomic information as well.

#rsync -av /Users/dickiei/Downloads/vt_types_fasta_from_05-06-2019.txt iad31@abacus:/share/data/people/ecoMyc/blastdb
#cd /share/data/people/ecoMyc/blastdb
#makeblastdb -in  /share/data/people/ecoMyc/blastdb/vt_types_fasta_from_05-06-2019.txt -dbtype nucl


######################
##One-off data preparation
## On Abacus make directories:
# mkdir /share/data/people/iad31/RamanaAMF2020/
# mkdir /share/data/people/iad31/RamanaAMF2020/dataP/
## From computer, send data files to abacus:
# rsync -av /Users/dickiei/Documents/MGS00318_John_Ramana_Delivery/MGS00318_1/processed/ iad31@abacus:/share/data/people/iad31/RamanaAMF2020/dataP/

## From the abacus side, unzip the files -- navigate to directory then gunzip
# cd /share/data/people/iad31/RamanaAMF2020/dataP/
# gunzip *.gz

## Remove the _001 from names
# for file in *; do [ -f "$file" ] && ( mv "$file" "$(echo $file | sed -e 's/_001//g')" ); done

## remove the underscore, which causes sequencing naming problems later, also removed "processed" for simplicity
# for file in *; do [ -f "$file" ] && ( mv "$file" "$(echo $file | sed -e 's/processed_//g')" ); done

######################
## R based pipeline
## The following can be run in qlogin mode

# qlogin -pe local 20 -q all.q
# R


##The following lines need to be modified for each taxa, along with the last block of lines
codeDir <- "/share/data/people/ecoMyc/PipelineFunctions/"
projectDir <- "/share/data/people/iad31/RamanaAMF2020/"
dataDir <- "dataP/"

dir <- paste(projectDir, dataDir, sep="")
source(paste(codeDir, "PipelineFunctions_abacusV.txt", sep=""))

systemP <- function(x)
    {
    system(x)
    print(x)
    }
runVparse <- function(dir, forwardPrimer, reversePrimer, maxee = 1.0, minseqlength = 200, minsize = 2, paths=NULL, skipMerge = FALSE)
    {
    if(is.null(paths))
        {
        paths <- setPaths()
        }
    systemP(paste("cd", dir))
    if(!skipMerge)
	{
	systemP(paste(paths$usearch, " -fastq_mergepairs ", dir,"*_R1.fastq  -fastqout ", dir,"merged_reads.fq -fastqout_notmerged_fwd ", dir,"notmerged_fwd.fastq -relabel @",sep=""))
    	}
    ## Cutadapt was not working for some reason and after testing was not necessary, so I've commented out the following two lines
    #systemP(paste(paths$cutadapt, " -g ", forwardPrimer, " -o  ", dir,"trimmed_merged_reads.fq  ", dir,"merged_reads.fq", sep="")) 
    #systemP(paste(paths$cutadapt, " -a ", strrevcomp(reversePrimer), " -o ", dir,"trimmed_merged_reads2.fq  ", dir,"trimmed_merged_reads.fq", sep="")) 
    systemP(paste(paths$vsearch," -fastq_filter ",dir,"merged_reads.fq -fastq_maxee ",maxee," -relabel Filt -fastaout  ",dir,"filtered_merged_reads.fa", sep=""))   
    systemP(paste(paths$vsearch," -derep_fulllength  ",dir,"filtered_merged_reads.fa -relabel Uniq --minseqlength ", minseqlength," -sizeout -output  ",dir,"uniques.fa -output ",dir,"uniques.fa", sep= ""))
    ## New unoise functions added to get zero-radius OTUs while still removing chimeras: (see: https://github.com/torognes/vsearch/pull/283 and https://drive5.com/usearch/manual/faq_uparse_or_unoise.html)
    systemP(paste(paths$vsearch, " --cluster_unoise ", dir ,"uniques.fa --centroids ", dir ,"zotus_chim.fa", sep =""))
    systemP(paste(paths$vsearch, " --sortbysize ", dir ,"zotus_chim.fa --output ", dir ,"zotus_sorted.fa", sep=""))
    systemP(paste(paths$vsearch, " --uchime_denovo ", dir ,"zotus_sorted.fa --abskew 16 --nonchimeras ", dir ,"zotus.fa", sep=""))
    
    systemP(paste(paths$usearch, " -cluster_otus ", dir,"uniques.fa -minsize ", minsize," -otus ", dir, "otus.fa -relabel Otu -uparseout ",dir,"uparse.out.txt", sep=""))
    systemP(paste(paths$vsearch," --fastq_filter  ",dir,"merged_reads.fq  --fastaout ",dir,"merged_reads.fa",sep=""))
    systemP(paste(paths$vsearch," --usearch_global ",dir,"merged_reads.fa -db ", dir, "otus.fa -strand plus -id 0.97 --blast6out ", dir, "sequence_match_to_Otus.txt",sep=""))
    
    systemP(paste(paths$vsearch," --usearch_global ",dir,"merged_reads.fa -db ", dir, "zotus.fa -strand plus -id 1 --blast6out ", dir, "sequence_match_to_Zotus.txt",sep=""))
    }


##The following lines need to be modified for each taxa, along with the last block of lines
dir <- paste(projectDir, dataDir, sep="")
database <- "vt_types_fasta_from_05-06-2019.txt"
fPrimer <- "GTGYCAGCMGCCGCGGTAA"  
rPrimer <- "GGACTACNVGGGTWTCTAAT"
minseqlength <- 100
maxee <- 2
dir

## No modifications should be needed in this block of code

##filter reverse reades: 
#for file in *_R2.fastq; do [ -f "$file" ] && ( usearch -fastq_filter $file -fastq_trunclen 150 -fastqout $file); done 
#usearch -fastq_filter JRL93_S60_L001_R2.fastq -fastq_trunclen 150 -fastqout $file JRL93_S60_L001_R2_trunc.fastq
#usearch -fastq_filter *_R2.fastq -fastq_trunclen 150 
#runVparse(dir, fPrimer, rPrimer, minseqlength = minseqlength, maxee = maxee, skipMerge=FALSE)


##RESULTS WERE HORRIBLE -- Very poor merging (< 1%) and many non-AMF.  Need to run with only forward and compare.

#Combine forward reads:

setwd(dir)
files <- list.files()
files <- files[grep("R1.fastq",files)]
for(file in files)
    {
    tq <- readFastq(file)
    tq <- ShortReadQ(sread(tq), quality(tq), BStringSet(paste(rep(file,length(tq)),":",1:length(tq), sep="")))
    if(file==files[1])
        {
        writeFastq(tq, "mergedForward1.fastq", mode="w", compress=FALSE)
        } else {
        writeFastq(tq, "mergedForward1.fastq", mode="a", compress=FALSE)
        }
    print(paste(file, "completed"))
    }
    

runVparseForward <- function(dir, forwardPrimer, reversePrimer, maxee = 1.0, minseqlength = 200, minsize = 2, paths=NULL, skipMerge = FALSE)
    {
    if(is.null(paths))
        {
        paths <- setPaths()
        }
    systemP(paste("cd", dir))
    
    systemP(paste(paths$vsearch," -fastq_filter ",dir,"mergedForward1.fastq -fastq_maxee ",maxee," -relabel Filt -fastaout  ",dir,"filtered_merged_reads.fa", sep=""))   
    systemP(paste(paths$vsearch," -derep_fulllength  ",dir,"filtered_merged_reads.fa -relabel Uniq --minseqlength ", minseqlength," -sizeout -output  ",dir,"uniques.fa -output ",dir,"uniques.fa", sep= ""))
    ## New unoise functions added to get zero-radius OTUs while still removing chimeras: (see: https://github.com/torognes/vsearch/pull/283 and https://drive5.com/usearch/manual/faq_uparse_or_unoise.html)
    systemP(paste(paths$vsearch, " --cluster_unoise ", dir ,"uniques.fa --centroids ", dir ,"zotus_chim.fa", sep =""))
    systemP(paste(paths$vsearch, " --sortbysize ", dir ,"zotus_chim.fa --output ", dir ,"zotus_sorted.fa", sep=""))
    systemP(paste(paths$vsearch, " --uchime_denovo ", dir ,"zotus_sorted.fa --abskew 16 --nonchimeras ", dir ,"zotus.fa", sep=""))
    
    systemP(paste(paths$usearch, " -cluster_otus ", dir,"uniques.fa -minsize ", minsize," -otus ", dir, "otus.fa -relabel Otu -uparseout ",dir,"uparse.out.txt", sep=""))
    systemP(paste(paths$vsearch," --fastq_filter  ",dir,"mergedForward1.fastq  --fastaout ",dir,"merged_reads.fa",sep=""))
    systemP(paste(paths$vsearch," --usearch_global ",dir,"merged_reads.fa -db ", dir, "otus.fa -strand plus -id 0.97 --blast6out ", dir, "sequence_match_to_Otus.txt",sep=""))
    
    systemP(paste(paths$vsearch," --usearch_global ",dir,"merged_reads.fa -db ", dir, "zotus.fa -strand plus -id 1 --blast6out ", dir, "sequence_match_to_Zotus.txt",sep=""))
    }

runVparseForward(dir, fPrimer, rPrimer, minseqlength = minseqlength, maxee = maxee, skipMerge=FALSE)

    
joinedReads <- readFastq(paste(dir, "mergedForward1.fastq", sep=""))
uniquesfa <- readFAST(paste(dir, "uniques.fa", sep=""))
OtusFasta <- readFAST(paste(dir, "otus.fa", sep=""))
ZotusFasta <- readFAST(paste(dir, "zotus.fa", sep=""))


otus <- runBlastOnOtus(dir, OtusFasta, database, cores=20)

## Did not run on resticted database. Note that otuMatches are to a number -- will need to resolve.


uniquesMatches <- match(as.character(sread(joinedReads)), uniquesfa$sequences)
seqData <- data.frame(ID = as.character(id(joinedReads)), uniqueSeq = uniquesfa$descs[uniquesMatches])

otuMatching <- read.table(paste(dir, "sequence_match_to_Otus.txt", sep=""))
zotuMatching <- read.table(paste(dir, "sequence_match_to_Zotus.txt", sep=""))

##plot(otus$identity,otus$identityRestricted, cex=sqrt(as.numeric(otus$lengthRestricted)/300),pch=16,main="Effect of restricting database on identification, sizes proportional to restricted length")

seqData$otu <- otuMatching$V2[match(seqData$ID, as.character(otuMatching$V1))]
seqData$zotu <- otuMatching$V2[match(seqData$ID, zotuMatching$V1)]


seqData$ID <- as.character(seqData$ID)
seqData$cleanID <- unlist(lapply(strsplit(seqData$ID, "_"), function(x) x[1]))
 ## Final reporting code : will need to replace the taxa name in this block.
fungiSeqData <- seqData
fungiOtus <- otus
fungiResM <- table(fungiSeqData$cleanID[!is.na(fungiSeqData$otu)], fungiSeqData$otu[!is.na(fungiSeqData$otu)])
save(file=paste(projectDir, "RamanaDemo_uparse_", Sys.Date(), "_small_3", sep=""), list=c("fungiSeqData", "fungiOtus", "fungiResM"))
save(file=paste(projectDir, "RamanaDemo_uparse_", Sys.Date(), "_full_3", sep=""), list=ls())

## From home computer
rsync -av iad31@abacus:/share/data/people/iad31/RamanaAMF2020/RamanaDemo_uparse_2020-01-29_full_3 /Users/dickiei/
rsync -av iad31@abacus:/share/data/people/iad31/RamanaAMF2020/RamanaDemo_uparse_2020-01-29_small_3 /Users/dickiei/

#########################################################################
### In R on home computer
load('/Users/dickiei/RamanaDemo_uparse_2020-01-29_small_3')

fungiOtus$length <- as.numeric(fungiOtus$length)
fungiOtus$identity <- as.numeric(fungiOtus$identity )
seqData$ID2 <- paste(unlist(lapply(strsplit(as.character(id(joinedReads)), "_"), function(x) x[1])), unlist(lapply(strsplit(as.character(id(joinedReads)), ":"), function(x) x[2])), sep=".")
fungiOtus$abundance <- colSums(fungiResM)[match(fungiOtus$otu, colnames(fungiResM))]
fungiOtus$frequency <- colSums(fungiResM>0)[match(fungiOtus$otu, colnames(fungiResM))]
fungiResMprop <- sweep(fungiResM, 1, rowSums(fungiResM), "/")
fungiOtus$dominance <- colSums(fungiResMprop)[match(fungiOtus$otu, colnames(fungiResMprop))]/fungiOtus$frequency


plot(fungiOtus$identity~fungiOtus$length, cex=sqrt(fungiOtus$abundance)/50)
plot(fungiOtus$dominance~fungiOtus$frequency, log="xy")


fungiOtus[fungiOtus$length > 200 & fungiOtus$identity < 80 & fungiOtus$abundance > 10000,]
##Mite, snail, centipede, and collembola

fungiOtus[fungiOtus$length > 200 & fungiOtus$identity < 93 & fungiOtus$identity > 90 & fungiOtus$abundance > 10000,]
##Interesting -- Diversispora, but very low identity match

OtuInfo <- read.csv("~/Documents/AMF_taxonomy_VTX.csv")
fungiOtus$VT <- unlist(sapply(strsplit(fungiOtus$otuMatch, "_"), "[", 2))
fungiOtus <- merge(fungiOtus, OtuInfo[,!colnames(OtuInfo) == "DNA"], all.x=TRUE)

fungiOtusClean <- fungiOtus[fungiOtus$length > 200 & fungiOtus$identity > 90,]

table(fungiOtusClean$Order[])
sum(fungiOtusClean$abundance)

fungiResMclean <- fungiResM[,colnames(fungiResM)%in%fungiOtusClean$otu]
plot(rowSums(fungiResMclean), rowSums(fungiResMclean>0))

plot(fungiOtusClean$identity~fungiOtusClean$length, cex=sqrt(fungiOtusClean$abundance)/50)
plot(fungiOtusClean$dominance~fungiOtusClean$frequency, log="xy")

fungiOtusClean$Family <- as.factor(as.character(fungiOtusClean$Family))
plot(fungiOtusClean$dominance~fungiOtusClean$frequency, log="xy", pch=16, 
    col=rainbow(10)[as.numeric(fungiOtusClean$Family)], bty="n")
legend(1,2,levels(fungiOtusClean$Family), pch=16, col=rainbow(10), xpd=NA, bty="n")