##requirements:
## instalation of vsearch, Usearch, cutadapt, blast
## vsearch here: https://github.com/torognes/vsearch
## usearch here: http://drive5.com/usearch/manual/install.html
## cutadapt here: http://cutadapt.readthedocs.io/en/stable/installation.html
## blast here: ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.4.0/ncbi-blast-2.4.0+.dmg
## blast here: https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download
##
## I'm pretty sure blast will install makeblastdb as part of the package, but we will need the whole thing.

library(ShortRead)
library("parallel")
library("foreach")
library("doParallel")
library(vegan)

############################################
## General purpose DNA manipulation functions
############################################

writeFast <- function(seqs, desc, file="Temp.fasta")
    #Much faster function to write a fasta-like format without any line breaks within sequences
    {
    fasta <- c()
    seqs <- as.character(seqs)
    desc <- as.character(desc)
    if(substr(desc[1],1,1)!=">")
        {
        desc <- paste(">", desc, sep="")
        }
    fasta[seq(1,length(seqs)*2, by=2)] <- desc
    fasta[seq(2,length(seqs)*2, by=2)] <- seqs
    write(fasta, file=file)
    }
    
readFAST <- function(file)
    { ##the function readFASTA in Biostrings is slow as based on a loop. This is faster...but not by tons (10% faster)
    in.file <- readLines(file,  encoding="latin1")
    in.file <- iconv(in.file, from="latin1", to="ASCII", ".")
    desc.lines <- grep(">", in.file)

    seq.lines <- data.frame(desc.lines+1, stop.line = c(desc.lines[-1]-1, length(in.file)))
    sequence.line.list <- apply(seq.lines, 1, function(x) seq(x[1], x[2]))
    if(class( sequence.line.list) == "matrix")
  	{
    	sequence.line.list <- split(t(sequence.line.list), f=in.file[desc.lines])
    	descs <- names(sequence.line.list)
  	} else {
    	descs <- gsub(">", "", in.file[desc.lines])
  	}
    sequences <- lapply(sequence.line.list, function(x) paste(in.file[x], collapse=""))
    list(descs = descs, sequences = sequences)
            }
            
iupacGrep <- function(pattern,x,...)
    #Convenience package, this uses regular expressions to run grep with IUPAC codes
    {
    pattern <- gsub("R", "[A|G]", pattern)
    pattern <- gsub("Y", "[C|T]", pattern)
    pattern <- gsub("S", "[G|C]", pattern)
    pattern <- gsub("W", "[A|T]", pattern)
    pattern <- gsub("K", "[G|T]", pattern)
    pattern <- gsub("M", "[A|C]", pattern)
    pattern <- gsub("B", "[C|G|T]", pattern)
    pattern <- gsub("D", "[A|G|T]", pattern)
    pattern <- gsub("H", "[A|C|T]", pattern)
    pattern <- gsub("V", "[A|C|G]", pattern)
    pattern <- gsub("N", "[A|C|T|G]", pattern)
    pattern <- gsub("I", "[A|C|T|G]", pattern)  ##Note that Deoxyinosine is not a true universal match, see http://sg.idtdna.com/site/Catalog/Modifications/Category/7
    grep(pattern, x, ...)
    }

iupacSplit <- function(x, pattern,...)
    #Convenience package, this uses regular expressions to run grep with IUPAC codes
    {
    pattern <- gsub("R", "[A|G]", pattern)
    pattern <- gsub("Y", "[C|T]", pattern)
    pattern <- gsub("S", "[G|C]", pattern)
    pattern <- gsub("W", "[A|T]", pattern)
    pattern <- gsub("K", "[G|T]", pattern)
    pattern <- gsub("M", "[A|C]", pattern)
    pattern <- gsub("B", "[C|G|T]", pattern)
    pattern <- gsub("D", "[A|G|T]", pattern)
    pattern <- gsub("H", "[A|C|T]", pattern)
    pattern <- gsub("V", "[A|C|G]", pattern)
    pattern <- gsub("N", "[A|C|T|G]", pattern)
    pattern <- gsub("I", "[A|C|T|G]", pattern)  ##Note that Deoxyinosine is not a true universal match, see http://sg.idtdna.com/site/Catalog/Modifications/Category/7
    strsplit( x, pattern, ...)
    }
    
iupacAgrep <- function(pattern, x, ...)
    #Convenience package, this uses regular expressions to run grep with IUPAC codes
    {
    pattern <- gsub("R", "[A|G]", pattern)
    pattern <- gsub("Y", "[C|T]", pattern)
    pattern <- gsub("S", "[G|C]", pattern)
    pattern <- gsub("W", "[A|T]", pattern)
    pattern <- gsub("K", "[G|T]", pattern)
    pattern <- gsub("M", "[A|C]", pattern)
    pattern <- gsub("B", "[C|G|T]", pattern)
    pattern <- gsub("D", "[A|G|T]", pattern)
    pattern <- gsub("H", "[A|C|T]", pattern)
    pattern <- gsub("V", "[A|C|G]", pattern)
    pattern <- gsub("N", "[A|C|T|G]", pattern)
    pattern <- gsub("I", "[A|C|T|G]", pattern)  ##Note that Deoxyinosine is not a true universal match, see http://sg.idtdna.com/site/Catalog/Modifications/Category/7
    agrep(pattern, x, ...)
    }
    
#DNA functions to reverse or reverse-complement strings
strrev <- function(x) paste(rev(strsplit(x, "")[[1]]), collapse = "") #Reverse any string of letters

strcomp <- function(x) paste(c("A","T","G","C","Y","R","S","W","K","M","B","D","H","V","N", ".")[match(strsplit(x, "")[[1]], 
                                c("T","A","C","G","R","Y","S","W","M","K","V","H","D","B","N", "." ))], collapse = "")

strrevcomp <- function(x) paste(c("A","T","G","C","Y","R","S","W","K","M","B","D","H","V","N", ".")[match(rev(strsplit(x, "")[[1]]), 
                                c("T","A","C","G","R","Y","S","W","M","K","V","H","D","B","N", "." ))], collapse = "")

CollapseHomopolys <- function(x)
    ##Largely obsolete function to get rid of homopolymers for 454 analysis
    {
    x <- gsub("T{3,}", paste(rep("T",3), collapse=""), x)
    x <- gsub("C{3,}", paste(rep("C",3), collapse=""), x)
    x <- gsub("G{3,}", paste(rep("G",3), collapse=""), x)
    x <- gsub("A{3,}", paste(rep("A",3), collapse=""), x)
    x
    }

extractFASTA <- function(file, ids, outfile="output.fsa")
    ##Not sure if this is used. It lacked good documentation, but
    ##I think it is primarily for dealing with odd characters in FASTA file headers
    {
    in.file <- readLines(file,  encoding="latin1")
    in.file <- iconv(in.file, from="latin1", to="ASCII", ".")
    desc.lines <- grep(">", in.file)
    seq.lines <- data.frame(desc.lines+1, stop.line = c(desc.lines[-1]-1, length(in.file)))
    sequence.line.list <- apply(seq.lines, 1, function(x) seq(x[1], x[2]))	
    write.lines<-c()
    for(get.i in 1:length(ids))
        {
        write.lines <- append(write.lines, paste(">", ids[get.i]))
        write.lines <- append(write.lines, in.file[sequence.line.list[[grep(ids[get.i], in.file[desc.lines])]]])
        }
    write(write.lines, file=outfile)
    }

### Pipeline functions

setPaths <- function(user = Sys.info()[7])
    {
    ##Convenience function used to set the paths to software based on the user.
    ##Define the code necessary to run usearch, vsearch, cutadapt, and blast plus the paths for blast database directory.
    ##Pretty specific to my computers and computers of my colleagues who use my pipeline. Note that the
    ##final "else" set is for the University of Canterbury Abacus cluster.
     if(user=="dickiei")
        {
        paths <- list(usearch="usearch",
                      vsearch="~/vsearch-2.3.0-osx-x86_64/bin/vsearch",
                      cutadapt = "cutadapt", 
                      blast = "/Users/dickiei/BLAST/bin/blastn",
                      dbPath = "/Users/dickiei/BLAST/db/")
        system("export PATH=$PATH:~/vsearch-2.3.0-osx-x86_64/bin/")  ##This seemed necessary on my MacPro, probably due to
                                                                     ##mistake in installation of vsearch.
        }
    if(user=="andi")
        {
        paths <- list(usearch="~/NGSprograms/usearchv9/usearch9.0.2132_i86linux32",
                      vsearch="~/NGSprograms/vsearch-2.3.0/bin/vsearch",
                      cutadapt = "~/.local/bin/cutadapt", 
                      blast = "~/BLAST/ncbi-blast-2.5.0+/bin/blastn",
                      dbPath = "~/BLAST/ncbi-blast-2.5.0+/db/")
        }
     else 
        {
        paths <- list(usearch="/share/apps/usearch/usearch",
                      vsearch="/share/apps/vsearch/bin/vsearch",
                      cutadapt="/usr/bin/cutadapt",
                      blast="/usr/bin/blastn",
                      dbPath="/data/people/ecoMyc/blastdb/")
        }
     return(paths)
     }
     
movefiles <- function(from, to, rootdir = "")
    ## Move files from one directory to another
    ## from is a vector of the directories where raw sequences are currently held
    ## to the directory where sequences should be put.
    ## root dir is pasted onto the front of from and to.
    {
    for(i in from)
        {
        system(paste("cd ", rootdir, i, sep=""))
        system(paste("gzip -d ", rootdir, i,"/*.gz", sep=""))
        system(paste("mv ", rootdir, i,"/*.fastq ", rootdir, to, sep=""))
        }
    }

systemP <- function(x)
    { ## Runs any command x but first prints the command to screen. Useful if using paste to build commands.
    print(x)
    system(x)
    }
    
runVparse <- function(dir, forwardPrimer, reversePrimer, maxee = 1.0, minseqlength = 200, minsize = 2, paths=NULL, skipMerge = FALSE)
    {  
    ## The main pipeline function in my code, this will run: usearch -fastq_mergepairs, vsearch -fastq_filter, vsearch -derep_fulllength
    ## and then form both Zero-radius OTUs (zotus) and OTUs.
    if(is.null(paths))
        {
        paths <- setPaths()
        }
    systemP(paste("cd", dir))
    if(!skipMerge)
	{
	systemP(paste(paths$usearch, " -fastq_mergepairs ", dir,"*_R1.fastq  -fastqout ", dir,"merged_reads.fq -fastqout_notmerged_fwd ", dir,"notmerged_fwd.fastq -relabel @",sep=""))
    	}
    systemP(paste(paths$vsearch," -fastq_filter ",dir,"merged_reads.fq -fastq_maxee ",maxee," -relabel Filt -fastaout  ",dir,"filtered_merged_reads.fa", sep=""))   
    systemP(paste(paths$vsearch," -derep_fulllength  ",dir,"filtered_merged_reads.fa -relabel Uniq --minseqlength ", minseqlength," -sizeout -output  ",dir,"uniques.fa -output ",dir,"uniques.fa", sep= ""))
    ## New unoise functions added to get zero-radius OTUs while still removing chimeras: (see: https://github.com/torognes/vsearch/pull/283 and https://drive5.com/usearch/manual/faq_uparse_or_unoise.html)
    systemP(paste(paths$vsearch, " --cluster_unoise ", dir ,"uniques.fa --centroids ", dir ,"zotus_chim.fa", sep =""))
    systemP(paste(paths$vsearch, " --sortbysize ", dir ,"zotus_chim.fa --output ", dir ,"zotus_sorted.fa", sep=""))
    systemP(paste(paths$vsearch, " --uchime_denovo ", dir ,"zotus_sorted.fa --abskew 16 --nonchimeras ", dir ,"zotus.fa", sep=""))
    ## cluster Otus
    systemP(paste(paths$usearch, " -cluster_otus ", dir,"uniques.fa -minsize ", minsize," -otus ", dir, "otus.fa -relabel Otu -uparseout ",dir,"uparse.out.txt", sep=""))
    systemP(paste(paths$vsearch," --fastq_filter  ",dir,"merged_reads.fq  --fastaout ",dir,"merged_reads.fa",sep=""))
    systemP(paste(paths$vsearch," --usearch_global ",dir,"merged_reads.fa -db ", dir, "otus.fa -strand plus -id 0.97 --blast6out ", dir, "sequence_match_to_Otus.txt",sep=""))
    ## match sequences to Otus
    systemP(paste(paths$vsearch," --usearch_global ",dir,"merged_reads.fa -db ", dir, "zotus.fa -strand plus -id 1 --blast6out ", dir, "sequence_match_to_Zotus.txt",sep=""))
    }
    
runBlastOnOtus <- function(dir, OtusFasta, database, cores=NULL)
    {
    #Take OTUs and identify their best match using BLASTn.
    if(is.null(cores))
        {
        cores <- detectCores()-1
        }
    if(length(OtusFasta$descs) < 50)
        {
        cores <- 1
        }
    blockstarts <- floor(seq(1, length(OtusFasta$descs), length.out=cores+1))[-(cores+1)]
    blockends <- c(blockstarts[-1]-1,length(OtusFasta$descs))
    print(paste("Running", cores, "cores."))
   
    for(core.i in 1:cores)
        {
        writeFast(unlist(OtusFasta$sequences[blockstarts[core.i]:blockends[core.i]]),OtusFasta$descs[blockstarts[core.i]:blockends[core.i]], file=paste(dir, "tempFasta",core.i,".fasta", sep=""))
        }
    print(paste("file divided into N =",cores,"parts as tempFastaN.fasta in",dir))
    cl <- makeCluster(cores)
    registerDoParallel(cl, cores=cores)
    blastpath <- setPaths()$blast
    dbpath <-  setPaths()$dbPath
    tempOtuBlast <- foreach(core.i = 1:cores, .packages = c("stats"), .combine=rbind) %dopar% 
        {
        system(paste(blastpath, " -query ", dir, "tempFasta",core.i,
        ".fasta -db ", dbpath, database," -task blastn -dust no -outfmt 7 -max_target_seqs 1 -out ", dir,"OTU_Blast_matching",
        core.i, ".txt", sep=""))
        }
    closeAllConnections()
    print("Blast run complete, wrapping up dataFrame")
    df <- collectBlastRes(dir, OtusFasta, database, cores=cores)
    return(df)
    }

collectBlastRes <- function(dir, OtusFasta, database, cores=NULL)
    ## If the  BLAST run completes but computer crashes before results saved, this can be run independently to recover results
    {
    temp <- c()
     if(is.null(cores))
        {
        cores <- detectCores()-1
        }
    if(length(OtusFasta$descs) < 50)
        {
        cores <- 1
        }
       
    for(core.i in 1:cores)
        {  
        temp <- append(temp, readLines(paste(dir, "OTU_Blast_matching", core.i, ".txt", sep="")))
        }
    matches <- grep("# 1 hits found", temp)+1
    df <- (data.frame(matrix(unlist(strsplit(temp[matches], "\t")), byrow=TRUE, ncol=length(unlist(strsplit(temp[matches][1], "\t")))), stringsAsFactors=FALSE)[,c(1:4,11,12)])
    names(df) <- c("otu", "otuMatch", "identity", "length", "evalue","bitscore")
    df$sequence <- OtusFasta$sequences[match(df$otu, OtusFasta$descs)] 
    return(df)
    }

expandNames <- function(otus, nameColumn = "otuMatchRestricted")
    {
    otuNames <- otus[,colnames(otus)==nameColumn]
    otus$nearestKingdom <-  as.factor(unlist(sapply(strsplit(otuNames, "[;|]"), 
        function(x) gsub("k__", "", ifelse(length(x[grep("k__", x)])>0, x[grep("k__", x)], NA)))))
    otus$nearestPhylum <-  as.factor(unlist(sapply(strsplit(otuNames, "[;|]"), 
        function(x) gsub("p__", "", ifelse(length(x[grep("p__", x)])>0, x[grep("p__", x)], NA)))))
    otus$nearestClass <-  as.factor(unlist(sapply(strsplit(otuNames, "[;|]"), 
        function(x) gsub("c__", "", ifelse(length(x[grep("c__", x)])>0, x[grep("c__", x)], NA)))))
    otus$nearestOrder <-  as.factor(unlist(sapply(strsplit(otuNames, "[;|]"), 
        function(x) gsub("o__", "", ifelse(length(x[grep("o__", x)])>0, x[grep("o__", x)], NA)))))
    otus$nearestFamily <-  as.factor(unlist(sapply(strsplit(otuNames, "[;|]"), 
        function(x) gsub("f__", "", ifelse(length(x[grep("f__", x)])>0, x[grep("f__", x)], NA)))))
    otus$nearestGenus <-  as.factor(unlist(sapply(strsplit(otuNames, "[;|]"), 
        function(x) gsub("g__", "", ifelse(length(x[grep("g__", x)])>0, x[grep("g__", x)], NA)))))
    otus$nearestSpecies <-  as.factor(unlist(sapply(strsplit(otuNames, "[;|]"), 
        function(x) gsub("s__", "", ifelse(length(x[grep("s__", x)])>0, x[grep("s__", x)], NA)))))
    return(otus)
    }

expandBacteriaNames <- function(otus, database, nameColumn = "otuMatchRestricted")
    {
    otuNames <- otus[,colnames(otus)==nameColumn]
    Bacteria <- readFAST(paste(setPaths()$dbPath,database, sep=""))
    ##Note that the above is a BIG file.
    bacteriaDescList <- strsplit(Bacteria$descs, " ")
    remove(list="Bacteria")  ##To save working memory, I remove the bacterial data after using it.
    bacteriaData <- cbind(OTU = sapply(bacteriaDescList, function(x) x[1]), 
                            genBank = sapply(bacteriaDescList, function(x) x[2]),
                            kingdom = sapply(bacteriaDescList, function(x) gsub("k__", "", x[grep("k__", x)])),
                            phylum = sapply(bacteriaDescList, function(x) gsub("p__", "", x[grep("p__", x)])),
                            class = sapply(bacteriaDescList, function(x) gsub("c__", "", x[grep("c__", x)])),
                            order = sapply(bacteriaDescList, function(x) gsub("o__", "", x[grep("o__", x)])),
                            family = sapply(bacteriaDescList, function(x) gsub("f__", "", x[grep("f__", x)])),
                            genus = sapply(bacteriaDescList, function(x) gsub("g__", "", x[grep("g__", x)])),
                            species = sapply(bacteriaDescList, function(x) gsub("s__", "", x[grep("s__", x)])))

    bacteriaData <- apply(bacteriaData, 2, function(x) gsub(";", "", x))
    bacteriaData <- apply(bacteriaData, 2, function(x) gsub("character(0)", NA, x))
    bacteriaData <- as.data.frame(bacteriaData)
    otus$nearestKingdom <- bacteriaData$kingdom[match(otuNames, bacteriaData$OTU)]
    otus$nearestPhylum <- bacteriaData$phylum[match(otuNames, bacteriaData$OTU)]
    otus$nearestClass <- bacteriaData$class[match(otuNames, bacteriaData$OTU)]
    otus$nearestOrder <- bacteriaData$order[match(otuNames, bacteriaData$OTU)]
    otus$nearestFamily <- bacteriaData$family[match(otuNames, bacteriaData$OTU)]
    otus$nearestGenus <- bacteriaData$genus[match(otuNames, bacteriaData$OTU)]
    otus$nearestSpecies <- bacteriaData$species[match(otuNames, bacteriaData$OTU)]
    return(otus)
}


subsetDB <- function(db, out, remove="f__(unid)|(Incert)", require=NULL)
    ##This last function is just moved here from another set of code, as uses the function "makeblastdb" which we will
    ## need installed on cluster...
    {
    fullDB <- readFAST(db)
    if(!is.null(require))
        {
        fullDB$sequences <- fullDB$sequences[grep(require,fullDB$descs)]
        fullDB$descs <- fullDB$descs[grep(require,fullDB$descs)]
        }
    writeFast(fullDB$sequences[!grepl(remove, fullDB$descs)], fullDB$descs[!grepl(remove, fullDB$descs)], file=out)
    system(paste("makeblastdb -in ", out," -dbtype nucl "))    
    }
    
