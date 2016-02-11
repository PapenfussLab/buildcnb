#!/usr/bin/Rscript --vanilla
args <- commandArgs(TRUE)

# For testing
RSTUDIO_DEBUG <- FALSE
if(RSTUDIO_DEBUG){ args <- "/home/jmarkham/mm/cnv/buildcnb/super_batch_1_2/buildcnb_input.dcf" }

isCommandLine <- !RSTUDIO_DEBUG && args[1]!="RStudio" &&  !grepl("R$",args)

if(length(args)!=1 && isCommandLine) {
  print("Builds files for cnb. Usage: makeCNTables.R /path/to/config_file.dcf")
  quit(save = "no", status = 0)
} 

library(futile.logger)
library(RCurl)
library(httr)
library(Rsubread)
library(Biostrings)
library(rtracklayer)
library(ggplot2)

if(!exists("gv")) {
  session <- NULL
  gv <- list()
  gv$log <- ""
  isolate <- function (f) { return (f)}
} # For stand-alone

# Some globals
wg_coarse_bin_size <- 1000000
wg_medium_bin_size <- 50000
wg_fine_bin_size <- 5000
loess_scales_badly_limit  <- 100000 # when loess gets slow
loess_min_usable_counts  <- 50 # For fitting
loess_min_usable_gc  <- 0.25 
loess_span <- 0.4
loess_bin_width <- 0.01
bed2bamtypes <- data.frame(bamtypes =c ("wg","wg","wg","targeted","targeted"), 
                           bedtypes = c("wg_coarse","wg_fine","wg_medium","targeted","off_target"),
                           stringsAsFactors=FALSE)
rownames(bed2bamtypes) <- bed2bamtypes$bedtypes
  
ip <- (installed.packages())
Rsubread_version <- ip["Rsubread","Version"]
if(DEBUG) { gv$log <- paste0(isolate(gv$log),flog.info("Found RSubread version: %s",Rsubread_version)) }

# -------------------------------------------------------------------
makeCNTables <- function(configFileName) {
# -------------------------------------------------------------------
  # setwd("~/mm/cnv/buildcnb") # for debugging
  # print("makeCNTables()")

  # This is the useful bit that makes all the tables that the CN browser needs
  if(!exists("gv")) {
    flog.info("Running makeCNTables() stand-alone with config file: %s",configFileName)
  } # For stand-alone

  df_cfg <- ftry(read.dcf(configFileName))
  if(length(df_cfg)==1) {
    gv$log <- paste0(isolate(gv$log),flog.error("Unable to read config file: %s",configFileName))
    return(FALSE)
  }
  gv$log <- paste0(isolate(gv$log),flog.info("Successfully read config file: %s",configFileName))
  df_cfg <- as.data.frame(df_cfg,stringsAsFactors=FALSE)
  rownames(df_cfg) <- df_cfg$name

  no_clobber <- as.logical(df_cfg["no_clobber","value"])
  output_path <- df_cfg["output_path","value"]
  reference_path <- df_cfg["reference_path","value"]
  chrom_info_file <- df_cfg["chrom_info_file","value"]
  number_of_samples <- as.integer(df_cfg["number_of_samples","value"])

  df_cfg <- df_cfg[!(row.names(df_cfg) %in% paste0("sample_",(number_of_samples+1):1000)), ] # remove extras
  
  bedfiles <- grep("bedfile",row.names(df_cfg),value=TRUE) 
  idx <-  grepl("bams",df_cfg$section) | ((df_cfg$section=="files") & !grepl("path",row.names(df_cfg)))
  downloadable_files <- row.names(df_cfg)[idx]
  annotation_files  <- c("wg_annotation_file","targeted_annotation_file","cytobands_file","chrom_info_file")
  ret <- sapply(downloadable_files,optionally_download,df_cfg)
  # Pull out URLs
  idx <- (grepl("bams",df_cfg$section) | grepl("files",df_cfg$section) ) & isURL(df_cfg$value)
  df_cfg$value[idx] <- sub("^.*/","",df_cfg$value[idx])
  # Put in paths in bams
  idx <- grepl("bams",df_cfg$section) & !hasPath(df_cfg$value)
  df_cfg$value[idx] <- file.path(df_cfg["input_path","value"],df_cfg$value[idx])

  # Put in paths in refs - except for input_file
  idx <- grepl("files",df_cfg$section) & !grepl("path",row.names(df_cfg)) & !hasPath(df_cfg$value) & !grepl("input_file",row.names(df_cfg))
  df_cfg$value[idx] <- file.path(df_cfg["reference_path","value"],df_cfg$value[idx])

  optionally_make_bins(wg_coarse_bin_size,df_cfg["wg_coarse_bedfile","value"],df_cfg)
  optionally_make_bins(wg_medium_bin_size,df_cfg["wg_medium_bedfile","value"],df_cfg)
  optionally_make_bins(wg_fine_bin_size,df_cfg["wg_fine_bedfile","value"],df_cfg)

  ret <- sapply(bedfiles,optionally_make_gc,df_cfg)

  optionally_make_aliases(df_cfg)
  cnb_table_entry <-  do.call("rbind",lapply(bedfiles,optionally_make_counts,df_cfg))
#   gc_cor_file <- file.path(output_path <- df_cfg["output_path","value"],"gc_cor_all.tsv")
#   write.table(gc_cor_df,gc_cor_file,row.names=FALSE,quote = FALSE,sep="\t")
  optionally_make_infile(cnb_table_entry,annotation_files,df_cfg) 
  # make_all_gc_plots(df_cfg)
  # TODO: OPTIONALLY : Pre-make best reference, blacklists from inter and intra-sample variability
}

# -------------------------------------------------------------------
optionally_make_infile <- function(cnb_table_entry,annotation_files,df_cfg) {
  no_clobber <- as.logical(df_cfg["no_clobber","value"])
  # output_path <- df_cfg["output_path","value"]
  # reference_path <- df_cfg["reference_path","value"]
  # input_file <- file.path(output_path,input_file)
  input_file <- df_cfg["input_file","value"]
  if((!file.exists(input_file) || no_clobber==FALSE)) {
    input_df <- df_cfg[annotation_files,]
    input_df$type <-input_df$name
    input_df$label <-"ref"
    # cnb_table_entry$value <- file.path(df_cfg["output_path","value"],basename(cnb_table_entry$value))
    df <- rbind(input_df,cnb_table_entry)
    output_df <- data.frame(file = df$type, 
                            sample = df$label, 
                            path = df$value, 
                            is_targeted = grepl("target",df$section) ,
                            is_wgs = grepl("wg",df$section))
    
    # Needs to be coerced in order not to break write.table
    # df <- data.frame(lapply(df, as.character), stringsAsFactors=FALSE)
    gv$log <- paste0(isolate(gv$log),flog.info("Making: %s",input_file))
    ret <- ftry(write.table(output_df,input_file,row.names=FALSE,sep="\t",quote = FALSE))
  } else {
    gv$log <- paste0(isolate(gv$log),flog.info("Already have: %s",input_file))
  }
}


# -------------------------------------------------------------------
bedrow2gcrow <- function(s,df) { 
  row <- df[s,] 
  rownames(row) <- sub("bed","gc",rownames(row))
  row$name <- sub("bed","gc",row$name)
  row$label <- sub("bed","gc",row$label)  
  row$value <- sub(".bed","_gc.bed",row$value)
  row$value <- gsub("^.*/","",row$value) # strip path
  return(row)
}

# -------------------------------------------------------------------
isURL <- function(s) { return(grepl("^https://",s) | grepl("^http://",s) | grepl("^ftp://",s)) }
hasPath <- function(s) { return(grepl("/",s)) }
# -------------------------------------------------------------------
  
# -------------------------------------------------------------------
remove_files_by_suffix  <- function(suffix) {
  files <- list.files(pattern=paste0("*.",suffix)) 
  if(length(files)>0) {
    gv$log <- paste0(isolate(gv$log),flog.info("Removing temporary file: %s",files))  
    unlink(files)
  }
}

# -------------------------------------------------------------------
optionally_make_aliases <- function(df_cfg) {
  chrAliasesFile <- file.path(df_cfg["chr_alias","value"])
  no_clobber <- as.logical(df_cfg["no_clobber","value"])
  gv$log <- paste0(isolate(gv$log),flog.info("Making: %s",chrAliasesFile))
  if((!file.exists(chrAliasesFile) || no_clobber==FALSE)) {
    nochr <- c(paste0("",c(1:22)),"X","Y")
    chr<- c(paste0("chr",c(1:22)),"chrX","chrY")
    chrAliases <- data.frame(x=c(nochr,chr), X=c(chr,nochr))
    ret <- ftry(write.csv(chrAliases,chrAliasesFile,row.names=FALSE,quote = FALSE))
  }
  gv$log <- paste0(isolate(gv$log),flog.info("Have: %s",chrAliasesFile))
}

# -------------------------------------------------------------------
optionally_make_counts  <-  function(bedfile,df_cfg) {
  # If needed Make read counts with featureCounts
  bed_type <- sub("_bedfile","",df_cfg[bedfile,"name"])
  output_path <- df_cfg["output_path","value"]
  bam_type <-  bed2bamtypes[bed_type,"bamtypes"]
  df_bamfiles <- df_cfg[grepl("bams",df_cfg$section) & grepl(bam_type,df_cfg$section),]
  num_samples <- nrow(df_bamfiles)
  if(num_samples>0) {
    df_cnb <-  do.call("rbind",lapply(1:num_samples,write_gc_corrected_counts,df_bamfiles,bedfile,df_cfg))
    return(df_cnb)
  }
  return(NULL)
}

# -------------------------------------------------------------------
write_gc_corrected_counts <- function(idx,df_bamfiles,bedfile,df_cfg) {
  no_clobber <- as.logical(df_cfg["no_clobber","value"])
  output_path <- df_cfg["output_path","value"]
  number_of_cores <- as.integer(df_cfg["number_of_cores","value"])
  bamname <- sub(".bam$","",basename(df_bamfiles[idx,"value"]))
  bed_type <- sub("_bedfile","",df_cfg[bedfile,"name"])
  bed_suffix <- paste0("_",bed_type,".tsv")
  outfile <- sub(".bam$",bed_suffix,basename(df_bamfiles[idx,"value"]))
  outfile <- file.path(output_path,outfile)
  countsfile <- sub(".tsv$","_fc_counts.tsv",outfile)
  annotationfile <- sub(".tsv$","_fc_annotation.tsv",outfile)
  bedpath <- df_cfg[bedfile,"value"]
  chrAliasesFile <- file.path(df_cfg["chr_alias","value"])
  
  if(file.exists(countsfile) && no_clobber) {
    gv$log <- paste0(isolate(gv$log),flog.info("Skipping making counts for %s / %s",bedfile,bamname))
  } else {
    gv$log <- paste0(isolate(gv$log),flog.info("Making counts for: %s / %s",bedfile,bamname))
    gv$log <- paste0(isolate(gv$log),flog.info("This bit takes a long time - of order 1 min/10^6 reads/bedfile or 1 hour/5GB/bedfile."))
    df_bedfile <- read.table(bedpath,header=FALSE,stringsAsFactors=FALSE)
    colnames(df_bedfile) <- c('chr','start','end','id','score')
    df_id <- df_bedfile[,4:5]
    gr_bedfile <- import.bed(bedpath)
    mcols(gr_bedfile) <- df_id
    annotation_file <- createAnnotationFile(gr_bedfile)
    # system("vmstat -S M -s")
    # change into destination directory to avoid temp file collisions
    cwd <- getwd()
    setwd(output_path)
    if(Rsubread_version=="1.16.1") {
      fc <- featureCounts(files=df_bamfiles[idx,"value"],
                          isPairedEnd=TRUE,
                          annot.ext=annotation_file,
                          countMultiMappingReads=TRUE,
                          nthreads = number_of_cores,
                          useMetaFeatures = TRUE,
                          allowMultiOverlap = TRUE,
                          minReadOverlap = 60,
                          chrAliases = chrAliasesFile)
    }
    else {
      fc <- featureCounts(files=df_bamfiles[idx,"value"],
                          isPairedEnd=TRUE,
                          annot.ext=annotation_file,
                          countMultiMappingReads=TRUE,
                          nthreads = number_of_cores,
                          useMetaFeatures = TRUE,
                          allowMultiOverlap = TRUE,
                          minOverlap = 60,
                          fraction = TRUE,     
                          chrAliases = chrAliasesFile)
    }
    # system("vmstat -S M -s")
    setwd(cwd)
    save_list_elements <- function(s) { write.table(fc[[s]],sub(".tsv$",paste0("_fc_",s,".tsv"),outfile),quote = FALSE,sep="\t") }
    lapply(names(fc), save_list_elements) # save everything that featureCounts() returns
    # remove temporary files
    lapply(c("fixbam","sam","tmp"), remove_files_by_suffix)
  }
  
  if(file.exists(outfile) && no_clobber) {
    gv$log <- paste0(isolate(gv$log),flog.info("Skipping GC correction for for %s / %s",bedfile,bamname))
    cnb_table_entry <- list(file = bed_type, path = outfile)
    df_cnb <- df_bamfiles[idx,]
    df_cnb$type <- bed_type
    df_cnb$value <- outfile
    return(df_cnb)
  } else {
    fc <- list()
    fc$counts <- read.table(countsfile,header=TRUE,stringsAsFactors=FALSE)
    fc$annotation <- read.table(annotationfile,header=TRUE,stringsAsFactors=FALSE)
    gv$log <- paste0(isolate(gv$log),flog.info("Making GC correction for: %s",outfile))
    df_bedfile <- read.table(bedpath,header=FALSE,stringsAsFactors=FALSE)
    colnames(df_bedfile) <- c('chr','start','end','id','score')
    if(nrow(df_bedfile)>loess_scales_badly_limit) {
      gv$log <- paste0(isolate(gv$log),flog.info("Bedfile %s has %d entries so this could take some time",bedpath,nrow(df_bedfile)))
    }
    gc_file <- bedrow2gcrow(bedfile,df_cfg)$value
    reference_path <- df_cfg["reference_path","value"]
    gc_path <- file.path(reference_path,gc_file)
    df_gc <- read.table(gc_path,stringsAsFactors=FALSE,header=TRUE)
    gr_gc  <- GRanges(seqnames=df_gc$chr,ranges=IRanges(start=df_gc$start,end=df_gc$end),score=df_gc$gc)
    counts <- fc$counts
    colnames(counts) <- "counts"
    annotation <- fc$annotation
    df_counts <- cbind(counts,annotation)
    rownames(df_counts) <- paste0(df_counts$Chr,":",df_counts$Start,"-",df_counts$End)
    idx <- paste0(df_gc$chr,":",df_gc$start,"-",df_gc$end)
    if(nrow(df_gc) != nrow(counts)) {
      gv$log <- paste0(isolate(gv$log),flog.error("Number of raw counts does not match GC correction bins: %d != %d",nrow(df_gc),length(counts)))  
    }
    df_gc$counts <- df_counts[idx,"counts"]
    
    binSize <- df_gc$end - df_gc$start + 1
    binCorrection <- median(binSize) / binSize
    
    df_gc$counts  <- binCorrection * df_gc$counts

    idx <- df_gc$counts>loess_min_usable_counts & df_gc$gc>loess_min_usable_gc
    
    if(!is.na(df_cfg["gc_fitting_blacklist","value"])) {
      gv$log <- paste0(isolate(gv$log),flog.info("Reading GC blacklist"))    
      gc_fitting_blacklist  <- import.bed(df_cfg["gc_fitting_blacklist","value"])
      gv$log <- paste0(isolate(gv$log),flog.info("Finding overlaps"))    
      idx <- is.na(findOverlaps(gr_gc,gc_fitting_blacklist,select="first")) & idx 
    }
    gv$log <- paste0(isolate(gv$log),flog.info("Done"))    
    df_gc_fit <- df_gc[idx,]
    
    # loess_predict <- loess(counts ~ gc, data = df_gc_fit,  span = loess_span) # subset = !blacklist, span=TODO?
    loess_predict <- tryCatch(
      loess(counts ~ gc, data = df_gc_fit,  span = loess_span),
      error=function(cond) { 
        message("Loess failed:") 
        gv$log <- paste0(isolate(gv$log),flog.info("Loess failed: %s",cond))
        return(NULL)
      }
    )
    
    if(is.null(loess_predict)) {
      correction <- double(length=nrow(df_gc)) + 1
    } else {
      # Even predict fails
      correction <- tryCatch(
        predict(loess_predict,df_gc),
        error=function(cond) { 
          message("Loess predict failed:") 
          gv$log <- paste0(isolate(gv$log),flog.info("Loess predict failed: %s",cond))
          return(NULL)
        }
      )
      if(is.null(correction)) {
        loess_predict <- NULL
        correction <- double(length=nrow(df_gc)) + 1  
      }
    }
    df_correctable <- df_gc[!is.na(correction),]
    correction_not_na <- correction[!is.na(correction)]
    min_gc_correction  <- correction_not_na[which.min(df_correctable$gc)]
    max_gc_correction <- correction_not_na[which.max(df_correctable$gc)]
    correction[is.na(correction) & (df_gc$gc < min(df_correctable$gc))] <- min_gc_correction
    correction[is.na(correction) & (df_gc$gc > max(df_correctable$gc))] <- max_gc_correction
    correction <- max(correction) / correction
    gc_corrected <- correction * df_gc$counts
    # Combine bedfile with counts 
    df <- df_gc
    df$cor <- gc_corrected
    colnames(df)[colnames(df)=="counts"] <- "raw"
    df$cor[is.na(df$cor)] <- 0 # shouldn't be but is sometimes
    if(!is.na(df_cfg["display_blacklist","value"])) {
      gv$log <- paste0(isolate(gv$log),flog.info("Reading display blacklist"))    
      display_blacklist  <- import.bed(df_cfg["display_blacklist","value"])
      idx <- is.na(findOverlaps(gr_gc,display_blacklist,select="first"))
      df <- df[idx,]
    }
    gv$log <- paste0(isolate(gv$log),flog.info("Writing GC corrected counts: %s",outfile))
    ftry(write.table(df,file=outfile,row.names=FALSE,sep="\t",quote = FALSE))
    bins <- seq(0, 1, by = loess_bin_width)
    fapply_median <- function(n,bins) { median(df_gc_fit$counts[df_gc_fit$gc>=(bins[n] - loess_bin_width) & df_gc_fit$gc<(bins[n] + loess_bin_width)]) }
    bin_medians <- unlist(lapply(1:length(bins),fapply_median,bins))
    loessCorrection_df <- data.frame(bam =  bamname,bed = bedfile, gc = bins, median = bin_medians)
    if(is.null(loess_predict)) {
      loessCorrection_df$counts <- bin_medians
    } else {
      loessCorrection_df$counts <- predict(loess_predict,loessCorrection_df$gc)
    }
    loessfile <- sub(".tsv$",paste0("_fc_loess_corrections.tsv"),outfile)
    write.table(loessCorrection_df,loessfile,row.names=FALSE,quote = FALSE,sep="\t")
    df_cnb <- df_bamfiles[idx,]
    df_cnb$type <- bed_type
    df_cnb$value <- outfile
    return(df_cnb)
  }
}

# -------------------------------------------------------------------
optionally_download <- function(svalue,df_cfg) {
# -------------------------------------------------------------------
  no_clobber <- as.logical(df_cfg["no_clobber","value"])
  s <- df_cfg[svalue,"value"]
  dest_path <- df_cfg["reference_path","value"]
  if(grepl("bams",df_cfg[svalue,"section"])) { dest_path <- df_cfg["input_path","value"] }
  if(isURL(s)) {
    url <- parse_url(s)
    filename <- basename(url$path)
    destfile <- file.path(dest_path,filename)
    if((!file.exists(destfile) || no_clobber==FALSE)) {
      gv$log <- paste0(isolate(gv$log),flog.info("Downloading %s",destfile))          
      ret <- ftry(download.file(s,destfile))
      if(ret==0) {
        gv$log <- paste0(isolate(gv$log),flog.info("Successfully downloaded %s",destfile))          
      } else {
        gv$log <- paste0(isolate(gv$log),flog.error("Unable to download file: %s from %s",destfile,s))
      }
    }  
  }
}


# -------------------------------------------------------------------
optionally_make_bins <- function(bin_size,bedfile,df_cfg) {
  reference_path <- df_cfg["reference_path","value"]
  chrom_info_file <- df_cfg["chrom_info_file","value"]
  no_clobber <- as.logical(df_cfg["no_clobber","value"])
  destfile <- bedfile
  chromInfo <- read.table(chrom_info_file,stringsAsFactors=FALSE)
  colnames(chromInfo) <- c("chr","end","file")
  rownames(chromInfo) <- chromInfo$chr
  if((!file.exists(destfile) || no_clobber==FALSE)) {
    gv$log <- paste0(isolate(gv$log),flog.info("Making: %s",destfile))
    chr_names <- chromInfo$chr
    chr_names <- chr_names[!grepl("[_M]",chr_names)]
    make_intervals <- function(chr,chromInfo) { 
      options(scipen=20)
      gv$log <- paste0(isolate(gv$log),flog.info("Making bins for %s:%s",chr,destfile))    
      endList <- as.integer(seq(bin_size,chromInfo[chr,"end"],bin_size))
      startList <- as.integer(endList-bin_size)
      id <- paste0(chr,":",startList,"-",endList)
      score <- endList-startList
      strand <- "+"
      df <- data.frame(chr = chr, start = startList, end = endList,id = id, score = score, strand = strand)
      return((df))
    }
    df_intervals <-  do.call("rbind",lapply(chr_names,make_intervals,chromInfo))
    options(scipen=20)
    ftry(write.table(df_intervals,destfile,col.names=FALSE,row.names=FALSE,quote = FALSE,sep="\t"))
    gv$log <- paste0(isolate(gv$log),flog.info("Making bins for %s: Done",destfile))    
  }
  gv$log <- paste0(isolate(gv$log),flog.info("Have: %s",destfile))
}

# -------------------------------------------------------------------
optionally_make_gc <- function(svalue,df_cfg) {
  reference_path <- df_cfg["reference_path","value"]
  bedfile <- df_cfg[svalue,"value"]
  chrom_file <- df_cfg["chrom_file","value"]
  no_clobber <- as.logical(df_cfg["no_clobber","value"])
  s <- bedrow2gcrow(svalue,df_cfg)$value
  destfile <- file.path(reference_path,s)
  
  if((!file.exists(destfile) || no_clobber==FALSE)) {
    gv$log <- paste0(isolate(gv$log),flog.info("Making GC corrections. Reading: %s",bedfile))
    intervals <- import.bed(bedfile)
    gv$log <- paste0(isolate(gv$log),flog.info("Making GC corrections. Reading: %s",chrom_file))
    reference <- readDNAStringSet(chrom_file)
    
    gv$log <- paste0(isolate(gv$log),flog.info("Making: %s",destfile))
    # for each chr, something like this
    find_gc_content <- function(idx) { 
      seq <- reference[idx]
      # Assumes that naming fields in fasta file are chr2 with nothing else.
      gv$log <- paste0(isolate(gv$log),flog.info("Making GC content for %s:%s",names(seq),destfile))    
      chr_int <- intervals[seqnames(intervals)==names(seq)] 
      if(length(chr_int)==0) {
        gv$log <- paste0(
          isolate(gv$log),
          flog.warn("Seqname %s in %s not found in reference %s",names(seq),chrom_file,bedfile)
        )    
      }
      str <- Views(seq[[1]],start(chr_int),end(chr_int))
      gc <- letterFrequency(str,letters="CG",as.prob=TRUE)
      return(data.frame(chr=seqnames(chr_int),start=start(chr_int),end=end(chr_int),gc=gc))
    }
    df_gc <-  do.call("rbind",lapply(1:length(reference),find_gc_content))
    if(nrow(df_gc)==0) {
      gv$log <- paste0(
        isolate(gv$log),
        flog.error("Seqname %s and reference %s have no sequence names in common. Mismatched chr naming convention?",names(seq),bedfile)
      )
    }
    colnames(df_gc) <- c("chr","start","end","gc")
    write.table(df_gc,destfile,row.names=FALSE,quote = FALSE,sep="\t")
    gv$log <- paste0(isolate(gv$log),flog.info("Making GC content for %s: Done",destfile))  
    return(df_gc)
  }
  gv$log <- paste0(isolate(gv$log),flog.info("Have: %s",destfile))
}

make_all_gc_plots <- function(df_cfg) {
  output_path <- df_cfg["output_path","value"]
  suffixes <- c("_targeted.tsv","_off_target.tsv","_medium.tsv","_coarse.tsv","_fine.tsv")
  filelist <- unlist(lapply(suffixes,function(x) { dir(output_path,pattern=x) }))
  lapply(file.path(output_path,filelist),make_gc_plot)
}

make_gc_plot <- function(gc_file) {
  gv$log <- paste0(isolate(gv$log),flog.info("Plotting:%s",gc_file))  
  # Optionally not ocerwrite because this takes a while
  plot_file <- sub(".tsv",".png",gc_file)
  df_gc <- read.table(gc_file,header=TRUE,stringsAsFactors=FALSE)
  idx <- df_gc$cor>0 & df_gc$raw>0 & df_gc$gc>0.2
  df_gc <- df_gc[idx,]
  df_gc$cor <- log10(df_gc$cor)
  df_gc$raw <- log10(df_gc$raw)

  loess_predict <- loess(raw ~ gc, data = df_gc,  span = loess_span) # subset = !blacklist, span=TODO?
  correction <- predict(loess_predict,df_gc)
  df_correctable <- df_gc[!is.na(correction),]
  correction_not_na <- correction[!is.na(correction)]
  min_gc_correction  <- correction_not_na[which.min(df_correctable$gc)]
  max_gc_correction <- correction_not_na[which.max(df_correctable$gc)]
  correction[is.na(correction) & (df_gc$gc < min(df_correctable$gc))] <- min_gc_correction
  correction[is.na(correction) & (df_gc$gc > max(df_correctable$gc))] <- max_gc_correction
  correction <- max(correction) / correction
  gc_corrected <- correction * df_gc$counts
  # Combine bedfile with counts 
  bins <- seq(0, 1, by = loess_bin_width)
  fapply_median <- function(n,bins) { median(df_gc$raw[df_gc$gc>=(bins[n] - loess_bin_width) & df_gc$gc<(bins[n] + loess_bin_width)]) }
  bin_medians <- unlist(lapply(1:length(bins),fapply_median,bins))
  loessCorrection_df <- data.frame(gc = bins, median = bin_medians)
  loessCorrection_df$cor <- predict(loess_predict,loessCorrection_df$gc)
  loessCorrection_df <- loessCorrection_df[!(is.na(loessCorrection_df$cor)|is.na(loessCorrection_df$median)),]
  
  gc <- loessCorrection_df$gc
  bin_medians <- loessCorrection_df$median
  bin_cor <- loessCorrection_df$cor
  # TODO: Why do things look pre-smoothed?
  
  p <- ggplot(NULL)  + 
    #geom_point(aes(x=-1,y=0,size=1)) + 
    scale_size_identity(aes(size=0), guide=FALSE) + 
    geom_point(data = df_gc,aes(x = gc, y = raw,size=1,colour="blue")) + 
    geom_point(data = df_gc,aes(x = gc, y = cor,size=1,colour="red")) + 
    geom_line(data = loessCorrection_df, aes(x = gc, y = cor,size=1,colour="red")) + 
    geom_line(data = loessCorrection_df, aes(x = gc, y = median,size=1,colour="blue")) + 
    guides(colour = "none", size = "none") 
  p
  ggsave(filename=plot_file, plot=p)
}

if(RSTUDIO_DEBUG) { 
  makeCNTables(args) 
  return(0)
}

if(isCommandLine) {
  makeCNTables(args) 
  quit(save = "no", status = 0)
}

