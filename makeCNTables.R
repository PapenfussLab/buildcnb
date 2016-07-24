#!/usr/bin/Rscript --vanilla
args <- commandArgs(TRUE)
# For testing
# setwd("~/mm/cnv/buildcnb")
# source('~/mm/cnv/buildcnb/makeCNTables.R', echo=TRUE)
# RSTUDIO_DEBUG <- TRUE
RSTUDIO_DEBUG <- FALSE
if(RSTUDIO_DEBUG){ args <- "/home/jmarkham/mm/cnv/buildcnb/buildcnb_input_160426_NS500817_0060_AH3FYWBGXY.txt"}
if(RSTUDIO_DEBUG){ args <- "buildcnb_test_160330_NS500817_0055_AH27FYBGXY.dcf" }

isCommandLine <- !RSTUDIO_DEBUG && length(args)>0 && args[1]!="RStudio" &&  !grepl("R$",args) 

if(length(args)!=1 && isCommandLine) {
  print("Builds files for cnb. Usage: makeCNTables.R /path/to/config_file.dcf")
  quit(save = "no", status = 0)
} 


libnames <- c(
  "plyr",
  "futile.logger",
  "RCurl",
  "httr",
  "Rsubread",
  "Biostrings",
  "rtracklayer",
  "ggplot2",
  "VariantAnnotation",
  "matrixStats",
  "DNAcopy",
  "Hmisc",
  "modeest"
  )
for (libname in libnames) { 
  library(libname,character.only=TRUE,warn.conflicts=FALSE,quietly=TRUE) 
}
if(!exists("gv")) {
  session <- NULL
  gv <- list()
  gv$log <- ""
  isolate <- function (f) { return (f)}
} # For stand-alone

# Some globals
WG_COARSE_BIN_SIZE <- 1000000
WG_MEDIUM_BIN_SIZE <- 50000
WG_SMALL_BIN_SIZE <- 5000
CBS_ALPHA = 0.01
CBS_NPERM = 10000
LOESS_SCALES_BADLY_LIMIT  <- 100000 # when loess gets slow
LOESS_MIN_USABLE_COUNTS  <- 50 # For fitting
LOESS_MIN_USABLE_GC  <- 0.25 
LOESS_SPAN <- 0.4
LOESS_BIN_WIDTH <- 0.01
MINIMUM_OVERLAP <- 120
NORM_CUTOFF <- 10   # Counts in a bin less that this are thrown away fornormalisation
MEDIAN_NORM_CUTOFF <- 10   # Counts in a bin with a dist whose median less that this are thrown away for  normalisation
MODE_NORM_CUTOFF <- 100   # Counts in a bin with a dist whose mode less that this are thrown away for  normalisation
COUNTS_CUTOFF <- 10   # Counts less that this are thrown away for CN estimation
SNP_MIN_MISMATCHES <- 1    # Before a SNP is called 
SNP_MIN_READS <- 20        # Before a SNP is called - these two also give the min B-allele freq
SNP_MAX_READS <- 1000   # For PCR artifacts
SEGMENT_VERBOSITY  <- 0
FUSION_SECTION <- "targeted_bams"
BAF_SECTION <- "targeted_bams"
VCF_SECTION <- "targeted_bams"
DEL_SECTION <- "wg_bams"

CFG <- list()
DF_CFG <- data.frame()

DF_BEDTYPE2BAMTYPE <- data.frame(bamtypes =c ("wg","wg","wg","targeted","targeted"), 
                           bedtypes = c("wg_coarse","wg_fine","wg_medium","targeted","off_target"),
                           bedfilename = c("wg_coarse_bedfile","wg_fine_bedfile","wg_medium_bedfile","targeted_bedfile","off_target_bedfile"),
                           stringsAsFactors=FALSE)
ANNOTATION_FILES  <- c("wg_annotation_file","targeted_annotation_file","cytobands_file","chrom_info_file") # TODO: From input file section, WG or not
BEDTYPE2BAMTYPE <- DF_BEDTYPE2BAMTYPE$bamtypes
names(BEDTYPE2BAMTYPE) <- DF_BEDTYPE2BAMTYPE$bedtypes

BEDFILENAME2BEDTYPE <- DF_BEDTYPE2BAMTYPE$bedtypes
names(BEDFILENAME2BEDTYPE) <- DF_BEDTYPE2BAMTYPE$bedfilename

ip <- (installed.packages())
Rsubread_version <- ip["Rsubread","Version"]
if(RSTUDIO_DEBUG) { gv$log <- paste0(isolate(gv$log),flog.info("Found RSubread version: %s",Rsubread_version)) }


# -------------------------------------------------------------------
makeCNTables <- function(configFileName) {
  # setwd("~/mm/cnv/buildcnb") # for debugging
  # print("makeCNTables()")

  # This is the useful bit that makes all the tables that the CN browser needs
  if(!exists("gv")) {
    flog.info("Running makeCNTables() stand-alone with config file: %s",configFileName)
  } # For stand-alone

  df_cfg <- ftry(read.dcf(configFileName))
  if(length(df_cfg)==1) {
    gv$log <<- paste0(isolate(gv$log),flog.error("Unable to read config file: %s",configFileName))
    return(FALSE)
  }
  gv$log <<- paste0(isolate(gv$log),flog.info("Successfully read config file: %s",configFileName))
  df_cfg <- as.data.frame(df_cfg,stringsAsFactors=FALSE)
  rownames(df_cfg) <- df_cfg$name
  number_of_samples <- as.integer(df_cfg["number_of_samples","value"])
  df_cfg <- df_cfg[!(row.names(df_cfg) %in% paste0("sample_",(number_of_samples+1):1000)), ] # remove extras
  idx <-  grepl("bams",df_cfg$section) | ((df_cfg$section=="files") & !grepl("path",row.names(df_cfg)))
  downloadable_files <- row.names(df_cfg)[idx]
  ret <- sapply(downloadable_files,download_files,df_cfg)
  df_cfg  <- massage_paths(df_cfg) # swap out URLs for local paths, adjust other locations as necessary
  #  CONFIG and DF_CONFIG are const globals
  DF_CFG <<- df_cfg
  CFG <<-  as.list(df_cfg$value)
  names(CFG) <<- df_cfg$name
  # bedfiles <- grep("bedfile",row.names(df_cfg),value=TRUE) 
  bedfiles <- grep("wg_coarse_bedfile|wg_medium_bedfile|wg_fine_bedfile|targeted_bedfile",row.names(df_cfg),value=TRUE) 
  make_bins(WG_COARSE_BIN_SIZE,df_cfg["wg_coarse_bedfile","value"])
  make_bins(WG_MEDIUM_BIN_SIZE,df_cfg["wg_medium_bedfile","value"])
  make_bins(WG_SMALL_BIN_SIZE,df_cfg["wg_fine_bedfile","value"])
  ret <- sapply(bedfiles,make_gc)
  make_aliases()
  make_bams()
  make_vcfs()
  df_bafs_tables <- make_bafs()
  df_pre_built_references <- do.call("rbind.fill",lapply(bedfiles,makeAllReferences))
  df_cnb_table_entry <-  do.call("rbind.fill",lapply(bedfiles,make_counts))
  df_deletions_table <- do.call("rbind.fill",lapply(bedfiles[!grepl("fine",bedfiles) & !grepl("coarse",bedfiles) & !grepl("target",bedfiles)],make_deletions)) # TODO: Clean up
  df_translocations <- make_translocations()
#   gc_cor_file <- file.path(output_path <- df_cfg["output_path","value"],"gc_cor_all.tsv")
#   write.table(gc_cor_df,gc_cor_file,row.names=FALSE,quote = FALSE,sep="\t")
  make_infile(rbind.fill(df_cnb_table_entry,df_deletions_table,df_bafs_tables,df_translocations,df_pre_built_references)) 
  # make_all_gc_plots(df_cfg)
  output_log_file <- gsub(".t.[tv]$",".log",CFG$input_file)
  write.table(gv$log,output_log_file,row.names = FALSE,col.names=FALSE,quote=FALSE)
}
# -------------------------------------------------------------------

# -------------------------------------------------------------------
make_infile <- function(df_cnb_table_entry) {
  no_clobber <- CFG$no_clobber
  input_file <- CFG$input_file
  if((!file.exists(input_file) || no_clobber==FALSE)) {
    input_df <- as.data.frame(DF_CFG[ANNOTATION_FILES,])
    input_df$file <-"reference"
    input_df$sample <- input_df$name
    input_df$resolution <- ifelse(grepl("targeted",input_df$name),"targeted","wg")
    input_df$is_targeted <- grepl("targeted",input_df$name)
    input_df$is_wgs <- !grepl("targeted",input_df$name)
    colnames(input_df)[colnames(input_df)=="value"] <- "path"
    output_df <- rbind.fill(input_df,df_cnb_table_entry)
    output_df$section <- NULL
    output_df$value <- NULL
    output_df$label <- NULL
    gv$log <<- paste0(isolate(gv$log),flog.info("Making: %s",input_file))
    ret <- ftry(write.table(output_df,input_file,row.names=FALSE,sep="\t",quote = FALSE))
  } else {
    gv$log <<- paste0(isolate(gv$log),flog.info("Already have: %s",input_file))
  }
}
# -------------------------------------------------------------------

# -------------------------------------------------------------------
bedrow2gcrow <- function(s) { 
  row <- DF_CFG[s,] 
  rownames(row) <- sub("bed","gc",rownames(row))
  row$name <- sub("bed","gc",row$name)
  row$label <- sub("bed","gc",row$label)  
  row$value <- sub(".bed","_gc.bed",row$value)
  # row$value <- gsub("^.*/","",row$value) # strip path
  return(row)
}
# -------------------------------------------------------------------

# -------------------------------------------------------------------
isURL <- function(s) { return(grepl("^https://",s) | grepl("^http://",s) | grepl("^ftp://",s)) }
# -------------------------------------------------------------------

# -------------------------------------------------------------------
hasPath <- function(s) { return(grepl("/",s)) }
# -------------------------------------------------------------------

# -------------------------------------------------------------------
massage_paths  <- function (df_cfg) {
  # Pull out URLs
  idx <- (grepl("bams",df_cfg$section) | grepl("files",df_cfg$section) ) & isURL(df_cfg$value)
  df_cfg$value[idx] <- sub("^.*/","",df_cfg$value[idx])
  # Put paths in bams
  idx <- grepl("bams",df_cfg$section) & !hasPath(df_cfg$value)
  samples <- sub(".bam$","",df_cfg$value[idx]) # maybe use label column
  bam_root <- makepath(panel = df_cfg["panel_path","value"] ,run = df_cfg["run_path","value"] ,sample = samples ,filetype = "bam",df_cfg = df_cfg)
  df_cfg$value[idx] <- paste0(bam_root,".bam")
  
  # Put in paths in refs - except for input_file
  idx_add_path <- grepl("files",df_cfg$section) & 
    !grepl("path",row.names(df_cfg)) & 
    !hasPath(df_cfg$value) & 
    !grepl("input_file",row.names(df_cfg))
  # A bit yuck. Anything without "list" and "target" in the name must be for ALL panels
  idx_all_ref <- !grepl("list",row.names(df_cfg)) & !grepl("target",row.names(df_cfg))
  idx <- idx_add_path & idx_all_ref
  reference_path <- makepath(panel = NULL ,run = NULL ,sample = NULL ,filetype = "reference",df_cfg = df_cfg) 
  df_cfg$value[idx] <- file.path(reference_path,df_cfg$value[idx])
  # Anything with "list" in the name must be for this particular panel panels
  idx_panel_ref <- grepl("list",row.names(df_cfg)) | grepl("target",row.names(df_cfg))
  idx <- idx_add_path & idx_panel_ref
  reference_path <- makepath(panel = df_cfg["panel_path","value"] ,run = NULL ,sample = NULL ,filetype = "reference",df_cfg = df_cfg) 
  df_cfg$value[idx] <- file.path(reference_path,df_cfg$value[idx])
  return (df_cfg)  
}
# -------------------------------------------------------------------

# -------------------------------------------------------------------
makepath <- function(panel = NULL ,run = NULL ,sample = NULL ,filetype = NULL,df_cfg = DF_CFG) { 
  base_path <- df_cfg["base_path","value"]
  filetype_path <- df_cfg[paste0(filetype,"_path"),"value"]
  if(sum(grepl(filetype,c("cnv","vcf","sv","baf","fastq","bam")))>0) {
    ret <- file.path(base_path,panel,run,sample,filetype_path,sample)
  } else if (filetype == "reference") {
    if(is.null(panel)) { ret <- file.path(base_path,filetype_path) }
    else if(is.null(run)) { ret <- file.path(base_path,panel,filetype_path) }
    else if(is.null(sample)) { ret <- file.path(base_path,panel,run,filetype_path) }
    # else if(is.null(sample)) { ret <- file.path(base_path,panel,run,df_cfg$value["run_stats_path"],filetype_path) }
    else { ret <- file.path(base_path,panel,run,filetype_path,sample) }
  }
  return(ret)
}
# -------------------------------------------------------------------


# -------------------------------------------------------------------
download_files <- function(svalue,df_cfg) {
  no_clobber <- as.logical(df_cfg$value["no_clobber"])
  s <- df_cfg[svalue,"value"]
  dest_path <- CFG$reference_path
  if(grepl("bams",df_cfg[svalue,"section"])) { dest_path <- df_cfg$value["input_path"] }
  if(isURL(s)) {
    url <- parse_url(s)
    filename <- basename(url$path)
    destfile <- file.path(dest_path,filename)
    if((!file.exists(destfile) || no_clobber==FALSE)) {
      gv$log <<- paste0(isolate(gv$log),flog.info("Downloading %s",destfile))          
      ret <- ftry(download.file(s,destfile))
      if(ret==0) {
        gv$log <<- paste0(isolate(gv$log),flog.info("Successfully downloaded %s",destfile))          
      } else {
        gv$log <<- paste0(isolate(gv$log),flog.error("Unable to download file: %s from %s",destfile,s))
      }
    }  
  }
}
# -------------------------------------------------------------------

# -------------------------------------------------------------------
remove_files_by_suffix  <- function(suffix) {
  files <- list.files(pattern=paste0("*.",suffix)) 
  if(length(files)>0) {
    gv$log <<- paste0(isolate(gv$log),flog.info("Removing temporary file: %s",files))  
    unlink(files)
  }
}
# -------------------------------------------------------------------


# -------------------------------------------------------------------
make_aliases <- function() {
  chrAliasesFile <- file.path(CFG$chr_alias)
  no_clobber <- as.logical(CFG$no_clobber)
  gv$log <<- paste0(isolate(gv$log),flog.info("Making: %s",chrAliasesFile))
  if((!file.exists(chrAliasesFile) || no_clobber==FALSE)) {
    nochr <- c(paste0("",c(1:22)),"X","Y")
    chr<- c(paste0("chr",c(1:22)),"chrX","chrY")
    chrAliases <- data.frame(x=c(nochr,chr), X=c(chr,nochr))
    ret <- ftry(write.csv(chrAliases,chrAliasesFile,row.names=FALSE,quote = FALSE))
  }
  gv$log <<- paste0(isolate(gv$log),flog.info("Have: %s",chrAliasesFile))
}
# -------------------------------------------------------------------

# -------------------------------------------------------------------
make_counts  <-  function(bedfile) {
  # If needed Make read counts with featureCounts
  bed_type <- BEDFILENAME2BEDTYPE[bedfile]
  bamtype <-  BEDTYPE2BAMTYPE[bed_type]
  refs <- readAllReferences(bed_type)
  if(length(refs)==0) {
    gv$log <<- paste0(isolate(gv$log),flog.error("Unable to find references.")) 
    return (NULL)
  }
  df_bamfiles <- DF_CFG[grepl("bams",DF_CFG$section) & grepl(bamtype,DF_CFG$section),]
  num_samples <- nrow(df_bamfiles)
  if(num_samples>0) {
    df_cnb <-  do.call("rbind.fill",lapply(1:num_samples,write_corrected_counts,df_bamfiles,bedfile,refs))
    return(df_cnb)
  }
  return(NULL)
}
# -------------------------------------------------------------------

# -------------------------------------------------------------------
write_corrected_counts <- function(idx,df_bamfiles,bedfile,references) {
  no_clobber <- as.logical(CFG$no_clobber)
  number_of_cores <- as.integer(CFG$number_of_cores)
  bamname <- sub(".bam$","",basename(df_bamfiles[idx,"value"]))
  bam_dir <- dirname(df_bamfiles[idx,"value"])
  bed_type <- BEDFILENAME2BEDTYPE[bedfile]
  bam_type <-  BEDTYPE2BAMTYPE[bed_type]
  correct_gc <- as.logical(DF_CFG[paste0(bam_type,"_gc_correct"),"value"])
  cnv_suffix <- paste0("_",bed_type,".tsv")
  counts_root <- makepath(panel = CFG$panel,run = CFG$run_path,sample = bamname,filetype = "cnv")
  outpath <- paste0(counts_root,cnv_suffix)
  cnb_info_path <- sub(".tsv$","_info.tsv",outpath)
  countspath <- sub(".tsv$","_fc_counts.tsv",outpath)
  annotationpath <- sub(".tsv$","_fc_annotation.tsv",outpath)
  bedpath <- CFG[[bedfile]]
  chrAliasesFile <- CFG$chr_alias
  
  if(file.exists(countspath) && no_clobber) {
    gv$log <<- paste0(isolate(gv$log),flog.info("Skipping making counts for: %s / %s. Already have %s",bedfile,bamname,countspath))
  } else {
    gv$log <<- paste0(isolate(gv$log),flog.info("Making counts for: %s / %s. In %s",bedfile,bamname,countspath))
    gv$log <<- paste0(isolate(gv$log),flog.info("This bit takes a long time - of order 1 min/10^6 reads/bedfile or 1 hour/5GB/bedfile."))
    df_bedfile <- read.table(bedpath,header=FALSE,stringsAsFactors=FALSE)
    colnames(df_bedfile) <- c('chr','start','end','id','score')
    df_id <- df_bedfile[,4:5]
    gr_bedfile <- import.bed(bedpath)
    mcols(gr_bedfile) <- df_id
    annotation_file <- createAnnotationFile(gr_bedfile)
    # system("vmstat -S M -s")
    # change into destination directory to avoid temp file collisions
    cwd <- getwd()
    chrAliasesFile <- file.path(cwd,chrAliasesFile)
    setwd(bam_dir)
    if(Rsubread_version=="1.16.1") {
      fc <- featureCounts(files=basename(df_bamfiles[idx,"value"]),
                          isPairedEnd=TRUE,
                          annot.ext=annotation_file,
                          countMultiMappingReads=TRUE,
                          nthreads = number_of_cores,
                          useMetaFeatures = TRUE,
                          allowMultiOverlap = TRUE,
                          minReadOverlap = MINIMUM_OVERLAP,
                          chrAliases = chrAliasesFile)
    }
    else {
      fc <- featureCounts(files=basename(df_bamfiles[idx,"value"]),
                          isPairedEnd=TRUE,
                          annot.ext=annotation_file,
                          countMultiMappingReads=TRUE,
                          nthreads = number_of_cores,
                          useMetaFeatures = TRUE,
                          allowMultiOverlap = TRUE,
                          largestOverlap = TRUE,
                          minOverlap = MINIMUM_OVERLAP,
#                           reportReads = TRUE # Probably not. Writes out too much
                          fraction = TRUE,     
                          chrAliases = chrAliasesFile)
    }
    # system("vmstat -S M -s")
    gv$log <<- paste0(isolate(gv$log),flog.info("Removing temporary files from %s",bam_dir))
    lapply(c("fixbam","sam","tmp","bin"), remove_files_by_suffix)
    setwd(cwd)
    save_list_elements <- function(s) { write.table(fc[[s]],sub(".tsv$",paste0("_fc_",s,".tsv"),outpath),quote = FALSE,sep="\t") }
    lapply(names(fc), save_list_elements) # save everything that featureCounts() returns
    # remove temporary files
    lapply(c("fixbam","sam","tmp","bin"), remove_files_by_suffix)
  }

  # Write entry for input file for cnb
  df_cnb <- data.frame(
    resolution = bed_type,
    path = outpath,
    file = "counts",
    name = df_bamfiles[idx,"name"],
    is_wgs = (bam_type=="wg"),
    is_targeted = (bam_type=="targeted"),
    sample = bamname,
    sex = "unknown",
    tumour_normal = "unknown",
    stringsAsFactors = FALSE)

  df_reference_file_list <- read.table(CFG$reference_file_list,header = TRUE,stringsAsFactors = FALSE)
  if(sum(grepl("Sex",colnames(df_reference_file_list)))>0) {
    df_cnb$sex <- df_reference_file_list$Sex[df_reference_file_list$SampleName==df_cnb$label][1]  
  }
  if(sum(grepl("Tumour",colnames(df_reference_file_list)))>0) {
    df_cnb$tumour_normal <- df_reference_file_list$TumourNormal[df_reference_file_list$SampleName==df_cnb$label][1]
  }

  if(file.exists(outpath) && no_clobber) {
    gv$log <<- paste0(isolate(gv$log),flog.info("Skipping GC correction, segmentation, reference correction for %s / %s",bedfile,bamname))
    if(file.exists(cnb_info_path)) {
      df_cnb_saved <- read.table(cnb_info_path,header=TRUE,stringsAsFactors=FALSE)
      df_cnb_saved[colnames(df_cnb)] <- df_cnb
      df_cnb <- df_cnb_saved
    } else {
      gv$log <<- paste0(isolate(gv$log),flog.warn("Unable to find %s but pressing on regardless.",cnb_info_path))
    }
    return(df_cnb)
  } else {
    fc <- list()
    fc$counts <- read.table(countspath,header=TRUE,stringsAsFactors=FALSE)
    fc$annotation <- read.table(annotationpath,header=TRUE,stringsAsFactors=FALSE)
    gv$log <<- paste0(isolate(gv$log),flog.info("Making GC correction for: %s",outpath))
    df_bedfile <- read.table(bedpath,header=FALSE,stringsAsFactors=FALSE)
    colnames(df_bedfile) <- c('chr','start','end','id','score')
    if(nrow(df_bedfile)>LOESS_SCALES_BADLY_LIMIT) {
      gv$log <<- paste0(isolate(gv$log),flog.info("Bedfile %s has %d entries so this could take some time",bedpath,nrow(df_bedfile)))
    }
    gc_path <- bedrow2gcrow(bedfile)$value
    # reference_path <- makepath(filetype = "reference") 
    df_gc <- read.table(gc_path,stringsAsFactors=FALSE,header=TRUE)
    gr_gc  <- GRanges(seqnames=df_gc$chr,ranges=IRanges(start=df_gc$start,end=df_gc$end),score=df_gc$gc)
    counts <- fc$counts
    colnames(counts) <- "counts"
    annotation <- fc$annotation
    df_counts <- cbind(counts,annotation)
    rownames(df_counts) <- paste0(df_counts$Chr,":",df_counts$Start,"-",df_counts$End)
    idx_chr_start_end <- paste0(df_gc$chr,":",df_gc$start,"-",df_gc$end)
    if(nrow(df_gc) != nrow(counts)) {
      gv$log <<- paste0(isolate(gv$log),flog.error("Number of raw counts does not match GC correction bins: %d != %d",nrow(df_gc),length(counts)))  
    }
    df_gc$counts <- df_counts[idx_chr_start_end,"counts"]
    # find usable bins
    binSize <- df_gc$end - df_gc$start + 1
    binCorrection <- median(binSize) / binSize
    df_gc$counts  <- binCorrection * df_gc$counts
    idx_gc_fit <- df_gc$counts>LOESS_MIN_USABLE_COUNTS & df_gc$gc>LOESS_MIN_USABLE_GC
    # Mask points in black list
    if(!is.na(CFG$gc_fitting_blacklist)) {
      gv$log <<- paste0(isolate(gv$log),flog.info("Reading GC blacklist"))    
      gv$log <<- paste0(isolate(gv$log),flog.info(CFG$gc_fitting_blacklist))    
      gc_fitting_blacklist  <- import.bed(CFG$gc_fitting_blacklist)
      gv$log <<- paste0(isolate(gv$log),flog.info("Finding overlaps"))    
      idx_gc_fit <- is.na(findOverlaps(gr_gc,gc_fitting_blacklist,select="first")) & idx_gc_fit 
    }
    gv$log <<- paste0(isolate(gv$log),flog.info("Done"))    
    df_gc_fit <- df_gc[idx_gc_fit,]
    
    loess_predict <- tryCatch(
      loess(counts ~ gc, data = df_gc_fit,  span = LOESS_SPAN),
      error=function(cond) { 
        message("Loess failed:") 
        gv$log <<- paste0(isolate(gv$log),flog.info("Loess failed: %s",cond))
        return(NULL)
      }
    )
    
    if(correct_gc==FALSE || is.null(loess_predict)) {
      correction <- double(length=nrow(df_gc)) + 1
    } else {
      # Even predict fails
      correction <- tryCatch(
        predict(loess_predict,df_gc),
        error=function(cond) { 
          message("Loess predict failed:") 
          gv$log <<- paste0(isolate(gv$log),flog.info("Loess predict failed: %s",cond))
          return(NULL)
        }
      )
      if(is.null(correction)) {
        loess_predict <- NULL
        correction <- double(length=nrow(df_gc)) + 1  
      }
    }
    # For the GC range outside where Loess worked, set it to the boundary values
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
    if(!is.na(CFG$display_blacklist)) {
      gv$log <<- paste0(isolate(gv$log),flog.info("Reading display blacklist"))    
      gv$log <<- paste0(isolate(gv$log),flog.info(CFG$display_blacklist))    
      display_blacklist  <- import.bed(CFG$display_blacklist)
      idx_not_black_listed <- is.na(findOverlaps(gr_gc,display_blacklist,select="first"))
      df$blackListed <- !idx_not_black_listed
    }
    
    gv$log <<- paste0(isolate(gv$log),flog.info("::Finding best reference for %s",outpath))
    bestReference <- findBestReference(references,df)
    gv$log <<- paste0(isolate(gv$log),flog.info("Writing GC corrected counts: %s",outpath))
    ftry(write.table(bestReference$corrected,file=outpath,row.names=FALSE,sep="\t",quote = FALSE))
    
    # Write out the loess correction across the range of GC values
    bins <- seq(0, 1, by = LOESS_BIN_WIDTH)
    fapply_median <- function(n,bins) { median(df_gc_fit$counts[df_gc_fit$gc>=(bins[n] - LOESS_BIN_WIDTH) & df_gc_fit$gc<(bins[n] + LOESS_BIN_WIDTH)]) }
    bin_medians <- unlist(lapply(1:length(bins),fapply_median,bins))
    loessCorrection_df <- data.frame(bam =  bamname,bed = bedfile, gc = bins, median = bin_medians)
    if(is.null(loess_predict)) {
      loessCorrection_df$counts <- bin_medians
    } else {
      loessCorrection_df$counts <- predict(loess_predict,loessCorrection_df$gc)
    }
    loesspath <- sub(".tsv$",paste0("_fc_loess_corrections.tsv"),outpath)
    write.table(loessCorrection_df,loesspath,row.names=FALSE,quote = FALSE,sep="\t")
    
    # Finish entry for input file for cnb
    df_reference_file_list <- read.table(CFG$reference_file_list,header = TRUE,stringsAsFactors = FALSE)
    if(sum(grepl("Sex",colnames(df_reference_file_list)))==0) {
      if(grepl("Female",bestReference$refname)) {
        df_cnb$sex <- "Female"
      } else if (grepl("Male",bestReference$refname)) {
        df_cnb$sex <- "Male"
      } 
    }
    df_cnb$refpoolsize <- median(bestReference$corrected$ref_N)
    df_cnb$refname <- bestReference$refname
    df_cnb$smoothness_diff <- bestReference$comparisons$smoothness_diff
    df_cnb$quantised_cn_diff <- bestReference$comparisons$quantised_cn_diff
    df_cnb$segments <- bestReference$comparisons$segments
    df_cnb$binSize <- median(df_gc$end - df_gc$start + 1)
    write.table(df_cnb,file = cnb_info_path,row.names=FALSE,sep="\t",quote = FALSE)
    return(df_cnb)
  }
}
# -------------------------------------------------------------------


# -------------------------------------------------------------------
# reference_file_list is the Master file for all samples. Contains:
# RunName  SampleName PanelName Date TumourNormal Sex SampleType CountsPath BafPath BamPath MedianFragmentLength TotalReads	PercentReadsMapped	PercentMappedReadsDuplicates
# Reference selected by vaues in these columns:
#   section: params
#   name: reference_variables
#   label: Reference Variables
#   value: TumourNormal=T|N;Sex=M|F;SampleType=Fresh|FFPE;MedianFragmentLength=short:120-150|medium:150-190|long:190-300;PercentReadsMapped=98pc:98-100;Date=Post2014:20140101-20200101
#   section: params
#   name: min_reference_sample
#   label: Minimum Reference Samples
#   value: 10
# Loop over all combinations of criteria and when there are enough samples, build a reference
# -------------------------------------------------------------------
makeAllReferences = function (bedfile) {
  no_clobber <- as.logical(CFG$no_clobber)
  bed_type <- BEDFILENAME2BEDTYPE[bedfile]
  bam_type <- BEDTYPE2BAMTYPE[bedfile]
  reference_file_list_file_name <- CFG$reference_file_list
  gc_path <- bedrow2gcrow(bedfile)$value
  no_clobber <- as.logical(CFG$no_clobber)
  output_path <- makepath(panel = CFG$panel_path,run = CFG$run_path,filetype = "reference")
  df_gc <- read.table(gc_path,stringsAsFactors=FALSE,header=TRUE)
  
  df_reference_file_list <- read.table(reference_file_list_file_name,header = TRUE,stringsAsFactors = FALSE)
  counts_ending <- paste0("_",bed_type,"_fc_counts.tsv")
  df_reference_file_list$CountsFile <- paste0(df_reference_file_list$SampleName,counts_ending)
  df_reference_file_list$CountsPath <- paste0 (
    makepath(panel = CFG$panel_path,
             run = df_reference_file_list$RunName,
             sample = df_reference_file_list$SampleName,
             filetype = "cnv"),
    counts_ending
  )
    
  num_samples <- nrow(df_reference_file_list)
  min_reference_samples <- CFG$min_reference_samples
  gv$log <<- paste0(isolate(gv$log),flog.info("Successfully read reference_file_list file: %s with %d samples",
                                             reference_file_list_file_name,num_samples))

  criteria <-  DF_CFG[grepl("params",DF_CFG$section) & grepl("reference_variables",DF_CFG$name),"value"]
  # Parse the reference_variables field
  cs <- strsplit(strsplit(criteria,";")[[1]],"=")
  df_criteria_col <- data.frame(SampleColumn = character,
                                Val = character,
                                isRange = logical,
                                Min = double,
                                Max = double,
                                stringsAsFactors = FALSE)
  ret <- ""
  for (i in 1:length(cs)) {
    colName <- cs[[i]][1]
    colCriteria <- strsplit(cs[[i]][2],"[|]")[[1]]
    for (j in 1:length(colCriteria) ) {
      colVal <- colCriteria[j]
      colValList <- strsplit(colVal,":")[[1]]
      if(length(colValList)==1) {
        df_criteria_col_entry <- data.frame(
          SampleColumn=colName,Val=colVal,isRange=FALSE,Min=NA,Max=NA,stringsAsFactors=FALSE)
      } else {
        colVal <- colValList[1]
        colMinMax <- strsplit(colValList[2],"[-]")[[1]]
        df_criteria_col_entry <- data.frame(
          SampleColumn=colName,Val=colVal,isRange=TRUE,Min=as.double(colMinMax[1]),Max=as.double(colMinMax[2]),stringsAsFactors=FALSE)
      }
      df_criteria_col <- rbind(df_criteria_col,df_criteria_col_entry)
    }
    # When there is more than one category, also do one with All in it.
    if(length(colCriteria)>1) {
      df_criteria_col_entry <- data.frame(
        SampleColumn=colName,Val="All",isRange=FALSE,Min=NA,Max=NA,stringsAsFactors=FALSE)
      df_criteria_col <- rbind(df_criteria_col,df_criteria_col_entry)
    }
  }
  # Make all the combinations of the chosen variables
  uniqueCols <- unique(df_criteria_col$SampleColumn)
  numCategories <- length(uniqueCols)
  df_combinations <- list()
  for (i in 1:numCategories) {
    df_combinations[[uniqueCols[i]]] <- which(df_criteria_col$SampleColumn==uniqueCols[i])
  }
  df_all_combinations <- expand.grid(df_combinations)
  numCombinations <- nrow(df_all_combinations)
  gv$log <<- paste0(isolate(gv$log),flog.info("Need %d references",numCombinations))
  df_files <- df_all_combinations
  df_files$file <- "counts_reference"
  df_files$resolution <- bed_type
  df_files$is_wgs <- (bam_type=="wg")
  df_files$is_targeted <- !df_files$is_wgs
  
  for (i in 1:numCombinations) {
    reference_file_name <- ""
    criteria_rows <- df_all_combinations[i,]
    idx <- !logical(num_samples)
    for (j in criteria_rows) {
      df_row <- df_criteria_col[j,]
      if(df_row$Val!="All") {
        if(df_row$isRange) {
          idx_from_criterion <- (df_reference_file_list[,df_row$SampleColumn] >= df_row$Min) & (df_reference_file_list[,df_row$SampleColumn] <= df_row$Max)
        } else {
          idx_from_criterion <- (df_reference_file_list[,df_row$SampleColumn] ==  df_row$Val)
        }
        idx <- idx & idx_from_criterion
      }
      df_files[i,df_row$SampleColumn] <- df_row$Val
      reference_file_name <- paste0(reference_file_name,df_row$Val)
    }
    df_files[i,"sample"] <- reference_file_name
    reference_file_name <- paste0("pre_built_ref_",reference_file_name,"_",bed_type)
    df_files[i,"name"] <- reference_file_name
    reference_file_name <- paste0(reference_file_name,".tsv")
    reference_file_path <- file.path(output_path,reference_file_name)
    
    if(sum(idx) >= as.integer(min_reference_samples)) {
      gv$log <<- paste0(isolate(gv$log),flog.info("Need reference %s with %d samples",reference_file_name,sum(idx)))
      if((!file.exists(reference_file_path) || no_clobber==FALSE)) {
        gv$log <<- paste0(isolate(gv$log),flog.info("Making: %s",reference_file_path))
        makeReference(df_reference_file_list[idx,],reference_file_path,df_gc)
      } else {
        gv$log <<- paste0(isolate(gv$log),flog.info("Already have: %s",reference_file_path))
      }
      ret <- append(ret,"reference_file_name")
      df_files[i,"path"] <- reference_file_path
    }
  }
  return(df_files[!is.na(df_files$path) & !is.null(df_files$path),])
}
# -------------------------------------------------------------------

# -------------------------------------------------------------------
# Pre-build using reference and put in with ref_median, ref_mad, ref_var, N, Ns 
# Read in files specified in df_reference_file_list
# Compute stats
# Save for later
# Assumed to contain featureCounts output.
# From command line using saf file
# Geneid  Chr  Start	End	Strand	Length	/researchers/gisela.mirarnau/AgilentSureSelect_SUPER/160302_NS500817_0052_AH2F3TBGXY/15K0577/Bam/15K0577.bam
# chr1:12807-13219	chr1	12807	13219	+	413	121
# chr1:131415-132302	chr1	131415	132302	+	888	298
# chr1:133454-133812	chr1	133454	133812	+	359	73
# From R version:
# X.pathology.NGS.Samples.Molpath.160201_M00139_0305_000000000.AM03V.14M4928.1.14M4928.1.bam
# chr1:0-5000	0
# chr1:5000-10000	0
# -------------------------------------------------------------------
makeReference = function (df_reference_file_list,reference_file_name,df_gc) {
  file_names <- df_reference_file_list$CountsPath
  references_exist <- file.exists(file_names)
  num_existing_references <- sum(references_exist)
  if(num_existing_references==0) {
    gv$log <<- paste0(isolate(gv$log),flog.error("Unable to make reference file: %s",reference_file_name)) 
    return(NULL)
  } else {
    if(num_existing_references<nrow(df_reference_file_list)) {
      gv$log <<- paste0(isolate(gv$log),flog.warn("Unable to find all the inputs for reference file: %s",reference_file_name)) 
      gv$log <<- paste0(isolate(gv$log),flog.warn("Missing: %s",file_names[!references_exist])) 
      file_names <- file_names[references_exist]
      df_reference_file_list <- df_reference_file_list[references_exist,]
    }
  }
  df_counts <-  do.call("cbind",lapply(file_names,
                                       function (file_name) {
                                         if(file.exists(file_name)) { 
                                            gv$log <<- paste0(isolate(gv$log),flog.info("Reading %s",file_name))
                                            read.table(file_name,header=TRUE,stringsAsFactors=FALSE)
                                         } 
                                         else if(file.exists(paste0(file_name,".gz"))) { 
                                           gv$log <<- paste0(isolate(gv$log),flog.info("Reading %s.gz",file_name))
                                           read.table(paste0(file_name,".gz"),header=TRUE,stringsAsFactors=FALSE)
                                         } else {
                                          gv$log <<- paste0(isolate(gv$log),flog.error("Unable to read %s or %s.gz",file_name,file_name))
                                          return(NULL)                                             
                                         }
                                       }
                                       ))
  # Select the columns with the counts in them
  col_idx <- grep(paste0(df_reference_file_list$SampleName,sep='',collapse='|'),colnames(df_counts))
  df_counts <- df_counts[,col_idx]
  
  
  gv$log <<- paste0(isolate(gv$log),flog.info("Calculating references"))
  if(nrow(df_gc) != nrow(df_counts)) {
    gv$log <<- paste0(isolate(gv$log),flog.error("GC file differs in rows to counts table: %d vs %s",nrow(df_gc),nrow(df_counts)))
  }
  # Align the rows of counts file with gc file
  idx_chr_start_end <- paste0(df_gc$chr,":",as.integer((df_gc$start-1)),"-",df_gc$end)
  df_counts <- df_counts[idx_chr_start_end,]
  # For a reference we want everything to be 2N
  # For Men X -> 2X, Y -> 2Y, For women set Y=0 
  x_idx <- grepl("X",df_gc$chr)
  y_idx <- grepl("Y",df_gc$chr)
  auto_idx <- (!(x_idx | y_idx))
  have_sex <- FALSE
  
  if(sum("Sex"==colnames(df_reference_file_list))==1) {
    female_idx <- grepl("Female",df_reference_file_list$Sex)
    df_reference_file_list$Sex[female_idx] <- "Female"
    gv$log <<- paste0(isolate(gv$log),flog.info("From table: Males:Females = %d:%d",sum(!female_idx),sum(female_idx)))
    have_sex <- TRUE
  } 
  gv$log <<- paste0(isolate(gv$log),flog.info("Guessing Sex using X_counts/Autosome_counts"))
  auto_means <- colMeans(df_counts[auto_idx,],na.rm = TRUE)
  x_means <- colMeans(df_counts[x_idx,],na.rm = TRUE)
  y_means <- colMeans(df_counts[y_idx,])
  x_ratio <- x_means / auto_means
  y_ratio <- y_means / auto_means
  sex <- kmeans(x_ratio,2)
  female_idx <- sex$cluster==2 
  gv$log <<- paste0(isolate(gv$log),flog.info("From counts: Males:Females = %d:%d",sum(!female_idx),sum(female_idx)))
  if(!have_sex) { df_reference_file_list$Sex <- ifelse(female_idx,"Female","Male") }
  # df_counts[y_idx,female_idx] <- 0
  df_counts[x_idx,!female_idx] <- df_counts[x_idx,!female_idx] * 2
  df_counts[y_idx,!female_idx] <- df_counts[y_idx,!female_idx] * 2
  
  # Normalise the counts with the lower end of the dist truncated
  df_counts_with_na <- df_counts
  df_counts_with_na[df_counts < NORM_CUTOFF] <- NA
  
  sampleMedians <- colMedians(as.matrix(df_counts_with_na),na.rm=TRUE)
  ensembleMedian <- median(sampleMedians)
  # Counts median normalised 
  nc <- as.matrix(df_counts) %*% diag(ensembleMedian/sampleMedians) 
  nc[df_counts < NORM_CUTOFF] <- NA
  row_medians <- rowMedians(nc,na.rm=TRUE)
  df_counts <- data.frame(chr = df_gc$chr, 
                    start = df_gc$start, 
                    end = df_gc$end, 
                    gc = df_gc$gc, 
                    raw = row_medians,
                    cor = row_medians, 
                    ref_median = rowMedians(nc,na.rm=TRUE), 
                    ref_mad = rowMads(nc,na.rm=TRUE), 
                    ref_var = rowVars(nc,na.rm=TRUE), 
                    ref_sd = rowSds(nc,na.rm=TRUE), 
                    ref_N = rowSums(!is.na(nc)), 
                    N = 2 * row_medians/ensembleMedian, 
                    Ns = 2 * row_medians/ensembleMedian,
                    stringsAsFactors=FALSE)
  gv$log <<- paste0(isolate(gv$log),flog.info("Segmenting"))
  cbc_list <- do_cbc(df_counts$N,df_counts)   
  df_counts$Ns <-  cbc_list$Ns
  gv$log <<- paste0(isolate(gv$log),flog.info("Writing: %s",reference_file_name))
  write.table(df_counts,reference_file_name,row.names=FALSE,sep="\t",quote = FALSE)
} 
# -------------------------------------------------------------------

# -------------------------------------------------------------------
segments.pval  <- function (x, ngrid = 100, tol = 1e-06, alpha = 0.05, search.range = 100, nperm = 1000) 
{
  if (!inherits(x, "DNAcopy")) 
    stop("First arg must be the result of segment")
  xdat <- x$data
  xout <- x$output
  nsample <- ncol(xdat) - 2
  sampleid <- colnames(xdat)[-(1:2)]
  chrom0 <- xdat$chrom
  maploc0 <- xdat$maploc
  uchrom <- unique(chrom0)
  nchrom <- length(uchrom)
  bstat <- pval <- lcl <- ucl <- rep(NA, nrow(xout))
  ll <- 0
  iisamp <- 2
  for (isamp in sampleid) {
    iisamp <- iisamp + 1
    genomdat <- xdat[, iisamp]
    ina <- which(is.finite(genomdat))
    genomdat <- genomdat[ina]
    chrom <- chrom0[ina]
    maploc <- maploc0[ina]
    for (ichrom in uchrom) {
      kk <- sum(1 * (xout$ID == isamp & xout$chrom == ichrom))
      if (kk > 1) {
        gendat <- genomdat[chrom == ichrom]
        seglen <- xout$num.mark[xout$ID == isamp & xout$chrom == 
                                  ichrom]
        segmean <- xout$seg.mean[xout$ID == isamp & xout$chrom == 
                                   ichrom]
        xresid <- gendat - rep(segmean, seglen)
        ibstat <- ipval <- ilcl <- iucl <- rep(NA, kk)
        lo <- 1
        hi <- sum(seglen[1:2])
        for (i in 1:(kk - 1)) {
          gendati <- gendat[lo:hi]
          xresidi <- xresid[lo:hi]
          gendati <- (gendati - mean(gendati))/sd(xresidi)
          n <- length(gendati)
          zzz <- .Fortran("bsegp", as.integer(n), as.double(gendati), 
                          ostat = double(1), pval = double(1), as.integer(ngrid), 
                          as.double(tol), PACKAGE = "DNAcopy")
          ibstat[i] <- zzz$ostat
          ipval[i] <- zzz$pval
          k <- seglen[i]
          sr <- c(max(2, k - search.range), min(n - 2, 
                                                k + search.range))
          sumxk <- sum(gendati[1:k])
          var.factor <- n/((1:n) * (n:1 - 1))
          var.factor[n] <- 0
          zzz <- .Fortran("bsegci", as.integer(n), as.integer(k), 
                          as.double(sumxk), as.double(gendati), px = double(n), 
                          sr = as.integer(sr), vfact = as.double(var.factor), 
                          as.integer(nperm), bsloc = integer(nperm), 
                          PACKAGE = "DNAcopy")
          bsloc <- zzz$bsloc
          bsci <- quantile(bsloc, c(alpha/2, 1 - alpha/2), 
                           type = 1)
          ilcl[i] <- bsci[1]
          iucl[i] <- bsci[2]
          lo <- lo + seglen[i]
          if (i < kk - 1) 
            hi <- hi + seglen[i + 2]
        }
        ibstat[kk] <- ipval[kk] <- ilcl[kk] <- iucl[kk] <- NA
      }
      else {
        seglen <- ibstat <- ipval <- ilcl <- iucl <- NA
      }
      if(kk>0) {
        bstat[ll + (1:kk)] <- ibstat
        pval[ll + (1:kk)] <- ipval
        lcl[ll + (1:kk)] <- maploc[chrom == ichrom][cumsum(seglen) + 
                                                      (ilcl - seglen)]
        ucl[ll + (1:kk)] <- maploc[chrom == ichrom][cumsum(seglen) + 
                                                      (iucl - seglen)]
      }
      ll <- ll + kk
    }
  }
  cbind(xout, bstat, pval, lcl, ucl)
}
# -------------------------------------------------------------------

# -------------------------------------------------------------------
readAllReferences <- function (bed_type) {
  output_path <- makepath(panel = CFG$panel_path,run = CFG$run_path,filetype = "reference")
  pattern <- paste0("^pre_built_ref.*",bed_type,".tsv$")
  file_list <- list.files(output_path,pattern)
  ref_names <- gsub("pre_built_ref_","",(gsub(paste0("_",bed_type,".tsv"),"",file_list)))
  file_paths <- file.path(output_path,file_list)
  refs <- lapply(file_paths,read.table,header=TRUE,stringsAsFactors=FALSE)
  names(refs) <- ref_names
  return(refs)
}
# -------------------------------------------------------------------

# -------------------------------------------------------------------
findBestReference <- function (refs,df) {
  comp <- lapply(names(refs),applyReference,refs,df)
  refNames <-  names(refs)
  names(comp) <- refNames
  
  nRefs <- length(comp)
  sm  <- double(nRefs)
  cn  <- double(nRefs)
  seg  <- double(nRefs)
  for(i in 1:nRefs) {
    sm[i] <- comp[[i]]$comparisons$smoothness_diff
    cn[i] <- comp[[i]]$comparisons$quantised_cn_diff
    seg[i] <- comp[[i]]$comparisons$segments
  }
  plot(-log10(sm),-log10(cn))
  plot(-log10(sm),seg)
  plot(-log10(cn),seg)
  df_best <- comp[[which.min(sm)]]
  df_best$refname <- refNames[which.min(sm)]
  # df_best$corrected 
  # df_best$comparisons$smoothness_diff 
  # df_best$comparisons$quantised_cn_diff 
  # df_best$comparisons$segments
  return(df_best)
}
# -------------------------------------------------------------------

# -------------------------------------------------------------------
applyReference  <- function (ref_name,refs,df) {
  df_reference <- refs[[ref_name]]
  gv$log <<- paste0(isolate(gv$log),flog.info("Doing applyReference(%s).",ref_name))
  
  # Median normalise truncated raw counts
  idx <- !is.na(df$raw) & df$raw>=MEDIAN_NORM_CUTOFF & !df$blackListed
  nc  <-  df$raw / median(df$raw[idx])
  idx <- !is.na(df_reference$ref_median) & df_reference$ref_median>=MEDIAN_NORM_CUTOFF & !df$blackListed 
  nref  <-  df_reference$ref_median / median(df_reference$ref_median[idx])
  ratio <- nc / nref
  idx <- !is.nan(ratio) & is.finite(ratio) & df$raw>=COUNTS_CUTOFF & df_reference$raw>=COUNTS_CUTOFF & !df$blackListed 
  cor <- ratio
  cor[!idx] <- NA
  # Strip out divide by zero stuff for evaluating metrics
  ratio <- ratio[idx]
  nmax <- length(ratio)
  fom <- list()
  cn_mode <- mlv(cor[idx & df$raw>MODE_NORM_CUTOFF & df_reference$raw>MODE_NORM_CUTOFF], method = "lientz", bw = 0.2)$M
  cn_median <- median(cor[idx])
  warn_norm <- (abs(cn_mode-cn_median)/(cn_mode+cn_median) > 0.1)
  if(warn_norm) { gv$log <<- paste0(isolate(gv$log),flog.info("##############################################"))}
  gv$log <<- paste0(isolate(gv$log),flog.info("::Normalising %s: Mode = %e Median = %e",ref_name,cn_mode,cn_median))
  if(warn_norm) { gv$log <<- paste0(isolate(gv$log),flog.info("##############################################"))}
  cor <- cor / cn_mode
  N = cor * 2
  if(nmax>1) {
    gv$log <<- paste0(isolate(gv$log),flog.info("Segmenting %s: Found %d bins.",ref_name,nmax))
    fom$smoothness_diff <- mad( ratio[2:nmax] - ratio[1:(nmax-1)]) 
    fom$quantised_cn_diff <- mad(ratio-1)
    cbc_list <- do_cbc(N,df_reference)    
    fom$segments <- nrow(cbc_list$df_seg)
    df$Ns <- cbc_list$Ns
    df$cor <- cor 
    df$N <- N
    df$ref_median <- df_reference$ref_median
    df$ref_mad <- df_reference$ref_mad
    df$ref_var <- df_reference$ref_var
    df$ref_sd <-df_reference$ref_sd
    df$ref_N <- df_reference$ref_N
  }
  else {
    # Flag error
    gv$log <<- paste0(isolate(gv$log),flog.error("No valid corrected bins found in applyReference()."))
  }
  gv$log <<- paste0(isolate(gv$log),flog.info("::Applying %s: dn/dbase = %f mad = %f segments = %f.",
                                              ref_name,fom$smoothness_diff,fom$quantised_cn_diff,fom$segments))
  ret <- list(comparisons = fom, corrected = df)
  gv$log <<- paste0(isolate(gv$log),flog.info("Returning corrected counts dim = %d x %d",dim(df)[1],dim(df)[2]))
  return(ret)
}
# -------------------------------------------------------------------


# -------------------------------------------------------------------
do_cbc <- function (N,df_reference) {
  Ns <- N
  mad_on_med <- median(df_reference$ref_mad/df_reference$ref_median,na.rm=TRUE)
  ref_mad <- df_reference$ref_mad
  idx <- !is.na(ref_mad) & ref_mad==0
  ref_mad[idx] <- mad_on_med * df_reference$ref_median[idx] 
  idx <- !is.na(df_reference$ref_median) & df_reference$ref_median!=0 & !is.na(ref_mad)
  weights <- 1/ref_mad^2
  weights <- weights[idx]/sum(weights[idx])
  CNA.object <- CNA(log2(N[idx]),df_reference$chr[idx],df_reference$start[idx],data.type="logratio",presorted=TRUE)
  smoothed.CNA.object <- smooth.CNA(CNA.object)
  segment.smoothed.CNA.object <- segment(x = smoothed.CNA.object,weights = weights, verbose = SEGMENT_VERBOSITY,alpha = CBS_ALPHA, nperm = CBS_NPERM )
  numSegments <- nrow(segment.smoothed.CNA.object$output)
  
  idx_na <- is.na(segment.smoothed.CNA.object$segRows$startRow)
  segment.smoothed.CNA.object$segRows$startRow[idx_na] <- segment.smoothed.CNA.object$segRows$endRow[idx_na] # repair bug in CNA
  if(sum(is.nan(segment.smoothed.CNA.object$segRows$startRow) | 
           is.na(segment.smoothed.CNA.object$segRows$startRow) |
           is.nan(segment.smoothed.CNA.object$segRows$endRow) |
           is.na(segment.smoothed.CNA.object$segRows$endRow) 
  )) {
    gv$log <<- paste0(isolate(gv$log),flog.error("Invalid smoothed segmented value."))  
  } 
  # Convert to old co-ordinates
  idx_conv <- 1:length(idx)
  idx_conv <- idx_conv[idx]
  for(i in 1:numSegments) {
    startSeg <- idx_conv[segment.smoothed.CNA.object$segRows$startRow[i]]
    endSeg   <- idx_conv[segment.smoothed.CNA.object$segRows$endRow[i]]
    Ns[startSeg:endSeg] <- 2^(segment.smoothed.CNA.object$output$seg.mean[i])
  }
  df_seg <- segments.pval(segment.smoothed.CNA.object,alpha = CBS_ALPHA, nperm = CBS_NPERM )
  colnames(df_seg) <- c("Gene","Chr","Start","End","Bins","log2N","bstat","pval","lcl","ucl")
  ret <- list(df_seg = df_seg, Ns = Ns, cna_output = segment.smoothed.CNA.object)
  return(ret)
}
# -------------------------------------------------------------------

# -------------------------------------------------------------------
make_bams <- function() {
  no_clobber <- as.logical(CFG$no_clobber)
  df_bamfiles <- DF_CFG[grepl("bams",DF_CFG$section),]
  chrom_file <- CFG$chrom_file
  chrom_index_file <- CFG$chrom_index_file
  number_of_cores <- CFG$number_of_cores
  num_index_files <- length(list.files(path = dirname(chrom_index_file), pattern = paste0(basename(chrom_index_file),".*")))
  if(num_index_files == 0) {
    gv$log <<- paste0(isolate(gv$log),flog.info("Making index for %s",chrom_file))
    buildindex(chrom_index_file,chrom_file)
  } else {
    gv$log <<- paste0(isolate(gv$log),flog.info("Already have: %s.*",chrom_index_file))
  }
  for (bampath in df_bamfiles$value) {
    sample <- gsub(".bam$","",basename(bampath))
    bam_root <- makepath(panel = CFG$panel_path ,run = CFG$run_path ,sample = sample ,filetype = "bam")
    fastq_root <- dirname(makepath(panel = CFG$panel_path ,run = CFG$run_path ,sample = sample ,filetype = "fastq"))
    bampath <- paste0(bam_root,".bam")
    chrom_index_file <- CFG$chrom_index_file
    
    if((!file.exists(bampath) || no_clobber==FALSE)) {
      gv$log <<- paste0(isolate(gv$log),flog.info("Aligning: making %s",bampath))
      fastqpaths <- list.files(path = fastq_root,"*fastq.gz",full.names = TRUE)
      if(length(fastqpaths) != 2) {
        gv$log <<- paste0(isolate(gv$log),flog.error("Unable to read fastq file matching %s (%d matches)",fastq_root,length(fastqpath)))
        return(FALSE)
      }
      cwd <- getwd()
      fastqpaths <- file.path(cwd,fastqpaths)
      chrom_index_file <- file.path(cwd,chrom_index_file)
      setwd(dirname(bam_root))
      align(chrom_index_file,
            fastqpaths[1],fastqpaths[2],
            type="dna",
            input_format="gzFASTQ",
            output_format="BAM",
            output_file = basename(bampath),
            nthreads = number_of_cores,
            detectSV=TRUE
            )
      lapply(c("fixbam","sam","tmp","bin"), remove_files_by_suffix)
      setwd(cwd)
    } else {
      gv$log <<- paste0(isolate(gv$log),flog.info("Already have: %s",bampath))
    }
  }
}
# -------------------------------------------------------------------

# -------------------------------------------------------------------
make_vcfs <- function() {
  no_clobber <- as.logical(CFG$no_clobber)
  df_files <- DF_CFG[grepl(VCF_SECTION,DF_CFG$section),]
  df_files$file <- "vcf"
  for (j in 1:nrow(df_files)) {
    bampath <- df_files[j,"value"]
    sample <- gsub(".bam$","",basename(bampath))
    bam_root <- makepath(panel = CFG$panel_path ,run = CFG$run_path ,sample = sample ,filetype = "bam")
    bampath <- paste0(bam_root,".bam")
    vcf_root <- makepath(panel = CFG$panel_path ,run = CFG$run_path ,sample = sample ,filetype = "vcf")
    vcfpath <- paste0(vcf_root,".vcf")
    vcffile <- basename(vcfpath)
    refGenomeFile <- CFG$chrom_file
    snp_annotation_file <- CFG$snp_annotation_file
    if((!file.exists(vcfpath) || no_clobber==FALSE)) {
      gv$log <<- paste0(isolate(gv$log),flog.info("Making: %s",vcfpath))
      cwd <- getwd()
      bampath <- file.path(cwd,bampath)
      vcfpath <- file.path(cwd,vcfpath)
      refGenomeFile <- file.path(cwd,refGenomeFile)
      snp_annotation_file <- file.path(cwd,snp_annotation_file)
      setwd(dirname(bam_root))
      exactSNP(bampath,
               isBAM=TRUE,
               refGenomeFile=refGenomeFile,
               SNPAnnotationFile=snp_annotation_file,
               outputFile=vcffile,
               # qvalueCutoff=12,
               # minAllelicFraction=0,
               minAllelicBases=SNP_MIN_MISMATCHES,
               minReads=SNP_MIN_READS,
               maxReads=SNP_MAX_READS,
               # minBaseQuality=13,
               # nTrimmedBases=3,
               nthreads=1)
      gv$log <<- paste0(isolate(gv$log),flog.info("Copying %s to %s",vcffile,vcfpath))
      file.copy(from = vcffile, to = vcfpath)
      # file.remove(vcffile)
      lapply(c("fixbam","sam","tmp","bin"), remove_files_by_suffix)
      setwd(cwd)
    } else {
      gv$log <<- paste0(isolate(gv$log),flog.info("Already have: %s",vcfpath))
    }
    df_files[j,"path"] <- vcfpath
    df_files[j,"sample"] <- sample
  }
}
# -------------------------------------------------------------------

# -------------------------------------------------------------------
make_bafs <- function() {
  no_clobber <- as.logical(CFG$no_clobber)
  df_files <- DF_CFG[grepl(BAF_SECTION,DF_CFG$section),]
  df_files$file <- "baf"
  df_files$resolution <- BAF_SECTION
  df_files$is_targeted <- BAF_SECTION=="targeted"
  df_files$is_wgs <- BAF_SECTION=="targeted"
  
  for (j in 1:nrow(df_files)) {
    bampath <- df_files[j,"value"]
    sample <- sub(".bam$","",basename(bampath))
    baf_root <- makepath(panel = CFG$panel_path ,run = CFG$run_path ,sample = sample ,filetype = "baf")
    bafpath <- paste0(baf_root,"_baf.tsv")
    baffile <- basename(bafpath)
    vcf_root <- makepath(panel = CFG$panel_path ,run = CFG$run_path ,sample = sample ,filetype = "vcf")
    vcfpath <- paste0(vcf_root,".vcf")
    vcffile <- basename(vcfpath)
    df_files[j,"path"] <- bafpath
    df_files[j,"sample"] <- sample
    
    if((!file.exists(bafpath) || no_clobber==FALSE)) {
      vcf <- readVcf(vcfpath,genome = "hg19")
      if(nrow(vcf)>0) {
        idx <- isSNV(vcf, singleAltOnly=FALSE)
        snvinfo <- info(vcf[idx])
        df_baf <- data.frame(chr = paste0("chr",t(as.data.frame(strsplit(rownames(snvinfo),":")))[,1]), 
                             start = start(ranges(GenomicRanges::rowRanges(vcf[idx]))), 
                             end = end(ranges(GenomicRanges::rowRanges(vcf[idx]))), 
                             freq = as.double(snvinfo$MM)/as.double(snvinfo$DP), 
                             ref_reads = as.double(snvinfo$DP), 
                             alt_reads = as.double(snvinfo$MM),
                             stringsAsFactors = FALSE)
      }
      else {
        gv$log <<- paste0(isolate(gv$log),flog.warn("Empty VCF file. Writing a fake row."))    
        df_baf <- data.frame(chr = "MT", start = 1, end = 1, freq = 0, ref_reads = 0, alt_reads = 0,stringsAsFactors = FALSE)   
      }        
      idx <- grepl("GL",df_baf$chr) | grepl("MT",df_baf$chr) | df_baf$alt_reads >= SNP_MAX_READS | df_baf$ref_reads >= SNP_MAX_READS 
      df_baf <- df_baf[!idx,]
      gv$log <<- paste0(isolate(gv$log),flog.info("Writing %s with %d entries",bafpath,sum(idx)))    
      write.table(df_baf,bafpath,row.names=FALSE,quote = FALSE,sep="\t")
    } else {
      gv$log <<- paste0(isolate(gv$log),flog.info("Already have: %s",bafpath))
    }
  }
  return(df_files)
}
# -------------------------------------------------------------------

# -------------------------------------------------------------------
make_deletions <- function(bedfile) {
  df_files <- DF_CFG[grepl(DEL_SECTION,DF_CFG$section),]
  whitelist_path = CFG$common_deletion_list
  no_clobber <- as.logical(CFG$no_clobber)
  bed_type <- BEDFILENAME2BEDTYPE[bedfile]
  bam_type <- BEDTYPE2BAMTYPE[bedfile]
  df_whitelist <- read.table(whitelist_path,header=TRUE,stringsAsFactors=FALSE)
  df_files$file <- "deletions"
  df_files$resolution <- bed_type
  df_files$is_wgs <- (bam_type=="wg")
  df_files$is_targeted <- !df_files$is_wgs
  
  for (j in 1:nrow(df_files)) {
    bampath <- df_files[j,"value"]
    sample <- sub(".bam$","",basename(bampath))
    del_root <- makepath(panel = CFG$panel_path ,run = CFG$run_path ,sample = sample ,filetype = "sv")
    del_ending <- paste0("_",bed_type,"_del.tsv")
    delpath <- paste0(del_root,del_ending)
    delfile <- basename(delpath)
    cnv_suffix <- paste0("_",bed_type,".tsv")
    counts_root <- makepath(panel = CFG$panel,run = CFG$run_path,sample = sample,filetype = "cnv")
    corrected_counts_path <- paste0(counts_root,cnv_suffix)
    df_files[j,"sample"] <- sample
    df_files[j,"path"] <- delpath
    if((!file.exists(delpath) || no_clobber==FALSE)) {
      df_counts <- read.table(corrected_counts_path,header=TRUE,stringsAsFactors=FALSE)
      gv$log <<- paste0(isolate(gv$log),flog.info("Have %s with %d rows",corrected_counts_path,nrow(df_counts)))
      cbc_list <- do_cbc(df_counts$N,df_counts)
      df_del <- cbc_list$df_seg
      gv$log <<- paste0(isolate(gv$log),flog.info("Found %d segments",nrow(df_del)))     
      df_del$whiteListed <- FALSE
      df_del$Arm <- ""
      df_del$N <- 2^df_del$log2N
      df_del$weighted_mean <- 0
      df_del$weighted_var <- 0
      df_del$weighted_sem <- 0
      idx_ok <- !is.na(df_counts$ref_mad) & 
        !is.na(df_counts$ref_median) & 
        df_counts$ref_mad!=0 & 
        df_counts$ref_median!=0 & 
        df_counts$blackListed==FALSE
      
      for (n in 1:nrow(df_del)) {
        row <- df_del[n,]
        idx_counts <- df_counts$start >= row$Start &
          df_counts$end <= row$End &
          df_counts$chr == row$Chr & 
          idx_ok
        numCounts <- sum(idx_counts)
        if(numCounts > 1) {
          weights <- 1 / (df_counts$ref_mad[idx_counts])
          counts <- df_counts$N[idx_counts]
          weights <- weights / sum(weights) # normalise
          df_del$weighted_mean[n] <- sum(counts * weights,na.rm=TRUE)
          df_del$weighted_var[n] <- weighted.var (counts,weights,na.rm=TRUE)
          df_del$weighted_sem[n] <- weighted.var.se(counts,weights,na.rm=TRUE)
          # If the pval is NA then it's a whole-chromosome event
          if(!is.na(row$pval) || (df_del$weighted_mean[n] - 2)^2 > 3 * df_del$weighted_sem[n]) {
            idx <- row$Chr == df_whitelist$Chr &
              row$Start >= df_whitelist$Start &
              row$End <= df_whitelist$End &
              0 != df_whitelist$GainLoss
              # row$Start <= df_whitelist$RequiredStart &
              # row$End >= df_whitelist$RequiredEnd &
              # "pq" != df_whitelist$Arm) $
            if(sum(idx)==1) {
              df_del$whiteListed[n] <- TRUE
              df_del$Arm[n] <- df_whitelist$Arm[idx][1]
              df_del$Gene[n] <- df_whitelist$Gene[idx][1]
              gain_loss <- ifelse(df_whitelist$GainLoss[idx][1]>0,"duplication","deletion")
              gv$log <<- paste0(isolate(gv$log),
                               flog.info("Found whitelisted entry: %s%s",
                                         row$Chr,df_del$Arm[n]))
            }
          }
        }
      }
      gv$log <<- paste0(isolate(gv$log),flog.info("Writing %s with %d entries",delpath,nrow(df_del)))    
      # df_del <- df_del[df_del$whiteListed,]
      # df_del$whiteListed <- NULL
      colnames(df_del)[colnames(df_del)=="Chr"] <- "chr"
      colnames(df_del)[colnames(df_del)=="Start"] <- "start"
      colnames(df_del)[colnames(df_del)=="End"] <- "end"
      colnames(df_del)[colnames(df_del)=="Bins"] <- "bins"
      colnames(df_del)[colnames(df_del)=="Gene"] <- "gene"
      write.table(df_del,delpath,row.names=FALSE,quote = FALSE,sep="\t")
    } else {
      gv$log <<- paste0(isolate(gv$log),flog.info("Already have: %s",delpath))
    }
  }
  return(df_files)
}
# -------------------------------------------------------------------

# -------------------------------------------------------------------
# Computes the variance of a weighted mean following 
# Cochran W.G (1977) Sampling Techniques (3rd Edn). Wiley, New York
# Copied from http://stats.stackexchange.com/questions/25895/computing-standard-error-in-weighted-mean-estimation 
weighted.var.se <- function(x, w, na.rm=FALSE) {
  if (na.rm) { w <- w[i <- !is.na(x)]; x <- x[i] }
  n = length(w)
  xWbar = wtd.mean(x,w,na.rm=na.rm)
  wbar = mean(w)
  out = n/((n-1)*sum(w)^2)*(sum((w*x-wbar*xWbar)^2)-2*xWbar*sum((w-wbar)*(w*x-wbar*xWbar))+xWbar^2*sum((w-wbar)^2))
  return(out)
}
# -------------------------------------------------------------------

# -------------------------------------------------------------------
weighted.var <- function(x, w, na.rm = FALSE) {
  if (na.rm) {
    w <- w[i <- !is.na(x)]
    x <- x[i]
  }
  sum.w <- sum(w)
  sum.w2 <- sum(w^2)
  mean.w <- sum(x * w) / sum(w)
  (sum.w / (sum.w^2 - sum.w2)) * sum(w * (x - mean.w)^2, na.rm = na.rm)
}
# -------------------------------------------------------------------

# -------------------------------------------------------------------
# Not doing this for the moment
make_deletion_backgrounds_distributions <- function(ref_cn_tables) {
  # nbins from counts
  # permutations = nbins * factorial(nreference) or nbins * nreference
  # df_pr_null <- matrix(data = 0, nrow = permutations,ncol = 45,dimnames = c("breakpoint","Chr_Dir"),colnames = -22-22)
  # Loop over chromosomes, directions, replicates and breakpoints
  # Compute FET for each breakpoint
  # Does it vary for chromosome?
  # Does it vary for position?
  # Does it vary for resolution?
  # When you include true positives what does the clustering say?
  # Look at statistic in CBC paper
}
# -------------------------------------------------------------------

# -------------------------------------------------------------------
find_putative_deletions <- function(start, end, cn_table, ref_cn_tables) {
  df_pr <- matrix(data = 0, nrow = nbins,ncol = 45,dimnames = c("breakpoint","Chr_Dir"),colnames = -22:22)  
}
# -------------------------------------------------------------------

# -------------------------------------------------------------------
make_translocations <- function() {
  #Chr  Location	Chr	Location	SameStrand	nSupport  
  # 12M1411.bam.breakpoints.txt 
  # 12M1411.bam.fusion.txt
  # ID  Partner	Gene	Chr	Arm	Band	Start	End	Strand	MinSupport
  df_files <- DF_CFG[grepl(FUSION_SECTION,DF_CFG$section),]
  no_clobber <- as.logical(CFG$no_clobber)
  fusions_white_list <- CFG$fusions_white_list
  df_whitelist <- read.table(fusions_white_list,header=TRUE,stringsAsFactors=FALSE)
  df_whitelist$Chr <- sub("^chr","",df_whitelist$Chr)
  df_files$file <- "fusion"
  df_files$resolution <- FUSION_SECTION
  df_files$is_targeted <- FUSION_SECTION=="targeted"
  df_files$is_wgs <- FUSION_SECTION=="targeted"
  
  for (j in 1:nrow(df_files)) {
    bampath <- df_files[j,"value"]
    sample <- sub(".bam$","",basename(bampath))
    fusions_suffix <- paste0("_fusions.tsv")
    fusions_root <- makepath(panel = CFG$panel,run = CFG$run_path,sample = sample,filetype = "sv")
    fusion_list <- paste0(fusions_root,fusions_suffix)
    putative_fusion_list <- sub(".bam$",".fusion.txt",bampath)
    df_files[j,"path"] <- fusion_list
    df_files[j,"sample"] <- sample
    if(!file.exists(fusion_list) || no_clobber==FALSE) {
      if(file.exists(putative_fusion_list)) {
        gv$log <<- paste0(isolate(gv$log),flog.info("Reading %s",putative_fusion_list))
        df_putative_fusions <- read.table(putative_fusion_list,header=TRUE,stringsAsFactors=FALSE)
        gv$log <<- paste0(isolate(gv$log),flog.info("Making: %s",fusion_list))
        colnames(df_putative_fusions) <- c("chr",  "start",  "chr2",	"end",	"sameStrand",	"nSupport")
        rownames(df_putative_fusions) <- df_putative_fusions$ID
        df_putative_fusions$idx1 <- 0
        df_putative_fusions$idx2 <- 0
        df_putative_fusions$whiteListed <- FALSE
        df_putative_fusions$name <- ""
        for (n in 1:nrow(df_putative_fusions)) {
          putative_fusion <- df_putative_fusions[n,]
          idx1 <- (putative_fusion$chr == df_whitelist$Chr & 
                     putative_fusion$start >= df_whitelist$Start & 
                     putative_fusion$start <= df_whitelist$End & 
                     putative_fusion$nSupport >= df_whitelist$MinSupport)
          idx2 <- (putative_fusion$chr2 == df_whitelist$Chr & 
                        putative_fusion$end >= df_whitelist$Start & 
                        putative_fusion$end <= df_whitelist$End & 
                        putative_fusion$nSupport >= df_whitelist$MinSupport)
          t1 <-which(idx1)[1] 
          t2 <-which(idx2)[1] 
          # df_putative_fusions[n,"start"] <- ifelse(sum(idx1)>0,t1,0)
          # df_putative_fusions[n,"end"] <- ifelse(sum(idx2)>0,t2,0)
          if( sum(idx1) > 0 && sum(idx2) > 0 ) {
            partners <- as.vector(df_whitelist[t1,"ID"])
            if(sum(partners == as.vector(df_whitelist[t2,"ID"]))) {
              df_putative_fusions[n,"whiteListed"] <- TRUE
              df_putative_fusions[n,"name"] <- paste0(
                "t(",
                  putative_fusion$chr,";",putative_fusion$chr2,
                ") ",
                "(",
                  df_whitelist[t1,"Arm"],df_whitelist[t1,"Band"],";",
                  df_whitelist[t2,"Arm"],df_whitelist[t2,"Band"],
                ") ",
                df_whitelist[t1,"Arm"],df_whitelist[t1,"Gene"]," / ",
                df_whitelist[t2,"Arm"],df_whitelist[t2,"Gene"]
              ) # paste
            }
          }
        } #for
        # df_fusions <- df_putative_fusions[df_putative_fusions$WhiteListed,]
        df_fusions <- df_putative_fusions
        if(sum(grepl("^chr",df_fusions$chr)) ==0) { df_fusions$chr  <- paste0("chr",df_fusions$chr ) }
        if(sum(grepl("^chr",df_fusions$chr2))==0) { df_fusions$chr2 <- paste0("chr",df_fusions$chr2) }
        write.table(df_fusions,fusion_list,row.names=FALSE,col.names=TRUE,quote = FALSE,sep="\t")
      } else {
        gv$log <<- paste0(isolate(gv$log),flog.warn("Can't find: %s. Need to run gridss or Rsubread::Align().",putative_fusion_list))
      }
    } else {
      gv$log <<- paste0(isolate(gv$log),flog.info("Already have: %s",fusion_list))
    }
  } # for
  return(df_files)
}
# -------------------------------------------------------------------

# -------------------------------------------------------------------
make_bins <- function(bin_size,bedfile) {
  if(is.na(bedfile) || length(bedfile)<1) {
    browser()
    gv$log <<- paste0(isolate(gv$log),flog.warn("Skipping bedfile for bin size %d",bin_size))  
    df <- data.frame(chr = character(0), start = integer(0), end = integer(0),id = integer(0), score = integer(0), strand = character(0))
    return(df)
  }
  chrom_info_file <- CFG$chrom_info_file
  no_clobber <- as.logical(CFG$no_clobber)
  destfile <- bedfile
  chromInfo <- read.table(chrom_info_file,stringsAsFactors=FALSE)
  colnames(chromInfo) <- c("chr","end","file")
  rownames(chromInfo) <- chromInfo$chr
  if((!file.exists(destfile) || no_clobber==FALSE)) {
    gv$log <<- paste0(isolate(gv$log),flog.info("Making: %s",destfile))
    chr_names <- chromInfo$chr
    chr_names <- chr_names[!grepl("[_M]",chr_names)]
    make_intervals <- function(chr,chromInfo) { 
      options(scipen=20)
      gv$log <<- paste0(isolate(gv$log),flog.info("Making bins for %s:%s",chr,destfile))    
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
    gv$log <<- paste0(isolate(gv$log),flog.info("Making bins for %s: Done",destfile))    
  }
  gv$log <<- paste0(isolate(gv$log),flog.info("Have: %s",destfile))
}
# -------------------------------------------------------------------

# -------------------------------------------------------------------
make_gc <- function(svalue) {
  bedfile <- CFG[[svalue]]
  chrom_file <- CFG$chrom_file
  no_clobber <- as.logical(CFG$no_clobber)
  destfile <- bedrow2gcrow(svalue)$value
  
  if((!file.exists(destfile) || no_clobber==FALSE)) {
    gv$log <<- paste0(isolate(gv$log),flog.info("Making GC corrections. Reading: %s",bedfile))
    intervals <- import.bed(bedfile)
    gv$log <<- paste0(isolate(gv$log),flog.info("Making GC corrections. Reading: %s",chrom_file))
    reference <- readDNAStringSet(chrom_file)
    
    gv$log <<- paste0(isolate(gv$log),flog.info("Making: %s",destfile))
    # for each chr, something like this
    find_gc_content <- function(idx) { 
      seq <- reference[idx]
      # Assumes that naming fields in fasta file are chr2 with nothing else.
      gv$log <<- paste0(isolate(gv$log),flog.info("Making GC content for %s:%s",names(seq),destfile))    
      chr_int <- intervals[seqnames(intervals)==names(seq)] 
      if(length(chr_int)==0) {
        gv$log <<- paste0(
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
      gv$log <<- paste0(
        isolate(gv$log),
        flog.error("Seqname %s and reference %s have no sequence names in common. Mismatched chr naming convention?",names(seq),bedfile)
      )
    }
    colnames(df_gc) <- c("chr","start","end","gc")
    write.table(df_gc,destfile,row.names=FALSE,quote = FALSE,sep="\t")
    gv$log <<- paste0(isolate(gv$log),flog.info("Making GC content for %s: Done",destfile))  
    return(df_gc)
  }
  gv$log <<- paste0(isolate(gv$log),flog.info("Have: %s",destfile))
}
# -------------------------------------------------------------------

# -------------------------------------------------------------------
make_all_gc_plots <- function(output_path) {
  suffixes <- c("_targeted.tsv","_off_target.tsv","_medium.tsv","_coarse.tsv","_fine.tsv")
  filelist <- unlist(lapply(suffixes,function(x) { dir(output_path,pattern=x) }))
  lapply(file.path(output_path,filelist),make_gc_plot)
}
# -------------------------------------------------------------------

# -------------------------------------------------------------------
make_gc_plot <- function(gc_file) {
  gv$log <<- paste0(isolate(gv$log),flog.info("Plotting:%s",gc_file))  
  # Optionally not ocerwrite because this takes a while
  plot_file <- sub(".tsv",".png",gc_file)
  df_gc <- read.table(gc_file,header=TRUE,stringsAsFactors=FALSE)
  idx <- df_gc$cor>0 & df_gc$raw>0 & df_gc$gc>0.2
  df_gc <- df_gc[idx,]
  df_gc$cor <- log10(df_gc$cor)
  df_gc$raw <- log10(df_gc$raw)

  loess_predict <- loess(raw ~ gc, data = df_gc,  span = LOESS_SPAN) # subset = !blacklist
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
  bins <- seq(0, 1, by = LOESS_BIN_WIDTH)
  fapply_median <- function(n,bins) { median(df_gc$raw[df_gc$gc>=(bins[n] - LOESS_BIN_WIDTH) & df_gc$gc<(bins[n] + LOESS_BIN_WIDTH)]) }
  bin_medians <- unlist(lapply(1:length(bins),fapply_median,bins))
  loessCorrection_df <- data.frame(gc = bins, median = bin_medians)
  loessCorrection_df$cor <- predict(loess_predict,loessCorrection_df$gc)
  loessCorrection_df <- loessCorrection_df[!(is.na(loessCorrection_df$cor)|is.na(loessCorrection_df$median)),]
  
  gc <- loessCorrection_df$gc
  bin_medians <- loessCorrection_df$median
  bin_cor <- loessCorrection_df$cor
  
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
# -------------------------------------------------------------------

# -------------------------------------------------------------------
if(RSTUDIO_DEBUG) { 
  makeCNTables(args) 
  return(0)
}
# -------------------------------------------------------------------
if(isCommandLine) {
  makeCNTables(args) 
  quit(save = "no", status = 0)
}
# -------------------------------------------------------------------

