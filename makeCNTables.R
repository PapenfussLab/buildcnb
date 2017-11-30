#!/usr/bin/env Rscript
#!/usr/bin/Rscript --vanilla # ideally, but how to?
args <- commandArgs(TRUE)

######################################### 
# Global config defaults
######################################### 
CFG_DEFAULTS <- list(
  rsubread_version = "1.24.0", # Currently preferred version
  require_subread_version = TRUE, # break if the right one isn't used
  use_own_subread = TRUE, # look in ./lib/ for subread package
  rstudio_debug = FALSE,
  debug = FALSE,
  is_command_line = TRUE,
  spawn_external_aligner = TRUE, # execute external aligner using system2() ?
  wg_coarse_bin_size = 1000000,
  wg_medium_bin_size = 50000,
  wg_fine_bin_size = 5000,
  targeted_bin_size = 119, # supplied in the bed file, but for me it's about the bait size
  cbs_alpha = 0.01,
  cbs_nperm = 10000,
  loess_scales_badly_limit  = 100000, # when loess gets slow
  loess_min_usable_counts  = 50, # For fitting
  loess_min_usable_gc  = 0.25, 
  loess_span = 0.4,
  loess_bin_width = 0.01,
  minimum_overlap = 1, 
  norm_cutoff = 10,   # Counts in a bin less that this are thrown away for normalisation
  median_norm_cutoff = 10,   # Counts in a bin with a dist whose median less that this are thrown away for ALL normalisation
  mode_norm_cutoff = 100,   # Counts in a bin with a dist whose mode less that this are thrown away for mode-based normalisation
  max_proportion_zero_entries = 0.5, # If more than this proportion of counts is less that the cutoff, something is wrong.
  median_data_points_cutoff = 3,
  median_magic = 0.2, # Magic number.   mode_norm_cutoff = efault is for targeted coverage in the hundreds. For lower coverage, use 0.2 * median
  mode_data_points_cutoff = 3,
  maximum_female_y_on_auto_ratio = 0.125,
  minimum_female_x_on_auto_ratio = 0.75,
  minimum_samples_for_sex_determination_by_clustering = 16,
  counts_cutoff = 10,   # Counts less that this are thrown away for CN estimation
  snp_min_mismatches = 1,    # Before a SNP is called 
  snp_min_reads = 20,        # Before a SNP is called - these two also give the min B-allele freq
  snp_max_reads = 1000,   # For PCR artifacts # TODO: Change for RNA-seq
  segment_verbosity  = 0,
  fusion_section = "targeted",
  minimumFusionSuport = 10, # minimum reads spanning two breakpoints to be reported
  baf_section = "targeted",
  vcf_section = "targeted",
  del_section = "wg_medium",
  del_zscore_section = "targeted",
  use_paired_end = TRUE,
  minimum_overlap = 70,
  normalisation_statistic = "median",
  baf_merge_distance = 1e4,
  baf_max_rows = 4e5, 
  baf_max_merge_iterations = 100,
  enforce_matched_reference = TRUE,
  minimum_hist_bins = 100,
  cluster_type = "FORK",
  cluster_cpus = 6, # how many cluster doParallel asks for. Giving problems on slurm.
  prebuild_references = TRUE,
  store_by_panel = TRUE,
  make_what = "alias,bin,gc_table,reference,bam,count,references_for_bedtype,corrected_count,vcf,baf,deletion,deletions_zscore,translocation,infile",
  clobber_locked_files = "FALSE",
  reference_variables = "TumourNormal=Tumour|Normal;Sex=Male|Female;SampleType=Fresh|FFPE",
  min_reference_samples = 4,
  wg_gc_correct = FALSE,
  targeted_gc_correct = TRUE,
  off_target_gc_correct = FALSE,
  no_clobber = TRUE,
  number_of_cores = 1,
  wg_label = "Whole Genome",
  targeted_label = "Targeted",
  dcf_file = "default.dcf",
  pq_boundary = list (chr1 = 125000000, chr10 =  40200000,  chr11 =  53700000,  chr12 =  35800000,  chr13 =	17900000,  chr14 =	17600000,  chr15 =	19000000,  chr16 =	36600000,  chr17 =	24000000,  chr18 =	17200000,  chr19 =	26500000,  chr2 =	93300000,  chr20 =	27500000,  chr21 =	13200000,  chr22 =	14700000,  chr3 =	91000000,  chr4 =	50400000,  chr5 =	48400000,  chr6 =	61000000,  chr7 =	59900000,  chr8 =	45600000,  chr9 =	49000000,  chrX =	60600000,  chrY =	12500000)
)

# Try and figure out if we are being debugged in Rstudio
# debugSource('~/mm/cnv/buildcnb/makeCNTables.R', echo=TRUE)
CFG_DEFAULTS$rstudio_debug <- (length(args)==0)
if(CFG_DEFAULTS$rstudio_debug){ args <- CFG_DEFAULTS$dcf_file }
CFG_DEFAULTS$is_command_line <- !CFG_DEFAULTS$rstudio_debug && length(args)>0 && args[1]!="RStudio" &&  !grepl("R$",args) 

# Wrong use from command line
if(length(args)!=1 && CFG_DEFAULTS$is_command_line) {
  print("Builds files for cnb. Usage: makeCNTables.R /path/to/config_dcf")
  quit(save = "no", status = 0)
} 

# Load a version of Rsubread - possibly modified
if(CFG_DEFAULTS$use_own_subread==TRUE) {
  library(devtools)
  load_all(file.path(getwd(),paste0("lib/Rsubread_",CFG_DEFAULTS$rsubread_version)))
  .libPaths(c("/pathology/NGS/Gaffa/gaffa/lib/",.libPaths()))
} else {
  library(Rsubread)
}
# Silently load other libraries
libnames <- c(
  "matrixStats",
  "plyr",
  "futile.logger",
  "RCurl",
  "httr",
  "Biostrings",
  "rtracklayer",
  "ggplot2",
  "ggrepel",
  "VariantAnnotation",
  "DNAcopy",
  "Hmisc",
  "modeest",
  "R.utils",
  "corrplot",
  "doParallel",
  "foreach"
)
for (libname in libnames) { 
  suppressPackageStartupMessages(library(libname,character.only=TRUE,warn.conflicts=FALSE,quietly=TRUE))
}

# As part of a shiny app, these functions and globals exist already
if(!exists("gv")) {
  session <- NULL
  gv <- list()
  gv$log <- ""
  gv$broken <- FALSE
  isolate <- function (f) { return (f)}
} # For stand-alone

# Initialise global constants
CFG <- list()
DF_CFG <- data.frame()
SAMPLES <- list()
DF_BEDTYPE2BAMTYPE <- data.frame(bamtypes =c ("wg","wg","wg","targeted","targeted"), 
                           bedtypes = c("wg_coarse","wg_fine","wg_medium","targeted","off_target"),
                           bedfilename = c("wg_coarse_bedfile","wg_fine_bedfile","wg_medium_bedfile","targeted_bedfile","off_target_bedfile"),
                           binsize = c(CFG_DEFAULTS$wg_coarse_bin_size,CFG_DEFAULTS$wg_fine_bin_size,CFG_DEFAULTS$wg_medium_bin_size,CFG_DEFAULTS$targeted_bin_size,CFG_DEFAULTS$targeted_bin_size),
                           stringsAsFactors=FALSE)
DF_FILETYPES <- data.frame(
  filetype = c("alias","bin","gc_table","reference","bam","count","references_for_bedtype","corrected_count","vcf","baf","deletion","deletions_zscore","translocation","infile"),
  per_sample  = c(FALSE,FALSE,FALSE,     FALSE,      TRUE, TRUE,    FALSE,        TRUE,             TRUE,  TRUE,  TRUE,    TRUE,              TRUE,           FALSE),
  per_bedtype = c(FALSE,TRUE, TRUE,      FALSE,      FALSE,TRUE,    TRUE,         TRUE,             FALSE, FALSE, FALSE,   TRUE,              FALSE,          FALSE),
  stringsAsFactors=FALSE)

ANNOTATION_FILES  <- c("wg_annotation_file","targeted_annotation_file","cytobands_file","chrom_info_file")
BEDTYPE2BAMTYPE <- DF_BEDTYPE2BAMTYPE$bamtypes
BEDTYPES <- DF_BEDTYPE2BAMTYPE$bedtypes
names(BEDTYPE2BAMTYPE) <- BEDTYPES
BINSIZES <- DF_BEDTYPE2BAMTYPE$binsize
names(BINSIZES) <- BEDTYPES
BEDFILENAMES <- DF_BEDTYPE2BAMTYPE$bedfilename
names(BEDFILENAMES) <- BEDTYPES
BEDFILENAME2BEDTYPE <- DF_BEDTYPE2BAMTYPE$bedtypes
names(BEDFILENAME2BEDTYPE) <- BEDFILENAMES
# Report package use
df_ip <- (installed.packages())
Rsubread_base_version <- df_ip["Rsubread","Version"]
Rsubread_version <- sessionInfo()$otherPkgs$Rsubread$Version
gv$log <- paste0(isolate(gv$log),flog.info("Found RSubread version: %s Installed base version: %s",Rsubread_version,Rsubread_base_version))

if(CFG_DEFAULTS$require_subread_version==TRUE && CFG_DEFAULTS$rsubread_version!=Rsubread_version) { 
  gv$log <- paste0(isolate(gv$log),flog.error("Found RSubread version: %s but require %s",Rsubread_version,CFG_DEFAULTS$rsubread_version))
  return(2)
}


# -------------------------------------------------------------------
makeCNTables <- function(configFileName) {
  # This makes all the tables that the CN browser needs
  if(!exists("gv")) {
    flog.info("Running makeCNTables() stand-alone with config file: %s",configFileName)
  } # For stand-alone
  df_cfg <- ftry(read.dcf(configFileName))
  if(length(df_cfg)==1) {
    gv$log <<- paste0(isolate(gv$log),flog.error("Unable to read config file: %s",configFileName))
    return(1)
  }
  gv$log <<- paste0(isolate(gv$log),flog.info("Found RSubread version: %s Installed base version: %s",Rsubread_version,Rsubread_base_version))
  gv$log <<- paste0(isolate(gv$log),flog.info("Successfully read config file: %s",configFileName))
  df_ip <- installed.packages()
  df_cfg <- as.data.frame(df_cfg,stringsAsFactors=FALSE)
  rownames(df_cfg) <- df_cfg$name
  number_of_samples <- as.integer(df_cfg["number_of_samples","value"])
  df_cfg <- df_cfg[!(row.names(df_cfg) %in% paste0("sample_",(number_of_samples+1):1000)), ] # remove extras
  idx <-  grepl("bams",df_cfg$section) | ((df_cfg$section=="files") & !grepl("path",row.names(df_cfg)))
  downloadable_files <- row.names(df_cfg)[idx]
  ret <- sapply(downloadable_files,download_files,df_cfg)
  df_cfg  <- massage_paths(df_cfg) # swap out URLs for local paths, adjust other locations as necessary
  df_cfg  <- find_samples(df_cfg)
  #  SAMPLES, CONFIG and DF_CONFIG are const globals that are set here
  DF_CFG <<- df_cfg
  CFG <<-  as.list(df_cfg$value)
  names(CFG) <<- df_cfg$name
  CFG[CFG=="TRUE"] <<- TRUE
  CFG[CFG=="FALSE"] <<- FALSE 
  CFG <<- add_defaults(CFG) # Add entries from CFG_DEFAULTS that are not in CFG
  SAMPLES <<- makesamplelist(df_cfg)
  if(length(SAMPLES)==0) { return (FALSE)}
  logfile <- gsub(".t.[tv]$",".log",CFG$input_file)  
  used_packages_file <- gsub(".t.[tv]$",".used_packages.log",CFG$input_file)  
  # No more modification of globals after this point
  gv$log <<- paste0(isolate(gv$log),"\n",sessionInfo())
  gv$log <<- paste0(isolate(gv$log),flog.info("Writing installed packages list to %s",used_packages_file))
  write.table(installed.packages(),used_packages_file,sep="\t",quote = FALSE)
  #########################################################################################
  # Main loop. Outer loop over file type (and possible resolution), inner loop over samples
  #########################################################################################
  df_table <- NULL
  make_what <- strsplit(CFG$make_what,",")[[1]]
  df_filetypes <- DF_FILETYPES[DF_FILETYPES$filetype %in% make_what,]
  for(j in 1:nrow(df_filetypes)) {
    row <-df_filetypes[j,] 
    if(row$per_sample) { samples <- SAMPLES } else { samples = "" }
    if(row$per_bedtype) { bedtypes <- BEDTYPES[!grepl("off_target",BEDTYPES)] } else { bedtypes = "" }
    if(row$filetype == "infile") { df_arg <- df_table } else { df_arg = "" }
    df_something <- make_something(row$filetype,bedtypes,samples,df_arg)
    if(gv$broken) {
      gv$log <<-  paste0(isolate(gv$log),flog.info("Broke making %s",row$filetype))
      return (1);
    }	
    df_table <- rbind.fill(df_table,df_something)
  }
  write.table(gv$log,logfile,row.names = FALSE,col.names=FALSE,quote=FALSE)
  return(0)
}
# -------------------------------------------------------------------

# -------------------------------------------------------------------
make_infiles <- function(df_cnb_table_entry) {
  no_clobber <- CFG$no_clobber
  input_file <- CFG$input_file
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
  output_df$name <- NULL
  
  if(!file_no_clobber(input_file) && file_apply_lock(input_file)) {
    gv$log <<- paste0(isolate(gv$log),flog.info("Making: %s",input_file))
    ftry(write.table(output_df,input_file,row.names=FALSE,sep="\t",quote = FALSE))
    file_remove_lock(input_file)
  } else {
    gv$log <<- paste0(isolate(gv$log),flog.info("Already have: %s",input_file))
  }
  
  return(output_df)
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
add_defaults  <- function (cfg) {
  for (f in names(CFG_DEFAULTS)) {
    if(is.null(cfg[[f]])) {
      gv$log <<- paste0(isolate(gv$log),flog.info("Adding default config entry %s = %s.",f,CFG_DEFAULTS[f]))
      cfg[f] <- CFG_DEFAULTS[f]
    }
  }
  return(cfg)
}
# -------------------------------------------------------------------

# -------------------------------------------------------------------
massage_paths  <- function (df_cfg) {
  # Pull out URLs
  idx <- (grepl("bams",df_cfg$section) | 
            grepl("fastqs",df_cfg$section) | 
            grepl("files",df_cfg$section) ) & isURL(df_cfg$value)
  df_cfg$value[idx] <- sub("^.*/","",df_cfg$value[idx])
  
  # Put in paths in refs - except for input_file
  idx_add_path <- grepl("files",df_cfg$section) & 
    !grepl("path",row.names(df_cfg)) & 
    !hasPath(df_cfg$value) & 
    !grepl("input_file",row.names(df_cfg))
  # A bit yuck. Anything without "list" and "target" in the name must be for ALL panels
  idx_all_ref <- !grepl("list",row.names(df_cfg)) & !grepl("target",row.names(df_cfg))
  idx <- idx_add_path & idx_all_ref
  if(sum(idx)>0) {
    reference_path <- df_cfg["reference_path","value"]
    df_cfg$value[idx] <- file.path(reference_path,df_cfg$value[idx])
  }
  # Anything with "list" in the name must be for this particular panel
  idx_panel_ref <- grepl("list",row.names(df_cfg)) | grepl("target",row.names(df_cfg))
  idx <- idx_add_path & idx_panel_ref
  # print(df_cfg[idx,]) 
  if(sum(idx)>0) {
    reference_path <- makepath(panel = df_cfg["panel_path","value"] ,run = NULL ,sample = NULL ,filetype = "reference",df_cfg = df_cfg)
    gv$log <<- paste0(isolate(gv$log),flog.info("Adding %s to %s",reference_path,df_cfg$value[idx]))
    df_cfg$value[idx] <- file.path(reference_path,df_cfg$value[idx]) 
  }
  return (df_cfg)  
}
# -------------------------------------------------------------------

# -------------------------------------------------------------------
makepath <- function(panel = NULL,
                     run = NULL,
                     sample = NULL,
                     filetype = NULL,
                     resolution = NULL,
                     extension = NULL,
                     df_cfg,
                     create_dir = TRUE) { 
  if(filetype == "fastq" ) {
    idx_fastqs <- grepl("fastqs",df_cfg[,"section"]) & df_cfg[,"label"] == sample
    if(sum(idx_fastqs)==1) {
      ret <- file.path(df_cfg[idx_fastqs,"value"],paste0(sample,extension))
      return(ret)
    }
  }
  if(filetype == "bam") {
    name_idx <- grepl(sample,df_cfg[,"value"])
    bam_idx <- grepl("bams",df_cfg[,"section"])
    idx <- name_idx & bam_idx
    if(sum(idx)>1 && (!is.null(resolution))) { 
      bam_type <-  BEDTYPE2BAMTYPE[resolution]
      bam_type_idx <- grepl(bam_type,df_cfg[,"section"])
      idx <- idx & bam_type_idx 
    }
    if(sum(idx==1)) { 
      ret <- df_cfg[idx,"value"] 
      return(ret)
    }
  }
  if(!is.null(resolution)) {
    extension <- paste0("_",resolution,extension)
  }
  base_path <- df_cfg["base_path","value"]
  filetype_path <- df_cfg[paste0(filetype,"_path"),"value"]
  per_sample_root_path <- ifelse(sum(names(CFG)=="per_sample_root_path")==1,CFG$per_sample_root_path,"")
  per_sample_root_path <- ifelse(per_sample_root_path==".","",per_sample_root_path)
  if(sum(grepl(filetype,c("cnv","vcf","sv","baf","fastq","bam")))>0) {
    # store by run unless asked to do otherwise
    # optionally put an EXTRA directory under sample to store the data
    ret <- ifelse(CFG$store_by_panel==TRUE,
                  file.path(base_path,panel,run,sample,per_sample_root_path,filetype_path,paste0(sample,extension)),
                  file.path(base_path,run,sample,per_sample_root_path,filetype_path,paste0(sample,extension)))
  } else if (filetype == "reference") {
    reference_path <- df_cfg["reference_path","value"]
    if(is.null(panel)) { ret <- reference_path }
    else if(is.null(run)) { ret <- file.path(reference_path,panel) }
    else if(is.null(sample)) { ret <- file.path(reference_path,panel,run) }
    else { ret <- file.path(reference_path,panel,run,paste0(sample,extension)) }
  }
  # assume that no extension means a directory
  path_to_file <- ifelse(is.null(extension),ret,dirname(ret))
  if(create_dir==TRUE && dir.exists(path_to_file)==FALSE) {
    gv$log <<- paste0(isolate(gv$log),flog.info("Making %s for %s",path_to_file,ret))
    dir.create(path_to_file, recursive = TRUE)
  }
  # gv$log <<- paste0(isolate(gv$log),flog.info("%s+%s+%s+%s+%s=%s",panel,run,sample,per_sample_root_path,filetype,ret))
  return(ret)
}
# -------------------------------------------------------------------


# -------------------------------------------------------------------
find_samples <- function(df_cfg) { 
  if(sum(grepl("regex",df_cfg$section))==0) {
    return (df_cfg)
  }
  gv$log <<- paste0(isolate(gv$log),flog.info("Finding samples via directory listing"))
  
  df_fastqs <- df_cfg[grepl("fastq_regex",df_cfg$section),]
  df_bams <- df_cfg[grepl("bam_regex",df_cfg$section),]
  if(nrow(df_fastqs)>0 || nrow(df_bams)>0) {
    dirlist <- grep(df_cfg["regex_path_match","value"],
                        dir(df_cfg["recursive_dir_path","value"],pattern = df_cfg["recursive_dir_pattern","value"], recursive = TRUE,full.names = TRUE),
                        value=TRUE)
    
    have_fastq <- (nrow(df_fastqs)>0)
    section <- ifelse(have_fastq,"fastqs","bams")
    values <- dirlist
    if(have_fastq) { values <- dirname(values)}
    labels <- sub(df_cfg["regex_gunge_removal","value"],"",basename(dirlist))
    gv$log <<- paste0(isolate(gv$log),flog.info("Looking in ls -R %s/%s",df_cfg["recursive_dir_path","value"],df_cfg["recursive_dir_pattern","value"]))
    gv$log <<- paste0(isolate(gv$log),flog.info("Found %d bam entries",nrow(df_bams)))
    gv$log <<- paste0(isolate(gv$log),flog.info("Making sample list from %s",section))
    gv$log <<- paste0(isolate(gv$log),flog.info("Found %d entries",length(dirlist)))
    if(length(dirlist)==0) {
      gv$log <<- paste0(isolate(gv$log),flog.error("Unable to find any matching samples."))
      return (df_cfg)
    }
    df_cfg_rows <- data.frame(section = section,
                            label = labels,
                            value = values,
                            stringsAsFactors = FALSE)
    df_cfg_rows <- df_cfg_rows[!duplicated(labels),]
    df_cfg_rows$name = paste0("sample_",1:nrow(df_cfg_rows))
  }    
  rownames(df_cfg_rows)  <- df_cfg_rows$name
  number_of_samples <- as.integer(df_cfg["number_of_samples","value"])
  if(number_of_samples != nrow(df_cfg_rows)) {
    gv$log <<- paste0(isolate(gv$log),flog.info("After duplicate removal: %d entries",nrow(df_cfg_rows)))
    gv$log <<- paste0(isolate(gv$log),flog.warn("This differs from the number_of_samples entry (%d)",number_of_samples))
  }
  return(rbind(df_cfg,df_cfg_rows))
}

# -------------------------------------------------------------------
makesamplelist <- function(df_cfg) { 
  gv$log <<- paste0(isolate(gv$log),flog.info("Making sample list"))
  ret <- NULL
  df_samples <- df_cfg[grepl("samples",df_cfg$section),]
  df_bams <- df_cfg[grepl("bams",df_cfg$section),]
  df_fastqs <- df_cfg[grepl("fastqs",df_cfg$section),]
  ret <- unique(c(df_samples$label,df_bams$label,df_fastqs$label))
  num_samples <- length(ret)
  gv$log <<- paste0(isolate(gv$log),flog.info("Found %d samples.",num_samples))
  return(ret)
}
# -------------------------------------------------------------------  
  
# -------------------------------------------------------------------
download_files <- function(svalue,df_cfg) {
  s <- df_cfg[svalue,"value"]
  label <- df_cfg[svalue,"label"]
  dest_path <- CFG$reference_path
  if(grepl("bams",df_cfg[svalue,"section"])) { 
    dest_path <- makepath(panel = df_cfg["panel_path","value"] ,run = df_cfg["run_path","value"] ,sample = label ,filetype = "bam",df_cfg = df_cfg)
  } else if(grepl("fastqs",df_cfg[svalue,"section"])) {
    dest_path <- makepath(panel = df_cfg["panel_path","value"] ,run = df_cfg["run_path","value"] ,sample = label ,filetype = "fastq",df_cfg = df_cfg)
  }
  if(isURL(s)) {
    url <- parse_url(s)
    filename <- basename(url$path)
    destfile <- file.path(dest_path,filename)
    if(!file_no_clobber(destfile) && file_apply_lock(destfile)) {
      gv$log <<- paste0(isolate(gv$log),flog.info("Downloading %s",destfile))          
      ret <- ftry(download.file(s,destfile))
      if(ret==0) {
        gv$log <<- paste0(isolate(gv$log),flog.info("Successfully downloaded %s",destfile))          
      } else {
        gv$log <<- paste0(isolate(gv$log),flog.error("Unable to download file: %s from %s",destfile,s))
      }
      file_remove_lock(destfile)
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
file_no_clobber  <- function(filename) {
  if(!file.exists(filename)) {
    gv$log <<- paste0(isolate(gv$log),flog.info("File does not exist: %s",filename))
    return(FALSE)
  }
  if (!file_is_locked(filename) && CFG$no_clobber==FALSE) {
    gv$log <<- paste0(isolate(gv$log),flog.info("Clobbering: %s",filename))
    return(FALSE)
  }
  if (file_is_locked(filename) && CFG$clobber_locked_files) {
    gv$log <<- paste0(isolate(gv$log),flog.info("Removing previous version of locked file: %s",filename))
    unlink(filename)
    file_remove_lock(filename)
    return(FALSE)
  }
  gv$log <<- paste0(isolate(gv$log),flog.info("File already exists: %s",filename))
  return(TRUE)
}
# -------------------------------------------------------------------

# -------------------------------------------------------------------
file_apply_lock  <- function(filename) {
  lockfile <- paste0(filename,".lock")
  gv$log <<- paste0(isolate(gv$log),flog.info("Locking file: %s",filename))
  ret <- tryCatch({
    if(!file.exists(lockfile)) {
      write(flog.info("Created lock: %s",filename),file = lockfile)
    } else {
      gv$log <<- paste0(isolate(gv$log),flog.warn("Already locked: %s",filename))
      write(flog.info("Relocked"),file = lockfile,append = TRUE)
    }
    return(TRUE)
  }, error = function(err) {
    gv$log <<- paste0(isolate(gv$log),flog.error("Unable to write lock: %s",filename))
    gv$broken <<- TRUE
    return(FALSE)
  })
  return(ret)
}
# -------------------------------------------------------------------
# -------------------------------------------------------------------
file_remove_lock  <- function(filename) {
  lockfile <- paste0(filename,".lock")
  gv$log <<- paste0(isolate(gv$log),flog.info("Unlocking file: %s",filename))
  if(file.exists(lockfile)) {
    unlink(lockfile)
  } else {
    gv$log <<- paste0(isolate(gv$log),flog.warn("No lock for file: %s",filename))
  }
}
# -------------------------------------------------------------------

# -------------------------------------------------------------------
file_is_locked  <- function(filename) { return(file.exists(paste0(filename,".lock")))}
# -------------------------------------------------------------------

# -------------------------------------------------------------------
make_something <- function(filetype,bedtypes = "",samples = "", df_arg = NULL) {
  gv$log <<- paste0(isolate(gv$log),flog.info("------------------------------------------------------------"))
  gv$log <<- paste0(isolate(gv$log),flog.info("make_%ss : n_bedtypes = %d n_samples = %d",filetype,length(bedtypes),length(samples)))
  gv$log <<- paste0(isolate(gv$log),flog.info("------------------------------------------------------------"))
  remaining_filetypes <- paste0(filetype,",",sub(paste0("^.*",filetype,","),"",paste(DF_FILETYPES$filetype,collapse=',')))
  gv$log <<- paste0(isolate(gv$log),flog.info("Recover from here with make_what=%s and clobber_locked_files=TRUE",remaining_filetypes))
  
  df <- NULL
  for (sample in samples) {
    for (bedtype in bedtypes) {
      args  <- list()
      if(sample != "") { args$sample = sample }
      if(bedtype != "") { args$bedtype = bedtype }
      if(is.data.frame(df_arg)) { args$df = df_arg }
      gv$log <<- paste0(isolate(gv$log),flog.info("sample = %s bedtype = %s",sample,bedtype))
      df_something <- do.call(paste0("make_",filetype,"s"),args) 
      df <- rbind.fill(df,df_something)
    }
  }
  return(df)
}
# -------------------------------------------------------------------


# -------------------------------------------------------------------
make_aliass <- function() {
  chrAliasesFile <- file.path(CFG$chr_alias)
  gv$log <<- paste0(isolate(gv$log),flog.info("Making: %s",chrAliasesFile))
  if(!file_no_clobber(chrAliasesFile) && file_apply_lock(chrAliasesFile)) {
    nochr <- c(paste0("",c(1:22)),"X","Y")
    chr<- c(paste0("chr",c(1:22)),"chrX","chrY")
    chrAliases <- data.frame(x=c(nochr,chr), X=c(chr,nochr))
    ret <- ftry(write.csv(chrAliases,chrAliasesFile,row.names=FALSE,quote = FALSE))
    file_remove_lock(chrAliasesFile)
  }
  gv$log <<- paste0(isolate(gv$log),flog.info("Have: %s",chrAliasesFile))
  return(NULL)
}
# -------------------------------------------------------------------

# -------------------------------------------------------------------
make_counts  <-  function(sample,bedtype) {
  bedfile <- BEDFILENAMES[bedtype]
  gv$log <<- paste0(isolate(gv$log),flog.info("Making counts for %d samples with %s and %s",length(SAMPLES),bedfile,bedtype))
  write_corrected_counts(sample,bedfile,"counts")
  return(NULL)
}
# -------------------------------------------------------------------
# -------------------------------------------------------------------
make_corrected_counts  <-  function(sample,bedtype) {
  df_cnb <- NULL
  bedfile <- BEDFILENAMES[bedtype]
  gv$log <<- paste0(isolate(gv$log),flog.info("Making counts for %d samples with %s and %s",length(SAMPLES),bedfile,bedtype))
  df_cnb <-  write_corrected_counts(sample,bedfile,"corrected_counts")
  gv$log <<- paste0(isolate(gv$log),flog.info("Processed %d samples.",ifelse(is.null(df_cnb),0,nrow(df_cnb))))
  return(df_cnb)
}
# -------------------------------------------------------------------

# -------------------------------------------------------------------
write_corrected_counts <- function(sample,bedfile,do_what) {
  bedtype <- BEDFILENAME2BEDTYPE[bedfile] # targeted, off_target, wg_coarse, wg_medium, wg_fine
  bam_type <-  BEDTYPE2BAMTYPE[bedtype] # wg, targeted, off_target
  number_of_cores <- as.integer(CFG$number_of_cores)
  bampath <- makepath(panel = CFG$panel_path ,run = CFG$run_path ,sample = sample ,filetype = "bam",extension = ".bam",df_cfg = DF_CFG)
  bamfile <- basename(bampath)
  correct_gc <- as.logical(DF_CFG[paste0(bam_type,"_gc_correct"),"value"])
  outpath <- makepath(panel = CFG$panel,run = CFG$run_path,sample = sample,filetype = "cnv",resolution = bedtype, extension = ".tsv",df_cfg = DF_CFG)
  cnb_info_path <- sub(".tsv$","_info.tsv",outpath)
  countspath <- sub(".tsv$","_fc_counts.tsv",outpath)
  annotationpath <- sub(".tsv$","_fc_annotation.tsv",outpath)
  bedpath <- CFG[[bedfile]]
  chrAliasesFile <- CFG$chr_alias
  
  # Write entry for input file for cnb
  df_cnb <- data.frame(
    resolution = bedtype,
    path = outpath,
    file = "counts",
    name = sample,
    is_wgs = ifelse(bam_type=="wg",TRUE,FALSE),
    is_targeted = ifelse(bam_type=="targeted",TRUE,FALSE),
    run = CFG$run_path,
    panel = CFG$panel_path,
    sample = sample,
    sex = "unknown",
    tumour_normal = "unknown",
    sample_type = "unknown",
    stringsAsFactors = FALSE)
  
  if(do_what=="counts") {
    df_reference_file_list <- NULL
  } else { # For corrected counts, attempt to find out sex etc
    df_reference_file_list <- read.table(CFG$reference_file_list,header = TRUE,stringsAsFactors = FALSE,sep = "\t")
    sample_idx <- df_reference_file_list$SampleName==sample & df_reference_file_list$RunName==CFG$run_path
    if(sum(sample_idx)==1) { 
      gv$log <<- paste0(isolate(gv$log),flog.info("Found entry for %s in %s.",sample,CFG$reference_file_list))
      if(sum(grepl("Sex",colnames(df_reference_file_list)))>0) {
        df_cnb$sex <- df_reference_file_list$Sex[sample_idx][1]
      }
      if(sum(grepl("Tumour",colnames(df_reference_file_list)))>0) {
        df_cnb$tumour_normal <- df_reference_file_list$TumourNormal[sample_idx][1]
      }
      if(sum(grepl("SampleType",colnames(df_reference_file_list)))>0) {
        df_cnb$sample_type <- df_reference_file_list$SampleType[sample_idx][1]
      }
    }
  }   

  if(do_what=="counts" && (!file_no_clobber(annotationpath)) && file_apply_lock(annotationpath)) {    
    gv$log <<- paste0(isolate(gv$log),flog.info("Making %s",annotationpath))
    df_bedfile <- read.table(bedpath,header=FALSE,stringsAsFactors=FALSE)
    colnames(df_bedfile) <- c('chr','start','end','id','score')
    df_id <- df_bedfile[,4:5]
    gr_bedfile <- import.bed(bedpath)
    mcols(gr_bedfile) <- df_id
    annotation_file <- createAnnotationFile(gr_bedfile)
    write.table(annotation_file,annotationpath,row.names=FALSE,sep="\t",quote = FALSE)
    file_remove_lock(annotationpath)  
  }
  

  if(do_what=="counts" && (!file_no_clobber(countspath)) && file_apply_lock(countspath)) {    
    gv$log <<- paste0(isolate(gv$log),flog.info("Making counts for: %s / %s. In %s",bedfile,sample,countspath))
    gv$log <<- paste0(isolate(gv$log),flog.info("This bit takes a long time - of order 1 min/10^6 reads/bedfile or 1 hour/5GB/bedfile."))
    # system("vmstat -S M -s")
    # change into destination directory to avoid temp file collisions
    cwd <- getwd()
    chrAliasesFile <- ifelse(isAbsolutePath(chrAliasesFile),chrAliasesFile,file.path(cwd,chrAliasesFile))
    annotationpath <- ifelse(isAbsolutePath(annotationpath),annotationpath,file.path(cwd,annotationpath))
    bampath <- ifelse(isAbsolutePath(bampath),bampath,file.path(cwd,bampath))
    setwd(dirname(bampath))
    if(Rsubread_version=="1.16.1") {
      fc <- featureCounts(files=bamfile,
                          isPairedEnd=CFG$use_paired_end,
                          annot.ext=annotationpath,
                          nthreads = number_of_cores,
                          useMetaFeatures = TRUE,
                          ignoreDup = FALSE, # may or may not be what you want
                          countMultiMappingReads = TRUE,
                          minReadOverlap = CFG$minimum_overlap,
                          requireBothEndsMapped = FALSE,
                          allowMultiOverlap = TRUE,
                          fraction = TRUE,
                          chrAliases = chrAliasesFile)
    }
    else {
      fc <- featureCounts(files=bamfile,
                          isPairedEnd=CFG$use_paired_end,
                          annot.ext=annotationpath,
                          nthreads = number_of_cores,
                          useMetaFeatures = TRUE,
                          largestOverlap = TRUE,
                          minOverlap = 70,
                          fraction = TRUE,     
                          ignoreDup = FALSE, # may or may not be what you want
                          countMultiMappingReads = TRUE,
                          requireBothEndsMapped = FALSE,
                          allowMultiOverlap = TRUE,
                          checkFragLength=FALSE,
                          chrAliases = chrAliasesFile)
    }
    gv$log <<- paste0(isolate(gv$log),flog.info("Removing temporary files from %s",dirname(bampath)))
    lapply(c("fixbam","sam","tmp","bin"), remove_files_by_suffix)
    setwd(cwd)
    save_list_elements <- function(s) { write.table(fc[[s]],sub(".tsv$",paste0("_fc_",s,".tsv"),outpath),quote = FALSE,sep="\t") }
    lapply(names(fc), save_list_elements) # save everything that featureCounts() returns
    # remove temporary files
    lapply(c("fixbam","sam","tmp","bin"), remove_files_by_suffix)
    file_remove_lock(countspath)  
  }
  
  if(do_what=="counts" || (do_what=="corrected_counts" && file_no_clobber(outpath))) {
    gv$log <<- paste0(isolate(gv$log),flog.info("Skipping GC correction, segmentation, reference correction for %s / %s",bedfile,sample))
    if(file.exists(cnb_info_path)) {
      df_cnb_saved <- read.table(cnb_info_path,header=TRUE,stringsAsFactors=FALSE)
      df_cnb_saved[colnames(df_cnb)] <- df_cnb
      df_cnb <- df_cnb_saved
    } else {
      gv$log <<- paste0(isolate(gv$log),flog.warn("Unable to find %s but pressing on regardless.",cnb_info_path))
    }
    return(df_cnb)
  } else if(do_what=="corrected_counts"){
    file_apply_lock(outpath)
    fc <- list()
    fc$counts <- read.table(countspath,header=TRUE,stringsAsFactors=FALSE)
    fc$annotation <- read.table(annotationpath,header=TRUE,stringsAsFactors=FALSE)
    gv$log <<- paste0(isolate(gv$log),flog.info("Making GC correction for: %s",outpath))
    df_bedfile <- read.table(bedpath,header=FALSE,stringsAsFactors=FALSE)
    colnames(df_bedfile) <- c('chr','start','end','id','score')
    if(nrow(df_bedfile)>CFG$loess_scales_badly_limit) {
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
    idx_gc_fit <- df_gc$counts>CFG$loess_min_usable_counts & df_gc$gc>CFG$loess_min_usable_gc
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
      loess(counts ~ gc, data = df_gc_fit,  span = CFG$loess_span),
      error=function(cond) { 
        message("Loess failed:") 
        gv$log <<- paste0(isolate(gv$log),flog.warn("Loess failed: %s",cond))
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
          gv$log <<- paste0(isolate(gv$log),flog.warn("Loess predict failed: %s",cond))
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
    references <- readAllReferences(bedtype)
    if(length(references)==0) {
      gv$log <<- paste0(isolate(gv$log),flog.error("Unable to find references."))
      gv$broken <<- TRUE 
      return (NULL)
    }
    if(CFG$enforce_matched_reference == TRUE) {
      sex_idx <- (grepl("Female",names(references)) & grepl("Female",df_cnb$sex)) | (grepl("Male",names(references)) & grepl("Male",df_cnb$sex))
      tumour_normal_idx <- (grepl("Tumour",names(references)) & grepl("Tumour",df_cnb$tumour_normal)) | (grepl("Normal",names(references)) & grepl("Normal",df_cnb$tumour_normal))
      sample_type_idx <- (grepl("Fresh",names(references)) & grepl("Fresh",df_cnb$sample_type)) | (grepl("FFPE",names(references)) & grepl("FFPE",df_cnb$sample_type))
      idx <- sex_idx & tumour_normal_idx & sample_type_idx
      if(sum(idx)==0) {
        idx <- sex_idx & tumour_normal_idx
      }
      if(sum(idx)==0) {
        idx <- sex_idx
      }
      if(sum(idx)!=0) {
        gv$log <<- paste0(isolate(gv$log),flog.info("Enforcing a choice of %d from %d references",sum(idx),length(references)))
        references <- references[idx]
      }
    }
    bestReference <- findBestReference(references,df)
    gv$log <<- paste0(isolate(gv$log),flog.info("Writing GC corrected counts: %s",outpath))
    ftry(write.table(bestReference$corrected,outpath,row.names=FALSE,sep="\t",quote = FALSE))
    
    # Write out the loess correction across the range of GC values
    bins <- seq(0, 1, by = CFG$loess_bin_width)
    fapply_median <- function(n,bins) { median(df_gc_fit$counts[df_gc_fit$gc>=(bins[n] - CFG$loess_bin_width) & df_gc_fit$gc<(bins[n] + CFG$loess_bin_width)]) }
    bin_medians <- unlist(lapply(1:length(bins),fapply_median,bins))
    loessCorrection_df <- data.frame(gc = bins, median = bin_medians)
    loessCorrection_df$bam <-  sample
    loessCorrection_df$bed <- bedfile
    
    if(is.null(loess_predict)) {
      loessCorrection_df$counts <- bin_medians
    } else {
      loessCorrection_df$counts <- predict(loess_predict,loessCorrection_df$gc)
    }
    loesspath <- sub(".tsv$",paste0("_fc_loess_corrections.tsv"),outpath)
    write.table(loessCorrection_df,loesspath,row.names=FALSE,quote = FALSE,sep="\t")
    
    # Finish entry for input file for cnb
    df_reference_file_list <- read.table(CFG$reference_file_list,header = TRUE,stringsAsFactors = FALSE,sep = "\t")
    if(sum(grepl("Sex",colnames(df_reference_file_list)))==0 && df_cnb$sex!="Male" && df_cnb$sex!="Female") {
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
    file_remove_lock(outpath) 
    return(df_cnb)
  }
}
# -------------------------------------------------------------------

# -------------------------------------------------------------------
make_references = function () {
  # make reference for aligner - once only
  chrom_file <- CFG$chrom_file
  chrom_index_file <- CFG$chrom_index_file
  number_of_cores <- CFG$number_of_cores
  num_index_files <- length(list.files(path = dirname(chrom_index_file), pattern = paste0(basename(chrom_index_file),".*")))
  if(num_index_files == 0) {
    gv$log <<- paste0(isolate(gv$log),flog.info("Making index for %s",chrom_file))
    buildindex(chrom_index_file,chrom_file)
  } else {
    gv$log <<- paste0(isolate(gv$log),flog.info("Doing alignment: index file = %s cores = %s",CFG$chrom_index_file,CFG$number_of_cores))
    gv$log <<- paste0(isolate(gv$log),flog.info("Already have: %d files matching %s.*",num_index_files,chrom_index_file))
  }
  return(NULL)    
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
make_references_for_bedtypes = function (bedtype) {
  # make resolution-specific references for others
  bam_type <- BEDTYPE2BAMTYPE[bedtype]
  reference_file_list_file_name <- CFG$reference_file_list
  gc_path <- bedrow2gcrow(BEDFILENAMES[bedtype])$value
  output_path <- makepath(panel = CFG$panel_path,run = CFG$run_path,filetype = "reference",df_cfg = DF_CFG)
  df_gc <- read.table(gc_path,stringsAsFactors=FALSE,header=TRUE)
  df_reference_file_list <- read.table(reference_file_list_file_name,header = TRUE,stringsAsFactors = FALSE, sep = "\t")
  counts_ending <- paste0("_",bedtype,"_fc_counts.tsv")
  df_reference_file_list$CountsFile <- paste0(df_reference_file_list$SampleName,counts_ending)
  for(j in 1:nrow(df_reference_file_list)) {
    df_reference_file_list$CountsPath[j] <- 
      makepath(panel = CFG$panel_path,
               run = df_reference_file_list$RunName[j],
               sample = df_reference_file_list$SampleName[j],
               filetype = "cnv",
               resolution = bedtype,
               extension = "_fc_counts.tsv"
               ,df_cfg = DF_CFG)
  }
    
  num_samples <- nrow(df_reference_file_list)
  min_reference_samples <- CFG$min_reference_samples
  gv$log <<- paste0(isolate(gv$log),flog.info("Successfully read reference_file_list file: %s with %d samples",
                                             reference_file_list_file_name,num_samples))
  if(grepl(CFG$del_zscore_section,bedtype)==TRUE) {
    make_deletions_stats(TRUE,df_reference_file_list$CountsPath,bedtype)
  }
  
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
  df_files$resolution <- bedtype
  df_files$is_wgs <- (bam_type=="wg")
  df_files$is_targeted <- !df_files$is_wgs
  
  for (i in 1:numCombinations) {
    reference_file_name <- ""
    criteria_rows <- df_all_combinations[i,]
    idx <- !logical(num_samples)
    for (j in criteria_rows) {
      df_row <- df_criteria_col[j,]
      if(df_row$Val!="All") {
        if(sum(df_row$SampleColumn==colnames(df_reference_file_list))!=1) {
          gv$log <<- paste0(isolate(gv$log),flog.error("Unable to find the column %s in %s",df_row$SampleColumn,reference_file_list_file_name))
        }
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
    reference_file_name <- paste0("pre_built_ref_",reference_file_name,"_",bedtype)
    df_files[i,"name"] <- reference_file_name
    reference_file_name <- paste0(reference_file_name,".tsv")
    reference_file_path <- file.path(output_path,reference_file_name)
    idx[is.na(idx)] <- FALSE    
    if(sum(idx) >= as.integer(min_reference_samples)) {
      gv$log <<- paste0(isolate(gv$log),flog.info("Need reference %s with %d samples",reference_file_name,sum(idx)))
      if(!file_no_clobber(reference_file_path) && file_apply_lock(reference_file_path)) {
        gv$log <<- paste0(isolate(gv$log),flog.info("Making: %s",reference_file_path))
        makeReference(df_reference_file_list[idx,],reference_file_path,df_gc)
        file_remove_lock(reference_file_path)
      } 
      ret <- append(ret,"reference_file_name")
      df_files[i,"path"] <- reference_file_path
    }
  }
  if(DEBUG==TRUE) {
    return(df_files[!is.na(df_files$path) & !is.null(df_files$path),])
  } else {
    return(NULL)
  }
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
    gv$log <<- paste0(isolate(gv$log),flog.warn("Unable to make reference file due to lack of counts files: %s",reference_file_name))
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
                                          gv$broken <<- TRUE
                                          return(NULL)                                             
                                         }
                                       }
                                       ))
  # Select the columns with the counts in them
  regex_cols <- paste0(df_reference_file_list$SampleName,sep='',collapse='|')
  regex_cols <- gsub("-",".",regex_cols)  # match the sample names that may have had funny things stripped out
  col_idx <- grep(regex_cols,colnames(df_counts))
  #print(sum(col_idx))
  #print(as.integer(CFG$min_reference_samples))
  #print(colnames(df_counts))
  #print(regex_cols)
  if(sum(col_idx) < as.integer(CFG$min_reference_samples)) {
    gv$log <<- paste0(isolate(gv$log),flog.info("Need %d samples, but only have %d.",as.integer(CFG$min_reference_samples),sum(col_idx)))
    return(NULL)
  }  
  df_counts <- df_counts[,col_idx]
  gv$log <<- paste0(isolate(gv$log),flog.info("Calculating references"))
  if(nrow(df_gc) != nrow(df_counts)) {
    gv$log <<- paste0(isolate(gv$log),flog.error("GC file differs in rows to counts table: %d vs %d",nrow(df_gc),nrow(df_counts)))
    return(NULL)
  }

  if(sum(colSums(df_counts)==0)>0) {
    gv$log <<- paste0(isolate(gv$log),flog.error("%d of the counts table cols is zero:",sum(colSums(df_counts)>0)))
    print(colSums(df_counts)==0)
    return(NULL)
  }
  if(nrow(df_reference_file_list) != ncol(df_counts)) {
    gv$log <<- paste0(isolate(gv$log),flog.error("Reference list differs in rows to counts table cols: %d vs %d",nrow(df_reference_file_list),ncol(df_counts)))
    print(df_reference_file_list$SampleName)
    print(colnames(df_counts))
    print(paste0(df_reference_file_list$SampleName,sep='',collapse='|'))
    return(NULL)
  }

  # Align the rows of counts file with gc file
  idx_chr_start_end <- paste0(df_gc$chr,":",as.integer((df_gc$start-1)),"-",df_gc$end)
  df_counts <- df_counts[idx_chr_start_end,]
  # For a reference we want everything to be 2N
  # For Men X -> 2X, Y -> 2Y, For women set Y=0 
  x_idx <- grepl("X",df_gc$chr)
  y_idx <- grepl("Y",df_gc$chr)
  auto_idx <- (!(x_idx | y_idx))
  threshold_idx <- (rowSums(df_counts >= CFG$norm_cutoff & (!is.na(df_counts)))==ncol(df_counts))
  have_sex <- FALSE
  if(sum(threshold_idx)==0) {
    gv$log <<- paste0(isolate(gv$log),flog.info("No bin coud be found with enough counts across all samples."))
    gv$log <<- paste0(isolate(gv$log),flog.info("Printing number of samples with no counts for each bin:"))
    print(rowSums(df_counts==0,na.rm=TRUE))
    gv$log <<- paste0(isolate(gv$log),flog.info("Number of good columns by row:"))
    print(rowSums(df_counts >= CFG$norm_cutoff & (!is.na(df_counts))))
    return(NULL)
  }
  if(sum("Sex"==colnames(df_reference_file_list))==1) {
    female_idx <- grepl("Female",df_reference_file_list$Sex)
    male_idx <- grepl("Male",df_reference_file_list$Sex)
    df_reference_file_list$Sex[female_idx] <- "Female"
    gv$log <<- paste0(isolate(gv$log),flog.info("From table: Males:Females = %d:%d",sum(!female_idx),sum(female_idx)))
    have_sex <- TRUE
  } 
  gv$log <<- paste0(isolate(gv$log),flog.info("Guessing Sex using X_counts/Autosome_counts"))
  auto_medians <- colMedians(as.matrix(df_counts[auto_idx & threshold_idx,]),na.rm = TRUE)
  if(sum(auto_medians==0)>0) {gv$log <<- paste0(isolate(gv$log),flog.warn("Autosome counts are zero for NA."))}
  auto_medians <- ifelse(is.na(auto_medians) | auto_medians==0,1,auto_medians)
  if(sum(x_idx & threshold_idx)>0) {
    x_medians <- colMedians(as.matrix(df_counts[x_idx & threshold_idx,]),na.rm = TRUE)
  } else {
    x_medians <- 0
  }
  if(sum(y_idx & threshold_idx)>0) {
    y_medians <- colMedians(as.matrix(df_counts[y_idx & threshold_idx,]),na.rm = TRUE)
  } else {
    y_medians <- 0    
  }
  x_ratio <- x_medians / auto_medians
  y_ratio <- y_medians / auto_medians
  
  male_idx_counts <- (y_ratio > CFG$maximum_female_y_on_auto_ratio) | (x_ratio < CFG$minimum_female_x_on_auto_ratio)
  female_idx_counts <- !male_idx
  if(!have_sex) { 
    gv$log <<- paste0(isolate(gv$log),flog.warn("Setting sex from thresholding median Y/autosome coverage."))
    df_reference_file_list$Sex <- ifelse(female_idx_counts,"Female","Male") 
  }
  gv$log <<- paste0(isolate(gv$log),
                    flog.info("From counts: Males:Females:NA = %d:%d:%d",
                              sum(!female_idx_counts,na.rm=TRUE),sum(female_idx_counts,na.rm=TRUE),sum(is.na(female_idx_counts))
                              )
                    )
  
  print(data.frame(Sample = df_reference_file_list$SampleName,
                   X = x_medians,
                   Y = y_medians, 
                   AUTO = auto_medians, 
                   X_ON_AUTO = x_ratio,
                   Y_ON_AUTO = y_ratio,
                   COUNTS = ifelse(female_idx_counts,"Female","Male"),
                   TABLE = df_reference_file_list$Sex
  ))
  
  if(length(x_ratio)<=CFG$minimum_samples_for_sex_determination_by_clustering) {
      gv$log <<- paste0(isolate(gv$log),flog.info("Guessing sex from clustering."))
    sex <- kmeans(x_ratio,2)
    female_idx_from_clustering <- sex$cluster==2
    if(mean(x_ratio[female_idx_from_clustering],na.rm=TRUE)) { 
      gv$log <<- paste0(isolate(gv$log),flog.info("Cluster sex-labels swapped"))
      female_idx_from_clustering <- !female_idx_from_clustering 
    } 
    if(sum(female_idx_from_clustering) > 0 &&  mean(y_ratio[female_idx_from_clustering],na.rm=TRUE) > CFG$maximum_female_y_on_auto_ratio) {
      gv$log <<- paste0(isolate(gv$log),flog.info("Both clusters same sex?"))
      # female_idx[female_idx] <- FALSE 
    } 
    if(sum(female_idx_from_clustering) > 0 && mean(x_ratio[female_idx_from_clustering],na.rm=TRUE) < CFG$minimum_female_x_on_auto_ratio) {
      gv$log <<- paste0(isolate(gv$log),flog.info("Both clusters same sex?"))
      # female_idx[!female_idx] <- TRUE
    } 
    gv$log <<- paste0(isolate(gv$log),flog.info("From clustering: Males:Females = %d:%d",sum(!female_idx_from_clustering),sum(female_idx_from_clustering)))
    # if(!have_sex) { df_reference_file_list$Sex <- ifelse(female_idx,"Female","Male") }
  }
  if(sum(x_idx)>0 && sum(y_idx)>0) {
    if(sum(!female_idx)>0) {
      gv$log <<- paste0(isolate(gv$log),flog.info("Adjusting X and Y counts for %d samples",sum(!female_idx)))
      df_counts[x_idx,!female_idx] <- df_counts[x_idx,!female_idx] * 2
      df_counts[y_idx,!female_idx] <- df_counts[y_idx,!female_idx] * 2
    }
    if(sum(female_idx)>0) {
      gv$log <<- paste0(isolate(gv$log),flog.info("Setting Y counts = NA for %d samples",sum(female_idx)))
      df_counts[y_idx,female_idx] <- NA
    }
  }
  
  # Normalise the counts with the lower end of the dist truncated
  df_counts_with_na <- df_counts
  df_counts_with_na[df_counts < CFG$norm_cutoff] <- NA
  
  sampleMedians <- colMedians(as.matrix(df_counts_with_na),na.rm=TRUE)
  ensembleMedian <- median(sampleMedians)
  # Counts median normalised 
  nc <- as.matrix(df_counts) %*% diag(ensembleMedian/sampleMedians) 
  nc[df_counts < CFG$norm_cutoff] <- NA
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
  return(df_counts)
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
readAllReferences <- function (bedtype) {
  gv$log <<- paste0(isolate(gv$log),flog.info("Reading references."))
  output_path <- makepath(panel = CFG$panel_path,run = CFG$run_path,filetype = "reference",df_cfg = DF_CFG)
  pattern <- paste0("^pre_built_ref.*",bedtype,".tsv$")
  file_list <- list.files(output_path,pattern)
  ref_names <- gsub("pre_built_ref_","",(gsub(paste0("_",bedtype,".tsv"),"",file_list)))
  file_paths <- file.path(output_path,file_list)
  refs <- lapply(file_paths,
    function (filename) {
      gv$log <<- paste0(isolate(gv$log),flog.info("Reading counts from %s.",filename))
      read.table(filename,header=TRUE,stringsAsFactors=FALSE)
    })
  names(refs) <- ref_names
  gv$log <<- paste0(isolate(gv$log),flog.info("Read %d references matching %s/%s",length(ref_names),output_path,pattern))
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
  # plot(-log10(sm),-log10(cn))
  # plot(-log10(sm),seg)
  # plot(-log10(cn),seg)
  df_best <- comp[[which.min(sm)]]
  df_best$refname <- refNames[which.min(sm)]
  # df_best$corrected 
  # df_best$comparisons$smoothness_diff 
  # df_best$comparisons$quantised_cn_diff 
  # df_best$comparisons$segments
  gv$log <<- paste0(isolate(gv$log),flog.info("Best reference is: %s",df_best$refname))
  return(df_best)
}
# -------------------------------------------------------------------

# -------------------------------------------------------------------
applyReference  <- function (ref_name,refs,df) {
  gv$log <<- paste0(isolate(gv$log),flog.info("Doing applyReference(%s).",ref_name))
  df_reference <- refs[[ref_name]]
  
  if(nrow(df) != nrow(df_reference)) {
    gv$log <<- paste0(isolate(gv$log),flog.warn("Reference and sample have different numbers of counts. Sample  = %d, Reference = %d",
                                                nrow(df),nrow(df_reference)))
    # gv$broken <<- TRUE
    # return(NULL)
  } 
  
  idx <- (!is.na(df$raw)) & (df$raw>=CFG$median_norm_cutoff) & (!df$blackListed)
  if(sum(df$raw==0) > length(idx) * CFG$max_proportion_zero_entries) {
    gv$log <<- paste0(isolate(gv$log),flog.warn("Failed to apply reference properly or bad sample? Too many counts entries are zero. Entries = %d, Zero Entries = %d, Usable Entries = %d",
                                                length(idx),sum(idx),sum(df$raw==0)))
    # gv$broken <<- TRUE
    # return(NULL)
  } 
  
  normalise_on_mode  <- ifelse(CFG$normalisation_statistic=="mode",TRUE,FALSE)
  
  # Median normalise truncated raw counts
  idx <- (!is.na(df$raw)) & (df$raw>=CFG$median_norm_cutoff) & (!df$blackListed)
  if(sum(df$raw==0) > length(idx) * CFG$max_proportion_zero_entries) {
    gv$log <<- paste0(isolate(gv$log),flog.warn("Failed to apply reference properly or bad sample? Too many counts entries are zero. Something is wrong. Entries = %d, Zero Entries = %d, Usable Entries = %d",
         length(idx),sum(idx),sum(df$raw==0)))
    # gv$broken <<- TRUE
    # return(NULL)
  } 
  nc  <-  df$raw / median(df$raw[idx])
  idx <- !is.na(df_reference$ref_median) & df_reference$ref_median>=CFG$median_norm_cutoff & !df$blackListed 
  nref  <-  df_reference$ref_median / median(df_reference$ref_median[idx])
  ratio <- nc / nref
  idx <- !is.nan(ratio) & is.finite(ratio) & df$raw>=CFG$counts_cutoff & df_reference$raw>=CFG$counts_cutoff & !df$blackListed 
  cor <- ratio
  cor[!idx] <- NA
  # Strip out divide by zero stuff for evaluating metrics
  ratio <- ratio[idx]
  nmax <- length(ratio)
  fom <- list()
  if(sum(idx) < CFG$median_data_points_cutoff) {
    gv$log <<- paste0(isolate(gv$log),flog.error("Unable to normalise. Only have %d points out of %d. Need at least %d to calculate median.",sum(idx),length(idx),CFG$median_data_points_cutoff))
    gv$broken <<- TRUE
    return(df)
  }
  cn_median <- median(cor[idx])
  h <- hist(cor[idx],breaks=max(length(cor[idx])/100,CFG$minimum_hist_bins),plot=FALSE)
  # rmax is index to the mode
  rzero <- which.max(h$counts)
  cs <- cumsum(h$density)
  # rzero is the index to 1 percentile (hopefully past noise)
  nupper <- nupper <- h$mids[which(cs>0.75)[1]]
  hist_mode <- h$mids[rzero]

  mode_norm_cutoff <- min(CFG$median_magic * cn_median,CFG$mode_norm_cutoff)
  cn_mode <- 0	
  idx_mode  <- idx & df$raw>mode_norm_cutoff & df_reference$raw>mode_norm_cutoff
  if(sum(idx_mode) < CFG$mode_data_points_cutoff) {
    gv$log <<- paste0(isolate(gv$log),flog.error("Unable to normalise. Only have %d points out of %d. Need at least %d bigger than %f to calculate mode.",sum(idx_mode),length(idx_mode),CFG$mode_data_points_cutoff))
    gv$log <<- paste0(isolate(gv$log),flog.error("sample/ref= NA: %d / %d   zero_entries: %d / %d median: %f / %f mean: %f / %f max: %f / %f",
      sum(is.na(df$raw)),sum(is.na(df_reference$raw)),sum(df$raw==0,na.rm=TRUE),sum(df_reference$raw==0,na.rm=TRUE),median(df$raw,na.rm=TRUE),median(df_reference$raw,na.rm=TRUE),mean(df$raw,na.rm=TRUE),mean(df_reference$raw,na.rm=TRUE),max(df$raw,na.rm=TRUE),max(df_reference$raw,na.rm=TRUE)))
    gv$broken <<- TRUE
    if(normalise_on_mode==TRUE) { return(df) }
  } else {
    cn_mode <- tryCatch(
      mlv(cor[idx_mode], method = "lientz", bw = 0.2)$M,
      error=function(cond) { 
        message("mlv failed:") 
        gv$log <<- paste0(isolate(gv$log),flog.warn("mlv failed: %s",cond))
        return(0)
      }
    )
  }
  if(normalise_on_mode==TRUE && cn_mode==0) {
    gv$log <<- paste0(isolate(gv$log),flog.error("Unable to normalise.cn_mode=0"))
    gv$broken <<- TRUE
    return(df)
  }

  warn_norm <- (abs(cn_mode-cn_median)/(cn_mode+cn_median) > 0.1)
  if(warn_norm) { gv$log <<- paste0(isolate(gv$log),flog.info("##############################################"))}
  gv$log <<- paste0(isolate(gv$log),flog.info("::Normalising %s: Mode = %e Median = %e UpperQuartile = %e ModeFromBins = %e",ref_name,cn_mode,cn_median,nupper,hist_mode))
  if(warn_norm) { gv$log <<- paste0(isolate(gv$log),flog.info("##############################################"))}
  cor <- cor / cn_mode
  N = cor * 2
  if(nmax>1) {
    gv$log <<- paste0(isolate(gv$log),flog.info("Applying reference %s: Found %d valid points.",ref_name,nmax))
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
    df$nupper <- nupper
    df$hist_mode <- hist_mode
  }
  else {
    # Flag error
    gv$log <<- paste0(isolate(gv$log),flog.error("No valid corrected bins found in applyReference()."))
    gv$broken <<- TRUE
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
  weights <- 1/ref_mad
  weights <- weights[idx]/sum(weights[idx])
  CNA.object <- CNA(log2(N[idx]),df_reference$chr[idx],df_reference$start[idx],data.type="logratio",presorted=TRUE)
  # smoothed.CNA.object <- smooth.CNA(CNA.object)
  segment.CNA.object <- segment(x = CNA.object, verbose = CFG$segment_verbosity,alpha = CFG$cbs_alpha, nperm = CFG$cbs_nperm )
  numSegments <- nrow(segment.CNA.object$output)
  
  idx_na <- is.na(segment.CNA.object$segRows$startRow)
  segment.CNA.object$segRows$startRow[idx_na] <- segment.CNA.object$segRows$endRow[idx_na] # repair bug in CNA
  if(sum(is.nan(segment.CNA.object$segRows$startRow) | 
           is.na(segment.CNA.object$segRows$startRow) |
           is.nan(segment.CNA.object$segRows$endRow) |
           is.na(segment.CNA.object$segRows$endRow) 
  )) {
    gv$log <<- paste0(isolate(gv$log),flog.error("Invalid smoothed segmented value.")) 
    gv$broken <<- TRUE 
  } 
  # Convert to old co-ordinates
  idx_conv <- 1:length(idx)
  idx_conv <- idx_conv[idx]
  for(i in 1:numSegments) {
    startSeg <- idx_conv[segment.CNA.object$segRows$startRow[i]]
    endSeg   <- idx_conv[segment.CNA.object$segRows$endRow[i]]
    Ns[startSeg:endSeg] <- 2^(segment.CNA.object$output$seg.mean[i])
  }
  df_seg <- segments.pval(segment.CNA.object,alpha = CFG$cbs_alpha, nperm = CFG$cbs_nperm )
  colnames(df_seg) <- c("Gene","Chr","Start","End","Bins","log2N","bstat","pval","lcl","ucl")
  ret <- list(df_seg = df_seg, Ns = Ns, cna_output = segment.CNA.object)
  return(ret)
}
# -------------------------------------------------------------------

# -------------------------------------------------------------------
make_bams <- function(sample) {
  ret <- 0
  chrom_file <- CFG$chrom_file
  chrom_index_file <- CFG$chrom_index_file
  number_of_cores <- CFG$number_of_cores
  num_index_files <- length(list.files(path = dirname(chrom_index_file), pattern = paste0(basename(chrom_index_file),".*")))
  if(num_index_files == 0) {
    gv$log <<- paste0(isolate(gv$log),flog.warning("Making index for %s. This should already have been done.",chrom_file))
    buildindex(chrom_index_file,chrom_file)
  } else {
    gv$log <<- paste0(isolate(gv$log),flog.info("Doing alignment: index file = %s cores = %s",CFG$chrom_index_file,CFG$number_of_cores))
    gv$log <<- paste0(isolate(gv$log),flog.info("Already have: %d files matching %s.*",num_index_files,chrom_index_file))
  }
  bampath <- makepath(panel = CFG$panel_path ,run = CFG$run_path ,sample = sample, filetype = "bam",extension = ".bam",df_cfg = DF_CFG)
  fastq_dir <- dirname(makepath(panel = CFG$panel_path ,run = CFG$run_path ,sample = sample ,filetype = "fastq",extension = ".fastq.gz", df_cfg = DF_CFG))
  if(!file_no_clobber(bampath) && file_apply_lock(bampath)) {  
    gv$log <<- paste0(isolate(gv$log),flog.info("Aligning:  %s",bampath))
    fastqpaths <- list.files(path = fastq_dir,"*fastq.gz",full.names = TRUE)
    if(length(fastqpaths) > 2) {
      gv$log <<- paste0(isolate(gv$log),flog.info("Need to join files matching %s (%d matches)",fastq_dir,length(fastqpaths)))
      fastqpaths[1] <-  paste0(sample,"_R1.fastq.gz")
      fastqpaths[2] <-  paste0(sample,"_R2.fastq.gz")
    } else if (length(fastqpaths) < 1) {
       gv$log <<- paste0(isolate(gv$log),flog.error("Not enough matching fastq files in %s",fastq_dir))
       gv$broken <<- TRUE
       return(NULL)
    }
    cwd <- getwd()
    chrom_index_file <- ifelse(isAbsolutePath(chrom_index_file),chrom_index_file,file.path(cwd,chrom_index_file))
    setwd(dirname(bampath))
    if(length(fastqpaths) > 2) {
        system2('cat',paste0(fastq_dir,"/*R1*fastq.gz"),stdout = fastqpaths[1])
        system2('cat',paste0(fastq_dir,"/*R2*fastq.gz"),stdout = fastqpaths[2])
    }

    if(!file.exists(fastqpaths[1]) || file.size(fastqpaths[1])<1) {
        gv$log <<- paste0(isolate(gv$log),flog.warn("fastq file not found or zero length: %s",fastqpaths[1]))
        gv$broken <<- TRUE
        return(NULL)
    } else {
      fastqpaths1 <- fastqpaths[1]
    }
    if(!file.exists(fastqpaths[2]) || file.size(fastqpaths[2])<1000) {
        gv$log <<- paste0(isolate(gv$log),flog.warn("fastq file not found or looks too short: %s",fastqpaths[2]))
	fastqpaths2 <- NULL
    } else {
        fastqpaths2 <- fastqpaths[2]
    }
    gv$log <<- paste0(isolate(gv$log),flog.info("Have %s and %s",fastqpaths1,fastqpaths2))
    
    if(CFG$spawn_external_aligner==TRUE) {
      if(is.null(fastqpaths2)) {
        ret <- system2('subread-align',c('-i',chrom_index_file,'-r',fastqpaths1,'-o',basename(bampath),
                                  '-t','1','--sv','-T',number_of_cores,'-B',16))
      } else {
        ret <- system2('subread-align',c('-i',chrom_index_file,'-r',fastqpaths1,'-R',fastqpaths2,'-o',basename(bampath),
                                  '-t','1','--sv','-T',number_of_cores,'-B',16))
      }
    } else {
      align(chrom_index_file,
            fastqpaths1,
            readfile2=fastqpaths2,
            type="dna",
            input_format="gzFASTQ",
            output_format="BAM",
            output_file = basename(bampath),
            nthreads = number_of_cores,
            indels = 2,
            unique = FALSE,
            nBestLocations = 16,
            detectSV=TRUE
      )
      if(file.info(bampath)$size <= 1) { ret = -1 }
    }
    lapply(c("fixbam","sam","tmp","bin","fastq.gz"), remove_files_by_suffix)
    if(ret!=0) {
         gv$log <<- paste0(isolate(gv$log),flog.error("subread-align failed."))
         gv$broken <<- TRUE
    } else {
         
         file_remove_lock(bampath)
         # gv$log <<- paste0(isolate(gv$log),flog.info("Created %s: %d bytes",bampath,file.info(bampath)$size))
    }
    setwd(cwd)
    return (NULL)
  } else {
    gv$log <<- paste0(isolate(gv$log),flog.info("Already have: %s",bampath))
    return (NULL)
  }
  return (NULL)
}
# -------------------------------------------------------------------

# -------------------------------------------------------------------
make_vcfs <- function(sample) {
  number_of_cores <- as.integer(CFG$number_of_cores)
  df_files <- data.frame(
    sample = sample,
    path = "",
    file = "vcf",
    resolution = CFG$vcf_section,
    is_targeted = (CFG$vcf_section=="targeted"),
    is_wgs = (CFG$vcf_section=="targeted"),
    stringsAsFactors = FALSE
    )
  bampath <- makepath(panel = CFG$panel_path ,run = CFG$run_path ,sample = sample ,filetype = "bam",extension = ".bam",df_cfg = DF_CFG)
  vcfpath <- makepath(panel = CFG$panel_path ,run = CFG$run_path ,sample = sample ,filetype = "vcf", extension = ".vcf",df_cfg = DF_CFG)
  vcffile <- basename(vcfpath)
  bamfile <- basename(bampath)
  sample_idx <- grepl(sample,df_files$sample)
  df_files[sample_idx,"path"] <- vcfpath
  refGenomeFile <- CFG$chrom_file
  snp_annotation_file <- CFG$snp_annotation_file
  if(!file_no_clobber(vcfpath) && file_apply_lock(vcfpath)) {
    gv$log <<- paste0(isolate(gv$log),flog.info("Making: %s",vcfpath))
    cwd <- getwd()
    bampath <- ifelse(isAbsolutePath(bampath),bampath,file.path(cwd,bampath))
    vcfpath <- ifelse(isAbsolutePath(vcfpath),vcfpath,file.path(cwd,vcfpath))
    refGenomeFile <- ifelse(isAbsolutePath(refGenomeFile),refGenomeFile,file.path(cwd,refGenomeFile))
    snp_annotation_file <- ifelse(isAbsolutePath(snp_annotation_file),snp_annotation_file,file.path(cwd,snp_annotation_file))
    setwd(dirname(bampath))
    exactSNP(bamfile,
             isBAM=TRUE,
             refGenomeFile=refGenomeFile,
             SNPAnnotationFile=snp_annotation_file, 
             outputFile=vcffile,
             # qvalueCutoff=12,
             # minAllelicFraction=0,
             minAllelicBases=CFG$snp_min_mismatches,
             minReads=CFG$snp_min_reads,
             maxReads=CFG$snp_max_reads,
             # minBaseQuality=13,
             # nTrimmedBases=3,
             nthreads=number_of_cores)
    gv$log <<- paste0(isolate(gv$log),flog.info("Copying %s to %s",vcffile,vcfpath))
    file.copy(from = vcffile, to = vcfpath)
    # file.remove(vcffile)
    lapply(c("fixbam","sam","tmp","bin"), remove_files_by_suffix)
    file_remove_lock(vcfpath)
    setwd(cwd)
  } else {
    gv$log <<- paste0(isolate(gv$log),flog.info("Already have: %s",vcfpath))
  }
  return(NULL)
}
# -------------------------------------------------------------------

# -------------------------------------------------------------------
make_bafs <- function(sample) {
  df_files <- data.frame(
  sample = sample,
  path = "",
  file = "baf",
  resolution = CFG$baf_section,
  is_targeted = (CFG$baf_section=="targeted"),
  is_wgs = (CFG$baf_section=="targeted"),
  stringsAsFactors = FALSE
  )
  sample_idx <- grepl(sample,df_files$sample)
  bafpath <- makepath(panel = CFG$panel_path ,run = CFG$run_path ,sample = sample ,filetype = "baf",extension = "_baf.tsv",df_cfg = DF_CFG)
  baffile <- basename(bafpath)
  vcfpath <- makepath(panel = CFG$panel_path ,run = CFG$run_path ,sample = sample ,filetype = "vcf",extension = ".vcf",df_cfg = DF_CFG)
  vcffile <- basename(vcfpath)
  df_files[sample_idx,"path"] <- bafpath
  
  if(!file_no_clobber(bafpath) && file_apply_lock(bafpath)) {
    gv$log <<- paste0(isolate(gv$log),flog.info("Reading VCF file %s",vcfpath))   
    vcf <- readVcf(vcfpath,genome = "hg19")
    if(nrow(vcf)>0) {
      idx <- isSNV(vcf, singleAltOnly=FALSE)
      gv$log <<- paste0(isolate(gv$log),flog.info("VCF file has %d entries and %d SNVs",nrow(vcf),sum(idx)))   
      snvinfo <- info(vcf[idx])
      df_baf <- data.frame(chr = paste0("chr",t(as.data.frame(strsplit(rownames(snvinfo),":")))[,1]), 
                           start = start(ranges(rowRanges(vcf[idx]))), 
                           end = end(ranges(rowRanges(vcf[idx]))), 
                           freq = as.double(snvinfo$MM)/as.double(snvinfo$DP), 
                           read_depth = as.double(snvinfo$DP), 
                           alt_reads = as.double(snvinfo$MM),
                           ref_reads = (as.double(snvinfo$DP) - as.double(snvinfo$MM)),
                           stringsAsFactors = FALSE)
    }
    else {
      gv$log <<- paste0(isolate(gv$log),flog.warn("Empty VCF file. Writing a fake row."))    
      df_baf <- data.frame(chr = "MT", start = 1, end = 1, freq = 0, ref_reads = 0, alt_reads = 0,stringsAsFactors = FALSE)   
    } 
    gv$log <<- paste0(isolate(gv$log),flog.info("Have %d entries from vcf before removal unwanted rows.",nrow(df_baf)))        
    idx <- grepl("GL",df_baf$chr) | 
          grepl("MT",df_baf$chr) | 
          df_baf$alt_reads >= CFG$snp_max_reads | 
          df_baf$ref_reads >= CFG$snp_max_reads 
#          df_baf$alt_reads <= CFG$snp_min_mismatches |
#          df_baf$ref_reads <= CFG$snp_min_mismatches 
    df_baf <- df_baf[!idx,]
    if(nrow(df_baf) > CFG$baf_max_rows) { gv$log <<- paste0(isolate(gv$log),flog.info("Have %d entries - need to merge some.",nrow(df_baf))) }
	    else { gv$log <<- paste0(isolate(gv$log),flog.info("Have %d entries from vcf after removal of GL, MT and support >  %d.",nrow(df_baf),CFG$snp_max_reads)) }
    niter = 1
    merge_distance <- CFG$baf_merge_distance
    while (nrow(df_baf) > CFG$baf_max_rows && niter < CFG$baf_max_merge_iterations) {
      df_baf$start_next <- 0
      df_baf$end_next <- 0
      df_baf$end_prev <- 0
      df_baf$chr_next <- 0
      df_baf$chr_prev <- 0
      df_baf$alt_reads_next <- 0
      df_baf$ref_reads_next <- 0
      df_baf$freq_next <- 0
      df_baf$merge_prev <- FALSE
      df_baf$start_next[1:(nrow(df_baf)-1)] <- df_baf$start[2:(nrow(df_baf))]
      df_baf$end_next[1:(nrow(df_baf)-1)] <- df_baf$end[2:(nrow(df_baf))]
      df_baf$end_prev[2:(nrow(df_baf))] <- df_baf$end[1:(nrow(df_baf)-1)]
      df_baf$chr_next[1:(nrow(df_baf)-1)] <- df_baf$chr[2:(nrow(df_baf))]
      df_baf$chr_prev[2:(nrow(df_baf))] <- df_baf$chr[1:(nrow(df_baf)-1)]
      df_baf$alt_reads_next[1:(nrow(df_baf)-1)] <- df_baf$alt_reads[2:(nrow(df_baf))]
      df_baf$ref_reads_next[1:(nrow(df_baf)-1)] <- df_baf$ref_reads[2:(nrow(df_baf))]
      df_baf$freq_next[1:(nrow(df_baf)-1)] <- df_baf$freq[2:(nrow(df_baf))]
      # merge if the next one is close enough and the previous one is far enough
      df_baf$merge_next <- ((df_baf$start_next-df_baf$end) <= merge_distance) & 
        ((df_baf$start - df_baf$end_prev) > merge_distance) &
        (df_baf$chr==df_baf$chr_next) & 
        (df_baf$chr==df_baf$chr_prev)
      df_baf$merge_prev[2:(nrow(df_baf))] <- df_baf$merge_next[1:(nrow(df_baf)-1)]
      # Some of the SNVs being merged onto will be replaced with the read counts of the other SNV
      # All will have their end co-ord replaced.
      replace_with_next <- df_baf$merge_next & ((df_baf$ref_reads + df_baf$alt_reads) < (df_baf$ref_reads_next + df_baf$alt_reads_next))
      df_baf$ref_reads[replace_with_next] <- df_baf$ref_reads_next[replace_with_next] 
      df_baf$alt_reads[replace_with_next] <- df_baf$alt_reads_next[replace_with_next]
      df_baf$freq[replace_with_next] <- df_baf$freq_next[replace_with_next]
      df_baf$end[replace_with_next] <- df_baf$end_next[replace_with_next]
      
      df_baf_new <- df_baf[!df_baf$merge_prev,]
      baf_len <- df_baf_new$end - df_baf_new$start
      baf_len <- baf_len[!is.na(baf_len)]
      merge_distance <- median(df_baf$start_next - df_baf$end,na.rm=T)
      flog.info("Merging BAF table: iter = %d,  rows = %d -> %d, median_dist = %f, merge_dist = %f, median_length = %f, med(len!=0) = %f, N(len!=0) = %d",
                niter,nrow(df_baf),nrow(df_baf_new),median(df_baf$start_next - df_baf$end,na.rm=T),merge_distance,median(baf_len,na.rm=T),median(baf_len[baf_len!=0],na.rm=T),sum(baf_len!=0,na.rm=T))
      df_baf <- df_baf_new
      df_baf$start_next <- NULL
      df_baf$end_next <- NULL
      df_baf$end_prev <- NULL
      df_baf$chr_next <- NULL
      df_baf$chr_prev <- NULL
      df_baf$start_next <- NULL
      df_baf$merge_next <- NULL
      df_baf$merge_prev <- NULL
      df_baf$alt_reads_next <- NULL
      df_baf$ref_reads_next <- NULL
      df_baf$freq_next <- NULL
      niter = niter + 1
    }
    
    if(nrow(df_baf)==0) {
      gv$log <<- paste0(isolate(gv$log),flog.warn("Empty VCF file. Writing a fake row."))    
      df_baf <- data.frame(chr = "MT", start = 1, end = 1, freq = 0, ref_reads = 0, alt_reads = 0,stringsAsFactors = FALSE)   
    }
    gv$log <<- paste0(isolate(gv$log),flog.info("Writing %s with %d entries",bafpath,nrow(df_baf)))    
    write.table(df_baf,bafpath,row.names=FALSE,quote = FALSE,sep="\t")
    file_remove_lock(bafpath)
  } else {
    gv$log <<- paste0(isolate(gv$log),flog.info("Already have: %s",bafpath))
  }
  return(df_files)
}
# -------------------------------------------------------------------

# -------------------------------------------------------------------
make_deletions <- function(sample) {
  whitelist_path = CFG$common_deletion_list
  df_whitelist <- read.table(whitelist_path,header=TRUE,stringsAsFactors=FALSE)
  df_files <- data.frame(
    sample = sample,
    path = "",
    file = "deletions",
    resolution = CFG$del_section,
    is_targeted = (CFG$del_section=="targeted"),
    is_wgs = (CFG$del_section=="targeted"),
    stringsAsFactors = FALSE
  )
  

  sample_idx <- grepl(sample,df_files$sample)
  delpath <- makepath(panel = CFG$panel_path ,run = CFG$run_path ,sample = sample ,filetype = "sv",resolution = CFG$del_section, extension = "_del.tsv",df_cfg = DF_CFG)
  delfile <- basename(delpath)
  corrected_counts_path <- makepath(panel = CFG$panel,run = CFG$run_path,sample = sample,filetype = "cnv",resolution = CFG$del_section, extension = ".tsv",df_cfg = DF_CFG)
  df_files[sample_idx,"path"] <- delpath
  if(!file_no_clobber(delpath) && file_apply_lock(delpath)) {
    gv$log <<- paste0(isolate(gv$log),flog.info("Need: %s",delpath))
    df_counts <- read.table(corrected_counts_path,header=TRUE,stringsAsFactors=FALSE)
    gv$log <<- paste0(isolate(gv$log),flog.info("Have %s with %d rows",corrected_counts_path,nrow(df_counts)))
    cbc_list <- do_cbc(df_counts$N,df_counts)
    df_del <- cbc_list$df_seg
    gv$log <<- paste0(isolate(gv$log),flog.info("Found %d segments",nrow(df_del)))     
    df_del$whiteListed <- FALSE
    df_del$Arm <- ""
    df_del$Gene <- "NA"
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
        # If the pval is NA then it's a whole-chromosome event. If other stuff is NA, it's dodgy.
        if(!(is.na(row$pval) || is.na(df_del$weighted_mean[n]) || is.na(df_del$weighted_sem[n]))) {
          if(!is.na(row$pval) || (df_del$weighted_mean[n] - 2)^2 > 3 * df_del$weighted_sem[n]) {
            idx <- (!is.na(df_whitelist$Gene)) &
              row$Chr == df_whitelist$Chr &
              ((row$Start <= df_whitelist$TestEnd & row$Start>= df_whitelist$TestStart) | 
                (row$End >= df_whitelist$TestStart & row$End <= df_whitelist$TestEnd) |
                (row$Start <= df_whitelist$TestStart & row$End >= df_whitelist$TestEnd)) &
              0 != df_whitelist$GainLoss
            # row$Start <= df_whitelist$RequiredStart &
            # row$End >= df_whitelist$RequiredEnd &
            # "pq" != df_whitelist$Arm) $
            if(sum(idx)>=1) {
              df_del$whiteListed[n] <- TRUE
              df_del$Arm[n] <- df_whitelist$Arm[idx][1]
              df_del$Gene[n] <- df_whitelist$Gene[idx][1]
              gain_loss <- ifelse(df_whitelist$GainLoss[idx][1]>0,"duplication","deletion")
              gv$log <<- paste0(isolate(gv$log),
                                flog.info("Found whitelisted entry: %s%s",
                                          row$Chr,df_del$Arm[n]))
            }
            if(sum(idx)>=2) {
              df_del$Gene[n] <- paste0(df_whitelist$Gene[idx],collapse = ",")
              gv$log <<- paste0(isolate(gv$log),flog.info("Multi-gene deletion: %s%s: %s",row$Chr,df_del$Arm[n],df_del$Gene[n]))
            }
          }
        } # if not dodgy
        else {
          gv$log <<- paste0(isolate(gv$log),flog.warn("Segmentation stats look odd in make_deletions(): Entry # %ds",n)) 
        }
      }
    } # for
    gv$log <<- paste0(isolate(gv$log),flog.info("Writing %s with %d entries",delpath,nrow(df_del)))    
    # df_del <- df_del[df_del$whiteListed,]
    # df_del$whiteListed <- NULL
    colnames(df_del)[colnames(df_del)=="Chr"] <- "chr"
    colnames(df_del)[colnames(df_del)=="Start"] <- "start"
    colnames(df_del)[colnames(df_del)=="End"] <- "end"
    colnames(df_del)[colnames(df_del)=="Bins"] <- "bins"
    colnames(df_del)[colnames(df_del)=="Gene"] <- "gene"
    write.table(df_del,delpath,row.names=FALSE,quote = FALSE,sep="\t")
    file_remove_lock(delpath)
  } else {
    gv$log <<- paste0(isolate(gv$log),flog.info("Already have: %s",delpath))
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
  xWbar = wtd.mean(x,w,na.na.rm.rm)
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
make_translocations <- function(sample) {
  #Chr  Location	Chr	Location	SameStrand	nSupport  
  # 12M1411.bam.breakpoints.txt 
  # 12M1411.bam.fusion.txt
  # ID  Partner	Gene	Chr	Arm	Band	Start	End	Strand	MinSupport
  df_files <- data.frame(
    sample = sample,
    path = "",
    file = "fusion",
    resolution = CFG$fusion_section,
    is_targeted = (CFG$fusion_section=="targeted"),
    is_wgs = (CFG$fusion_section=="targeted"),
    stringsAsFactors = FALSE
  )
  fusions_white_list <- CFG$fusions_white_list
  df_whitelist <- read.table(fusions_white_list,header=TRUE,stringsAsFactors=FALSE)
  df_whitelist$Chr <- sub("^chr","",df_whitelist$Chr)
  
  bampath <- makepath(panel = CFG$panel_path ,run = CFG$run_path ,sample = sample ,filetype = "bam",extension = ".bam",df_cfg = DF_CFG)
  fusions_path <- makepath(panel = CFG$panel,run = CFG$run_path,sample = sample,filetype = "sv",extension = "_fusions.tsv",df_cfg = DF_CFG)
  # putative_fusion_list <- sub(".bam$",".bam.fusion.txt",bampath)
  putative_fusion_list <- sub(".bam$",".bam.breakpoints.txt",bampath)
  sample_idx <- grepl(sample,df_files$sample)
  df_files[sample_idx,"path"] <- fusions_path

  if(!file_no_clobber(fusions_path) && file_apply_lock(fusions_path)) {
    if(file.exists(putative_fusion_list)) {
      gv$log <<- paste0(isolate(gv$log),flog.info("Reading %s",putative_fusion_list))
      df_putative_fusions <- read.table(putative_fusion_list,header=TRUE,comment.char = "#",stringsAsFactors=FALSE)
      gv$log <<- paste0(isolate(gv$log),flog.info("Found: %d rows %d cols",nrow(df_putative_fusions),ncol(df_putative_fusions)))
      colnames(df_putative_fusions) <- c("chr",  "start",  "chr2",  "end",      "sameStrand",   "nSupport","BreakPoint1_GoUp","BreakPoint2_GoUp")[1:ncol(df_putative_fusions)]
      df_putative_fusions <- df_putative_fusions[df_putative_fusions$nSupport>=CFG$minimumFusionSuport,]
      gv$log <<- paste0(isolate(gv$log),flog.info("After filtering out support < %d: %d rows %d cols",CFG$minimumFusionSuport,nrow(df_putative_fusions),ncol(df_putative_fusions)))
      gv$log <<- paste0(isolate(gv$log),flog.info("Making: %s",fusions_path))
      if(nrow(df_putative_fusions) > 0) {
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
      } else {
        df_fusions <- df_putative_fusions
      } # if no rows
    } else {
      # write null file
      gv$log <<- paste0(isolate(gv$log),flog.warn("Can't find: %s. Need to run gridss or Rsubread::Align().",putative_fusion_list))
      df_fusions <- data.frame(chr=character,start=integer,chr2=character,end=integer,sameStrand=logical,nSupport=integer)
    }
    write.table(df_fusions,fusions_path,row.names=FALSE,col.names=TRUE,quote = FALSE,sep="\t")
    file_remove_lock(fusions_path)
  } else {
    gv$log <<- paste0(isolate(gv$log),flog.info("Already have: %s",fusions_path))
  }
  return(df_files)
}
# -------------------------------------------------------------------

# -------------------------------------------------------------------
make_bins <- function(bedtype) {
  bin_size <- BINSIZES[bedtype]
  bedfile <- CFG[[BEDFILENAMES[bedtype]]]
  chrom_info_file <- CFG$chrom_info_file
  
  if(is.na(bedfile) || length(bedfile)<1) {
    gv$log <<- paste0(isolate(gv$log),flog.warn("Skipping bedfile for bin size %d",bin_size))  
    return(NULL)
  }
  destfile <- bedfile
  if(!file_no_clobber(destfile) && file_apply_lock(destfile)) {
    gv$log <<- paste0(isolate(gv$log),flog.info("Reading: %s",chrom_info_file))
    chromInfo <- read.table(chrom_info_file,stringsAsFactors=FALSE)
    colnames(chromInfo) <- c("chr","end","file")
    rownames(chromInfo) <- chromInfo$chr
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
    file_remove_lock(destfile)
    gv$log <<- paste0(isolate(gv$log),flog.info("Making bins for %s: Done",destfile))    
  }
  gv$log <<- paste0(isolate(gv$log),flog.info("Have: %s",destfile))
  return(NULL)
}
# -------------------------------------------------------------------

# -------------------------------------------------------------------
make_gc_tables <- function(bedtype) {
  chrom_file <- CFG$chrom_file
  bedfile <- CFG[[BEDFILENAMES[bedtype]]]
  gc_file <- bedrow2gcrow(BEDFILENAMES[bedtype])$value
  
  if(!file_no_clobber(gc_file) && file_apply_lock(gc_file)) {
    gv$log <<- paste0(isolate(gv$log),flog.info("Making GC corrections. Reading: %s",bedfile))
    intervals <- import.bed(bedfile)
    gv$log <<- paste0(isolate(gv$log),flog.info("Making GC corrections. Reading: %s",chrom_file))
    reference <- readDNAStringSet(chrom_file)
    
    gv$log <<- paste0(isolate(gv$log),flog.info("Making: %s",gc_file))
    # for each chr, something like this
    find_gc_content <- function(idx) { 
      seq <- reference[idx]
      seq_names <- names(seq)
      seq_names <-sub("^ *","",seq_names)
      seq_names <-sub(" .*$","",seq_names)

      print(seq_names)
      # Add chr if needed, because bedfile has one
      if(sum(!grepl("^chr",seq_names))>0) {
        seq_names <- paste0("chr",seq_names)
      }
      gv$log <<- paste0(isolate(gv$log),flog.info("Making GC content for %s:%s",seq_names,gc_file))    
      chr_int <- intervals[seqnames(intervals)==seq_names] 
      if(length(chr_int)==0) {
        gv$log <<- paste0(
          isolate(gv$log),
          flog.warn("Seqname %s in %s not found in reference %s",seq_names,chrom_file,bedfile)
        )
        return(NULL)
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
    write.table(df_gc,gc_file,row.names=FALSE,quote = FALSE,sep="\t")
    file_remove_lock(gc_file)
    gv$log <<- paste0(isolate(gv$log),flog.info("Making GC content for %s: Done",gc_file))  
    return(NULL)
  }
  gv$log <<- paste0(isolate(gv$log),flog.info("Have: %s",gc_file))
  return(NULL)
}
# -------------------------------------------------------------------

# -------------------------------------------------------------------
make_gc_plots <- function(output_path) {
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

  loess_predict <- loess(raw ~ gc, data = df_gc,  span = CFG$loess_span) # subset = !blacklist
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
  bins <- seq(0, 1, by = CFG$loess_bin_width)
  fapply_median <- function(n,bins) { median(df_gc$raw[df_gc$gc>=(bins[n] - CFG$loess_bin_width) & df_gc$gc<(bins[n] + CFG$loess_bin_width)]) }
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
  
  ggsave(filename=plot_file, plot=p)
}
# -------------------------------------------------------------------


# -------------------------------------------------------------------
read_segmentation_table <- function (filename) {
  print(paste0("Reading ",filename))
  if(!file.exists(filename)) { 
    print(paste0("Missing ",filename,": Skipping"))
    return (NULL)
  }
  df_seg <- read.table(filename,header=TRUE,stringsAsFactors=FALSE,sep = '\t')
  df_seg$chr_arm <- paste0(df_seg$chr,"_",df_seg$Arm)
  seg_mode <-  mlv(df_seg$N[df_seg$N > 1 & !is.na(df_seg$N)], method = "lientz", bw = 0.2)$M
  # regions <- unique(df_seg$chr_arm[df_seg$whiteListed==TRUE])
  regions <- unique(df_seg$chr_arm)
  df_output <- data.frame(sample = sub("_.*tsv","",basename(filename)),
                          run = gsub("/.*$","",gsub("^.*/160","160",gsub("^.*/161","161",filename))),
                          regions = regions, 
                          chr = "NA", 
                          arm = "NA",
                          start = 0, 
                          end = 0, 
                          bins = 0,  
                          weighted_mean = 0, 
                          mode_N = seg_mode,
                          stringsAsFactors = FALSE)
  rownames(df_output) <- regions
  for (region in regions) {
    idx <- df_seg$chr_arm==region & df_seg$whiteListed==TRUE
    # if no whitelisted, take the other
    if(sum(idx)>0) {
      df_region <- df_seg[idx,]
    } else {
      df_region <- df_seg[df_seg$chr_arm==region,]
    }
    row <- df_region[1,]
    df_output[region,"chr"] <- row$chr
    df_output[region,"arm"] <- ifelse(row$Arm=="",NA,row$Arm)
    df_output[region,"start"] <- min(df_region$start)
    df_output[region,"end"] <- max(df_region$end)
    df_output[region,"bins"] <- sum(df_region$bins)
    df_output[region,"weighted_mean"]  <-  weighted.mean(df_region$weighted_mean,df_region$bins/sum(df_region$bins),na.rm=TRUE)
    df_output[region,"minus_log_cbc_pval"] <- ifelse(sum(!is.na(df_region$pval))==0,NA,-log10(min(df_region$pval,na.rm=TRUE)))
    df_output[region,"minus_pseudo_log_likelihood_var"] <- (sum(df_region$bins * (df_region$weighted_mean - seg_mode)^2 / df_region$weighted_var))
    df_output[region,"minus_pseudo_log_likelihood_sem"] <- (sum((df_region$weighted_mean - seg_mode)^2 / df_region$weighted_sem))
  }
  return(df_output)
}
# -------------------------------------------------------------------
  

# -------------------------------------------------------------------
read_counts_table <- function(filename,whitelist_path,in_parfor = FALSE) {
  if(!in_parfor) { gv$log <<- paste0(isolate(gv$log),flog.info("Reading whitelist from %s.",whitelist_path)) }
  df_whitelist <- read.table(whitelist_path,header=TRUE,stringsAsFactors=FALSE)
  df_roi <- df_whitelist[!is.na(df_whitelist$Gene),c("Gene","Chr","Arm","TestStart","TestEnd")]
  rownames(df_roi) <- df_roi$Gene
  # if(!in_parfor) { gv$log <<- paste0(isolate(gv$log),flog.info("Looking at regions for %s.",paste0(df_roi$Gene,collapse = ",")))}
  if(!in_parfor) { gv$log <<- paste0(isolate(gv$log),flog.info("Reading counts from %s.",filename)) }
  if(file.exists(filename)) {
    df <- read.table(filename,header=TRUE,stringsAsFactors=FALSE,sep = '\t')  
  } else {
    if(!in_parfor) { gv$log <<- paste0(isolate(gv$log),flog.warn("Unable to open %s.",filename)) }
    return(NULL)
  }
  
  s <- sub(".*cnv/","",filename)
  s <- sub("_.*tsv","",s)
  if(ncol(df)==1) {
    colnames(df)[1] <- "raw"
    df$range <- rownames(df)
    df$chr <- sub(":.*$","",rownames(df))
    df$start <- as.integer(sub("-.*$","",sub("^.*:","",rownames(df))))
    df$end <- as.integer(sub("^.*-","",sub("^.*:","",rownames(df))))
  }
  df$filename <- s
  df$arm <- ifelse(df$start < CFG$pq_boundary[df$chr],"p","q")
  df$gene <- ""
  df$exon <- ""
  # Mark a region as being in a gene/exon if there is overlap - judged by the middle of the region being in the feature
  # For regions that overlap two of these features, assign the features to the second  feature
  df_mid <- 0.5 * (df$start+df$end)
  for (gene_exon in df_roi$Gene) {
    row <- df_roi[gene_exon,]
    idx <- ifelse(df$chr == row$Chr & (df_mid <= row$TestEnd) & (df_mid >= row$TestStart),TRUE,FALSE)
    df$exon[idx] <- gene_exon
    df$gene[idx] <- sub("_.*$" ,"",gene_exon)
  }
  return(df)
}
# -------------------------------------------------------------------


# -------------------------------------------------------------------
make_deletions_stats <- function (make_distributions,file_list,bedtype) {
  if(make_distributions==TRUE) {
    gv$log <<- paste0(isolate(gv$log),flog.info("Making NULL distributions for deletions."))
  }
  
  #outpath <- makepath(panel = CFG$panel_path,run = CFG$run_path,resolution = bedtype,filetype = "reference",df_cfg = DF_CFG)
  #reference_file_name <- paste0(reference_file_name,".tsv")
  samples_file = makepath(panel = CFG$panel_path,sample = "normal_deletions_samples",run = CFG$run_path,resolution = bedtype,filetype = "reference",extension = ".tsv",df_cfg = DF_CFG)  
  mean_file = makepath(panel = CFG$panel_path,sample = "normal_deletions_mean",run = CFG$run_path,resolution = bedtype,filetype = "reference",extension = ".tsv",df_cfg = DF_CFG)  
  med_file = makepath(panel = CFG$panel_path,sample = "normal_deletions_med",run = CFG$run_path,resolution = bedtype,filetype = "reference",extension = ".tsv",df_cfg = DF_CFG)  
  madcv_file = makepath(panel = CFG$panel_path,sample = "normal_deletions_madcv",run = CFG$run_path,resolution = bedtype,filetype = "reference",extension = ".tsv",df_cfg = DF_CFG)  
  cv_file = makepath(panel = CFG$panel_path,sample = "normal_deletions_cv",run = CFG$run_path,resolution = bedtype,filetype = "reference",extension = ".tsv",df_cfg = DF_CFG)  
  n_file = makepath(panel = CFG$panel_path,sample = "normal_deletions_n",run = CFG$run_path,resolution = bedtype,filetype = "reference",extension = ".tsv",df_cfg = DF_CFG)
  nraw_file = makepath(panel = CFG$panel_path,sample = "normal_deletions_nraw",run = CFG$run_path,resolution = bedtype,filetype = "reference",extension = ".tsv",df_cfg = DF_CFG)  
  whitelist_path = CFG$common_deletion_list
  # df_samples <- do.call(rbind, lapply(file_list,read_counts_table,whitelist_path,TRUE))
  if(make_distributions==TRUE) {
    if((!file_no_clobber(samples_file)) && file_apply_lock(samples_file)) { 
      gv$log <<- paste0(isolate(gv$log),flog.info("Making deletions references. Reading counts."))
      cores=detectCores()
      gv$log <<- paste0(isolate(gv$log),flog.info("Found %d cores. Making cluster with %d of them.",cores,CFG$cluster_cpus))
      cl <- makeCluster(CFG$cluster_cpus,type=CFG$cluster_type)
      gv$log <<- paste0(isolate(gv$log),flog.info("Registering cluster."))
      registerDoParallel(cl)
      gv$log <<- paste0(isolate(gv$log),flog.info("Splitting loop"))
      read_counts_table_function <- read_counts_table
      if(CFG$cluster_cpus==1) {
        df_samples <- do.call(rbind, lapply(file_list,read_counts_table,whitelist_path,FALSE))  
      } else {
        df_samples <- foreach (n=1:length(file_list), .combine=rbind, .verbose=TRUE) %dopar% { df_row <- read_counts_table_function(file_list[n],whitelist_path,TRUE); return(df_row) }
      }
      stopCluster(cl)  
      write.table(df_samples,samples_file,row.names=FALSE,sep="\t",quote = FALSE)
      file_remove_lock(samples_file)
    } else {
      gv$log <<- paste0(isolate(gv$log),flog.info("Already have deletions samples for references. Reading %s samples_file",samples_file))
      df_samples <- read.table(samples_file,sep="\t",header=TRUE,stringsAsFactors=FALSE)
    }
  } else {
    gv$log <<- paste0(isolate(gv$log),flog.info("Reading file list %s",file_list))
    df_samples <- read_counts_table(file_list,whitelist_path)
  }  
  
  chrs <- unique(df_samples$chr)
  genes <- unique(df_samples$gene[df_samples$gene!=""])
  exons <- unique(df_samples$exon[df_samples$exon!=""])
  arms <- unique(df_samples$arm)
  samples <- unique(df_samples$filename)
  chr_arms <- c(paste0(chrs,arms[1]),paste0(chrs,arms[2]))
  all_chrs_names <- c(chrs,paste0(chrs,arms[1]),paste0(chrs,arms[2]))
  subset_names <- c(chrs,paste0(chrs,arms[1]),paste0(chrs,arms[2]),genes,exons)
  reference_subset_names <- all_chrs_names
  
  n_chrs <- length(chrs)
  n_all_chrs_names <- length(all_chrs_names)
  n_chr_arms <- length(chr_arms)
  n_subsets <- length(subset_names)
  n_samples <- length(samples)
  n_genes <- length(genes)
  n_exons <- length(exons)
  n_counts <- nrow(df_samples)
  n_reference_subsets <- length(reference_subset_names)

  idx_subsets_samples <- vector("integer",n_counts)
  idx_subsets_exons <- vector("integer",n_counts)
  idx_subsets_genes <- vector("integer",n_counts)
  idx_subsets_chrs <- vector("integer",n_counts)
  idx_subsets_chr_arms <- vector("integer",n_counts)
  

  xraw <- df_samples$raw
  nraw <- matrix(0,n_samples,n_subsets)  
  
  if(make_distributions==FALSE || ((!file_no_clobber(nraw_file)) && file_apply_lock(nraw_file))) { 
    gv$log <<- paste0(isolate(gv$log),flog.info("Making indexes"))
    for (n in 1:n_samples) { idx_subsets_samples[df_samples$filename == samples[n]] <- n }
    for (n in 1:n_genes) { idx_subsets_genes[df_samples$gene == genes[n]] <- n }
    for (n in 1:n_exons) { idx_subsets_exons[df_samples$exon == exons[n]] <- n }
    for (n in 1:n_chrs) { idx_subsets_chrs[df_samples$chr == chrs[n]] <- n }
    for (n in 1:n_chr_arms) { idx_subsets_chr_arms[paste0(df_samples$chr,df_samples$arm) == chr_arms[n]] <- n }
    
    gv$log <<- paste0(isolate(gv$log),flog.info("Making %s.",nraw_file))
    cores=detectCores()
    gv$log <<- paste0(isolate(gv$log),flog.info("Found %d cores. Making cluster with %d them.",cores,CFG$cluster_cpus))
    cl <- makeCluster(CFG$cluster_cpus,type=CFG$cluster_type)
    gv$log <<- paste0(isolate(gv$log),flog.info("Registering cluster."))
    registerDoParallel(cl)
    gv$log <<- paste0(isolate(gv$log),flog.info("Splitting loop"))
    master_only_varnames <- c("nraw","df_samples","gv")
    nraw <- foreach (n=1:n_subsets, .combine=cbind, .verbose=TRUE,.noexport=master_only_varnames) %dopar% {
      # s1 <- subset_names[n]
      # gv$log <<- paste0(isolate(gv$log),flog.info("Combining bins n = %d : %s",n,s1))
      # print(paste0("Combining bins n = ",n," : ",s1))
      if(n <= n_chrs) { 
        idx1 <- (idx_subsets_chrs==n)
      } else if(n > n_chrs && n <= (n_chrs + n_chr_arms)) { 
        idx1 <- (idx_subsets_chr_arms==(n-n_chrs))
      } else if (n > (n_chrs + n_chr_arms) && n <= (n_chrs + n_chr_arms + n_genes)) { 
        idx1 <- (idx_subsets_genes==(n - n_chrs - n_chr_arms)) 
      } else if (n > (n_chrs + n_chr_arms + n_genes ) && n <= (n_chrs + n_chr_arms + n_genes + n_exons)) { 
        idx1 <- (idx_subsets_exons==(n - n_chrs - n_chr_arms - n_genes)) 
      } else {
        # Parallel index out of range
        # gv$log <<- paste0(isolate(gv$log),flog.error("Out of range index: %d",n))
      }
      nraw_col <- matrix(0,n_samples,1)   
      if(sum(idx1)>0) {
        for (m in 1:n_samples) {
          idx1s <- (idx1 & (idx_subsets_samples==m))
          if(sum(idx1s)>0) {
            nraw_col[m] <- sum(xraw[idx1s])
          }
        }
      }
      return(nraw_col)
    }
    stopCluster(cl)
    colnames(nraw) <- subset_names
    rownames(nraw) <- samples
    if(make_distributions==TRUE) { 
      gv$log <<- paste0(isolate(gv$log),flog.info("Done!!! Writing %s.",nraw_file))
      write.table(nraw,nraw_file,sep="\t",quote = FALSE)
      file_remove_lock(nraw_file)
      }
  } else if(make_distributions==TRUE) {
      gv$log <<- paste0(isolate(gv$log),flog.info("Reading %s.",nraw_file))
      nraw <- read.table(nraw_file,header=TRUE,stringsAsFactors=FALSE)
  }
    
  if(make_distributions==FALSE || ((!file_no_clobber(n_file))  && file_apply_lock(n_file))) { 
    gv$log <<- paste0(isolate(gv$log),flog.info("Making %s and other stats files.",n_file))
    cores=detectCores()
    gv$log <<- paste0(isolate(gv$log),flog.info("Found %d cores. Making cluster with %d them.",cores,CFG$cluster_cpus))
    cl <- makeCluster(CFG$cluster_cpus,type=CFG$cluster_type)
    gv$log <<- paste0(isolate(gv$log),flog.info("Registering cluster."))
    registerDoParallel(cl)
    gv$log <<- paste0(isolate(gv$log),flog.info("Splitting loop"))
    master_return <- foreach (idx_1=1:n_subsets, .combine=rbind, .verbose=TRUE) %dopar% {
    # for (idx_1 in 1:n_subsets) {
      s1 <- subset_names[idx_1]
      n1 <- nraw[,s1]
      col_ret <-  matrix(0,1,n_reference_subsets * 5)
      if(sum(n1)>0) {
        for (idx_2 in 1:n_reference_subsets) {
          s2 <- reference_subset_names[idx_2]
          n2 <- nraw[,s2]
          if(sum(n2)>0) {
            # gv$log <<- paste0(isolate(gv$log),flog.info("Calculating ratio for %s:%s",s1,s2))
            idx_n <- (n2>0)
            r12 <- n1[idx_n]/n2[idx_n]
            col_ret[idx_2 + 0 * n_reference_subsets] <-mean(r12)
            col_ret[idx_2 + 1 * n_reference_subsets] <-median(r12)
            col_ret[idx_2 + 2 * n_reference_subsets] <-sd(r12)
            col_ret[idx_2 + 3 * n_reference_subsets] <-mad(r12)
            col_ret[idx_2 + 4 * n_reference_subsets] <-length(r12)
            #xmean[s1,s2] <- mean(r12)  
            #xmed[s1,s2] <- median(r12)  
            #xsd[s1,s2] <- sd(r12)
            #xmad[s1,s2] <- mad(r12,xmed[s1,s2])
            #xn[s1,s2] <- length(r12)
          }
        }
      }
      return(col_ret)
    }
    stopCluster(cl)
    gv$log <<- paste0(isolate(gv$log),flog.info("Stopped Cluster. Extracting parallel matrices."))
    xmean <- master_return[,(1 + 0 * n_reference_subsets) : (1 * n_reference_subsets)]
    xmed  <- master_return[,(1 + 1 * n_reference_subsets) : (2 * n_reference_subsets)]
    xsd   <- master_return[,(1 + 2 * n_reference_subsets) : (3 * n_reference_subsets)]
    xmad  <- master_return[,(1 + 3 * n_reference_subsets) : (4 * n_reference_subsets)]
    xn    <- master_return[,(1 + 4 * n_reference_subsets) : (5 * n_reference_subsets)]
    gv$log <<- paste0(isolate(gv$log),flog.info("Done"))
    colnames(xmean) <- reference_subset_names
    rownames(xmean) <- subset_names
    colnames(xmed) <- reference_subset_names
    rownames(xmed) <- subset_names
    colnames(xsd) <- reference_subset_names
    rownames(xsd) <- subset_names
    colnames(xmad) <- reference_subset_names
    rownames(xmad) <- subset_names
    colnames(xn) <- reference_subset_names
    rownames(xn) <- subset_names
    if(make_distributions==TRUE) {
      gv$log <<- paste0(isolate(gv$log),flog.info("Writing ratio test references"))
      xmadcv <- xmad/xmed
      xcv <- xsd / xmean
      write.table(xcv,cv_file,sep="\t",quote = FALSE)
      write.table(xmean,mean_file,sep="\t",quote = FALSE)
      write.table(xmadcv,madcv_file,sep="\t",quote = FALSE)
      write.table(xmed,med_file,sep="\t",quote = FALSE)
      write.table(xn,n_file,sep="\t",quote = FALSE)
      file_remove_lock(n_file)
    }
  } 
  if (make_distributions==FALSE ){
    gv$log <<- paste0(isolate(gv$log),flog.info("Reading ratio test references"))
    dist_mean <- read.table(mean_file,sep="\t",header=TRUE,stringsAsFactors=FALSE)
    dist_cv <- read.table(cv_file,sep="\t",header=TRUE,stringsAsFactors=FALSE)
    gv$log <<- paste0(isolate(gv$log),flog.info("Doing ratio test"))
    # Reference can have extra stuff
    if(ncol(dist_mean) != ncol(xmean)) { dist_mean <- dist_mean[,colnames(xmean)]}
    if(ncol(dist_cv) != ncol(xmean)) { dist_cv <- dist_cv[,colnames(xmean)]}
    
    zscore <- (xmean - dist_mean)/dist_mean/dist_cv
    n <- (xmean/dist_mean) * 2.0
    row_idx <- rowSums(!is.na(zscore))>0  
    col_idx <- colSums(!is.na(zscore))>0  
    zscore_med <- rowMedians(as.matrix(zscore[row_idx,col_idx]),na.rm = TRUE)
    n_med <- rowMedians(as.matrix(n[row_idx,col_idx]),na.rm = TRUE)
    df <- rbind(zscore_med,n_med)
    colnames(df) <- rownames(n[row_idx,col_idx])
    return(data.frame(df))
  }
}
# -------------------------------------------------------------------

# -------------------------------------------------------------------
make_deletions_zscores <- function(sample,bedtype) {
  if(grepl(CFG$del_zscore_section, bedtype)==FALSE) {
    gv$log <<- paste0(isolate(gv$log),flog.info("Skipping %s",bedtype))
    return(NULL)
  }
  delpath <- makepath(panel = CFG$panel_path ,run = CFG$run_path ,sample = sample ,filetype = "sv",resolution = bedtype, extension = "_del_zscore.tsv",df_cfg = DF_CFG)
  df_files <- data.frame(
    sample = sample,
    path = delpath,
    file = "counts",
    resolution = "off_target", # DANGER: Temporary hack to put it in the browser
    is_targeted = grepl("target",CFG$del_zscore_section),
    is_wgs = !grepl("target",CFG$del_zscore_section),
    stringsAsFactors = FALSE
  )
  
  if(!file_no_clobber(delpath) && file_apply_lock(delpath)) {
    gv$log <<- paste0(isolate(gv$log),flog.info("Making deletions zscores: %s.",delpath))
    cnvpath <- makepath(panel = CFG$panel,run = CFG$run_path,sample = sample,filetype = "cnv",resolution = bedtype, extension = ".tsv",df_cfg = DF_CFG)
    countspath <- sub(".tsv$","_fc_counts.tsv",cnvpath)
    gv$log <<- paste0(isolate(gv$log),flog.info("Reading %s",countspath))
    df <- make_deletions_stats(FALSE,countspath,bedtype)
    df_cn <- data.frame(t(df))
    df_cn$id <- rownames(df_cn)
    df_whitelist <- read.table(CFG$common_deletion_list,header=TRUE,stringsAsFactors=FALSE)
    df_gene <- df_whitelist[grepl("_ex1$",df_whitelist$Gene),]
    df_gene$Gene <- sub("_ex1","",df_gene$Gene )
    chrs <- paste0("chr",c(1:22,"X","Y"),c(rep("",24),rep("p",24),rep("q",24)))
    df_chr <- data.frame(row.names = chrs,
                         chr = c(chrs[1:24],chrs[1:24],chrs[1:24]),
                         start = c(rep(1,24),rep(1,24),as.integer(unlist(CFG$pq_boundary[chrs[1:24]]))),
                         end = c(rep(1,24),rep(1,24),as.integer(unlist(CFG$pq_boundary[chrs[1:24]]))),    
                         stringsAsFactors = F)
    df_regions <- rbind(df_whitelist,df_gene)
    rownames(df_regions) <- df_regions$Gene
    df_regions <- df_regions[,c("Chr","TestStart","TestEnd")]
    colnames(df_regions) <- c("chr","start","end")
    df_regions <- rbind(df_chr,df_regions)
    df_regions$id <- rownames(df_regions)
    df <- merge(df_regions,df_cn,by = "id")
    df$N <- df$n_med
    df$Nsd <- (df$N-2) / df$zscore_med
    df$Nlow <- ifelse(df$N>df$Nsd,df$N-df$Nsd,0)
    df$Nhigh <- df$N+df$Nsd
    rownames(df) <- df$id
    df <- df[,c("chr","start","end","zscore_med","n_med","Nsd","Nlow","Nhigh","id")]
    colnames(df) <- c("chr","start","end","zscore","N","Nsd","Nlow","Nhigh","id")
    write.table(df,delpath,row.names=FALSE,sep="\t",quote = FALSE)
    file_remove_lock(delpath)
    gv$log <<- paste0(isolate(gv$log),flog.info("Making zscore plots"))
    plot_deletions_zscores(df,delpath)
    return (df_files)
  } else {
    gv$log <<- paste0(isolate(gv$log),flog.info("Already have deletion zscores for %s",delpath))
    return (df_files)
  }
}
# -------------------------------------------------------------------

# -------------------------------------------------------------------
plot_deletions_zscores <- function(dfz,delpath) {
  sample_name <- sub(".*/1","1",sub("_del_zscore.tsv","",delpath))
  dfz$bintype <- ifelse(grepl("_ex",rownames(dfz)),"exon","gene")
  # dfz$bintype[grepl("chr",rownames(dfz))] <- "chr"
  dfz$zsdr <-  dfz$zscore/(dfz$N-2)
  dfz$sd <-  (dfz$N-2)/dfz$zscore
  dfz$id <- rownames(dfz)
  z_sd <- sd(dfz$zscore)
  median_sd <- median(dfz$sd)
  
  dfz$outlier_z <- ((dfz$zscore > 2 * z_sd) | (dfz$zscore < -2 * z_sd))
  dfz$outlier_sd <-  dfz$sd > median_sd * 2
  dfz$outlier <- dfz$outlier_sd | dfz$outlier_z
  dfz$bintype <- "ok"
  dfz$bintype[dfz$outlier_sd] <- "sd"
  dfz$bintype[dfz$outlier_z] <- "z_score"
  dfout <- dfz[dfz$outlier,]

  dfzpath <- sub("_del_zscore.tsv","_plotted_del_zscore.tsv",delpath)
  write.table(dfz,dfzpath,row.names=TRUE,col.names = NA,sep="\t",quote = FALSE)
  
  p <- ggplot() + 
    geom_point(data=dfz,aes(x=zscore,y=N-2,colour=bintype),size = 1) +
    ggtitle(sample_name) +
    xlim(c(-8,8)) +
    ylim(c(-2,4))
    # geom_text_repel(data = dfout, aes(x=zscore_med,y=n_med-2, label =id), size = 3)
  pngfile <- sub("_del_zscore.tsv","_delta_n_vs_zscore.png",delpath)
  ggsave(pngfile)

  p <- ggplot() + 
    geom_point(data=dfz,aes(x=zscore,y=sd,colour=bintype),size = 1) +
    ggtitle(sample_name) +
    xlim(c(-8,8)) +
    ylim(c(0,2)) 
    # geom_text_repel(data = dfout, aes(x=zscore,y=sd, label = id), size = 3)
  pngfile <- sub("_del_zscore.tsv","_sd_vs_zscore.png",delpath)
  ggsave(pngfile)
  
}
# -------------------------------------------------------------------

# -------------------------------------------------------------------
if(CFG_DEFAULTS$rstudio_debug) { 
  ret <- makeCNTables(args) 
  return(ret)
}
# -------------------------------------------------------------------
if(CFG_DEFAULTS$is_command_line) {
  ret <- makeCNTables(args) 
  quit(save = "no", status = ret)
}
# -------------------------------------------------------------------
