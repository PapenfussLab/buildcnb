# library(Cairo)   # For nicer ggplot2 output when deployed on Linux
library(shiny)
library(futile.logger)
library(RCurl)
library(httr)

# TODO: Test with new cnb
# TODO: validateInputs() is called too much. Also, say what it found and validate - URL or file
# TODO: BUG: Problem with forward button AFTER back and no changes
# TODO: Offset everything from LH page margin
# TODO: Help links?

server <- function(input, output, session) {
  readSessionInfo <- reactive({
    init_table <- list()
    cdata <- session$clientData
    cnames <- names(cdata)
    if(nchar(cdata$url_search)>2) {
      query <- parseQueryString(cdata$url_search)
      print(paste(names(query), query, sep = " = ", collapse=", "))
      if(!is.null(query$buildfile)) {
        init_table$buildfile <- query$buildfile
      }
    }
    init_table
  })
  
  # print("isolate(readSessionInfo())")
  session_info <- isolate(readSessionInfo())
  progress <- shiny::Progress$new(session, min=0, max=1)
#   if(!is.null(session_info$buildfile)) {
#     progress$set(message = 'Preparing for build', detail = '...')
#     source(session_info$buildfile)
#     buildcnb(session = session)
#     return
#   }
  # gv <- init_globals(session_info,session)
  progress$set(message = 'Preparing for display', detail = '...')
  return  

  # -------------------------------------------------------------------
  # reactiveValues
  # -------------------------------------------------------------------
  gv <- reactiveValues(
    read_defaults_from_file = FALSE,
    no_clobber = TRUE,
    input_file = "",
    config_file = "",
    number_of_samples = 2,
    number_of_cores = 1,
    wg_label = "Targeted",
    targeted_label = "Whole Genome",
    current_section = 1,
    input_path = "",
    output_path = "",
    targeted_seq_bedfile = "",
    off_target_bedfile = "",
    wg_coarse_bedfile = "",
    wg_medium_bedfile = "",
    wg_fine_bedfile = "",
    wg_annotation_file = "",
    targeted_annotation_file = "",
    cytobands_file = "",
    chrom_info_file = "",
    chrom_file = "",
    min_samples = 2,
    max_samples = 100,
    min_cores = 1,
    max_cores = 32,
    name = vector("list",100),
    label = vector("list",100),
    wg_bamfile =  vector("list",100),
    targeted_bamfile =  vector("list",100),
    off_target_bamfile =  vector("list",100),
    log = "" # this is printed on the screen to indicate progress/errors
  )

  # Give it the local environment because it needs to modify gv
  source("makeCNTables.R",local=TRUE)

  # -------------------------------------------------------------------
  validateInputs <- function(table,current_section) {
  # -------------------------------------------------------------------
  # updateNumericInput(session, inputId, label = NULL, value = NULL, min = NULL, max = NULL, step = NULL)
  # updateCheckboxInput(session, inputId, label = NULL, value = NULL)
  # updateTextInput(session, inputId, label = NULL, value = NULL)
    ret <- TRUE
    if(current_section>=2) {
      ret <- ret && dir_exists_or_can_be_made(gv$input_path)
      ret <- ret && dir_exists_or_can_be_made(gv$output_path)
      ret <- ret && dir_exists_or_can_be_made(gv$reference_path)
      
      ret <- ret && url_or_file_exists(gv$targeted_bedfile,gv$reference_path)
      ret <- ret && url_or_file_exists(gv$off_target_bedfile,gv$reference_path)
      ret <- ret && url_or_file_exists(gv$wg_annotation_file,gv$reference_path)
      ret <- ret && url_or_file_exists(gv$targeted_annotation_file,gv$reference_path)
      ret <- ret && url_or_file_exists(gv$cytobands_file,gv$reference_path)
      ret <- ret && url_or_file_exists(gv$chrom_info_file,gv$reference_path)
      ret <- ret && url_or_file_exists(gv$chrom_file,gv$reference_path)
    }
  
    if(current_section>=3) {
      bam_files <- table[grepl("bams",table$section),]            
      for (i in 1:gv$number_of_samples){
        query <-  paste0("_",i,"$")
        idx <- grepl(query,table$name)
        for(filename in table$value[idx]) {
          ret <- ret && url_or_file_exists(filename,gv$input_path)
        }
      }
    }
  return(ret)
  }

  # -------------------------------------------------------------------
  url_or_file_exists <- function(s,spath) {
  # -------------------------------------------------------------------
    ret <- TRUE
    if(isURL(s)) {
      url <- parse_url(s)
      filename <- basename(url$path)
      if((file.exists(file.path(spath,filename)) && gv$no_clobber==TRUE) || url.exists(s)) {
        gv$log <- paste0(isolate(gv$log),flog.info("Found %s",s))          
        return(ret)
      }  
    } else if(file.exists(file.path(spath,s))) {
      gv$log <- paste0(isolate(gv$log),flog.info("Found %s",s))          
      return(ret)
    } else if(file.exists(s)) {
      gv$log <- paste0(isolate(gv$log),flog.info("Found %s",s))          
      return(ret)
    }   
    gv$log <- paste0(isolate(gv$log),flog.error("Unable to find %s",s)) 
    ret <- FALSE
    return(ret)
  }

  # -------------------------------------------------------------------
  dir_exists_or_can_be_made <- function(s) {
  # -------------------------------------------------------------------
    if(file.exists(s)) {
      gv$log <- paste0(isolate(gv$log),flog.info("Found %s",s))
      return(TRUE)
    }
    ret <- ftry(dir.create(s))
    if(ret != TRUE) {
      gv$log <- paste0(isolate(gv$log),flog.error("Unable to create directory for config file: %s",configDir))
      return(FALSE)
    }
    gv$log <- paste0(isolate(gv$log),flog.info("Made %s",s))
    return(TRUE)
  }

  # -------------------------------------------------------------------
  writeConfigFile <- function(table,filename) {
  # -------------------------------------------------------------------
    ret <- TRUE
    configDir <- dirname(filename)
    if(dir_exists_or_can_be_made(configDir)==FALSE) {
      return (FALSE)
    }
    bam_idx <- grepl("bams",table$section)
    keep_idx <- !bam_idx
    for (i in 1:gv$number_of_samples) {
      keep_idx <- keep_idx | (grepl(paste0("_",i,"$"),table$name) & bam_idx)
    }
    table <- table[keep_idx,]
    ret <- ftry(write.dcf(table,file = filename))
    if(!is.null(ret)) {
      gv$log <- paste0(isolate(gv$log),flog.error("Unable to write config file: %s",filename))
      return(FALSE)
    } 
    gv$log <- paste0(isolate(gv$log),flog.info("Wrote %s",filename))
    return(TRUE)
  }

  # -------------------------------------------------------------------
  tableDefaults <- function() {
  # -------------------------------------------------------------------
  # Set defaults values for fields
    number_of_samples <- gv$max_samples
    params <- NULL
    reference_files <- NULL
    bam_files <- NULL
  
    names <- c("no_clobber","number_of_cores","number_of_samples","wg_label","targeted_label")
    labels <- c("No Clobber","Number of CPU Cores","Number of Samples","Whole Geneome Label","Targeted Label")
    default_values <- c(TRUE,gv$min_cores,gv$min_samples,"WG","Targeted")
    params <- data.frame(section = "params", name = names,label = labels, value = default_values,stringsAsFactors = FALSE)
    
    names <- c("output_path",
             "input_path",
             "config_file",
             "reference_path",
             "input_file",
             "wg_coarse_bedfile",
             "wg_medium_bedfile",
             "wg_fine_bedfile",
             "targeted_bedfile",
             "off_target_bedfile",
             "wg_annotation_file",
             "targeted_annotation_file",
             "gc_fitting_blacklist",
             "display_blacklist",
             "cytobands_file",
             "chrom_info_file",
             "chrom_file")
    default_values <- c("batch_name",
                      "input_path",
                      "config_file.dcf",
                      "panel_name",
                      "input_file.txt",
                      "panel_name/hg19_wg_coarse.bed",
                      "panel_name/hg19_wg_medium.bed",
                      "panel_name/hg19_wg_fine.bed",
                      "panel_name/panel_targeted_bins.bed",
                      "panel_name/panel_off_target_bins.bed",
                      "all_panels/hg19_exons.tsv",
                      "panel_name/targeted_annotation.bed",
                      "panel_name/panel_excluded_from_gc_correction.bed",
                      "panel_name/panel_excluded_from_display.bed",
                      "all_panels/cytoBand.txt.gz",
                      "all_panels/chromInfo.txt.gz",
                      "all_panels/hg19.fa")
    labels <- c("Output path",
              "Input path",
              "Config file",
              "Reference path",
              "Input file",
              "WG_coarse_bedfile",
              "WG_medium_bedfile",
              "WG fine bedfile",
              "Targeted bedfile",
              "Off Target bedfile",
              "WG annotation bedfile",
              "Targeted annotation bedfile",
              "GC Fitting Blacklist",
              "Display Blacklist",
              "Cytobands file",
              "Chromosome info file",
              "Chromosome sequence file")
    reference_files <- data.frame(section = "files", name = names,label = labels, value = default_values, stringsAsFactors = FALSE)
  
    wg_bam_files <- data.frame(section = "wg_bams", 
                            name = paste0("wg_sample_",1:number_of_samples),
                            label = paste0("Sample Name ",1:number_of_samples), 
                            value = paste0("sample_wg_",1:number_of_samples,".bam"),
                            stringsAsFactors = FALSE)
    targeted_bam_files <- data.frame(section = "targeted_bams", 
                             name = paste0("targeted_sample_",1:number_of_samples),
                             label = paste0("Sample Name ",1:number_of_samples), 
                             value = paste0("sample_targeted_",1:number_of_samples,".bam"),
                             stringsAsFactors = FALSE)
    off_target_bam_files <- data.frame(section = "off_target_bams", 
                                   name = paste0("off_target_sample_",1:number_of_samples),
                                   label = paste0("Sample Name ",1:number_of_samples), 
                                   value = paste0("sample_off_target_",1:number_of_samples,".bam"),
                                   stringsAsFactors = FALSE)
  
    targeted_bam_files$label[1] <- "HG00096"
    targeted_bam_files$value[1] <- "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00096/exome_alignment/HG00096.chrom20.ILLUMINA.bwa.GBR.exome.20120522.bam"
    off_target_bam_files$label[1] <- "HG00096"
    off_target_bam_files$value[1] <- "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00096/exome_alignment/HG00096.chrom20.ILLUMINA.bwa.GBR.exome.20120522.bam"
    wg_bam_files$label[1] <- "HG00096"
    wg_bam_files$value[1] <- "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00096/alignment/HG00096.chrom20.ILLUMINA.bwa.GBR.low_coverage.20120522.bam"
    targeted_bam_files$label[2] <- "HG00097"
    targeted_bam_files$value[2] <- "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00097/exome_alignment/HG00097.chrom20.ILLUMINA.bwa.GBR.exome.20120522.bam"
    off_target_bam_files$label[2] <- "HG00097"
    off_target_bam_files$value[2] <- "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00097/exome_alignment/HG00097.chrom20.ILLUMINA.bwa.GBR.exome.20120522.bam"
    wg_bam_files$label[2] <- "HG00097"
    wg_bam_files$value[2] <- "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00097/alignment/HG00097.chrom20.ILLUMINA.bwa.GBR.low_coverage.20120522.bam"
    toTable <- rbind(params,reference_files,wg_bam_files,targeted_bam_files,off_target_bam_files)
    rownames(toTable) <- toTable$name
    # Put in defaults from the file if there is a slot for them
    if(gv$read_defaults_from_file) {
      fromTable <- gv$defaults_from_file
      idx <- fromTable$name %in% toTable$name
      toTable[fromTable$name[idx],]  <- fromTable[idx,]
    }
    return(toTable)
  }
  
  # -------------------------------------------------------------------
  tableFromInputs <- function(input) {
  # -------------------------------------------------------------------
  # Creates a table from currently visible fields
    table <- tableDefaults()
    bam_idx <- grepl("bams",table$section)
    not_bams <- table[!bam_idx,]            
    for (s in not_bams$name){
      v <- input[[s]]
      if(!is.null(v)) {
        table$value[table$name==s] <- v
      }
    }
    bams <- table[bam_idx,] 
  
    # number_of_samples <- sum(bam_idx)
    for (s in bams$name) {
      # s <- bams$name[i] 
      v <- input[[s]]
      if(!is.null(v)) {
        table[table$name==s,"value"] <- v
      }
    }
    return(table)
  }

  # -------------------------------------------------------------------
  inputsFromTable <- function(table) {
  # -------------------------------------------------------------------
  # Set default fields from table
  #
    idx_bams <- grepl("bams",table$section)
    bam_files <- table[idx_bams,]
    gv$name  <- bam_files$name
    gv$label  <- bam_files$label
    gv$value  <- bam_files$value
    #gv$wg_bamfile  <- bam_files$wg_bamfile
    #gv$targeted_bamfile  <- bam_files$targeted_bamfile
    not_bam_files <- table[!idx_bams,]
    for (s in not_bam_files$name){
      gv[[s]] <- not_bam_files[s,"value"]
    }
  # updateNumericInput(session, inputId, label = NULL, value = NULL, min = NULL, max = NULL, step = NULL)
  # updateCheckboxInput(session, inputId, label = NULL, value = NULL)
  # updateTextInput(session, inputId, label = NULL, value = NULL)
  }

# -------------------------------------------------------------------
# observers
# -------------------------------------------------------------------
#   observeEvent(input$no_clobber, {
#       gv$no_clobber <- input$no_clobber
#     })

  observeEvent(input$input_file_button, {
    gv$log <- ""
    # browser()
    filename <- input$input_file_button$datapath
    table <- ftry(read.dcf(filename,all=TRUE))
    if(!is.data.frame(table)) {
      gv$log <- paste0(isolate(gv$log),flog.error("Unable to read config file %s stored locally as %s",input$input_file_button$name,filename))
      return
    }
    gv$log <- paste0(isolate(gv$log),flog.info("Successfully read config file %s stored locally as %s",input$input_file_button$name,filename))
    gv$read_defaults_from_file <- TRUE
    rownames(table) <- table$name

    gv$defaults_from_file <- table
    inputsFromTable(table)
    isValid <- validateInputs(input,1)
  })

  observeEvent(input$number_of_samples, {
    if(gv$current_section==1) {
      num_samples <- as.integer(input$number_of_samples)
      # TODO: Check number
      gv$number_of_samples <- num_samples
    } else {
      updateNumericInput(session, inputId = "number_of_samples", value = gv$number_of_samples )
    }
  })

  observeEvent(input$section_next, {
    print(paste0("section_next: ",gv$current_section))
    gv$log <- ""
    table <- tableFromInputs(input)
    if(gv$current_section>=1) {
      params <- table[table$section=="params",]            
      for (i in 1:nrow(params)){
        s <- params$name[i]
        gv[[s]] <- params$value[i]
      }
    } 
    if(gv$current_section>=1) {
      files <- table[table$section=="files",]            
      for (i in 1:nrow(files)){
        s <- files$name[i]
        gv[[s]] <- files$value[i]
      }
    } 
    if(gv$current_section>=2) {
      bam_files <- table[grepl("bams",table$section),]            
      gv$name <- bam_files$name
      gv$label <- bam_files$label
      gv$value <- bam_files$value
    }
    isValid <- validateInputs(table,gv$current_section)
    if(isValid==TRUE) {
      gv$current_section <- gv$current_section + 1
    }
    print(paste0("section_next: ",gv$current_section))
  })
  observeEvent(input$section_prev, {
    print(paste0("section_prev: ",gv$current_section))
    if(gv$current_section > 1) {
      gv$current_section <- gv$current_section - 1
    }
  })
  # -------------------------------------------------------------------
  print("Render UI")
  # -------------------------------------------------------------------
  # Render UI
  # -------------------------------------------------------------------
  output$section_1 <- renderUI({ 
    ret <- NULL
    if(gv$current_section==1) {
      ret <- list(h4("Settings"),hr())
    } else {
      ret <- list(h4("Settings"))
    }
    ret <- append(ret,list(div(
      column(2,offset = 0.1,checkboxInput("no_clobber", label = "No Clobber", value = as.logical(gv$no_clobber))),
      column(2,numericInput("number_of_samples", "Number of Samples",  value = gv$number_of_samples,min = gv$min_samples, max = gv$max_samples)),
      column(2,numericInput("number_of_cores", "Number of CPU Cores",  value = gv$number_of_cores,min = gv$min_cores, max = gv$max_cores)),
      column(3,textInput("targeted_label", label = "Targeted Label", value = gv$targeted_label)),
      column(3,textInput("wg_label", label = "Whole Geneome Label", value = gv$wg_label))
    )))
    if(gv$current_section==1) {
      ret <- append(ret,list(
        div(
          actionButton("section_prev", "Prev",  icon = icon("arrow-left")),
          actionButton("section_next", "Next",  icon = icon("arrow-right"))
        ),
        hr()
      ))
    }
    ret
  })

  output$section_2 <- renderUI({  
    ret <- NULL
    if(gv$current_section>=2) {
      ret <- list(h4("References"))
      if(gv$current_section==2) {
        ret <- append(ret,list(hr()))
      }
    }
    if(gv$current_section>=2) {
      table <- tableDefaults()
      files <- table[table$section=="files",]
      for (i in 1:nrow(files)) {
        label <- files$label[i]
        name <- files$name[i]
        ret <- append(ret,list(div(column(4,offset = 0.1,strong(label)),column(8,textInput(name, label = "", value = gv[[name]])))))
      }
    }
    if(gv$current_section==2) {
      ret <- append(ret,list(
        div(
          actionButton("section_prev", "Prev",  icon = icon("arrow-left")),
          actionButton("section_next", "Next",  icon = icon("arrow-right"))
        ),
        hr()
      ))
    }
    ret
  })

  output$section_3 <- renderUI({  
    ret <- NULL
    if(gv$current_section>=3) {
      ret <- list(h4("Samples"))
      if(gv$current_section==3) {
        ret <- append(ret,list(hr()))
      }
      for (i in 1:gv$number_of_samples) {
        sample_name_id <- paste0("sample_",i)
        wg_bamfile_id <- paste0("wg_sample_",i)
        targeted_bamfile_id <- paste0("targeted_sample_",i)
        off_target_bamfile_id <- paste0("off_target_sample_",i)
        wg_bamfile_val <- gv$value[gv$name==wg_bamfile_id]
        targeted_bamfile_val<- gv$value[gv$name==targeted_bamfile_id]
        off_target_bamfile_val<- gv$value[gv$name==off_target_bamfile_id]
        # TODO: Deal with situation when no targeted bam file
        sample_label <- gv$label[gv$name==targeted_bamfile_id]
        row <- div(
          column(3,offset = 0.1,textInput(inputId = sample_name_id, label = paste0("Sample Name ",i), value = sample_label)),
          column(3,textInput(inputId = wg_bamfile_id, label = paste0("WG Bamfile ",i), value = wg_bamfile_val)),
          column(3,textInput(inputId = targeted_bamfile_id, label = paste0("Targeted Bamfile ",i), value = targeted_bamfile_val)),
          column(3,textInput(inputId = off_target_bamfile_id, label = paste0("Off Target Bamfile ",i), value = off_target_bamfile_val))
        )
        ret <- append(ret,list(row))
      }
    }
    if(gv$current_section==3) {
      ret <- append(ret,list(
        div(
          actionButton("section_prev", "Prev",  icon = icon("arrow-left")),
          actionButton("section_next", "Next",  icon = icon("arrow-right"))
        ),
        hr()
      ))
    }
    ret
  })

  output$section_4 <- renderUI({  
    ret <- NULL
    if(gv$current_section==4) {
      table <- tableFromInputs(input)
      gv$dcf_filename <- file.path(input$output_path,input$config_file)
      writeConfigFile(table,gv$dcf_filename)
      ret <- append(ret,list(
      hr(),
      div(
        actionButton("section_prev", "Prev",  icon = icon("arrow-left")),
        actionButton("section_next", "Run Build",  icon = icon("gears"))
      ),
      hr()
      ))
    }
    ret
  })
  output$section_5 <- renderUI({  
    ret <- NULL
    if(gv$current_section==5) {
      ret <- list(
        hr(),
        h5("Preparing all files. This could take a long time so please be patient - and don't click on anything!")
      )
      makeCNTables(gv$dcf_filename) 
    }
    ret
  })
  output$log_message <- renderText({  
     if(nchar(gv$log) > 0) {
       gv$log
     } else {
       NULL
     }
 })

# -------------------------------------------------------------------
}
