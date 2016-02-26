
ui <- fluidPage(
  # Some custom CSS
  tags$head(
    h3("Configure Copy Number Builder"),
    tags$style(HTML("
        /* Smaller font for preformatted text */
        pre, table.table {
          font-size: smaller;
        }

        body {
          min-height: 2000px;
        }

        .option-group {
          border: 1px solid #ccc;
          border-radius: 6px;
          padding: 0px 5px;
          margin: 5px -10px;
          background-color: #f5f5f5;
        }

        .option-header {
          color: #79d;
          text-transform: uppercase;
          margin-bottom: 5px;
        }
      "))
  ),
  
  fluidRow(
    div(
      fileInput("input_file_button", label = "Input File", multiple = FALSE, accept = NULL, width = NULL)
      # style = "font-size:60%"
    )
  ),
  fluidRow(
    #uiOutput("section_1_above"),
    #div(
      uiOutput("section_1")
      # checkboxInput("no_clobber", label = "No Clobber", value = TRUE, width = NULL),
      # textInput("number_of_samples", label = "Number of Samples", value = "1", width = NULL),
      # textInput("targeted_label", label = "Targeted Label", value = "Targeted", width = NULL),
      # textInput("wg_label", label = "Whole Geneome Label", value = "Whole Geneome", width = NULL),
      # actionButton("done_section_1", "Next",  icon = icon("arrow-right")),
      # style = "font-size:60%"
    #),
    #uiOutput("section_1_below")
  ),
  fluidRow(
    #div(
      uiOutput("section_2")
      # style = "font-size:60%"
    #)
  ),
  fluidRow(
    #div(
      uiOutput("section_3")
      # style = "font-size:60%"
    #)
  ),
  fluidRow(
    div(
      uiOutput("section_4")
      # style = "font-size:60%"
    )
  ),
  fluidRow(
    div(
      uiOutput("section_5")
      # style = "font-size:60%"
    )
  ),
  fluidRow(
    div(
      hr(),
      h4("Logging Output"),
      verbatimTextOutput("log_message")
      # style = "font-size:60%"
    )
  )
)
