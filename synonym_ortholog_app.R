## ============================================================
## synonym_ortholog_app.R
## Fast Gene Synonym / Ortholog Lookup (DepMap anchored)
## ============================================================

options(shiny.maxRequestSize = 200 * 1024^2)

library(shiny)
library(data.table)
library(DT)

## ============================================================
## UI
## ============================================================
ui <- fluidPage(
  
  titlePanel("Gene Synonym / Ortholog Lookup (DepMap anchored)"),
  
  sidebarLayout(
    sidebarPanel(
      
      fileInput(
        inputId = "ref_file",
        label   = "Upload reference file (depmap_gene_reference_clean_v8.tsv)",
        accept  = ".tsv"
      ),
      
      tags$hr(),
      
      textAreaInput(
        inputId = "genes",
        label   = "Paste gene names (one per line or space/comma separated)",
        value   = "",
        rows    = 10,
        placeholder = "tp53\ntrp53\nil4r\nil4ra\ncd124"
      ),
      
      actionButton(
        inputId = "search",
        label   = "Search",
        class   = "btn-primary"
      ),
      
      tags$hr(),
      
      strong("Status:"),
      verbatimTextOutput("status")
    ),
    
    mainPanel(
      DTOutput("results")
    )
  )
)

## ============================================================
## SERVER
## ============================================================
server <- function(input, output, session) {
  
  ref_data <- reactiveVal(NULL)
  indices  <- reactiveVal(NULL)
  status   <- reactiveVal("Waiting for reference file.")
  
  output$status <- renderText(status())
  
  ## ------------------------------------------------------------
  ## Helpers
  ## ------------------------------------------------------------
  explode_column <- function(dt, column) {
    
    tmp <- dt[
      ,
      .(
        depmap_symbol = depmap_symbol,
        key = tolower(unlist(strsplit(get(column), ",", fixed = TRUE)))
      ),
      by = seq_len(nrow(dt))
    ]
    
    tmp <- tmp[
      !is.na(key) & key != "",
      .(depmap_symbol = unique(depmap_symbol)),
      by = key
    ]
    
    tmp
  }
  
  ## ------------------------------------------------------------
  ## Load + index reference (FAST, vectorized)
  ## ------------------------------------------------------------
  observeEvent(input$ref_file, {
    
    req(input$ref_file$datapath)
    
    status("Reading reference file...")
    
    ref <- tryCatch(
      fread(input$ref_file$datapath, sep = "\t", showProgress = FALSE),
      error = function(e) NULL
    )
    
    if (is.null(ref)) {
      status("ERROR: Could not read reference file.")
      return()
    }
    
    required_cols <- c(
      "depmap_symbol",
      "all_names_low_risk",
      "all_names_mid_risk",
      "all_names_high_risk"
    )
    
    if (!all(required_cols %in% colnames(ref))) {
      status("ERROR: File is not a valid v8 reference (missing required columns).")
      return()
    }
    
    status(
      paste("Reference loaded (", nrow(ref),
            " genes). Building indices...", sep = "")
    )
    
    ## ---- FAST index construction ------------------------------
    idx_depmap <- setNames(
      as.list(ref$depmap_symbol),
      tolower(ref$depmap_symbol)
    )
    
    low_dt  <- explode_column(ref, "all_names_low_risk")
    mid_dt  <- explode_column(ref, "all_names_mid_risk")
    high_dt <- explode_column(ref, "all_names_high_risk")
    
    idx_low  <- split(low_dt$depmap_symbol,  low_dt$key)
    idx_mid  <- split(mid_dt$depmap_symbol,  mid_dt$key)
    idx_high <- split(high_dt$depmap_symbol, high_dt$key)
    
    indices(list(
      depmap = idx_depmap,
      low    = idx_low,
      mid    = idx_mid,
      high   = idx_high
    ))
    
    ref_data(ref)
    
    status(
      paste("Reference loaded (", nrow(ref),
            " genes). Ready to search.", sep = "")
    )
  })
  
  ## ------------------------------------------------------------
  ## Search
  ## ------------------------------------------------------------
  results <- eventReactive(input$search, {
    
    if (is.null(indices())) {
      showNotification(
        "Reference not loaded yet. Upload v8 file first.",
        type = "error"
      )
      return(data.table())
    }
    
    idx <- indices()
    
    queries <- tolower(unlist(strsplit(
      input$genes, "[\n\r\t ,]+"
    )))
    queries <- unique(queries[queries != ""])
    
    if (length(queries) == 0) {
      return(data.table())
    }
    
    rbindlist(lapply(queries, function(q) {
      
      depmap_hits <- idx$depmap[[q]]
      low_hits    <- idx$low[[q]]
      mid_hits    <- idx$mid[[q]]
      high_hits   <- idx$high[[q]]
      
      found_in <- c(
        if (!is.null(depmap_hits)) "depmap" else NULL,
        if (!is.null(low_hits))    "low_risk" else NULL,
        if (!is.null(mid_hits))    "mid_risk" else NULL,
        if (!is.null(high_hits))   "high_risk" else NULL
      )
      
      data.table(
        query_gene        = q,
        depmap_hit        = if (!is.null(depmap_hits)) paste(depmap_hits, collapse = ",") else "",
        low_risk_matches  = if (!is.null(low_hits))    paste(sort(low_hits),  collapse = ",") else "",
        mid_risk_matches  = if (!is.null(mid_hits))    paste(sort(mid_hits),  collapse = ",") else "",
        high_risk_matches = if (!is.null(high_hits))   paste(sort(high_hits), collapse = ",") else "",
        found_in          = paste(found_in, collapse = ",")
      )
    }))
  })
  
  ## ------------------------------------------------------------
  ## Render table
  ## ------------------------------------------------------------
  output$results <- renderDT({
    datatable(
      results(),
      options = list(
        pageLength = 25,
        scrollX = TRUE
      )
    )
  })
}

## ============================================================
## Run app
## ============================================================
shinyApp(ui = ui, server = server)
