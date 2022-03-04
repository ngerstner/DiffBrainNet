# module UI function
comparisonHubGenesUI <- function(id, label = "Upset Plot Hub Genes"){
  
  ns <- NS(id)
  
  #tabPanel("Comparison Brain Regions",
  tagList(
    sidebarPanel(
      pickerInput(
        inputId = ns("myPicker"),
        label = "Select brain regions",
        choices = c("Amygdala" = "AMY",
                    "Cerebellum" = "CER",
                    "Dorsal DG of Hippocampus" = "dDG",
                    "Dorsal CA1 of Hippocampus" = "dCA1",
                    "Prefrontal Cortex" = "PFC",
                    "PVN of Hypothalamus" = "PVN",
                    "Ventral DG of Hippocampus" = "vDG",
                    "Ventral CA1 of Hippocampus" = "vCA1"),
        options = list(
          `actions-box` = TRUE,
          size = 8,
          `selected-text-format` = "count > 3"
        ),
        multiple = TRUE
      ),
      radioButtons(
        ns("network_type"),
        label = "Select type of network to display:",
        choiceNames = list("Baseline", "Differential", "Treatment"),
        choiceValues = list("baseline", "differential", "treatment"),
        selected = "differential"
      ),
      sliderInput(
        inputId = ns("normBetween"),
        label = "Choose a threshold for the normalized nodebetweenness:",
        min = 0.5,
        max = 1.5,
        step = 0.05,
        value = 1.0
      ),
      sliderInput(
        ns("nintersects"),
        label = "Number of intersections",
        min = 2,
        max = 40,
        step = 1,
        value = 20
      ),
      downloadButton(ns("download"),"Download (filtered) table as csv"),
      width = 3
    ),
    mainPanel(
      br(),
      tags$style(".fa-brain {color:#2980B9}"),
      h3(p(em("Comparsion of hub genes across brain regions "),icon("brain",lib = "font-awesome"),style="color:black;text-align:center")),
      hr(),
      
      tags$style(type="text/css",
                 ".shiny-output-error { visibility: hidden; }",
                 ".shiny-output-error:before { visibility: hidden; }"
      ),
      strong("Compare"),span(" hub genes either at baseline, after Dexamethasone treatment or at the differential level"),
      strong(" between different brain regions"), span(". Please choose the brain regions you are interested 
             in on the left panel. You can also select a different threshold for the normalized
             nodebetweenness. The table on the bottom displays all genes in the dataset with the brain regions 
             they are hub genes in."),
      br(),br(),
      plotlyOutput(ns("plotly_upset"), height = "600px"),
      DT::dataTableOutput(ns("upset_table"))
      #h3("Intersection plots")
    ))
  
}

# module server function
comparisonHubGenesServer <- function(id) {
  moduleServer(
    id,
    
    function(input, output, session) {
      
      # Read data required for upset plot
      read_upset_data <- function(basepath, regions, between_cutoff = 1.0){
        
        # Read DE tables from all required regions
        list_genes_sig <- list()
        
        for (reg in regions){
          f <- Sys.glob(file.path(paste0(basepath, "tables/coExpression_kimono/",
                                         "03_AnalysisFuncoup/"), 
                                  paste0("*_", reg, "_funcoup_", input$network_type,
                                         "_nodebetweennessNorm_betacutoff0.01.csv")))
          res <- fread(file=f, sep=",")
          #na_indices <- which(is.na(res$padj))
          #res$padj[na_indices] <- 1
          res_sig <- res[res$nodebetweenness_norm >= between_cutoff,]
          list_genes_sig[[reg]] <- res_sig$ensembl_id
        }
        
        return(list_genes_sig)
      }
      
      # Read data required for data table below upset plot
      upset_datatable <- function(basepath, regions, between_cutoff = 1.0){
        
        # READ DATA 
        
        # Read DE tables from all required regions
        list_reg_sig <- list()
        
        for (reg in regions){
          f <- Sys.glob(file.path(paste0(basepath, "tables/coExpression_kimono/",
                                         "03_AnalysisFuncoup/"), 
                                  paste0("*_", reg, "_funcoup_", input$network_type,
                                         "_nodebetweennessNorm_betacutoff0.01.csv")))
          res <- fread(file=f, sep=",")
          #na_indices <- which(is.na(res$padj))
          #res$padj[na_indices] <- 1
          res_sig <- res[res$nodebetweenness_norm >= between_cutoff,]
          list_reg_sig[[reg]] <- res_sig
        }
        
        df <- bind_rows(list_reg_sig, .id="region")
        df <- df[,c("ensembl_id", "region")] %>% 
          group_by(ensembl_id) %>%
          dplyr::summarise(region = list(region)) %>%
          dplyr::select(ensembl_id, region)
        df$gene_symbol <- mapIds(org.Mm.eg.db, keys = df$ensembl_id,
                                 column = "SYMBOL", keytype = "ENSEMBL")
        df <- dplyr::select(df, gene_symbol, ensembl_id, region)
        
        return(df)
      }
      
      
      # reactive table 
      upset_data <- reactive({
        
        req(input$myPicker)
        data <- upset_datatable(basepath, input$myPicker, input$normBetween) %>%
          mutate(region = sapply(region, toString))
        data
        
      })
      
      
      # Data table below upset plot
      output$upset_table <- DT::renderDataTable({
        
        upset_data() %>%
          datatable(filter = list(position = 'top'))
      }, server = TRUE)
      
      # Download of (filtered) network table
      output$download <- downloadHandler(
        filename = "hubGenes.csv",
        content = function(file) {
          indices <- input$upset_table_rows_all
          print(upset_data())
          write.csv(upset_data()[indices, ], file)
        }
      )
      
      
      # UPSET PLOT 
      
      # TODO: mark here with a comment (SPDX) that following code snippet is 
      # under AGPL license https://github.com/pinin4fjords/shinyngs/blob/master/R/upset.R
      # --> whole repo needs to be AGPL but everything that is my own can be MIT
      
      # Accessor for user-selected sets
      
      getSelectedSetNames <- reactive({
        req(input$myPicker)
        input$myPicker
      })
      
      
      # Get the sets we're going to use based on nsets
      
      getSets <- reactive({
        selected_set_names <- getSelectedSetNames()
        req(length(selected_set_names) > 0)
        
        selected_sets <- read_upset_data(basepath, selected_set_names, input$normBetween)
      })
      
      
      # Accessor for the nintersections parameter
      
      getNintersections <- reactive({
        validate(need(!is.null(input$nintersects), "Waiting for nintersects"))
        input$nintersects
      })
      
      # Calculate intersections between sets
      
      calculateIntersections <- reactive({
        selected_sets <- getSets()
        
        withProgress(message = "Calculating set intersections", value = 0, {
          sets <- getSets()
          nsets <- length(sets)
          
          # Get all possible combinations of sets
          
          combinations <- function(items, pick) {
            x <- combn(items, pick)
            lapply(seq_len(ncol(x)), function(i)
              x[, i])
          }
          
          combos <- lapply(1:nsets, function(x) {
            combinations(1:length(selected_sets), x)
          })
          combos <- lapply(1:nsets, function(x) {
            combinations(1:nsets, x)
          })
          
          # Calculate the intersections of all these combinations
          
          withProgress(message = "Running intersect()", value = 0, {
            intersects <- lapply(combos, function(combonos) {
              lapply(combonos, function(combo) {
                Reduce(intersect, selected_sets[combo])
              })
            })
          })
          
          # For UpSet-ness, membership of higher-order intersections takes priority 
          # Otherwise just return the number of entries in each intersection
          
          intersects <- lapply(1:length(intersects), function(i) {
            intersectno <- intersects[[i]]
            if (i != length(intersects)){
              members_in_higher_levels <-
                unlist(intersects[(i + 1):length(intersects)])
            } else {
              members_in_higher_levels <- NULL
            }
            lapply(intersectno, function(intersect) {
              length(setdiff(intersect, members_in_higher_levels))
            })
          })
          
          combos <- unlist(combos, recursive = FALSE)
          intersects <- unlist(intersects)
          
          
          # Sort by intersect size
          
          combos <- combos[order(intersects, decreasing = TRUE)]
          intersects <-
            intersects[order(intersects, decreasing = TRUE)]
          list(combinations = combos, intersections = intersects)
          
        })
        
      })
      
      output$plotly_upset <- renderPlotly({
        grid <- upsetGrid()
        set_size_chart <- upsetSetSizeBarChart()
        intersect_size_chart <- upsetIntersectSizeBarChart()
        
        # add axis labels
        intersect_size_chart <-
          intersect_size_chart %>% layout(yaxis = list(title = "Intersections size"))
        
        set_size_chart <- 
          set_size_chart %>% layout(xaxis = list(title = "Set size"),
                                    yaxis = list(title = "Brain region"))
        
        # subplots
        s1 <-
          subplot(
            plotly_empty(type = "scatter", mode = "markers"),
            set_size_chart,
            nrows = 2,
            # widths = c(0.6, 0.4),
            titleX = TRUE,
            titleY = TRUE
          ) %>% layout(showlegend = FALSE)
        s2 <-
          subplot(intersect_size_chart,
                  grid,
                  nrows = 2,
                  shareX = TRUE, titleY = TRUE) %>% layout(showlegend = FALSE)
        
        subplot(s1, s2, widths = c(0.25, 0.75),
                shareX=F,
                shareY=F,
                titleX=T,
                titleY=T)
        
        
      })
      
      # Add some line returns to contrast names
      
      getSetNames <- reactive({
        selected_sets <- getSets()
        names(selected_sets)
      })
      
      # Make the grid of points indicating set membership in intersections
      
      upsetGrid <- reactive({
        selected_sets <- getSets()
        ints <- calculateIntersections()
        
        intersects <- ints$intersections
        combos <- ints$combinations
        
        # Reduce the maximum number of intersections if we don't have that many
        nintersections <- getNintersections()
        nintersections <- min(nintersections, length(combos))
        
        # Fetch the number of sets
        nsets <- length(getSelectedSetNames())
        setnames <- getSelectedSetNames()
        
        lines <-
          data.table::rbindlist(lapply(1:nintersections, function(combono) {
            data.frame(
              combo = combono,
              x = rep(combono, max(2, length(combos[[combono]]))),
              y = (nsets - combos[[combono]]) + 1,
              name = setnames[combos[[combono]]]
            )
          }))
        
        plot_ly(
          type = "scatter",
          mode = "markers",
          marker = list(color = "lightgrey", size = 8),
          customdata = "grid",
          source = "upset_hub"
        ) %>% add_trace(
          type = "scatter",
          x = rep(1:nintersections,
                  length(selected_sets)),
          y = unlist(lapply(1:length(selected_sets), function(x)
            rep(x - 0.5, nintersections))),
          hoverinfo = "none"
        ) %>% add_trace(
          type = "scatter",
          data = group_by(lines, combo),
          mode = "lines+markers",
          x = lines$x,
          y = lines$y - 0.5,
          line = list(color = "black", width = 3),
          marker = list(color = "black",
                        size = 10),
          hoverinfo = "text",
          text = ~ name
        ) %>% layout(
          xaxis = list(
            showticklabels = FALSE,
            showgrid = FALSE,
            zeroline = FALSE
          ),
          yaxis = list(
            showticklabels = FALSE,
            showgrid = TRUE,
            range = c(0, nsets),
            zeroline = FALSE,
            range = 1:nsets
          ),
          margin = list(t = 0, b = 40)
        )
        
      })
      
      # Make the bar chart illustrating set sizes
      
      upsetSetSizeBarChart <- reactive({
        setnames <- getSelectedSetNames()
        selected_sets <- getSets()
        
        plot_ly(
          x = unlist(lapply(selected_sets, length)),
          y = setnames,
          type = "bar",
          orientation = "h",
          marker = list(color = "black"),
          customdata = "setsize",
          source = "upset_hub"
        ) %>% layout(
          bargap = 0.4,
          yaxis = list(
            categoryarray = rev(setnames),
            categoryorder = "array"
          )
        )
      })
      
      # Make the bar chart illustrating intersect size
      
      upsetIntersectSizeBarChart <- reactive({
        ints <- calculateIntersections()
        intersects <- ints$intersections
        combos <- ints$combinations
        nintersections <- getNintersections()
        
        p <-
          plot_ly(showlegend = FALSE,
                  customdata = "intersectsize",
                  source = "upset_hub") %>% add_trace(
                    x = 1:nintersections,
                    y = unlist(intersects[1:nintersections]),
                    type = "bar",
                    marker = list(color = "black",
                                  hoverinfo = "none")
                  )
        
        bar_numbers <- TRUE
        
        if (bar_numbers) {
          p <-
            p %>% add_trace(
              type = "scatter",
              mode = "text",
              x = 1:nintersections,
              y = unlist(intersects[1:nintersections]) + (max(intersects) * 0.05),
              text = unlist(intersects[1:nintersections]),
              textfont = list(color = "black")
            )
        }
        
        p
      })
      
      observe({
        click <- event_data("plotly_click", source = "upset_hub")
        req(click)
        print(click)
      })
      
      return(NULL)
    }
    
    
    
    
  )
}