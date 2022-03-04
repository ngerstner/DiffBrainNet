# module UI function
upsetPlotUI <- function(id, label = "Upset Plot"){
  
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
             sliderTextInput(
               inputId = ns("padj_upset"),
               label = "Choose a threshold for FDR p-value:",
               choices = c(0.01, 0.05, 0.1),
               selected = 0.1,
               grid = TRUE
             ),
             sliderTextInput(
               inputId = ns("FC_upset"),
               label = "Choose a threshold for logFC:",
               choices = c(0.0, 0.5, 1.0, 2.0),
               selected = 0.0,
               grid = TRUE
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
             h3(p(em("Comparsion of DE genes across brain regions "),icon("brain",lib = "font-awesome"),style="color:black;text-align:center")),
             hr(),
             tags$style(type="text/css",
                        ".shiny-output-error { visibility: hidden; }",
                        ".shiny-output-error:before { visibility: hidden; }"
             ),
             strong("Compare"),span(" the differentially expressed (DE) genes after Dexamethasone treatment"),
             strong(" between different brain regions"), span(". Please choose the brain regions you are interested 
             in on the left panel. You can also select a different FDR p-value threshold and a different fold 
             change cutoff. The table on the bottom displays all genes in the dataset with the brain regions 
             they are differentially expressed in."),
             br(),br(),
             # plotlyOutput("upset_de")
             #plotOutput("upset_de"),
             plotlyOutput(ns("plotly_upset"), height = "600px"),
             DT::dataTableOutput(ns("upset_table"))
             #h3("Intersection plots")
           ))
  
}

# module server function
upsetPlotServer <- function(id) {
  moduleServer(
    id,
    
    function(input, output, session) {
      
      # Read data required for upset plot
      read_upset_data <- function(basepath, regions, padj_cutoff = 0.1, fc_cutoff = 0){
        
        # Read DE tables from all required regions
        list_genes_sig <- list()
        
        for (reg in regions){
          res <- fread(file=paste0(basepath, "tables/02_", reg, 
                                   "_deseq2_Dex_1_vs_0_lfcShrink.txt"),
                       sep="\t") 
          na_indices <- which(is.na(res$padj))
          res$padj[na_indices] <- 1
          res_sig <- res[res$padj <= padj_cutoff,]
          res_sig <- res_sig[abs(res_sig$log2FoldChange) >= fc_cutoff,]
          list_genes_sig[[reg]] <- res_sig$Ensembl_ID
        }
        
        return(list_genes_sig)
      }
      
      # Read data required for data table below upset plot
      upset_datatable <- function(basepath, regions, padj_cutoff = 0.1, fc_cutoff = 0){
        
        # READ DATA 
        
        # Read DE tables from all required regions
        list_reg_sig <- list()
        
        for (reg in regions){
          res <- fread(file=paste0(basepath, "tables/02_", reg, 
                                   "_deseq2_Dex_1_vs_0_lfcShrink.txt"),
                       sep="\t") 
          na_indices <- which(is.na(res$padj))
          res$padj[na_indices] <- 1
          res_sig <- res[res$padj <= padj_cutoff,]
          res_sig <- res_sig[abs(res_sig$log2FoldChange) >= fc_cutoff,]
          list_reg_sig[[reg]] <- res_sig
        }
        
        df <- bind_rows(list_reg_sig, .id="region")
        df <- df[,c("Ensembl_ID", "region")] %>% 
          group_by(Ensembl_ID) %>%
          dplyr::summarise(region = list(region)) %>%
          dplyr::select(Ensembl_ID, region) %>%
          dplyr::distinct(Ensembl_ID, .keep_all = TRUE)
        df$Gene_Symbol <- mapIds(org.Mm.eg.db, keys = df$Ensembl_ID,
                                 column = "SYMBOL", keytype = "ENSEMBL")
        df <- dplyr::select(df, Gene_Symbol, Ensembl_ID, region)
  
        # df <- df[,c("Ensembl_ID", "Gene_Symbol", "region", "padj", "log2FoldChange")] %>% 
        #   group_by(Ensembl_ID, Gene_Symbol) %>%
        #   dplyr::summarise(region = list(region), 
        #                    p_adjusted = list(padj),
        #                    log2FC = list(log2FoldChange)) %>%
        #   dplyr::select(Gene_Symbol, Ensembl_ID, region, p_adjusted, log2FC)
        
        return(df)
      }
      
      # reactive table 
      upset_data <- reactive({
        
        req(input$myPicker)
        data <- upset_datatable(basepath, input$myPicker, input$padj_upset, input$FC_upset) %>%
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
        filename = "DEGenes.csv",
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
        
        selected_sets <- read_upset_data(basepath, selected_set_names, input$padj_upset, input$FC_upset)
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
          source = "upset_de"
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
          source = "upset_de"
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
                  source = "upset_de") %>% add_trace(
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
        click <- event_data("plotly_click", source = "upset_de")
        req(click)
        print(click)
      })
      
      return(NULL)
    }
    
    
    
    
  )
}