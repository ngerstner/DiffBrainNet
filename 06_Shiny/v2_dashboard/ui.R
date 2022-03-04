library(semantic.dashboard)
library(shiny)
library(markdown)
library(plotly)

ui <- dashboardPage(dashboardHeader("Explore!"),
              ## Sidebar content
              dashboardSidebar(sidebarMenu(
                menuItem(tabName = "home", text = "Home", icon = icon("home")),
                menuItem(tabName = "de", text = "Diff. Gene Expression", icon = "heart"))
              ),
              dashboardBody(
                tabItems(
                  # First tab content
                  tabItem(tabName = "de",
                          fluidRow(
                            selectInput("region", "Choose a brain region:", 
                                            choices = c("AMY", "CER", "dDG", 
                                                        "dCA1", "PFC", 
                                                        "vDG", "vCA1")),
                            plotlyOutput("volcano_output"),
                            plotlyOutput("exp_plot")
                          ))
                )
              ))


library(data.table)
library(stringr)
library(ggplot2)


# Define server logic required to generate and plot 
server <- function(input, output) {
  
  output$volcano_plot <- renderPlotly({
    
    table <- fread(paste0("/Users/nathalie_gerstner/Documents/ownCloud/DexStim_RNAseq_Mouse/tables/02_",
                          input$region,"_deseq2_Dex_1_vs_0_lfcShrink.txt"))
    # table <- fread(paste0("/Users/nathalie_gerstner/Documents/ownCloud/DexStim_RNAseq_Mouse/tables/02_",
    #                       "AMY","_deseq2_Dex_1_vs_0_lfcShrink.txt"))
    table$sig <- as.factor(ifelse(table$padj <= 0.1, "p_adj <= 0.1", "p_adj > 0.1"))
    
    # plot volcano plot
    volcano_plot <- plot_ly(data = table, x = ~log2FoldChange, y = ~-log10(padj),
                            color = ~sig, colors = "Set1",
                            text = ~paste0("Gene: ", V1)) %>%
      layout(title = paste0("Volcano Plot ", input$region))
  })
  
  output$exp_plot <- renderPlotly({
    d <- event_data("plotly_click")
    print(d)
    if (is.null(d)) return(NULL) 
    
    table <- fread(paste0("/Users/nathalie_gerstner/Documents/ownCloud/DexStim_RNAseq_Mouse/tables/02_",
                          input$region,"_deseq2_Dex_1_vs_0_lfcShrink.txt"))
    ensembl_id <- table$V1[d$pointNumber + 1]
    
    exp_table <- fread(paste0("/Users/nathalie_gerstner/Documents/ownCloud/DexStim_RNAseq_Mouse/tables/02_",
                              input$region,"_deseq2_expression_vsd.txt"))
    # exp_table <- fread(paste0("/Users/nathalie_gerstner/Documents/ownCloud/DexStim_RNAseq_Mouse/tables/02_",
    #                           "AMY","_deseq2_expression_vsd.txt"))
    
    exp_gene <- transpose(exp_table[V1 == ensembl_id,2:ncol(exp_table)], keep.names = "col")
    exp_gene$col <- str_replace(exp_gene$col, ".*\\_","")
    exp_gene$mouse_id <- str_extract(exp_gene$col, pattern = "[0-9]+")
    exp_gene$col <- as.factor(str_replace(exp_gene$col, "[0-9]+",""))
    
    exp_plot <- ggplot(exp_gene, aes(x = col, y = V1)) +
      geom_boxplot() +
      geom_jitter(color="black", size=0.4, alpha=0.9) +
      theme_light() +
      theme(
        legend.position="none",
        plot.title = element_text(size=11)
      ) +
      ggtitle("A boxplot with jitter") +
      xlab("")
    ggplotly(exp_plot)
    
    # gSel <- tg %>% filter(dose %in% d$y) %>% group_by(supp) %>%  mutate(newLen=floor(len)) %>% 
    #   ggplot(aes(x=supp, fill=as.factor(newLen))) + geom_bar()
    # ggplotly(gSel)
    
    
  })
}

shinyApp(ui, server)
           