
ui <- fluidPage(theme = shinytheme("flatly"),
titlePanel(title=div(img(src="DiffBrainNet_logo.png", height = 50, width = 200), "Glucocorticoid receptor regulated gene expression in the mouse brain")),
navbarPage("Explore!",
           
           tabPanel(icon("home"),
                    fluidRow(column(tags$img(src="mousejavi_reversebrain.png",
                                             width="100%", height= "100%"),
                                             #width="450px",height="320px"), 
                                    width=floor(0.4*12)),
                             column(
                               p(strong("Abstract."), "Network analysis can identify the molecular connectivity 
                               that is underpinning function at baseline, after a stimulus or a disease state. We inferred 
                               regression-based prior-knowledge guided gene networks in 8 brain regions of the mouse brain: 
                               the prefrontal cortex, the amygdala, the paraventricular nucleus of the hypothalamus, the dorsal 
                               and ventral Cornu ammonis 1, the dorsal and ventral dentate gyrus and the cerebellar cortex. 
                               We constructed networks at baseline and treatment levels using KiMONo (Ogris et al.) and at 
                               the differential level using DiffGRN (Kim et al.). DiffGRN uses the regression coefficients 
                               of the baseline and treatment networks to calculate differential connections, thus providing 
                               us with the actual functional couplings of the genes after the stimulus that differ from the 
                               baseline ones.",br(),
                               "As a stimulus we used dexamethasone. Dexamethasone is a synthetic glucocorticoid that is used 
                       
                               to activate the glucocorticoid receptors. Glucocorticoid receptors, when coupled with glucocorticoids 
                               like dexamethasone, act as transcription factors modulating the transcriptional landscape. 
                               Glucocorticoid receptors chromatin binding is one of the main results of stress-axis activation, 
                               which in turn is an important component of the biology of stress-related psychiatric disorders. 
                               Thus, dexamethasone mimics some aspects of the altered transcriptomic landscape that is associated
                               with stress-related psychiatric disorders. We provide differential networks and differential 
                               expression analysis (DESeq2, Love et al.) that can be used to analyse the effects of dexamethasone 
                               both at the molecular connectivity and at the gene level in each brain region.", br(),
                               "DiffBrainNet is an analysis framework and a resource for studying the transcriptional landscape 
                               of 8 mouse brain regions at baseline, dexamethasone-treatment and differential levels. It can be 
                               used to pinpoint molecular pathways important for the basic function and response to glucocorticoids 
                               in a brain-region specific manner. DiffBrainNet can also support the identification and analysis 
                               of biological processes regulated by brain and psychiatric diseases risk genes at the baseline and 
                               differential levels.",
                                 style="text-align:justify;color:black;background-color:#2980B9;padding:15px;border-radius:10px"),
                               br(),
                               p(strong("Data and code availability. Reference to paper."), "Sed ut perspiciatis unde omnis iste 
                                 natus error sit voluptatem accusantium doloremque laudantium, totam rem aperiam, eaque ipsa quae 
                                 ab illo inventore veritatis et quasi architecto beatae vitae dicta sunt explicabo. Nemo enim 
                                 ipsam voluptatem quia voluptas sit aspernatur aut odit aut fugit, sed quia consequuntur magni 
                                 dolores eos qui ratione voluptatem sequi nesciunt. Neque porro quisquam est, qui dolorem ipsum 
                                 quia dolor sit amet, consectetur, adipisci velit, sed quia non numquam eius modi tempora incidunt 
                                 ut labore et dolore magnam aliquam quaerat voluptatem. Ut enim ad minima veniam, quis nostrum 
                                 exercitationem ullam corporis suscipit laboriosam, nisi ut aliquid ex ea commodi consequatur? 
                                 Quis autem vel eum iure reprehenderit qui in ea voluptate velit esse quam nihil molestiae 
                                 consequatur, vel illum qui dolorem eum fugiat quo voluptas nulla pariatur?",
                                 style="text-align:justify;color:black;background-color:#AED6F1;padding:15px;border-radius:10px"),
                               width = floor(0.45*12)
                             ),
                             column(
                               br(),
                               tags$img(src="mpilogo.png", width="250px", height="60px"),
                               br(),
                               br(),
                               p("For more information please visit the website of", em("the MPI of Psychiatry"),
                                 br(),
                                 a(href="https://www.psych.mpg.de/", "Here", target="_blank"),
                                 style="text-align:center;color:black"),
                               width=ceiling(0.15*12)
                             ))),
           
           navbarMenu("Diff. Gene Expression",
           tabPanel("Single Brain Region",
                    sidebarLayout(
                      sidebarPanel(
                        selectInput("region", "Choose a brain region:", 
                                    choices = c("Amygdala" = "AMY", 
                                                "Cerebellum" = "CER", 
                                                "Dorsal DG of Hippocampus" = "dDG", 
                                                "Dorsal CA1 of Hippocampus" = "dCA1", 
                                                "Prefrontal Cortex" = "PFC", 
                                                "PVN of Hypothalamus" = "PVN", 
                                                "Ventral DG of Hippocampus" = "vDG", 
                                                "Ventral CA1 of Hippocampus" = "vCA1")),
                        downloadButton("download1","Download (filtered) table as csv"),
                        width = 3
                      ),
                      mainPanel(
                        fluidPage(
                          br(),
                          tags$style(".fa-chart-bar {color:#2980B9}"),
                          h3(p(em("Differential Gene Expression "),icon("chart-bar",lib = "font-awesome"),style="color:black;text-align:center")),
                          hr(),
                          span("Explore the transcriptomic response to Dexamethasone treatment of a brain region of your 
                          interest on a "), strong("differential expression level"), span(". You can select one of 8
                          different brain regions on the left panel. The volcano plot shows the log2-transformed fold change 
                          and the -log10-transformed FDR corrected p-value of the genes detected in our dataset. Please "), 
                          strong("click on a point/gene"), span(" in the volcano plot to see its 
                          normalized expression levels before and after Dexamethasone treatment. You can filter and
                          download the table of the differentially expressed genes."),
                          br(),br(),
                          #  style="text-align:left;color:black"),
                          splitLayout(cellWidths = c("55%", "45%"),
                          plotlyOutput("volcano_plot"),
                          plotlyOutput("exp_plot")),
                          DT::dataTableOutput("de_table")
                        )
                      )
                    )
           ),
           tabPanel("Comparison Brain Regions",
                    # UI defined in upsetPlot module
                    upsetPlotUI("upsetDE")
           )
           ),
           
           navbarMenu(
             "Network Analysis",
             tabPanel(
               "Introduction diff. networks",
               fluidPage(
                 br(),
                 tags$style(".fa-project-diagram {color:#2980B9}"),
                 h3(p(
                   em("Introduction: Differential network analysis"),
                   icon("project-diagram", lib = "font-awesome"),
                   style = "color:black;text-align:center"
                 )),
                 hr(),
                 
                 tags$img(src="DiffNetworks.png",
                          width="60%", height= "60%", 
                          style="display: block; margin-left: auto; margin-right: auto;")
               )
                 
                 
               
             ),
             
             tabPanel(
               "Network overview",
               sidebarPanel(
                 selectInput(
                   "overview_region",
                   "Choose a brain region:",
                   choices = c(
                     "Amygdala" = "AMY",
                     "Cerebellum" = "CER",
                     "Dorsal DG of Hippocampus" = "dDG",
                     "Dorsal CA1 of Hippocampus" = "dCA1",
                     "Prefrontal Cortex" = "PFC",
                     "PVN of Hypothalamus" = "PVN",
                     "Ventral DG of Hippocampus" = "vDG",
                     "Ventral CA1 of Hippocampus" = "vCA1"
                   )
                 ),
                 radioButtons(
                   "overview_metric",
                   label = "Select network metric:",
                   choices = list("Nodedegree" = "nodedegrees",
                                  "Nodebetweenness" = "nodebetweenness"),
                   selected = "nodebetweenness"
                 ),
                 radioButtons(
                   "overview_network",
                   label = "Select network:",
                   choices = list("Baseline" = "baseline",
                                  "Differential" = "differential"),
                   selected = "differential"
                 ),
                 width = 3
               ),
               mainPanel(
                 fluidPage(
                   br(),
                   tags$style(".fa-project-diagram {color:#2980B9}"),
                   h3(p(
                     em("Network analysis of gene expression: Overview "),
                     icon("project-diagram", lib = "font-awesome"),
                     style = "color:black;text-align:center"
                   )),
                   hr(),
                   splitLayout(cellWidths = c("50%", "50%"),
                               h4(
                                 p("Baseline network", style = "color:black;text-align:center")
                               ),
                               h4(
                                 p("Differential network", style = "color:black;text-align:center")
                               )),
                   # plotlyOutput("histogram_network")
                   splitLayout(
                     cellWidths = c("50%", "50%"),
                     plotlyOutput("histogram_network"),
                     plotlyOutput("barplot_network")
                   )
                   # DT::dataTableOutput("network_table")
                 )
               )
             ),
             
             tabPanel(
               "Comparison hub genes",
               comparisonHubGenesUI("compHubGenes")
             ),
             
             tabPanel(
               "Network visualization",
               # UI from networkSingle module
               networkSingleUI("singleVisualization")
             ),
             
             tabPanel(
               "Multi-region network visualization",
               # UI from networkMulti module
               networkMultiUI("multiVisualization")
             )
           ), 
           
           
           navbarMenu("More",
                      tabPanel("Data download",
                               # DT::dataTableOutput("table")
                      ),
                      tabPanel("About",
                               # fluidRow(
                               #   column(6,
                               #          includeMarkdown("about.md")
                               #   ),
                               #   column(3,
                               #          img(class="img-polaroid",
                               #              src=paste0("http://upload.wikimedia.org/",
                               #                         "wikipedia/commons/9/92/",
                               #                         "1919_Ford_Model_T_Highboy_Coupe.jpg")),
                               #          tags$small(
                               #            "Source: Photographed at the Bay State Antique ",
                               #            "Automobile Club's July 10, 2005 show at the ",
                               #            "Endicott Estate in Dedham, MA by ",
                               #            a(href="http://commons.wikimedia.org/wiki/User:Sfoskett",
                               #              "User:Sfoskett")
                               #          )
                               #   )
                               # )
                      )
           )
)
)
