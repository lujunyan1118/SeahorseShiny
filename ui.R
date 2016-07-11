library(shiny)

load("shinyData.RData")

pca <- prcomp(na.omit(t(seaMain)))
pca <- pca$x

navbarPage("Tasks",
           tabPanel("Correlations",
                    titlePanel("Correlations between metabolic features and drug responses"),
                    sidebarLayout(
                      sidebarPanel(
                        radioButtons("screenSet","Choose drug screen dataset",c("Main screen" = "main", "EMBL screen" = "embl")),
                        checkboxInput("filterCor", "Pre-filter by correlation tests", value = FALSE),
                        conditionalPanel(condition = "input.filterCor == true",
                                         numericInput("coefCut","Coefficient cutoff",value=0, min = 0, max=1, step=0.1),
                                         numericInput("pCut","Raw P-value cutoff",value=1, min = 0, max=1,step=0.01)),
                        uiOutput("measureBox"),
                        uiOutput("drugBox"),
                        uiOutput("concBox"),
                        checkboxInput("colorGene", "Color samples by genetic background", value = FALSE),
                        conditionalPanel(condition = "input.colorGene == true",
                                         selectInput("geneBox","Select genetic feature",colnames(geneMain),size=5, selectize = FALSE),
                                         checkboxInput("showNA", "Show samples without genetic data", value = TRUE),
                                         checkboxInput("sperateRegression", "Perform regression separately", value = FALSE)
                        )
                      ),
                      mainPanel(plotOutput("scatterPlot")
                      )
                    )
           ),
           tabPanel("t-Test",
                    titlePanel("Associations between Seahorse measurements and genetic background"),
                    sidebarLayout(
                      sidebarPanel(
                          selectInput("seleSea1","Select a metabolic feature",rownames(seaMain),size=5, selectize = FALSE),
                          radioButtons("tMethod", "Select a t-test method", choices = c("Equal variance","Unequal variance")),
                          checkboxInput("filterP", "Filter results by raw or adjusted P values", value = FALSE),
                          conditionalPanel(condition = "input.filterP == true",
                                           numericInput("rawPcut","Raw p value threshold",value=0.05, min = 0, max=1, step=0.05),
                                           numericInput("adjPcut","Adjusted P-value threshold",value=0.1, min = 0, max=1,step=0.05)),
                          checkboxInput("treatIGHV","Consider IGHV status", value = FALSE ),
                          conditionalPanel(condition = "input.treatIGHV == true",
                                           radioButtons("methodIGHV","How to deal with IGHV status",
                                                        choices = c("Remove","Unmutated only","Mutated only"))),
                      width = 4),
                      mainPanel(column(DT::dataTableOutput("tTable"),width = 5),
                                column(plotOutput("tPlot",height = 500, width = 600),
                                       checkboxInput("showNA2","Show NA points", value=FALSE),width=7)
                    )
           )),
           tabPanel("t-SNE",
                    titlePanel("Visualization of the energy metabolism dataset"),
                    sidebarLayout(
                      sidebarPanel(
                        radioButtons("dataSet","Dataset",c("Seahorse","Main screen","EMBL screen","RNAseq")),
                        radioButtons("clusterMethod","Methods",c("t-SNE" = "ifTSNE", "PCA" = "ifPCA","k-means"= "ifKmeans")),
                        conditionalPanel(condition = "input.clusterMethod == 'ifTSNE'",
                          numericInput("perplexity","Perplexity",value=10),
                          numericInput("theta","Theta (0 for exact t-SNE)",value=0.5, min=0,max=1),
                          numericInput("maxIter","Maximum number of iterations", value = 1000)),
                        conditionalPanel(condition = "input.clusterMethod == 'ifKmeans'",
                          sliderInput("clusterNum", label = "Cluster number", min = 2, max = 20, value = 2)),
                        conditionalPanel(condition = "input.clusterMethod == 'ifPCA'",
                          selectInput("xaxis", label = NULL, choices = colnames(pca), selected = "PC1"),
                          selectInput("yaxis", label = NULL, choices = colnames(pca), selected = "PC2")),
                        uiOutput("geneBox2")
                      ),  
                      mainPanel(plotOutput("tsnePlot"))
                    )
           ),
           tabPanel("Scatter",
                    titlePanel("Plotting all data points colored by genetic backgrounds"),
                    sidebarLayout(
                      sidebarPanel(
                        radioButtons("dataSet1","Dataset",c("Seahorse","Main screen")),
                        conditionalPanel(condition = "input.dataSet1 == 'Main screen'",
                            radioButtons("patSet",label = "Patient set", choices = c("All patients","Only CLL patients")),
                            radioButtons("drugSet",label = "Drug set", choices = c("All drugs","High variant drugs","Drug set 1 (small)","Drug set 2 (large)"), selected = "Drug set 1 (small)"),
                            conditionalPanel(condition = "input.drugSet == 'High variant drugs'",
                                numericInput("varCut","Top percentage", min = 10, max=100, value = 10))
                        ),
                        uiOutput("geneBox3")
                      ),
                      mainPanel(plotOutput("scatterAll"))
                    )
           ),
           tabPanel("RNAseq",
                    titlePanel("Exploring the associations between Seahorse dataset and RNAseq dataset"),
                    sidebarLayout(
                      sidebarPanel(
                        selectInput("seleSea","Select metabolic feature",rownames(seaMain),size=5, selectize = FALSE),
                        radioButtons("corMethod", "Select correlation method", choices = c("pearson","spearman")),
                        checkboxInput("filterRNA", "Filter results by coefficients and p-values", value = TRUE),
                        conditionalPanel(condition = "input.filterRNA == true",
                                         numericInput("coefCut1","Coefficient cutoff",value=0.5, min = 0, max=1, step=0.1),
                                         numericInput("pCut1","Raw P-value cutoff",value=0.05, min = 0, max=1,step=0.01)),
                        checkboxInput("colorGene2", "Color by genetic background", value = FALSE),
                        conditionalPanel(condition = "input.colorGene2 == true",
                                         selectInput("geneBox4","Select genetic feature",colnames(geneMain),size=3, selectize = FALSE)),
                        checkboxInput("ifGSEA", "Perform GSEA according to correlations (may be time consuming)", value = FALSE),
                        conditionalPanel(condition = "input.ifGSEA == true",
                                         selectInput("setBox","Select geneset database",list.files(pattern = "\\.gmt$"), selectize = FALSE),
                                         numericInput("permNum","Permutation number",value=100,min = 10, max = 10000),
                                         actionButton("goGSEA", "Perform GSEA"),
                                         numericInput("sigLevel", label = "Significance", min = 0, max = 1, value = 0.01))
                      ,width=4),
                      mainPanel(downloadButton('downloadTable', 'Download current table'),
                        DT::dataTableOutput("tab1"),
                        column(plotOutput("plot1",height = 500, width = 500),width = 6),
                        column(DT::dataTableOutput("tab2"),width = 6)
                      )
                    )
           )
)
