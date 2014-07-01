library(shiny)

source("helper.R")

zlim <- max(round(abs(range(sapply(saved.data,function(x) x$table$logFC)))))
samples <- names(saved.data)

## Define UI for application that draws a histogram
shinyUI(fluidPage(
    
    ## Application title
    titlePanel("Interactive Gene Heatmap"),
    
    ## Sidebar with a slider input for the number of bins
    sidebarLayout(
        sidebarPanel(

            ## Selection samples to show
            checkboxGroupInput("samples","Samples to show:",samples,selected=samples),

            ## Selecting a Fold Change cutoff
            numericInput("FC",
                         "Keep only genes with absolute fold change of:",
                         value=2,
                         min = 0,
                         max = Inf,
                         step=0.5),

            ## Selecting a p-value cutoff
            numericInput("pval",
                         "Keep only genes with a p value lower then",
                         0.01,
                         min = 0,
                         max = 1,
                         step = 0.01),

            ## Selecting the heat map contrast
            hr(),
            h5("Ajusting the contrast"),

            sliderInput("zlim.low",
                        "Maximum Blue Caped At:",
                        min = -zlim,
                        max = 0,
                        value = -zlim,,
                        step=0.5
                        ),
            
            sliderInput("zlim.high",
                        "Maximum Yellow Caped At:",
                        min = 0,
                        max = zlim,
                        value = zlim,
                        step=0.5
                        )
            
            
            ),
        
        ## Show a heatmap
        mainPanel(
            plotOutput("plot",height="600px"),
            textOutput("text")
            )
        )
    ))

