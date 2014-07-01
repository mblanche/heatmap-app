library(shiny)

source("helper.R")

# Define server logic required to draw a heatmap
shinyServer(function(input, output) {

    ## A reactive expression computing what to cluster and display
    data <- reactive({
        samples <- input$samples

        ## Keep any gene that passes the FC and pval filter in one experiment
        filter <- rowSums(sapply(saved.data[samples], function(d){
                ## here, I am filtering on the pvalue and fold change the user chooses
            ## Creating two filters depending on reactive values from the UI
            fc.filter <- abs(d$table$logFC) >= log2(input$FC)
                pv.filter <- d$table$PV <  input$pval
            return(fc.filter & pv.filter)
        })) > 0
        
        ## get the fold change as a matrix
        FC.mat <- sapply(saved.data[samples],function(d) d$table$logFC)

        ## Filter the data
        return(FC.mat[filter,,drop=FALSE])

    })
    
    ## Generate the cluster, a second reactive expression
    ## here, I am making it reactive so that I could provide the user
    ## with different clustering algorithm choices (Pearson, Spearman, Euclidean, etc...)
    cluster <- reactive({
        ## filter the NA value first
        forHeatmap <- data()[!apply(is.na( data() ),1,any),]
        ## Returning the clustered data
        cluster <- hclust(dist( forHeatmap ))
    })

    ## Rendering our heatmap
    ## Have some logic to deal with user not clicking any samples
    observe({
        if(length(input$samples) == 0){
            ## Wipeout the ploting area
            output$plot <- renderPlot({ return(NULL) })
        } else {
            output$plot <- renderPlot({
                plotHeatMap(data(),cluster(),c(input$zlim.low,input$zlim.high) )
            })
        }
    })

})
