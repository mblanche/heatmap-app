library(shiny)
library(hwriter)

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
        FC.mat <- do.call(cbind,lapply(saved.data[samples],function(d){
            res <- data.frame(d$table$logFC)
            rownames(res) <- rownames(d$table)
            return(res)
        }))
        colnames(FC.mat) <- names(saved.data[samples])
        
        ## Filter the data
        return(FC.mat[filter,,drop=FALSE])
    })
    
    ## Generate the cluster, a second reactive expression
    ## here, I am making it reactive so that I could provide the user
    ## with different clustering algorithm choices (Pearson, Spearman, Euclidean, etc...)
    ## Using a observer on the data to make sure I have something to cluster.
    ## Could perhaps merge with lower evaluation
    cluster <- reactive({
        ## filter the NA value first
        forHeatmap <- data()[!apply(is.na( data() ),1,any),]
        ## Returning the clustered data
        cluster <- try(hclust(dist( forHeatmap )),silent=TRUE)
        if(class(cluster) == 'try-error'){
            return(NULL)
        } else {
            return(cluster)
        }
    })
            
    ## Create a set of reactive values to store tables of genes displayed in the heatmap
    values <- reactiveValues()
    ## Create a reactive context to assign the reactive values
    observe({
        if(nrow(data()) == 0){
            values <- NULL
        } else {
            for (s in input$samples) {
                ## The reactive values data() should trigger re-evalution of this bit on modification
                geneSymbol <- gene2name$external_gene_id[match(rownames(data()),gene2name$ensembl_gene_id)]
                
                linkOut <- 'http://flybase.org/reports/'
                
                links <- hwrite(geneSymbol, 
                                link = paste0(linkOut,rownames(data())),
                                table = FALSE)
                
                d <- data.frame('id'= links,
                                saved.data[[s]]$table[rownames(data()),c('logFC','PValue')])
                
                d$adj.p <- p.adjust(d$PV,'BH')
                colnames(d) <- c('Gene','log2(FC)','p-value','adj. p-value')
                values[[s]] <- d
            }
        }
    })
    
    ## Rendering our heatmap
    ## The ui main panel is render at the same time
    ## Have some logic to deal with user not clicking any samples
    observe({
        if(length(input$samples) == 0 || is.null(values) || is.null(cluster())){
            ## Wipeout the ploting area
            output$main <- renderUI({ return(NULL) })
            ## Wipeout the tabsets
            output$plotUI <- renderUI({ tabsetPanel("tabPanel") })
            ## Remove the message of number of genes selected
            output$message <- renderUI({ return(NULL) })
        } else {
            ## Print a message with the number of selected genes
            output$message <- renderUI({
                list(hr(),
                     h5(paste(nrow(data()),"genes are selected"))
                     )
            })
            ## Create a UI with tab panels, one for the plot and one for each slected samples
            output$plotUI <- renderUI({
                ## Dynamically render tabset based on the user selected samples
                do.call(tabsetPanel,
                        c(call("tabPanel","Plot",
                               call("plotOutput","plot",height='600px'),
                               call("downloadButton",'img','Save as png')
                                ),
                          lapply(input$samples,function(s){
                              call("tabPanel",s,
                                   call('textOutput',paste0("text_",s)),
                                   call('dataTableOutput',paste0("table_",s)),
                                   call("downloadButton",paste0("save_",s),'Save as csv')
                                   )
                          })
                           )
                        )
            })
            ## Render a heatmap and dendrogram in the ploting tab panel
             output$plot <- renderPlot({
                 plotHeatMap(data(),cluster(),c(input$zlim.low,input$zlim.high) )
             })
        }
    })
    
    
    ## Create a reactive context to populate the tab panels with content
    observe({
        lapply(names(values), function(s){
            ## Add a DataTable of gene selected in the heatmap
            output[[paste0('table_',s)]] <- renderDataTable(values[[s]],options=list(iDisplayLength=10))
            ## Add a download button to allow download of a csv file
            output[[paste0('save_',s)]] <- downloadHandler(
                filename = function() { paste0(s,".csv") },
                content = function(file) { write.csv( values[[s]],file=file) }
                )
            return(s)
        })
    })
    
    ## Function to download the heatmap as png
    output$img <- downloadHandler(
        filename = function() { "heatmap.png" },
        content = function(file) {
            pdf(file)
            plotHeatMap(data(),cluster(),c(input$zlim.low,input$zlim.high) )
            dev.off()
        }
        )
    
})



