library(shiny)
library(shinyIncubator)
library(hwriter)

source("helper.R")

# Define server logic required to draw a heatmap
shinyServer(function(input, output, session) {
    
    ## A reactive expression computing what to cluster and display
    data <- reactive({
        if (is.null(input$samples) || input$pval==0) {
            return(NULL)
        } else {
            samples <- input$samples ## Just to make the code smaller to read
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
        }
    })
    
    ## Generate the cluster, a second reactive expression
    ## here, I am making it reactive so that I could provide the user
    ## with different clustering algorithm choices (Pearson, Spearman, Euclidean, etc...)
    ## Using a observer on the data to make sure I have something to cluster.
    ## Could perhaps merge with lower evaluation
    cluster <- reactive({
        if(length(data()) == 0){
            return(NULL)
        } else {
            ## filter the NA value first
            forHeatmap <- data()[!apply(is.na( data() ),1,any),]
            ## Returning the clustered data
            hc <- try(as.dendrogram(hclust(dist( forHeatmap ))),silent=TRUE)
            if(class(cluster) == 'try-error'){
                return(NULL)
            } else {
                return(hc)
            }
        }
    })
    
    ## Computing a cluster colored based on selected height
    color.cluster <- reactive ({
        if (is.null(cluster())){
            return(NULL)
        } else {
            if (length(input$height) != 0){
                hc <- colBranches(cluster(),input$height,hc.cols)
            }
        }
    })

    ## Rendering a slider to select the height used to break the cluster
    output$heightSelector <- renderUI({
        if (is.null(cluster())){
            return(NULL)
        } else {
            h <- attributes(cluster())$height
            list(hr(),
                 h5("Creating clusters of genes"),
                 sliderInput("height",
                             "Break in clusters at height of:",
                             min = round(h*0.2,2),
                             max = h,
                             value = h,
                             format = '#.00'
                             )
                 )
        }
    })
    
    ## Create a set of reactive values to store tables of genes displayed in the heatmap
    values <- reactiveValues()
    ## Create a reactive context to assign the reactive values
    observe({
        if(is.null(data())){
            values <- NULL
        } else {
            for (s in input$samples) {
                ## The reactive values data() should trigger re-evalution of this bit on modification
                geneSymbol <- gene2name$external_gene_id[match(rownames(data()),gene2name$ensembl_gene_id)]
                linkOut <- 'http://flybase.org/reports/'
                links <- hwrite(geneSymbol, 
                                link = paste0(linkOut,rownames(data())),
                                target = s,
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
        if(length(input$samples) == 0 || is.null(values) || is.null(color.cluster())){
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
                          call("tabPanel","Clusters",
                               call("uiOutput","clusters")),
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

            ##withProgress(session, {
              ##  setProgress(message = "Recomputing the heatmap and cluster data",
                ##            detail = "This may take a few moments...")
            output$plot <- renderPlot({
                plotHeatMap(data(),
                            color.cluster(),
                            c(input$zlim.low,input$zlim.high),
                            input$height
                            )
            })
        }
    })
    
    ## Create a reactive context to populate the sample tab panels with content
    observe({
        lapply(names(values), function(s){
            ## Do some cleanup before saving
            d <- values[[s]]
            d$Gene <- sub("<a.+?>(.+)</a>","\\1",d$Gene)
            d <- cbind(FBid=rownames(d),d)
            ## Add a DataTable of gene selected in the heatmap
            output[[paste0('table_',s)]] <- renderDataTable(values[[s]],options=list(iDisplayLength=10))
            ## Add a download button to allow download of a csv file
            output[[paste0('save_',s)]] <- downloadHandler(
                filename = function() { paste0(s,".csv") },
                content = function(file) { write.csv(d,file=file,row.names=FALSE) }
                )
            return(s)
        })
    })
    
    ## Create a reactive context to populate the cluster tab panel
    observe({
        if (!is.null(color.cluster()) ){
            ## Cutting the clustering into sub-group
            cuts <- cut(color.cluster(),h=input$height) 
            sub.dendro <- rev(cuts$lower)
            groups <- lapply(sub.dendro,unlist)

            ## Creating a new UI for each group
            ui <- unlist(lapply(seq(groups),function(i){
                c(call("h3",paste("Cluster",i)),
                  call("plotOutput",paste0("subDendro_",i),height="100px",width="100px"),
                  call("dataTableOutput",paste0("cluster_",i)),
                  call("downloadButton",paste0("saveCluster_",i),'Save as csv'),
                  call("hr")
                  )
                
            }))
            output$clusters <- renderUI({ lapply(ui,eval) })

            ## Rendering the sub-dendro plot and a table of genes
            lapply(seq(groups),function(i){
                output[[paste0('subDendro_',i)]] <- renderPlot({
                    par(mar=rep(0,4))
                    plot(sub.dendro[[i]],
                         axes=FALSE,
                         yaxt='s',
                         yaxs='i',
                         xaxt='n',
                         xaxs='i',
                         horiz=TRUE,
                         leaflab='none')
                })

                ## Rendering the tables of gene for each cluster, linking out to flybase
                geneIds <- rownames(data())[groups[[i]]]
                geneSymbol <- gene2name$external_gene_id[match(geneIds,gene2name$ensembl_gene_id)]
                linkOut <- 'http://flybase.org/reports/'
                links <- hwrite(geneSymbol, 
                                link = paste0(linkOut,geneIds),
                                target = 'subCluster',
                                table = FALSE)

                d <- as.data.frame(sapply(input$samples,function(s) saved.data[[s]]$table[geneIds,'logFC']))
                names(d) <- paste(names(d),"Log2(FC)")
                
                output[[paste0('cluster_',i)]] <- renderDataTable(cbind(Genes=links,d),
                                                                  options=list(iDisplayLength=10))
                
                ## Creating functions to save the tables link to the buttons
                output[[paste0('saveCluster_',i)]] <- downloadHandler(
                    filename = function() { paste0("cluster_",i,"_h",round(input$height,1),".csv") },
                    content = function(file) {
                        write.csv(data.frame(gene.id=geneIds,symbol=geneSymbol),
                                  file=file,row.names=FALSE)
                    }
                    )
            })
        }
        
    })
    
    ## Function to download the heatmap as png
    output$img <- downloadHandler(
        filename = function() { "heatmap.png" },
        content = function(file) {
            png(file)
            plotHeatMap(data(),cluster(),c(input$zlim.low,input$zlim.high),input$height,noMarker=TRUE)
            dev.off()
        }
        )
})




