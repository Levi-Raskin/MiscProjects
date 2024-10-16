library(shiny)
library(phytools)
library(phangorn)
library(tictoc)
library(TreeTools)
library(ggtree)

ui <- fluidPage(

  sidebarLayout(
    sidebarPanel(
                 title = "phylogenetic inference with ambiguous data",
                 sliderInput(
                   inputId = "nTaxa",
                   label = "Number of taxa:",
                   min = 3,
                   max = 10,
                   value = 7
                 ),
                 sliderInput(
                   inputId = "nChar",
                   label = "Number of characters:",
                   min = 1,
                   max = 500,
                   value = 100,
                   ticks = F
                 ),
                 numericInput(
                   inputId = "rate",
                   label = "Rate of sequence evolution:",
                   min = 0,
                   max = 100,
                   value = 0.2
                 ),
                 sliderInput(
                   inputId = "percAmbig",
                   label = "Percent missing traits:",
                   min = 0,
                   max = 1,
                   value = 0,
                   ticks = F
                 ),
                 sliderInput(
                   inputId = "numTaxaAmbig",
                   label = "Number of taxa with ambiguous traits:",
                   min = 0,
                   max = 20,
                   value = 0,
                   ticks = F
                 ),
                 sliderInput(
                   inputId = "replicates",
                   label = "Number of replicates for ambiguity:",
                   min = 1,
                   max = 100,
                   value = 50,
                   ticks = F
                 ),
                 selectInput("method", 
                             label = "Inference method", 
                             choices = list(
                               "Parsimony"),
                             selected = NULL,
                             width = "100%"),
                 actionButton("reRunAnalysis", "Run analysis"),
                 actionButton("clearTable", "Clear table")
              ),
  
    mainPanel(
      splitLayout(cellWidths = c("25%", "75%"), plotOutput("tree"), tableOutput("resultsTable"))
    ),
    position = "left"
  )
)

server <- function(input, output) {
  
  observe({
    updateSliderInput(getDefaultReactiveDomain(), "numTaxaAmbig", max = input$nTaxa)
  })
  
  trueTree <- reactive({rcoal(input$nTaxa, tip.label = LETTERS[1:input$nTaxa])})
  output$tree <- renderPlot({
    ggtree(trueTree()) + 
      geom_tiplab(fontface = 4)
  })
  
  #action button script:
  results <- reactiveVal(data.frame("Method" = character(), 
                                    "Replicates" = integer(), 
                                    "Median distance" = numeric(), 
                                    "Proportion correct" = numeric(), 
                                    "Max distance" = numeric()
                                    ))
  
  observeEvent(input$reRunAnalysis,{
    #populate charMat with A, C, G, T (makes ambiguity easier)
    tic("Simulated data")
    datMat <- as.matrix(simSeq(trueTree(), l = input$nChar, type = "DNA", rate = input$rate))
    toc() #print statement to r console so that I can track speed
    
    #replace with ambiguous
    
    resVec <- rep(NA, input$replicates)
    if(input$method == "Parsimony"){
      if(input$numTaxaAmbig == 0 || input$percAmbig == 0){
        tic("BAB no ambig")
        datMat <- phyDat(datMat, type = "DNA")
        
        babTree <- try(bab(datMat), silent = T)
        
        while(class(babTree) =="try-error"){
          datMat <- as.matrix(simSeq(trueTree(), l = input$nChar, type = "DNA", rate = input$rate))
          datMat <- phyDat(datMat, type = "DNA")
          babTree <- try(bab(datMat), silent = T)
          print("in while loop")
        }
        if(class(babTree) == "phylo"){
          print(as.numeric(dist.topo(unroot(babTree), unroot(trueTree()))))
          resVec[1] <- as.numeric(dist.topo(unroot(babTree), unroot(trueTree())))
          
        }else{
          print(as.numeric(dist.topo(unroot(maxCladeCred(babTree)), unroot(trueTree()))))
          resVec[1] <- as.numeric(dist.topo(unroot(maxCladeCred(babTree)), unroot(trueTree())))
          
        }
        toc()
      }else{
        if (input$numTaxaAmbig > nrow(datMat)) {
          stop("Number of ambiguous taxa exceeds the number of taxa")
        }
        for(i in 1:input$replicates){
          whichTaxa <- sample(1:nrow(datMat), input$numTaxaAmbig, replace = F)
          ambigMat <- datMat
          ambigMatSubset <- ambigMat[whichTaxa,]
          
          nAmbigCells <- round(input$percAmbig * length(ambigMatSubset), 0)
          if (nAmbigCells > length(ambigMatSubset)) {
            stop("Number of ambiguous cells exceeds the number of available cells")
          }
          sampledCells <- sample(x = 1:length(ambigMatSubset), size = nAmbigCells)
          ambigMatSubset[sampledCells] <- "?"
          ambigMat[whichTaxa,] <- ambigMatSubset
          
          infMat <- phyDat(ambigMat, type = "DNA")
          
          tic("Tree inference with ambig")
          babTree <- try(bab(infMat), silent = T)
          
          while(class(babTree) =="try-error"){
            
            whichTaxa <- sample(1:nrow(datMat), input$numTaxaAmbig, replace = F)
            ambigMat <- datMat
            ambigMatSubset <- ambigMat[whichTaxa,]
            
            nAmbigCells <- round(input$percAmbig * length(ambigMatSubset), 0)
            if (nAmbigCells > length(ambigMatSubset)) {
              stop("Number of ambiguous cells exceeds the number of available cells")
            }
            sampledCells <- sample(x = 1:length(ambigMatSubset), size = nAmbigCells)
            ambigMatSubset[sampledCells] <- "?"
            ambigMat[whichTaxa,] <- ambigMatSubset
            
            infMat <- phyDat(ambigMat, type = "DNA")
            babTree <- try(bab(infMat), silent = T)
            print("in while loop")
          }
          
          if(length(babTree) == 1){
            resVec[i] <- as.numeric(dist.topo(unroot(babTree), unroot(trueTree())))
          }else{
            resVec[i] <- as.numeric(dist.topo(unroot(maxCladeCred(trueTree())), unroot(trueTree())))
          }
          toc()
          gc()
        }
        
      }
    }else{
      
    }
    
    if(input$numTaxaAmbig == 0 || input$percAmbig == 0){
      propAcc <- NA
    }else{
        propAcc <- sum(resVec == 0 ) / length(resVec)
      }
    
    # Append the new result to the existing results
    newResult <- data.frame(
      "Method" = input$method, 
      "Replicates" = length(na.omit(resVec)), 
      "Median distance" = median(resVec, na.rm = TRUE),
      "Proportion correct" = propAcc,
      "Max distance" = max(resVec, na.rm = T)
    )
    
    currentResults <- results()
    updatedResults <- rbind(currentResults, newResult)
    results(updatedResults)  # Update the reactive value
    
  })
  observeEvent(input$clearTable, {
    results(data.frame("Method" = character(), 
                       "Replicates" = integer(), 
                       "Median distance" = numeric(), 
                       "Proportion correct" = numeric(), 
                       "Max distance" = numeric())
            )
  })
  output$resultsTable <- renderTable({
    results()
  })
  
}

shinyApp(ui, server)