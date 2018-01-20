library(shiny)
library(gdsfmt)
library(SNPRelate)

# Define UI for app that draws a histogram ----
ui <- fluidPage(
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      h1("Settings"),
      
      selectInput(inputId = "format", 
                  label = "File format", 
                  choices = c("GDS", "VCF"), 
                  selected = "GDS"),
      
      textInput(inputId = "vcf", 
                label = "Path to genome file", 
                value = "~/Projekty/MBI/mbi-app/chr22.gds"),
      
      textInput(inputId = "pop", 
                label = "Path to file with population labels", 
                value = "~/Projekty/MBI/mbi-app/chr22.pop"),
      
      textInput(inputId = "ranges", 
                label = "Path to file with chromosome ranges to analyze", 
                value = ""),
      
      checkboxInput(inputId = "ld", 
                    label = "Use LD pruning (may take a longer time to complete)", 
                    value = FALSE),
      
      sliderInput(inputId = "maf",
                  label = "MAF threshold:",
                  min = 0,
                  max = 0.5,
                  value = 0.3),
      
      selectInput(inputId = "ev1", 
                  label = "First eigenvector (horizontal axis)", 
                  choices = c(1:20), 
                  selected = 1),
      
      selectInput(inputId = "ev2", 
                  label = "Second eigenvector (vertical axis)", 
                  choices = c(1:20), 
                  selected = 2),
      
      checkboxInput(inputId = "clust", 
                    label = "Generate dendrogram (hierarchical clustering) on the data", 
                    value = FALSE)
      
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      plotOutput(outputId = "pcaPlot"),
      plotOutput(outputId = "dendPlot")
    )
  )
)

# Define server logic required to draw a histogram ----
server <- function(input, output) {
  
  # Make PCA plot
  output$pcaPlot <- renderPlot({
    
    # Load VCF file
    if(input$format == "VCF") {
      snpgdsVCF2GDS(input$vcf, paste0(input$vcf,".gds"), method="biallelic.only")
      chr22 <- snpgdsOpen(paste0(input$vcf,".gds"))
    }
    # Load GDS file
    if(input$format == "GDS") {
      chr22 <- snpgdsOpen(input$vcf)
    }
    
    # Do LD-pruning if required
    if(input$ld) {
      snpset <- snpgdsLDpruning(
        chr22, 
        ld.threshold=0.2, 
        num.thread = 4,
        maf = input$maf)
      snpset.id <- unlist(snpset)
    }
    else {
      snpset.id = NULL
    }
    
    # Choose chromosome ranges of required
    if(input$ranges != "") {
      rng <- read.table(input$ranges, sep = '\t', header = TRUE)
      num_rng = length(rng$chromosomes)
      selected = c()
    }
    
    # Calculate PCA on samples
    pca <- snpgdsPCA(
      chr22, 
      snp.id = snpset.id, 
      num.thread=4, 
      maf = input$maf)
    # Close unnecessary GDS file
    snpgdsClose(chr22)
    # Read population codes
    # AFR = 1
    # AMR = 2
    # EAS = 3
    # EUR = 4
    # SAS = 5
    pop_code <- read.table(input$pop, 
                           sep = "", 
                           header = T, 
                           nrows = 2504)
    # Make a plot
    eigenvector1 <- strtoi(input$ev1)
    eigenvector2 <- strtoi(input$ev2)
    tab <- data.frame(sample.id = pca$sample.id,
                      pop = pop_code$super_pop,
                      EV1 = pca$eigenvect[,eigenvector1],    # the first eigenvector
                      EV2 = pca$eigenvect[,eigenvector2],    # the second eigenvector
                      stringsAsFactors = FALSE)
    plot(tab$EV2, 
         tab$EV1, 
         xlab="Eigenvector 1", 
         ylab="Eigenvector 2", 
         col=as.integer(tab$pop), 
         main = paste0("Populations PCA, MAF >= ", input$maf))
    
    legend("bottomright", 
           legend=levels(pop_code$super_pop), 
           text.col=1:nlevels(pop_code$super_pop))
    
  })
  
  output$dendPlot = renderPlot({
    
    if(input$clust) {
      # Load VCF file
      if(input$format == "VCF") {
        snpgdsVCF2GDS(input$vcf, paste0(input$vcf,".gds"), method="biallelic.only")
        chr22 <- snpgdsOpen(paste0(input$vcf,".gds"))
      }
      # Load GDS file
      if(input$format == "GDS") {
        chr22 <- snpgdsOpen(input$vcf)
      }
      
      # Do LD-pruning if required
      if(input$ld) {
        snpset <- snpgdsLDpruning(
          chr22, 
          ld.threshold=0.2, 
          num.thread = 4,
          maf = input$maf)
        snpset.id <- unlist(snpset)
      }
      else {
        snpset.id = NULL
      }
      
      # Read population codes
      # AFR = 1
      # AMR = 2
      # EAS = 3
      # EUR = 4
      # SAS = 5
      pop_code <- read.table(input$pop, 
                             sep = "", 
                             header = T, 
                             nrows = 2504)
      pop <- pop_code$super_pop
      # Do IBS anaysis
      ibs <- snpgdsIBS(chr22, 
                       snp.id = snpset.id, 
                       maf = input$maf,
                       num.thread=4)
      # Close unnecessary GDS file
      snpgdsClose(chr22)
      # Do hierarchical clustering
      ibs.hc <- snpgdsHCluster(ibs)
      rv <- snpgdsCutTree(ibs.hc, samp.group=as.factor(pop))
      
      plot(rv$dendrogram, 
           leaflab="none", 
           main=paste0("Populations dendrogram, MAF >= ", 
                       input$maf, ", ", 
                       length(levels(as.factor(pop))), " groups"))
      legend("bottomright", legend=levels(pop), text.col=1:nlevels(pop))
    }
  })
  
}

shinyApp(ui = ui, server = server)













