library(shiny)
library(gdsfmt)
library(SNPRelate)

# Define UI for app that draws a histogram ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("MBI-17Z"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      # Input: Slider for the number of bins ----
      sliderInput(inputId = "maf",
                  label = "MAF threshold:",
                  min = 0,
                  max = 0.45,
                  value = 0.05)
      
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: Histogram ----
      plotOutput(outputId = "distPlot")
      
    )
  )
)

# Define server logic required to draw a histogram ----
server <- function(input, output) {
  
  # Histogram of the Old Faithful Geyser Data ----
  # with requested number of bins
  # This expression that generates a histogram is wrapped in a call
  # to renderPlot to indicate that:
  #
  # 1. It is "reactive" and therefore should be automatically
  #    re-executed when inputs (input$bins) change
  # 2. Its output type is a plot
  output$distPlot <- renderPlot({
    
    chr22 <- snpgdsOpen("~/Projekty/MBI/mbi-app/chr22.gds")
    snpsetcsv <- read.csv("~/Projekty/MBI/chr22-pruned.csv")
    snpsetcsv.id = unlist(snpsetcsv$chr22)
    set_maf <- input$maf
    pca <- snpgdsPCA(chr22, snp.id = snpsetcsv.id, num.thread=4, maf = set_maf)
    snpgdsClose(chr22)
    pop_code <- read.table("~/Projekty/MBI/integrated_call_samples_v3.20130502.ALL.panel", sep = "", header = T, nrows = 2504)
    # AFR = 1
    # AMR = 2
    # EAS = 3
    # EUR = 4
    # SAS = 5
    tab <- data.frame(sample.id = pca$sample.id,
                      pop = pop_code$super_pop,
                      EV1 = pca$eigenvect[,1],    # the first eigenvector
                      EV2 = pca$eigenvect[,4],    # the second eigenvector
                      stringsAsFactors = FALSE)
    plot(tab$EV2, tab$EV1, xlab="eigenvector 2", ylab="eigenvector 1", col=as.integer(tab$pop))
    
    
  })
  
}

shinyApp(ui = ui, server = server)













