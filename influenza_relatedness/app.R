#Shiny influenza INTERACTIVE heatmaps
#Author: Lauren Tindale

library(shiny)
library(tidyverse)
library(plotly)
library(DT)

#read in heatmap data
data1 <- read.csv("Cohort_Analysis_Heatmap_Data_normalized_comparison.csv")
data1$Year <- factor(data1$Year)
data1$Position <- factor(data1$Position)
data1$GlycanShieldingScore[data1$GlycanShieldingScore == 0] <- NA

#reorder Position column
data1$Position <- fct_relevel(data1$Position, "131", "155", "156", "157", "158", "159", "160", "197", "128", "189", "190", "192", "193", "133", "135", "137", "138", "144", "145", "225", "226", "227")

#rename columns
colnames(data1)[which(colnames(data1) == 'MatchScore_Historic_vs_2a1b.131K')] <- '3C.2a1b_131K'
colnames(data1)[which(colnames(data1) == 'MatchScore_Historic_vs_3a')] <- '3C.3a'
colnames(data1)[which(colnames(data1) == 'MatchScore_Historic_vs_IVR.186')] <- '3C.2a1_IVR-186'
colnames(data1)[which(colnames(data1) == 'BlendedMatchScore_2a1b.131K')] <- '3C.2a1b_131K + 3C.2a1_IVR-186'
colnames(data1)[which(colnames(data1) == 'BlendedMatchScore_3a')] <- '3C.3a + 3C.2a1_IVR-186'

#read in bargraph data
bardata <- read.csv("Cohort_Analysis_Histogram_Data.csv", header = TRUE)

bsub <- na_if(bardata, 'Unknown (X)') #change "Unknown (X)" values to NA
bsub <- na.exclude(bsub) #exclude NA
colnames(bsub) <- sub("X", "AA", colnames(bsub)) #remove "X" from start of column names


########################################

ui <- navbarPage(
  title = 'Heatmaps',
  tabPanel(title = 'Influenza subtype',
           tags$h1('Relatedness between HA1 residues in historic relative to contemporary influenza A(H3N2) viruses'),
           p('Heatmap of historic-contemporary match scores. The historic AA is the the most frequent AA for each combination of year and AA position. The contemporary AA is the current variant consensus sequence at the same position.'),
           p('Scores are based on a normalized BLOSUM80 matrix where a match = 1.'),
           selectInput('typeInput1', 'Influenza subtype:', names(data1[c("3C.2a1b_131K", "3C.3a", "3C.2a1_IVR-186")])),
           radioButtons('AAoptions1', 'HA1 residues to show:', choices = list('All', 'Koel positions')),
           plotlyOutput('plot1', width ='850px', height = '400px')
  ),
  
  tabPanel(title = 'Blended scores',
           tags$h1('Relatedness between HA1 residues in historic relative to contemporary influenza A(H3N2) viruses'),
           p('Blended Score = (historic vs vaccine score) + (historic vs contemporary clade score) + (contemporary clade vs vaccine score)'),
           p('Historic-Contemporary-Vaccine Match and Mismatch were assigned discrete 1/0 values and added together to generate a blended 3-way theoretical score and to compare to BLOSUM80 matrix scoring in heatmaps.'),
           p('0 = historic, contemporary and vaccine all mismatch'),
           p('1 = only one combination matches (historic-contemporary or historic-vaccine or contemporary-vaccine)'),
           p('3 = historic, contemporary and vaccine all match'),
           selectInput('typeInput2', 'Influenza subtype:', names(data1[c("3C.2a1b_131K + 3C.2a1_IVR-186", "3C.3a + 3C.2a1_IVR-186")])),
           radioButtons('AAoptions2', 'HA1 residues to show:', choices = list('All', 'Koel positions')),
           plotlyOutput('plot2', width ='850px', height = '400px')
  ),
  
  tabPanel(title = 'Historic amino acids',
           tags$h1('Relatedness between HA1 residues in historic relative to contemporary influenza A(H3N2) viruses'),
           selectInput('typeInput3', 'Amino acid residue:', names(bsub[c("AA133", "AA135", "AA145", "AA155", "AA156", "AA158", "AA159", "AA189", "AA193")])),
           plotlyOutput('bargraph', width ='850px', height = '400px')
  
#  ),
  
#  tabPanel(title = 'Data',
#           tags$h1('Relatedness between HA1 residues in historic relative to contemporary influenza A(H3N2) viruses'),
#           DTOutput('results')
  )
)

########################################


server <- function(input, output) {
  
  #subset Koel positions
  Koel <- subset(data1, Position=='133' | Position=='135' | Position=='145' | Position=='155' | Position=='156' | Position=='158' | Position=='159' | Position=='189' | Position=='193')
  
  #create reactive variable - full dataset (data1) or Koel positions only
  AA1 <- reactive({
    req(input$AAoptions1)
    if(input$AAoptions1 == 'All'){
      data1
    } else{
      Koel
    }
  })
  
  AA2 <- reactive({
    req(input$AAoptions2)
    if(input$AAoptions2 == 'All'){
      data1
    } else{
      Koel
    }
  })
  
  bar_input <- reactive({
    input$typeInput3
  })
  
  
  output$plot1 <- renderPlotly({
    p1 <- ggplot(AA1(), aes(x = Year, y = Position, fill = AA1()[[input$typeInput1]],
                            text = paste(" Year: ", Year,
                                         "<br> Position: ", Position,
                                         "<br> Subtype Match Score: ", AA1()[[input$typeInput1]],
                                         "<br> Glycan Shielding Score: ", GlycanShieldingScore))) +
      geom_tile(colour = "white") +  #line colour between tiles
      xlab(label = "Year") +
      ylab(label = "HA1 Residue") +
      ggtitle(label = input$typeInput1) +
      labs(fill = "Relatedness") + #tile fill legend label
      theme(plot.title = element_text(hjust = 0.5)) + #centre main title
      theme(axis.text.x = element_text(angle =60, hjust = 0.6, size = 6.5),
            axis.ticks.x = element_blank(), #remove x axis ticks
            axis.ticks.y = element_blank()) + #remove y axis ticks
      scale_fill_distiller(palette = "RdYlBu", #choose colour of tile fill
                           direction=0, #flip colour direction of palette
                           limits=c(min(-0.75), max(1))) + #set min max of legend scale
      #theme(legend.position = "bottom") +
      geom_point(aes (size = GlycanShieldingScore), na.rm = TRUE) + #overlap point data
      scale_size(range = c(0,2.5)) + #choose range size of points
      theme(panel.background = element_rect(fill = "white"))
    
    ggplotly(p1, tooltip = 'text')
    
  })
  
  output$plot2 <- renderPlotly({
    p2 <- ggplot(AA2(), aes(x = Year, y = Position, fill = AA2()[[input$typeInput2]],
                            text = paste(" Year: ", Year,
                                         "<br> Position: ", Position,
                                         "<br> Subtype Match Score: ", AA2()[[input$typeInput2]],
                                         "<br> Glycan Shielding Score: ", GlycanShieldingScore))) +
      geom_tile(colour = "white") +  #line colour between tiles
      xlab(label = "Year") +
      ylab(label = "HA1 Residue") +
      ggtitle(label = input$typeInput2) +
      labs(fill = "Relatedness") + #tile fill legend label
      theme(plot.title = element_text(hjust = 0.5)) + #centre main title
      theme(axis.text.x = element_text(angle =60, hjust = 0.6, size = 6.5),
            axis.ticks.x = element_blank(), #remove x axis ticks
            axis.ticks.y = element_blank()) + #remove y axis ticks
      scale_fill_distiller(palette = "RdYlBu", #choose colour of tile fill
                           direction=1) +
      geom_point(aes (size = GlycanShieldingScore), na.rm = TRUE) + #overlap point data
      scale_size(range = c(0,2.5)) + #choose range size of points
      theme(panel.background = element_rect(fill = "white"))
    
    ggplotly(p2, tooltip = 'text')
    
  })
  
  output$bargraph <- renderPlotly({
    bargraph <- bsub %>% 
      ggplot(aes(as.factor(Year),fill = bar_input())) +
      geom_bar(position="fill") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 10)) +
      xlab("Year") +
      ylab("Percent of GISAID Sequences") +
      ggtitle(label = input$typeInput3) +
      scale_y_continuous(labels = scales::percent_format())
    
    ggplotly(bargraph)
    
  })  
  
  
  
  
  #remove intermediate data columns
#  subset <- select(data1, -'MatchScore_2a1b.131K_vs_IVR.186', -'MatchScore_3a_vs_IVR.186')
  
  #datatable
#  output$results <- DT::renderDT({
#    datatable(
#      subset,
#      colnames = c('Glycan Shielding Score' = 'GlycanShieldingScore'),
#      caption = 'Data table of glycan shielding scores and historic-contemporary match scores by year and AA position.',
#      filter = 'top',
#      options = list(
#        pageLength = 25
        #columnDefs = list(list(targets = -seq_len(6), searchable = FALSE)))
#      ))
#  })
}

shinyApp(ui = ui, server = server)