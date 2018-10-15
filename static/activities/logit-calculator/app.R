#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinythemes)

# Define UI for application that draws a histogram
ui <- navbarPage(
    "Calculate probability of Island Mouse occurrence",
    tabPanel("Calculator",
    
    fluidPage(theme = "flatly",
   
   # Application title
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
          fluidRow(
              column(6, 
                     numericInput("beta0",
                                  "Intercept",
                                  min = -100,
                                  max = 100,
                                  value = -25,
                                  step = 0.1))),
          fluidRow(
              column(6, 
                     numericInput("beta1",
                                  "Beta 1",
                                  min = -100,
                                  max = 100,
                                  value = 2.3,
                                  step = 0.1)),
              column(6,
                     numericInput("x1",
                                  "Average air temp.",
                                  min = -100,
                                  max = 100,
                                  value = 23.5,
                                  step = 0.1)
                     )
          ),
          fluidRow(
              column(6,
                     numericInput("beta2",
                                  "Beta 2",
                                  min = -100,
                                  max = 100,
                                  value = 1.2,
                                  step = 0.1)),
              column(6,
                     numericInput("x2",
                                  "Percent rock/sand/clay cover",
                                  min = -100,
                                  max = 100,
                                  value = 30,
                                  step = 0.1))
          ),
          fluidRow(
              column(6,
                     numericInput("beta3",
                                  "Beta 3",
                                  min = -100,
                                  max = 100,
                                  value = -4.3,
                                  step = 0.1)),
              column(6,
                     numericInput("x3",
                                  "Temperature range",
                                  min = -100,
                                  max = 100,
                                  value = 15,
                                  step = 0.1))
            ),
          fluidRow(
              column(6, offset = 3, 
                     actionButton("go",
                                  "Calculate!")
                     )
              
          )
          ),
                     
      
      # Show a plot of the generated distribution
      mainPanel(
         textOutput("p.occ")
      )
   )
)),

tabPanel("What's going on here?",
         includeHTML("equation.html")
         )
)

server <- function(input, output) {
   
    p <- eventReactive(input$go, {
        
        x1 = scale(input$x1)
        x2 = scale(input$x2)
        x3 = scale(input$x3)
        
        plogis(input$beta0 + input$beta1*input$x1 + input$beta2*input$x2 + input$beta3*input$x3)
    })
    
    
   output$p.occ <- renderText({
       paste("Probability of occurrence = ", round(p(), 3))
   })
   
}

# Run the application 
shinyApp(ui = ui, server = server)

