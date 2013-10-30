library(shiny)
library(shinyAce)



shinyServer(function(input, output) {


    
    options(warn=-1)
    
    
    anovakun <- reactive({
        
        dat <- read.csv(text=input$text, sep="\t")
        
    
# insert the source code of anovakun here by downloading it from:
# http://www11.atpages.jp/~riseki/pukiwikiplus/index.php?ANOVA%B7%AF
    
    
    # One-way ANOVA
        
        if (input$factor == "oneway") {

            type <- switch(input$one.design,
                        Between = "As",
                        Within = "sA")
            
            level1 <- input$factor1.level

            anovakun(dat, type, level1, mau=T, auto=T, holm = T, peta=T)
            
        } else { # Two-way ANOVA
            
            type <- switch(input$two.design,
                        Factor1Between_Factor2Between = "ABs",
                        Factor1Between_Factor2Within = "AsB",
                        Factor1Within_Factor2Within = "sAB")
            
            level1 <- input$factor1.level
            level2 <- input$factor2.level
            
            anovakun(dat, type, level1, level2, mau=T, auto=T, holm = T, peta=T)
            
        }

    })




    output$anovakun.out <- renderPrint({
        anovakun()
    })
        

})
