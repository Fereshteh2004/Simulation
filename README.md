# Load libraries
library(shiny)
library(shinythemes)
library(actuar)
library(VGAM)
library(moments)
library(evd)
library(LaplacesDemon)
# Define UI
ui = fluidPage(
  #Valid themes are: cerulean, cosmo, cyborg, darkly, flatly, journal, lumen, paper, readable, sandstone, simplex, slate, spacelab, superhero, united, yeti.
  theme=shinytheme("superhero"),
  
  titlePanel("Invers_Transform Methode for continuous distributions"),
  
  sidebarLayout(
    sidebarPanel(
      h3("Choose a continuous distribution"),
      selectInput("distribution", "Distribution:",
                  choices = list("Logistic" = "logistic", "Pareto" = "pareto", "Gumbel"="gumbel","Exponential"="exponential"
                                 ,"Lognormal"="lognormal","Beta"="beta","Laplace"="laplace","Negative Exponential"="nexp",
                                 "Weibull"="weibull"
                  )),
      
      # Parameters for Logistic distribution
      conditionalPanel(
        condition = "input.distribution == 'logistic'",
        h3("Logistic Distribution Parameters"),
        numericInput("location", "Location :", value = 5, step = 0.001),
        numericInput("scale", "Scale :", value = 6, min = 0.001, step = 0.001)
      ),
      
      # Parameters for Pareto distribution
      conditionalPanel(
        condition = "input.distribution == 'pareto'",
        h3("Pareto Distribution Parameters"),
        numericInput("scale", "Scale :", value = 5, min = 0.001, step = 0.001),
        numericInput("shape", "shape :", value = 6, min = 0.001, step = 0.001)
      ),
      # Parameters for Gumbel distribution
      conditionalPanel(
        condition = "input.distribution == 'gumbel'",
        h3("Gumbel Distribution Parameters"),
        numericInput("location", "location :", value = 5, step = 0.001),
        numericInput("scale", "Scale :", value = 6, min = 0.001, step = 0.001)
      ),
      # Parameters for Exponential distribution
      conditionalPanel(
        condition = "input.distribution == 'exponential'",
        h3("Exponential Distribution Parameters"),
        numericInput("rate", "Rate :", value = 6, min = 0.001, step = 0.001)
      ),
      # Parameters for Lognormal distribution
      conditionalPanel(
        condition = "input.distribution == 'lognormal'",
        h3("Lognormal Distribution Parameters"),
        numericInput("location", "Location (meanlog):", value = 0.05, step = 0.001),
        numericInput("scale", "Scale (sdlog):", value = 0.06, min = 0.001, step = 0.001)
      ),
      # Parameters for Beta distribution
      conditionalPanel(
        condition = "input.distribution == 'beta'",
        h3("Beta Distribution Parameters"),
        numericInput("aShape", "shape(alpha) :", value = 5, step = 0.001),
        numericInput("bShape", "shape(beta) :", value = 6, min = 0.001, step = 0.001)
      ),
      # Parameters for Laplace distribution
      conditionalPanel(
        condition = "input.distribution == 'laplace'",
        h3("Laplace Distribution Parameters"),
        numericInput("location", "Location :", value = 5, step = 0.001),
        numericInput("scale", "Scale :", value = 6, min = 0.001, step = 0.001)
      ),
      # Parameters for Negative Exponential distribution
      conditionalPanel(
        condition = "input.distribution == 'nexp'",
        h3("Negative Exponential Distribution Parameters"),
        numericInput("rate", "Rate :", value = 6, min = 0.001, step = 0.001)
      ),# Parameters for Weibull distribution
      conditionalPanel(
        condition = "input.distribution == 'weibull'",
        h3("Weibull Distribution Parameters"),
        numericInput("shape", "Shape :", value = 5, min = 0.001, step = 0.001),
        numericInput("scale", "Scale :", value = 6, min = 0.001, step = 0.001)
      ),
      
      
      numericInput("sample_size", "Sample Size:", value = 100, min = 1),
      actionButton("generate", "Generate")
    ),
    
    mainPanel(
      verbatimTextOutput("summary"),
      plotOutput("histPlot"),
      plotOutput("boxPlot"),
      plotOutput("qqPlot")
    )
  )
)

# Define Server
server = function(input, output) {
  
  ppareto=function(x,scale,shape){
    ifelse (x>=scale,1-(scale/x)^shape,0)
  }
  q_laplace= function(u, location, scale) {
    ifelse(u< 0.5,
           location+scale*log(2*u),
           location-scale*log(2*(1-u))
    )
  }
  
  # Function to generate random numbers based on selected distribution
  numbers = eventReactive(input$generate, {
    if (input$distribution == "logistic") {
      u = runif(input$sample_size)
      x = input$location - input$scale * log((1 / u) - 1)
      return(x)
    } else if (input$distribution == "pareto") {
      u = runif(input$sample_size)
      x = input$scale /((1 -u)^(1 / input$shape)) 
      return(x)
    }else if (input$distribution == "gumbel") {
      u = runif(input$sample_size)
      x = input$location-((input$scale)*log(-log(u)) )
      return(x)
      
    }else if (input$distribution == "exponential") {
      u = runif(input$sample_size)
      x = ((-1)*(1/input$rate) )*log(1-u)
      return(x)
    }else if (input$distribution == "lognormal") {
      u = runif(input$sample_size)
      x = qlnorm(u,meanlog = input$location,sdlog=input$scale)
      return(x)
    }else if (input$distribution == "beta") {
      u = runif(input$sample_size)
      x = qbeta(u,input$aShape,input$bShape)
      return(x)
    }else if (input$distribution == "laplace") {
      u = runif(input$sample_size)
      x = q_laplace(u,input$location,input$scale)
      return(x)
    }else if (input$distribution == "nexp") {
      u = runif(input$sample_size)
      x=-log(1-u)/input$rate
      return(x)
    }else if (input$distribution == "weibull") {
      u = runif(input$sample_size)
      x = input$scale*(-log(1-u))^(1/input$shape)
      return(x)
    }
  })
  
  output$summary = renderPrint({
    x = numbers()
    if (input$distribution == "pareto") {
      cat("Generated random numbers (Pareto):\n\n")
      print(x)
      
      mean_x = mean(x)
      pareto_mean = input$shape *input$scale / (input$shape - 1)
      cat("\n\n------------------------------------Comparison of means-------------------------------------")
      cat("\nMean of generated random numbers : ", mean_x, "\nReal mean : ", pareto_mean, "\n")
      t_test_result = t.test(x, mu = pareto_mean)
      cat("T-test result:\n")
      print(t_test_result)
      
      if (t_test_result$p.value > 0.05) {
        cat("The means are not significantly different (p-value > 0.05).\n")
      } else {
        cat("The means are significantly different (p-value <= 0.05).\n")
      }
      
      cat("\n\n------------------------------------Explanation for plots---------------------------------- ")
      cat("\nSkewness of the generated numbers:", skewness(x), "\n")
      if (skewness(x) < 0) {
        cat("\nThe histogram and boxplot of generated random numbers show negative skewness (the tail is on the left).\n")
      } else if (skewness(x) > 0) {
        cat("\nThe histogram and boxplot of generated random numbers show positive skewness (the tail is on the right).\n")
      } else {
        cat("\nThe histogram and boxplot of generated random numbers are symmetrical.\n")
      }
      
      cat("\n-----------------------------------Kolmogorov-Smirnov Test----------------------------------")
      
      ks_test = ks.test(x, "ppareto",input$scale ,input$shape )
      cat("\nKolmogorov-Smirnov test result:\n")
      print(ks_test)
      
      if (ks_test$p.value > 0.05) {
        cat("we can tell that the QQ plot is a straight line, so generated random numbers are likely to be following a Pareto distribution with scale parameter =", input$scale, "and shape parameter =", input$shape, "(p-value > 0.05).\n")
      } else {
        cat("The generated data does not follow a Pareto distribution (p-value <= 0.05).\n")
      }
      #====================================================================================================================================
    } else if (input$distribution == "logistic") {
      cat("Generated random numbers (Logistic):\n\n")
      print(x)
      
      mean_x = mean(x)
      logistic_mean = input$location
      cat("\n\n------------------------------------Comparison of means-------------------------------------")
      cat("\nMean of generated random numbers : ", mean_x, "\nReal mean : ", logistic_mean, "\n")
      t_test_result = t.test(x, mu = logistic_mean)
      cat("T-test result:\n")
      print(t_test_result)
      
      if (t_test_result$p.value > 0.05) {
        cat("The means are not significantly different (p-value > 0.05).\n")
      } else {
        cat("The means are significantly different (p-value <= 0.05).\n")
      }
      
      cat("\n\n------------------------------------Explanation for plots---------------------------------- ")
      cat("\nSkewness of the generated numbers:", skewness(x), "\n")
      if (skewness(x) < 0) {
        cat("\nThe histogram and boxplot of generated random numbers show negative skewness (the tail is on the left).\n")
      } else if (skewness(x) > 0) {
        cat("\nThe histogram and boxplot of generated random numbers show positive skewness (the tail is on the right).\n")
      } else {
        cat("\nThe histogram and boxplot of generated random numbers are symmetrical.\n")
      }
      
      cat("\n-----------------------------------Kolmogorov-Smirnov Test----------------------------------")
      ks_test = ks.test(x, "plogis", input$location, input$scale)
      cat("\nKolmogorov-Smirnov test result:\n")
      print(ks_test)
      
      if (ks_test$p.value > 0.05) {
        cat("we can tell that the QQ plot is a straight line, so generated random numbers are likely to be following a Logistic distribution with location parameter=", input$location, "and scale parameter =", input$scale, "(p-value > 0.05).\n")
      } else {
        cat("The generated data does not follow a Logistic distribution (p-value <= 0.05).\n")
      }
      
    } else if (input$distribution == "gumbel") {
      x = numbers()
      cat("Generated random numbers (Gumbel):\n\n")
      print(x)
      
      mean_x = mean(x)
      gumbel_mean = input$location+input$scale*0.577216
      cat("\n\n------------------------------------Comparison of means-------------------------------------")
      cat("\nMean of generated random numbers : ", mean_x, "\nReal mean : ", gumbel_mean, "\n")
      t_test_result = t.test(x, mu = gumbel_mean)
      cat("T-test result:\n")
      print(t_test_result)
      
      if (t_test_result$p.value > 0.05) {
        cat("The means are not significantly different (p-value > 0.05).\n")
      } else {
        cat("The means are significantly different (p-value <= 0.05).\n")
      }
      
      cat("\n\n------------------------------------Explanation for plots---------------------------------- ")
      cat("\nSkewness of the generated numbers:", skewness(x), "\n")
      if (skewness(x) < 0) {
        cat("\nThe histogram and boxplot of generated random numbers show negative skewness (the tail is on the left).\n")
      } else if (skewness(x) > 0) {
        cat("\nThe histogram and boxplot of generated random numbers show positive skewness (the tail is on the right).\n")
      } else {
        cat("\nThe histogram and boxplot of generated random numbers are symmetrical.\n")
      }
      
      cat("\n-----------------------------------Kolmogorov-Smirnov Test----------------------------------")
      ks_test = ks.test(x, "pgumbel", input$location, input$scale)
      cat("\nKolmogorov-Smirnov test result:\n")
      print(ks_test)
      
      if (ks_test$p.value > 0.05) {
        cat("we can tell that the QQ plot is a straight line, so generated random numbers are likely to be following a gumbel distribution with location parameter =", input$location, "and scale parameter =", input$scale, "(p-value > 0.05).\n")
      } else {
        cat("The generated data does not follow a gumbel distribution (p-value <= 0.05).\n")
      }
      
    } else if (input$distribution == "exponential") {
      cat("Generated random numbers (exponential):\n\n")
      print(x)
      
      mean_x = mean(x)
      exponential_mean = 1/input$rate
      cat("\n\n------------------------------------Comparison of means-------------------------------------")
      cat("\nMean of generated random numbers : ", mean_x, "\nReal mean : ", exponential_mean, "\n")
      t_test_result = t.test(x, mu = exponential_mean)
      cat("T-test result:\n")
      print(t_test_result)
      
      if (t_test_result$p.value > 0.05) {
        cat("The means are not significantly different (p-value > 0.05).\n")
      } else {
        cat("The means are significantly different (p-value <= 0.05).\n")
      }
      
      cat("\n\n------------------------------------Explanation for plots---------------------------------- ")
      cat("\nSkewness of the generated numbers:", skewness(x), "\n")
      if (skewness(x) < 0) {
        cat("\nThe histogram and boxplot of generated random numbers show negative skewness (the tail is on the left).\n")
      } else if (skewness(x) > 0) {
        cat("\nThe histogram and boxplot of generated random numbers show positive skewness (the tail is on the right).\n")
      } else {
        cat("\nThe histogram and boxplot of generated random numbers are symmetrical.\n")
      }
      
      cat("\n-----------------------------------Kolmogorov-Smirnov Test----------------------------------")
      ks_test = ks.test(x, "pexp", input$rate)
      cat("\nKolmogorov-Smirnov test result:\n")
      print(ks_test)
      
      if (ks_test$p.value > 0.05) {
        cat("we can tell that the QQ plot is a straight line, so generated random numbers are likely to be following a exponential distribution with location parameter=", input$location, "and scale parameter =", input$scale, "(p-value > 0.05).\n")
      } else {
        cat("The generated data does not follow a exponential distribution (p-value <= 0.05).\n")
      }
      
    }else if (input$distribution == "lognormal") {
      cat("Generated random numbers (Lognormal):\n\n")
      print(x)
      
      mean_x = mean(x)
      lognormal_mean = exp(input$location+ ((input$scale)^2)/2)
      cat("\n\n------------------------------------Comparison of means-------------------------------------")
      cat("\nMean of generated random numbers : ", mean_x, "\nReal mean : ", lognormal_mean, "\n")
      t_test_result = t.test(x, mu = lognormal_mean)
      cat("T-test result:\n")
      print(t_test_result)
      
      if (t_test_result$p.value > 0.05) {
        cat("The means are not significantly different (p-value > 0.05).\n")
      } else {
        cat("The means are significantly different (p-value <= 0.05).\n")
      }
      
      cat("\n\n------------------------------------Explanation for plots---------------------------------- ")
      cat("\nSkewness of the generated numbers:", skewness(x), "\n")
      if (skewness(x) < 0) {
        cat("\nThe histogram and boxplot of generated random numbers show negative skewness (the tail is on the left).\n")
      } else if (skewness(x) > 0) {
        cat("\nThe histogram and boxplot of generated random numbers show positive skewness (the tail is on the right).\n")
      } else {
        cat("\nThe histogram and boxplot of generated random numbers are symmetrical.\n")
      }
      
      cat("\n-----------------------------------Kolmogorov-Smirnov Test----------------------------------")
      ks_test = ks.test(x, "plnorm", input$location, input$scale)
      cat("\nKolmogorov-Smirnov test result:\n")
      print(ks_test)
      
      if (ks_test$p.value > 0.05) {
        cat("we can tell that the QQ plot is a straight line, so generated random numbers are likely to be following a Lognormal distribution with location parameter=", input$location, "and scale parameter =", input$scale, "(p-value > 0.05).\n")
      } else {
        cat("The generated data does not follow a Lognormal distribution (p-value <= 0.05).\n")
      }
      
    }else if (input$distribution == "beta") {
      cat("Generated random numbers (Beta):\n\n")
      print(x)
      
      mean_x = mean(x)
      beta_mean = input$aShape/(input$aShape+input$bShape)
      cat("\n\n------------------------------------Comparison of means-------------------------------------")
      cat("\nMean of generated random numbers : ", mean_x, "\nReal mean : ", beta_mean, "\n")
      t_test_result = t.test(x, mu = beta_mean)
      cat("T-test result:\n")
      print(t_test_result)
      
      if (t_test_result$p.value > 0.05) {
        cat("The means are not significantly different (p-value > 0.05).\n")
      } else {
        cat("The means are significantly different (p-value <= 0.05).\n")
      }
      
      cat("\n\n------------------------------------Explanation for plots---------------------------------- ")
      cat("\nSkewness of the generated numbers:", skewness(x), "\n")
      if (skewness(x) < 0) {
        cat("\nThe histogram and boxplot of generated random numbers show negative skewness (the tail is on the left).\n")
      } else if (skewness(x) > 0) {
        cat("\nThe histogram and boxplot of generated random numbers show positive skewness (the tail is on the right).\n")
      } else {
        cat("\nThe histogram and boxplot of generated random numbers are symmetrical.\n")
      }
      
      cat("\n-----------------------------------Kolmogorov-Smirnov Test----------------------------------")
      ks_test = ks.test(x, "pbeta", input$aShape, input$bShape)
      cat("\nKolmogorov-Smirnov test result:\n")
      print(ks_test)
      
      if (ks_test$p.value > 0.05) {
        cat("we can tell that the QQ plot is a straight line, so generated random numbers are likely to be following a Beta distribution with shape parameter (alpha) = ", input$aShape, "and shape parameter (beta)=", input$bShape, "(p-value > 0.05).\n")
      } else {
        cat("The generated data does not follow a Beta distribution (p-value <= 0.05).\n")
      }
      
    }else if (input$distribution == "laplace") {
      cat("Generated random numbers (Laplace):\n\n")
      print(x)
      
      mean_x = mean(x)
      laplace_mean = input$location
      cat("\n\n------------------------------------Comparison of means-------------------------------------")
      cat("\nMean of generated random numbers : ", mean_x, "\nReal mean : ", laplace_mean, "\n")
      t_test_result = t.test(x, mu = laplace_mean)
      cat("T-test result:\n")
      print(t_test_result)
      
      if (t_test_result$p.value > 0.05) {
        cat("The means are not significantly different (p-value > 0.05).\n")
      } else {
        cat("The means are significantly different (p-value <= 0.05).\n")
      }
      
      cat("\n\n------------------------------------Explanation for plots---------------------------------- ")
      cat("\nSkewness of the generated numbers:", skewness(x), "\n")
      if (skewness(x) < 0) {
        cat("\nThe histogram and boxplot of generated random numbers show negative skewness (the tail is on the left).\n")
      } else if (skewness(x) > 0) {
        cat("\nThe histogram and boxplot of generated random numbers show positive skewness (the tail is on the right).\n")
      } else {
        cat("\nThe histogram and boxplot of generated random numbers are symmetrical.\n")
      }
      
      cat("\n-----------------------------------Kolmogorov-Smirnov Test----------------------------------")
      ks_test = ks.test(x, "plaplace", input$location, input$scale)
      cat("\nKolmogorov-Smirnov test result:\n")
      print(ks_test)
      
      if (ks_test$p.value > 0.05) {
        cat("we can tell that the QQ plot is a straight line, so generated random numbers are likely to be following a Laplace distribution with shape parameter (alpha) = ", input$aShape, "and shape parameter (beta)=", input$bShape, "(p-value > 0.05).\n")
      } else {
        cat("The generated data does not follow a Laplace distribution (p-value <= 0.05).\n")
      }
      
    }else if (input$distribution == "nexp") {
      cat("Generated random numbers (Negative Exponential):\n\n")
      print(x)
      
      mean_x = mean(x)
      Nexponential_mean = 1/input$rate
      cat("\n\n------------------------------------Comparison of means-------------------------------------")
      cat("\nMean of generated random numbers : ", mean_x, "\nReal mean : ", Nexponential_mean, "\n")
      t_test_result = t.test(x, mu = Nexponential_mean)
      cat("T-test result:\n")
      print(t_test_result)
      
      if (t_test_result$p.value > 0.05) {
        cat("The means are not significantly different (p-value > 0.05).\n")
      } else {
        cat("The means are significantly different (p-value <= 0.05).\n")
      }
      
      cat("\n\n------------------------------------Explanation for plots---------------------------------- ")
      cat("\nSkewness of the generated numbers:", skewness(x), "\n")
      if (skewness(x) < 0) {
        cat("\nThe histogram and boxplot of generated random numbers show negative skewness (the tail is on the left).\n")
      } else if (skewness(x) > 0) {
        cat("\nThe histogram and boxplot of generated random numbers show positive skewness (the tail is on the right).\n")
      } else {
        cat("\nThe histogram and boxplot of generated random numbers are symmetrical.\n")
      }
      
      cat("\n-----------------------------------Kolmogorov-Smirnov Test----------------------------------")
      ks_test = ks.test(x, "pexp", input$rate)
      cat("\nKolmogorov-Smirnov test result:\n")
      print(ks_test)
      
      if (ks_test$p.value > 0.05) {
        cat("we can tell that the QQ plot is a straight line, so generated random numbers are likely to be following a Negative Exponential distribution with location parameter=", input$location, "and scale parameter =", input$scale, "(p-value > 0.05).\n")
      } else {
        cat("The generated data does not follow a Negative Exponential distribution (p-value <= 0.05).\n")
      }
      
    }else if (input$distribution == "weibull") {
      cat("Generated random numbers (Weibull):\n\n")
      print(x)
      
      mean_x = mean(x)
      weibull_mean = input$scale*gamma(1+1/input$shape)
      cat("\n\n------------------------------------Comparison of means-------------------------------------")
      cat("\nMean of generated random numbers : ", mean_x, "\nReal mean : ", weibull_mean, "\n")
      t_test_result = t.test(x, mu = weibull_mean)
      cat("T-test result:\n")
      print(t_test_result)
      
      if (t_test_result$p.value > 0.05) {
        cat("The means are not significantly different (p-value > 0.05).\n")
      } else {
        cat("The means are significantly different (p-value <= 0.05).\n")
      }
      
      cat("\n\n------------------------------------Explanation for plots---------------------------------- ")
      cat("\nSkewness of the generated numbers:", skewness(x), "\n")
      if (skewness(x) < 0) {
        cat("\nThe histogram and boxplot of generated random numbers show negative skewness (the tail is on the left).\n")
      } else if (skewness(x) > 0) {
        cat("\nThe histogram and boxplot of generated random numbers show positive skewness (the tail is on the right).\n")
      } else {
        cat("\nThe histogram and boxplot of generated random numbers are symmetrical.\n")
      }
      
      cat("\n-----------------------------------Kolmogorov-Smirnov Test----------------------------------")
      ks_test = ks.test(x, "pweibull", input$shape, input$scale)
      cat("\nKolmogorov-Smirnov test result:\n")
      print(ks_test)
      
      if (ks_test$p.value > 0.05) {
        cat("we can tell that the QQ plot is a straight line, so generated random numbers are likely to be following a Weibull distribution with shape parameter (alpha) = ", input$aShape, "and shape parameter (beta)=", input$bShape, "(p-value > 0.05).\n")
      } else {
        cat("The generated data does not follow a Weibull distribution (p-value <= 0.05).\n")
      }
      
    }
    
    
    
    
  })
  
  # Histogram plot
  output$histPlot = renderPlot({
    x = numbers()
    hist(x, col = "light blue", main = paste("Histogram of", input$distribution, "Distribution"))
  })
  
  # Box plot
  output$boxPlot = renderPlot({
    x = numbers()
    boxplot(x, col = "light green", main = paste("Boxplot of", input$distribution, "Distribution"))
  })
  
  # QQ plot
  output$qqPlot = renderPlot({
    x = numbers()
    if (input$distribution == "logistic") {
      qq_plot = qlogis(ppoints(input$sample_size), input$location, input$scale)
    } else if(input$distribution == "pareto"){
      qq_plot = qpareto(ppoints(input$sample_size), input$scale,input$shape )
    } else if(input$distribution == "gumbel"){
      qq_plot = qgumbel(ppoints(input$sample_size), input$location, input$shape)
    }else if(input$distribution == "exponential"){
      qq_plot = qexp(ppoints(input$sample_size), input$rate)
    }else if(input$distribution == "lognormal"){
      qq_plot = qlnorm(ppoints(input$sample_size), input$location, input$scale)
    }else if(input$distribution == "beta"){
      qq_plot = qbeta(ppoints(input$sample_size), input$aShape, input$bShape)
    }else if(input$distribution == "laplace"){
      qq_plot = qlaplace(ppoints(input$sample_size), input$location, input$scale)
    }else if(input$distribution == "nexp"){
      qq_plot = qexp(ppoints(input$sample_size), input$rate)
    }else if(input$distribution == "weibull"){
      qq_plot = qweibull(ppoints(input$sample_size), input$shape,input$scale)
    }
    
    plot(qq_plot, sort(x))
    abline(lm(sort(x) ~ qq_plot), col = "dark red")
  })
}

#to combine ui and server into an interactive app!
shinyApp(ui=ui,server=server)
