#distance intensity profile


library(shiny)
library(shinythemes)
library(shinyjs)


# Define UI for application that draws a histogram
ui <- fluidPage( theme = shinytheme("slate"),
                 useShinyjs(),
                 # Application title
                 titlePanel("Fluorescence profile along distance from nucleus"),
                 
                 fluidRow(column(3, # width is 3/12
                                 h2("1. Choose file"),
                                 fileInput("file1", "Choose .tif to analyze", multiple = FALSE, accept = NULL, width = NULL)),
                          
                          column(9, # with is 9/12
                                 textOutput("status_text"), # status text
                                 actionButton("reset_zoom", "Reset zoom"),
                                 fluidRow(column(3, # column divided into quarters
                                                 plotOutput("plot1", click = "plot_click", dblclick = "plot1_dblclick", brush = brushOpts(id = "brush_zoom", resetOnNew = TRUE))),
                                          column(3,
                                                 plotOutput("plot2", click = "plot_click", dblclick = "plot2_dblclick", brush = brushOpts(id = "brush_zoom", resetOnNew = TRUE))),
                                          column(3, 
                                                 plotOutput("plot3", click = "plot_click", dblclick = "plot3_dblclick", brush = brushOpts(id = "brush_zoom", resetOnNew = TRUE))),
                                          column(3,
                                                 plotOutput("plot4", click = "plot_click", dblclick = "plot4_dblclick", brush = brushOpts(id = "brush_zoom", resetOnNew = TRUE)))
                                 ))),
                 
                 hr(),
                 fluidRow(column(3, h2("2. Click into the plots to define the section lines"),
                                 helpText("Define section line 1. This should be the longes possible direction through the soma."),
                                 actionButton("place_nucleus1", "Place marker for nucleus 1"),
                                 actionButton("place_section1", "Place end point for section line 1"),
                                 hr(),
                                 helpText("Define section line 2. This should be about 45 degree clockwise from first section line."),
                                 actionButton("place_nucleus2", "Place marker for nucleus 2"),
                                 actionButton("place_section2", "Place end point for section line 2"),
                                 helpText("Define section line 3- This should be about 45 degrees counter-clockwise from frist section line."),
                                 actionButton("place_nucleus3", "Place marker for nucleus 3"),
                                 actionButton("place_section3", "Place end point for section line 3")),
                          column(6, fluidRow(column(3,radioButtons("which_big_plot", "Pick channel to enlarge:", c("Red" = 1, "Greed" = 2, "Blue" = 3), inline = T)),
                                             column(3, sliderInput("big_brightness", "Brightness (only affects display, not analyzed values):", min = -1, max = 1, value = 0, step = 0.1)),
                                             column(3, sliderInput("big_contrast", "Contrast (only affects display, not analyzed values):", min = 0, max = 5, value = 1, step = 0.5))),
                                 plotOutput("big_plot", click = "plot_click", dblclick = "plot1_dblclick",  width = "100%", height = 800))),
                 
                 hr(),
                 fluidRow(column(3,   h2("3. Analyze"),
                                 selectInput("dropdown", "Select condition", "no file loaded yet"),
                                 fluidRow(column(8, textInput("new_condition", "Or enter new condition")),
                                          column(4, actionButton("add_condition", "Add"))),
                                 tags$style(type='text/css', "#add_condition { width:100%; margin-top: 25px;}"),
                                 actionButton("save_analysis", "Save analysis"),
                                 textOutput("save_notification")),
                          column(6,
                                 plotOutput("mito_section1"),
                                 plotOutput("mito_section2"),
                                 plotOutput("mito_section3")))
                 
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  
  library(EBImage)
  library(tidyverse)
  library(ggplot2)
  
  all_points_chosen <- reactiveValues(ok = F)
  chosen_points <- reactiveValues(n1 = F, s1 = F, n2 = F, s2 = F, n3 = F, s3 = F)
  ranges <- reactiveValues(x = c(1, 512), y = c(1,512))
  big_image <- reactiveValues(image = NULL)
  save_notification_text <<-""
  
  output$status_text <- renderText({
    paste("Chosen file:", input$file1$name)
  })
  
  get_line_coordinates <- function(x_nuc, y_nuc, x_s, y_s) {
    # m = delta y / delta x ;  t = y_nuc - m*x_nuc ;  y = mx+t
    if (x_nuc == x_s) {x_s <<- x_s + 1} # correct end point one pixel to the right if at the same x position as nucleus
    m <- (y_s - y_nuc) / (x_s - x_nuc)
    t <- y_nuc - m*x_nuc
    fixed_coordinates <- seq(from = min(x_nuc, x_s), to = max(x_nuc, x_s), by = 1)
    fixed_coordinates <- rbind(fixed_coordinates, round((m*fixed_coordinates)+t))
    # this is regarding the line as a linear function, which can only have one y coordinate per x coordinate.
    # But for a steep slope, we need several. Extreme example: line from 0|0 to 1|100 - the current algorithm would have only two pixels colored
    # So, I will fill the holes in the line: Go from the first to the second pixel. If they are not adjacent, I add pixels with the same x as the
    # first, but with y moving towards the second pixel, until they are adjacent. Then go from second to third, ...
    filled_coordinates <- fixed_coordinates
    for (i in 1:(ncol(fixed_coordinates)-1)){
      y_dist <- abs(fixed_coordinates[2,i]-fixed_coordinates[2,(i+1)])
      if (y_dist > 1) {
        added_coordinates <- matrix(fixed_coordinates[1,i], 1, y_dist-1)
        added_coordinates <- rbind(added_coordinates, seq(from =(min(fixed_coordinates[2,i], fixed_coordinates[2,(i+1)])+1), to = (max(fixed_coordinates[2,i], fixed_coordinates[2,(i+1)])-1), by = 1 ))
        filled_coordinates <- cbind(filled_coordinates, added_coordinates)
      }
    }
    return(filled_coordinates)
  }
  
  observeEvent(input$file1, {
    #initialize global variables
    save_click_as <<- "waiting"
    all_points_chosen <- reactiveValues(ok = F)
    chosen_points <- reactiveValues(n1 = F, s1 = F, n2 = F, s2 = F, n3 = F, s3 = F)
    ranges <- reactiveValues(x = c(1, 512), y = c(1,512))
    big_image <- reactiveValues(image = NULL)
    save_notification_text <<-""
    # load conditions list
    old_conditions <<- NULL
    if (file.exists("condition_choices.rds")){
      old_conditions <<- readRDS(file = "condition_choices.rds")}
    updateSelectInput(getDefaultReactiveDomain(), "dropdown",
                      label = "Select condition",
                      choices = old_conditions)
    #load data
    rawimage <<- readImage(input$file1$datapath)
  })
  # Add entered condition as choice upon button press
  observeEvent(input$add_condition, {
    new_condition <- input$new_condition
    condition_choices <<- append(old_conditions, new_condition)
    saveRDS(condition_choices, file = "condition_choices.rds")
    updateSelectInput(getDefaultReactiveDomain(), "dropdown",
                      label = "Select condition",
                      choices = condition_choices,
                      selected = new_condition)
  })
  
  output$plot1 <- renderPlot({
    input$file1
    # split image into three color channels
    rgbr <- rgbImage(red = channel(rawimage,"r"))
    input$plot_click
    if (chosen_points$n1 == T) {
      if (chosen_points$s1 == T) {
        line1_coordinates <<- get_line_coordinates(x_n1, y_n1, x_s1, y_s1)
        for (i in 1:ncol(line1_coordinates)) {
          rgbr[(line1_coordinates[1,i]-1):(line1_coordinates[1,i]+1), (line1_coordinates[2,i]-1):(line1_coordinates[2,i]+1),] <- c(255,255,255)
        }
      }
    }
    if (chosen_points$n2 == T) {
      if (chosen_points$s2 == T) {
        line2_coordinates <<- get_line_coordinates(x_n2, y_n2, x_s2, y_s2)
        for (i in 1:ncol(line2_coordinates)) {
          rgbr[(line2_coordinates[1,i]-1):(line2_coordinates[1,i]+1), (line2_coordinates[2,i]-1):(line2_coordinates[2,i]+1),] <- c(255,255,255)
        }
      }
    }
    if (chosen_points$n3 == T) {
      if (chosen_points$s3 == T) {
        line3_coordinates <<- get_line_coordinates(x_n3, y_n3, x_s3, y_s3)
        for (i in 1:ncol(line3_coordinates)) {
          rgbr[(line3_coordinates[1,i]-1):(line3_coordinates[1,i]+1), (line3_coordinates[2,i]-1):(line3_coordinates[2,i]+1),] <- c(255,255,255)
        }
      }
    }
    plot(rgbr[ranges$x[1]:ranges$x[2], ranges$y[1]:ranges$y[2],])
  }, bg = "transparent")
  
  output$plot2 <- renderPlot({
    input$file1
    rgbg <- rgbImage(green = channel(rawimage,"g"))
    input$plot_click
    if (chosen_points$n1 == T) {
      if (chosen_points$s1 == T) {
        line1_coordinates <<- get_line_coordinates(x_n1, y_n1, x_s1, y_s1)
        for (i in 1:ncol(line1_coordinates)) {
          rgbg[(line1_coordinates[1,i]-1):(line1_coordinates[1,i]+1), (line1_coordinates[2,i]-1):(line1_coordinates[2,i]+1),] <- c(255,255,255)
        }
      }
    }
    if (chosen_points$n2 == T) {
      if (chosen_points$s2 == T) {
        line2_coordinates <<- get_line_coordinates(x_n2, y_n2, x_s2, y_s2)
        for (i in 1:ncol(line2_coordinates)) {
          rgbg[(line2_coordinates[1,i]-1):(line2_coordinates[1,i]+1), (line2_coordinates[2,i]-1):(line2_coordinates[2,i]+1),] <- c(255,255,255)
        }
      }
    }
    if (chosen_points$n3 == T) {
      if (chosen_points$s3 == T) {
        line3_coordinates <<- get_line_coordinates(x_n3, y_n3, x_s3, y_s3)
        for (i in 1:ncol(line3_coordinates)) {
          rgbg[(line3_coordinates[1,i]-1):(line3_coordinates[1,i]+1), (line3_coordinates[2,i]-1):(line3_coordinates[2,i]+1),] <- c(255,255,255)
        }
      }
    }
    plot(rgbg[ranges$x[1]:ranges$x[2], ranges$y[1]:ranges$y[2],])
  }, bg = "transparent")
  
  output$plot3 <- renderPlot({
    input$file1
    rgbb <- rgbImage(blue = channel(rawimage,"b"))
    input$plot_click
    if (chosen_points$n1 == T) {
      if (chosen_points$s1 == T) {
        line1_coordinates <<- get_line_coordinates(x_n1, y_n1, x_s1, y_s1)
        for (i in 1:ncol(line1_coordinates)) {
          rgbb[(line1_coordinates[1,i]-1):(line1_coordinates[1,i]+1), (line1_coordinates[2,i]-1):(line1_coordinates[2,i]+1),] <- c(255,255,255)
        }
      }
    }
    if (chosen_points$n2 == T) {
      if (chosen_points$s2 == T) {
        line2_coordinates <<- get_line_coordinates(x_n2, y_n2, x_s2, y_s2)
        for (i in 1:ncol(line2_coordinates)) {
          rgbb[(line2_coordinates[1,i]-1):(line2_coordinates[1,i]+1), (line2_coordinates[2,i]-1):(line2_coordinates[2,i]+1),] <- c(255,255,255)
        }
      }
    }
    if (chosen_points$n3 == T) {
      if (chosen_points$s3 == T) {
        line3_coordinates <<- get_line_coordinates(x_n3, y_n3, x_s3, y_s3)
        for (i in 1:ncol(line3_coordinates)) {
          rgbb[(line3_coordinates[1,i]-1):(line3_coordinates[1,i]+1), (line3_coordinates[2,i]-1):(line3_coordinates[2,i]+1),] <- c(255,255,255)
        }
      }
    }
    plot(rgbb[ranges$x[1]:ranges$x[2], ranges$y[1]:ranges$y[2],])
  }, bg = "transparent")
  
  output$plot4 <- renderPlot({
    input$file1
    input$plot_click
    rimage <- rawimage
    if (chosen_points$n1 == T) {
      if (chosen_points$s1 == T) {
        line1_coordinates <<- get_line_coordinates(x_n1, y_n1, x_s1, y_s1)
        for (i in 1:ncol(line1_coordinates)) {
          rimage[(line1_coordinates[1,i]-1):(line1_coordinates[1,i]+1), (line1_coordinates[2,i]-1):(line1_coordinates[2,i]+1),] <- c(255,255,255)
        }
      }
    }
    if (chosen_points$n2 == T) {
      if (chosen_points$s2 == T) {
        line2_coordinates <<- get_line_coordinates(x_n2, y_n2, x_s2, y_s2)
        for (i in 1:ncol(line2_coordinates)) {
          rimage[(line2_coordinates[1,i]-1):(line2_coordinates[1,i]+1), (line2_coordinates[2,i]-1):(line2_coordinates[2,i]+1),] <- c(255,255,255)
        }
      }
    }
    if (chosen_points$n3 == T) {
      if (chosen_points$s3 == T) {
        line3_coordinates <<- get_line_coordinates(x_n3, y_n3, x_s3, y_s3)
        for (i in 1:ncol(line3_coordinates)) {
          rimage[(line3_coordinates[1,i]-1):(line3_coordinates[1,i]+1), (line3_coordinates[2,i]-1):(line3_coordinates[2,i]+1),] <- c(255,255,255)
        }
      }
    }
    plot(rimage[ranges$x[1]:ranges$x[2], ranges$y[1]:ranges$y[2],])
  }, bg = "transparent")
  
  output$big_plot <- renderPlot({
    input$file1
    input$plot_click
    chan <- input$which_big_plot
    if (chan == 1) {
      channel_image <- rgbImage(red = channel(rawimage,"r"))
    } else if (chan == 2) {
      channel_image <- rgbImage(green = channel(rawimage,"g"))
    } else {
      channel_image <- rgbImage(blue = channel(rawimage,"b"))
    }
    channel_image <- (channel_image + input$big_brightness) * input$big_contrast
    input$plot_click
    if (chosen_points$n1 == T) {
      if (chosen_points$s1 == T) {
        line1_coordinates <<- get_line_coordinates(x_n1, y_n1, x_s1, y_s1)
        for (i in 1:ncol(line1_coordinates)) {
          channel_image[(line1_coordinates[1,i]-1):(line1_coordinates[1,i]+1), (line1_coordinates[2,i]-1):(line1_coordinates[2,i]+1),] <- c(255,255,255)
        }
      }
    }
    if (chosen_points$n2 == T) {
      if (chosen_points$s2 == T) {
        line2_coordinates <<- get_line_coordinates(x_n2, y_n2, x_s2, y_s2)
        for (i in 1:ncol(line2_coordinates)) {
          channel_image[(line2_coordinates[1,i]-1):(line2_coordinates[1,i]+1), (line2_coordinates[2,i]-1):(line2_coordinates[2,i]+1),] <- c(255,255,255)
        }
      }
    }
    if (chosen_points$n3 == T) {
      if (chosen_points$s3 == T) {
        line3_coordinates <<- get_line_coordinates(x_n3, y_n3, x_s3, y_s3)
        for (i in 1:ncol(line3_coordinates)) {
          channel_image[(line3_coordinates[1,i]-1):(line3_coordinates[1,i]+1), (line3_coordinates[2,i]-1):(line3_coordinates[2,i]+1),] <- c(255,255,255)
        }
      }
    }
    par(bg = "#272B30")
    plot(channel_image[ranges$x[1]:ranges$x[2], ranges$y[1]:ranges$y[2],])
  }, width = 800, height = 800)
  
  ## zoom
  observeEvent(input$brush_zoom, {
    brush <- input$brush_zoom
    if (!is.null(brush)) {
      center_x <- round(brush$xmin + ((brush$xmax - brush$xmin)/2))
      center_y <- round(brush$ymin + ((brush$ymax - brush$ymin)/2))
      radius <- round(max(((brush$xmax - brush$xmin)/2), ((brush$ymax - brush$ymin)/2)))
      ranges$x <- c(max(1,center_x - radius), min(512, center_x + radius))
      ranges$y <- c(max(1,center_y - radius), min(512, center_y + radius))
      
    } else {
      ranges$x <- c(1, 512)
      ranges$y <- c(1, 512)
    }
    runjs("document.getElementById('brush_zoom').remove()")
  })
  observeEvent(input$reset_zoom, {
    ranges$x <- c(1, 512)
    ranges$y <- c(1, 512)
  })
  
  
  observeEvent(input$place_nucleus1, {save_click_as <<- "nucleus1"})
  observeEvent(input$place_section1, {save_click_as <<- "section1"})
  observeEvent(input$place_nucleus2, {save_click_as <<- "nucleus2"})
  observeEvent(input$place_section2, {save_click_as <<- "section2"})
  observeEvent(input$place_nucleus3, {save_click_as <<- "nucleus3"})
  observeEvent(input$place_section3, {save_click_as <<- "section3"})
  
  observeEvent(input$plot_click, {
    x_coor <- round(input$plot_click$x + ranges$x[1])
    y_coor <- round(input$plot_click$y + ranges$y[1])
    if (save_click_as == "nucleus1") {
      x_n1 <<- x_coor
      y_n1 <<- y_coor
      chosen_points$n1 <- T
      save_click_as <<- "waiting"
    } else if (save_click_as == "section1") {
      x_s1 <<- x_coor
      y_s1 <<- y_coor
      chosen_points$s1 <- T
      save_click_as <<- "waiting"
    } else if (save_click_as == "nucleus2") {
      x_n2 <<- x_coor
      y_n2 <<- y_coor
      chosen_points$n2 <- T
      save_click_as <<- "waiting"
    } else if (save_click_as == "section2") {
      x_s2 <<- x_coor
      y_s2 <<- y_coor
      chosen_points$s2 <- T
      save_click_as <<- "waiting"
    } else if (save_click_as == "nucleus3") {
      x_n3 <<- x_coor
      y_n3 <<- y_coor
      chosen_points$n3 <- T
      save_click_as <<- "waiting"
    } else if (save_click_as == "section3") {
      x_s3 <<- x_coor
      y_s3 <<- y_coor
      chosen_points$s3 <- T
      save_click_as <<- "waiting"
    }
  })
  
  plot_section_line <- function(line_coordinates, x_n, y_n) {
    line <- data.frame(t(line_coordinates))
    colnames(line) <- c("x", "y")
    line <- line %>% mutate(distance = sqrt(((x-x_n1)^2) + ((y-y_n1)^2))) %>% arrange(distance) 
    distance_min <- min(line$distance)
    distance_max <- max(line$distance)
    line <- line %>% mutate(distance_norm = (distance-distance_min)/(distance_max - distance_min))
    mitovalues <- matrix(0, nrow(line), 1)
    for (i in 1:nrow(line)) {
      mitovalues[i,] <- imageData(rawimage)[line$x[i], line$y[i],1]
    }
    line <- line %>% mutate(mito_section_values = mitovalues)
    return(line)
    
  }
  
  # section plots
  output$mito_section1 <- renderPlot({
    input$plot_click
    if (chosen_points$n1 == T) {
      if (chosen_points$s1 == T) {
        line1 <<- plot_section_line(line1_coordinates, x_n1, y_n1)
        line1 %>% ggplot(aes(x = distance_norm, y = mito_section_values)) + geom_point() + theme_bw(base_size = 16) + labs(title = "Red (mito) intensity values along section line 1", x = "Distance from set nucleus point", y = "intensity")
      }
    }
  })
  
  output$mito_section2 <- renderPlot({
    input$plot_click
    if (chosen_points$n2 == T) {
      if (chosen_points$s2 == T) {
        line2 <<- plot_section_line(line2_coordinates, x_n2, y_n2)
        line2 %>% ggplot(aes(x = distance_norm, y = mito_section_values)) + geom_point() + theme_bw(base_size = 16) + labs(title = "Red (mito) intensity values along section line 2", x = "Distance from set nucleus point", y = "intensity")
      }
    }
  })
  
  output$mito_section3 <- renderPlot({
    input$plot_click
    if (chosen_points$n3 == T) {
      if (chosen_points$s3 == T) {
        line3 <<- plot_section_line(line3_coordinates, x_n3, y_n3)
        line3 %>% ggplot(aes(x = distance_norm, y = mito_section_values)) + geom_point() + theme_bw(base_size = 16) + labs(title = "Red (mito) intensity values along section line 3", x = "Distance from set nucleus point", y = "intensity")
      }
    }
  })
  
  observeEvent(input$save_analysis, {
    setwd(getSrcDirectory(function(){})[1])
    rec_name <- substr(input$file1$name, 1, (nchar(input$file1$name)-4))
    # load conditions lookup-table ( if exists), add this file and save
    condition_lookup <<- NULL
    if (file.exists("../../Results/Mitosections/condition_lookup.rds")){
      condition_lookup <<- readRDS(file = "../../Results/Mitosections/condition_lookup.rds")
      if (any(str_detect(condition_lookup[,2], rec_name))) {
        # if the current filename is already in the list
        condition_lookup[str_which(condition_lookup[,2], rec_name),] <- c(input$dropdown, rec_name)
      } else {
        condition_lookup <<- rbind(condition_lookup, c(input$dropdown, rec_name))
      }
    } else {
      condition_lookup <<- rbind(condition_lookup, c(input$dropdown, rec_name))  
    }
    saveRDS(condition_lookup, file = "../../Results/Mitosections/condition_lookup.rds")
    # generate thumbnail for saving
    x_mean <- mean(c(x_n1, x_n2, x_n3))
    y_mean <- mean(c(y_n1, y_n2, y_n3))
    max_dist_x <- max(abs(c((x_s1-x_mean), (x_s2-x_mean), (x_s3-x_mean))))
    max_dist_y <- max(abs(c((y_s1-y_mean), (y_s2-y_mean), (y_s3-y_mean))))
    max_dist <- max(c(max_dist_x, max_dist_y)) *1.5
    thumbnail <- rawimage[max(c(1,x_mean-max_dist)):min(c(512, x_mean+max_dist)), max(c(0, y_mean-max_dist)):min(c(512, y_mean+max_dist)),]
    # gather results for saving
    save_filename <- paste("../../Results/Mitosections/", rec_name, ".rds", sep="")
    results <- list()
    results$name <- rec_name
    results$condition <- input$dropdown
    results$images$rawimage <- rawimage
    results$images$thumbnail <- thumbnail
    results$sections$s1 <- line1
    results$sections$s2 <- line2
    results$sections$s3 <- line3
    results$timestamp <- c(format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
    # save
    saveRDS(results, file = save_filename)
    save_notification_text <<- "File saved. You can close this window or upload the next image to analyze."
  })
  
  output$save_notification <- renderText({
    input$save_analysis
    paste(save_notification_text)
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
