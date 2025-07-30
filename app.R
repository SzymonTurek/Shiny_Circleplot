library(shiny)
library(DT)  
library(dplyr)
library(circlize)



# App 1 UI and Server
ui1 <- fluidPage(
  titlePanel("Convert positions on contigs to positions on chromosomes"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("file1", "Choose CSV File",
                multiple = FALSE,
                accept = c(".csv")),
      checkboxInput("header", "Header", TRUE),
      radioButtons("sep", "Separator",
                   choices = c(Comma = ",",
                               Semicolon = ";",
                               Tab = "\t"),
                   selected = ","),
      
      radioButtons("disp", "Display",
                   choices = c(Head = "head",
                               All = "all"),
                   selected = "head"),
      hr(),
      downloadButton("downloadData", "Download converted CSV")
    ),
    
    mainPanel(
      DTOutput("contents")
    )
  )
)



server1 <- function(input, output, session) {
  
  df_reference <- read.csv("Contigstochromosomestable.csv", sep=";")
  df_reference <- df_reference[,c("Chr",  "Contig_name" , "Contig_length", "Contig_start", "Contig_end" , "ctg_start_sum" , "ctg_end_sum" )]
  output_df <- data.frame(  Searched_conting = character(), Searched_start_position =  numeric(), Searched_end_position = numeric(), Found_chromosome = character(), Found_start_position_on_chromosome = numeric(), Found_end_position_on_chromosome = numeric()
                            , stringsAsFactors = FALSE)
  
  
  data <- reactive({
    req(input$file1)
    
    tryCatch({
      df <- read.csv(input$file1$datapath,
                     header = input$header,
                     sep = input$sep,
                     stringsAsFactors = FALSE)
      df
    }, error = function(e) {
      showNotification("Error reading file. Please check the file format.", type = "error")
      NULL
    })
    
    
    
  })
  
  output$contents <- renderDT({
    df <- data()
    req(df)
    
    if (input$disp == "head") {
      datatable(head(df))
    } else {
      datatable(df)
    }
    
    
  })
  
  
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste0("Chromosome_positions_", Sys.Date(), ".csv")
    },
    content = function(file) {
      df <- data()
      req(df)
      
      
      for (i in 1:nrow(df)){
        contig <- df[i, "contig"]
        start_position <- df[i, "start_position"]
        end_position <- df[i, "end_position"]
        
        
        filtered_df <- df_reference[df_reference$Contig_name == contig, ]
        Chromosome <- filtered_df$Chr
        Chromosome <- paste0("CHR", Chromosome)
        position_start <-  filtered_df$ctg_start_sum + start_position
        
        position_end <- position_start + (end_position - start_position )
        new_row <- data.frame(Searched_conting = contig, Searched_start_position =  start_position, Searched_end_position = end_position, Found_chromosome = Chromosome, Found_start_position_on_chromosome = position_start, Found_end_position_on_chromosome = position_end)
        output_df <- rbind(output_df, new_row)
        
      }
      
      
      
      write.csv(output_df, file, row.names = FALSE)
      
    }
  )
}



# App 2 UI and Server

ui2 <- fluidPage(
  tags$style(type = "text/css", ".shiny-output-error { visibility: hidden; }"),
  tags$style(type = "text/css", ".shiny-output-error:before { visibility: hidden; }"),
  
  titlePanel("Create Circle Plot"),
  
  sidebarLayout(
    sidebarPanel(
      numericInput("num_files", "How many files do you want to upload?", value = 1, min = 1, max = 10),
      uiOutput("file_inputs"),

     
      hr(style = "border-top: 1px solid black;"),
      textInput("color1", "Color of first circle", value = "#B31B2A"),
      textInput("color2", "Color of second circle ", value = "#003228"),
      textInput("color3", "Color of third circle ", value = "#BEDCF0"),
      textInput("color4", "Color of fourth circle ", value = "#F01E32"),
      textInput("color5", "Color of fifth circle ", value = "#d66a27"),
      textInput("color6", "Color of sixth circle ", value = "#CD661D"),
      textInput("color7", "Color of seventh circle ", value = "#556B2F"),
      downloadButton("downloadPlot", "Download Plot")
      
      
    ),
    
    mainPanel(
      uiOutput("file_heads"),
      plotOutput("circosPlot", height = "900px")
    )
  )
)


server2 <- function(input, output, session) {
  
  output$file_inputs <- renderUI({
    req(input$num_files)
    n <- input$num_files
    input_list <- lapply(seq_len(n), function(i) {
      fileInput(inputId = paste0("file", i), label = paste("Upload File", i))
    })
    do.call(tagList, input_list)
  })
  
  output$file_heads <- renderUI({
    req(input$num_files)
    n <- input$num_files
    output_list <- lapply(seq_len(n), function(i) {
      output_name <- paste0("head", i)
      tableOutput(output_name)
    })
    do.call(tagList, output_list)
  })
  
  df_cuc = data.frame(
    name  = c("CHR1",  "CHR2", "CHR3",  "CHR4",  "CHR5", "CHR6",  "CHR7"),
    start = c(0, 0, 0, 0, 0, 0, 0),
    end   = c(32066787, 25497826, 41162263, 33084761, 40500461, 32941524, 22885705))
  
  

  
  observe({
    req(input$num_files)
    n <- input$num_files
    
    for (i in seq_len(n)) {
      local({
        j <- i  # local copy for correct closure
        output_name <- paste0("head", j)
        file_input_id <- paste0("file", j)
        obj_name <- paste0("csv_file_", j)
        file <- input[[file_input_id]]
        if (is.null(file)) return(NULL)
        

        df <- tryCatch({
          read.csv(file$datapath)#, nrows = 5)
        })
        df <- df[!grepl("ctg", df$Found_chromosome), ]
        # print(df)
        
        assign(obj_name, df, envir = .GlobalEnv)
        
        
      })
    }
    
    output$circosPlot <- renderPlot({
      
      if (input$num_files == 1) {
        req(get("csv_file_1", envir = .GlobalEnv))
        csv_f1 <- get("csv_file_1", envir = .GlobalEnv)
        
        circos.clear()
        circos.par("start.degree" = 83.5, "gap.degree" = c(2,2,2,2,2,2,12), "track.height" = 0.1, "track.margin" = c(0, 0.01))
        circos.genomicInitialize(df_cuc, sector.names = df_cuc$name, factors = df_cuc$name, axis.labels.cex = 0.4, plotType = "axis")
        
        
        circos.trackHist( factors = csv_f1[["Found_chromosome"]], x=csv_f1[["Found_start_position_on_chromosome"]], bin.size = 200000, col = input$color1, border = input$color1)
        circos.yaxis(side = "right", sector.index =CELL_META$sector.index,
                     track.index = get.cell.meta.data("track.index"), labels.cex = 0.4)
        
        circos.track(ylim = c(0, 1), 
                     bg.col = "#edece8", 
                     bg.border = "black", track.height = 0.05, panel.fun = function(x, y) {
                       circos.text(CELL_META$xcenter, CELL_META$ycenter, CELL_META$sector.index, facing = "inside", niceFacing = TRUE, cex = 0.5)})
        
        
        
        output$downloadPlot <- downloadHandler(
          filename = function() {
            paste("plot-", Sys.Date(), ".pdf", sep = "")
          },
          content = function(file) {
            pdf(file)
            circos.clear()
            circos.par("start.degree" = 83.5, "gap.degree" = c(2,2,2,2,2,2,12), "track.height" = 0.1, "track.margin" = c(0, 0.01))
            circos.genomicInitialize(df_cuc, sector.names = df_cuc$name, factors = df_cuc$name, axis.labels.cex = 0.4, plotType = "axis")
            
            
            circos.trackHist( factors = csv_f1[["Found_chromosome"]], x=csv_f1[["Found_start_position_on_chromosome"]], bin.size = 200000, col = input$color1, border = input$color1)
            circos.yaxis(side = "right", sector.index =CELL_META$sector.index,
                         track.index = get.cell.meta.data("track.index"), labels.cex = 0.4)
            
            circos.track(ylim = c(0, 1), 
                         bg.col = "#edece8", 
                         bg.border = "black", track.height = 0.05, panel.fun = function(x, y) {
                           circos.text(CELL_META$xcenter, CELL_META$ycenter, CELL_META$sector.index, facing = "inside", niceFacing = TRUE, cex = 0.5)})
            
            dev.off()
          }
        )
        
        
      }
      
      if (input$num_files == 2) {
        req(get("csv_file_2", envir = .GlobalEnv))
        req(get("csv_file_1", envir = .GlobalEnv))
        
        csv_f1 <- get("csv_file_1", envir = .GlobalEnv)
        csv_f2 <- get("csv_file_2", envir = .GlobalEnv)
        
        circos.clear()
        circos.par("start.degree" = 83.5, "gap.degree" = c(2,2,2,2,2,2,12), "track.height" = 0.1, "track.margin" = c(0, 0.01))
        circos.genomicInitialize(df_cuc, sector.names = df_cuc$name, factors = df_cuc$name, axis.labels.cex = 0.4, plotType = "axis")
        
        
        circos.trackHist( factors = csv_f1[["Found_chromosome"]], x=csv_f1[["Found_start_position_on_chromosome"]], bin.size = 200000, col = input$color1, border = input$color1)
        circos.yaxis(side = "right", sector.index =CELL_META$sector.index,
                     track.index = get.cell.meta.data("track.index"), labels.cex = 0.4)
        circos.trackHist( factors =csv_f2[["Found_chromosome"]], x=csv_f2[["Found_start_position_on_chromosome"]], bin.size = 200000, col = input$color2, border = input$color2)
        circos.yaxis(side = "right", sector.index =CELL_META$sector.index,
                     track.index = get.cell.meta.data("track.index"), labels.cex = 0.4)
        
        
        
        circos.track(ylim = c(0, 1), 
                     bg.col = "#edece8", 
                     bg.border = "black", track.height = 0.05, panel.fun = function(x, y) {
                       circos.text(CELL_META$xcenter, CELL_META$ycenter, CELL_META$sector.index, facing = "inside", niceFacing = TRUE, cex = 0.5)})
        
        
        
        
        output$downloadPlot <- downloadHandler(
          filename = function() {
            paste("plot-", Sys.Date(), ".pdf", sep = "")
          },
          content = function(file) {
            pdf(file)
            circos.clear()
            circos.par("start.degree" = 83.5, "gap.degree" = c(2,2,2,2,2,2,12), "track.height" = 0.1, "track.margin" = c(0, 0.01))
            circos.genomicInitialize(df_cuc, sector.names = df_cuc$name, factors = df_cuc$name, axis.labels.cex = 0.4, plotType = "axis")
            
            
            circos.trackHist( factors = csv_f1[["Found_chromosome"]], x=csv_f1[["Found_start_position_on_chromosome"]], bin.size = 200000, col = input$color1, border = input$color1)
            circos.yaxis(side = "right", sector.index =CELL_META$sector.index,
                         track.index = get.cell.meta.data("track.index"), labels.cex = 0.4)
            circos.trackHist( factors =csv_f2[["Found_chromosome"]], x=csv_f2[["Found_start_position_on_chromosome"]], bin.size = 200000, col = input$color2, border = input$color2)
            circos.yaxis(side = "right", sector.index =CELL_META$sector.index,
                         track.index = get.cell.meta.data("track.index"), labels.cex = 0.4)
            
            circos.track(ylim = c(0, 1), 
                         bg.col = "#edece8", 
                         bg.border = "black", track.height = 0.05, panel.fun = function(x, y) {
                           circos.text(CELL_META$xcenter, CELL_META$ycenter, CELL_META$sector.index, facing = "inside", niceFacing = TRUE, cex = 0.5)})
            
            dev.off()
          }
        )
        
        
      }
      
      if (input$num_files == 3) {
        req(get("csv_file_2", envir = .GlobalEnv))
        req(get("csv_file_1", envir = .GlobalEnv))
        req(get("csv_file_3", envir = .GlobalEnv))
        
        csv_f1 <- get("csv_file_1", envir = .GlobalEnv)
        csv_f2 <- get("csv_file_2", envir = .GlobalEnv)
        csv_f3 <- get("csv_file_3", envir = .GlobalEnv)
        
        circos.clear()
        circos.par("start.degree" = 83.5, "gap.degree" = c(2,2,2,2,2,2,12), "track.height" = 0.1, "track.margin" = c(0, 0.01))
        circos.genomicInitialize(df_cuc, sector.names = df_cuc$name, factors = df_cuc$name, axis.labels.cex = 0.4, plotType = "axis")
        
        
        circos.trackHist( factors = csv_f1[["Found_chromosome"]], x=csv_f1[["Found_start_position_on_chromosome"]], bin.size = 200000, col = input$color1, border = input$color1)
        circos.yaxis(side = "right", sector.index =CELL_META$sector.index,
                     track.index = get.cell.meta.data("track.index"), labels.cex = 0.4)
        circos.trackHist( factors =csv_f2[["Found_chromosome"]], x=csv_f2[["Found_start_position_on_chromosome"]], bin.size = 200000, col = input$color2, border = input$color2)
        circos.yaxis(side = "right", sector.index =CELL_META$sector.index,
                     track.index = get.cell.meta.data("track.index"), labels.cex = 0.4)
        circos.trackHist( factors =csv_f3[["Found_chromosome"]], x=csv_f3[["Found_start_position_on_chromosome"]], bin.size = 200000, col = input$color3, border = input$color3)
        circos.yaxis(side = "right", sector.index =CELL_META$sector.index,
                     track.index = get.cell.meta.data("track.index"), labels.cex = 0.4)
        circos.track(ylim = c(0, 1), 
                     bg.col = "#edece8", 
                     bg.border = "black", track.height = 0.05, panel.fun = function(x, y) {
                       circos.text(CELL_META$xcenter, CELL_META$ycenter, CELL_META$sector.index, facing = "inside", niceFacing = TRUE, cex = 0.5)})
        
        
        
        
        output$downloadPlot <- downloadHandler(
          filename = function() {
            paste("plot-", Sys.Date(), ".pdf", sep = "")
          },
          content = function(file) {
            pdf(file)
            circos.clear()
            circos.par("start.degree" = 83.5, "gap.degree" = c(2,2,2,2,2,2,12), "track.height" = 0.1, "track.margin" = c(0, 0.01))
            circos.genomicInitialize(df_cuc, sector.names = df_cuc$name, factors = df_cuc$name, axis.labels.cex = 0.4, plotType = "axis")
            
            
            circos.trackHist( factors = csv_f1[["Found_chromosome"]], x=csv_f1[["Found_start_position_on_chromosome"]], bin.size = 200000, col = input$color1, border = input$color1)
            circos.yaxis(side = "right", sector.index =CELL_META$sector.index,
                         track.index = get.cell.meta.data("track.index"), labels.cex = 0.4)
            circos.trackHist( factors =csv_f2[["Found_chromosome"]], x=csv_f2[["Found_start_position_on_chromosome"]], bin.size = 200000, col = input$color2, border = input$color2)
            circos.yaxis(side = "right", sector.index =CELL_META$sector.index,
                         track.index = get.cell.meta.data("track.index"), labels.cex = 0.4)
            circos.trackHist( factors =csv_f3[["Found_chromosome"]], x=csv_f3[["Found_start_position_on_chromosome"]], bin.size = 200000, col = input$color3, border = input$color3)
            circos.yaxis(side = "right", sector.index =CELL_META$sector.index,
                         track.index = get.cell.meta.data("track.index"), labels.cex = 0.4)
            
            circos.track(ylim = c(0, 1), 
                         bg.col = "#edece8", 
                         bg.border = "black", track.height = 0.05, panel.fun = function(x, y) {
                           circos.text(CELL_META$xcenter, CELL_META$ycenter, CELL_META$sector.index, facing = "inside", niceFacing = TRUE, cex = 0.5)})
            
            dev.off()
          }
        )
        
      }
      if (input$num_files == 4) {
        req(get("csv_file_2", envir = .GlobalEnv))
        req(get("csv_file_1", envir = .GlobalEnv))
        req(get("csv_file_3", envir = .GlobalEnv))
        req(get("csv_file_4", envir = .GlobalEnv))
        
        csv_f1 <- get("csv_file_1", envir = .GlobalEnv)
        csv_f2 <- get("csv_file_2", envir = .GlobalEnv)
        csv_f3 <- get("csv_file_3", envir = .GlobalEnv)
        csv_f4 <- get("csv_file_4", envir = .GlobalEnv)
        
        
        circos.clear()
        circos.par("start.degree" = 83.5, "gap.degree" = c(2,2,2,2,2,2,12), "track.height" = 0.1, "track.margin" = c(0, 0.01))
        circos.genomicInitialize(df_cuc, sector.names = df_cuc$name, factors = df_cuc$name, axis.labels.cex = 0.4, plotType = "axis")
        
        
        circos.trackHist( factors = csv_f1[["Found_chromosome"]], x=csv_f1[["Found_start_position_on_chromosome"]], bin.size = 200000, col = input$color1, border = input$color1)
        circos.yaxis(side = "right", sector.index =CELL_META$sector.index,
                     track.index = get.cell.meta.data("track.index"), labels.cex = 0.4)
        circos.trackHist( factors =csv_f2[["Found_chromosome"]], x=csv_f2[["Found_start_position_on_chromosome"]], bin.size = 200000, col = input$color2, border = input$color2)
        circos.yaxis(side = "right", sector.index =CELL_META$sector.index,
                     track.index = get.cell.meta.data("track.index"), labels.cex = 0.4)
        circos.trackHist( factors =csv_f3[["Found_chromosome"]], x=csv_f3[["Found_start_position_on_chromosome"]], bin.size = 200000, col = input$color3, border = input$color3)
        circos.yaxis(side = "right", sector.index =CELL_META$sector.index,
                     track.index = get.cell.meta.data("track.index"), labels.cex = 0.4)
        circos.trackHist( factors =csv_f4[["Found_chromosome"]], x=csv_f4[["Found_start_position_on_chromosome"]], bin.size = 200000, col = input$color4, border = input$color4)
        circos.yaxis(side = "right", sector.index =CELL_META$sector.index,
                     track.index = get.cell.meta.data("track.index"), labels.cex = 0.4)
        
        
        circos.track(ylim = c(0, 1), 
                     bg.col = "#edece8", 
                     bg.border = "black", track.height = 0.05, panel.fun = function(x, y) {
                       circos.text(CELL_META$xcenter, CELL_META$ycenter, CELL_META$sector.index, facing = "inside", niceFacing = TRUE, cex = 0.5)})
        
        
        output$downloadPlot <- downloadHandler(
          filename = function() {
            paste("plot-", Sys.Date(), ".pdf", sep = "")
          },
          content = function(file) {
            pdf(file)
            circos.clear()
            circos.par("start.degree" = 83.5, "gap.degree" = c(2,2,2,2,2,2,12), "track.height" = 0.1, "track.margin" = c(0, 0.01))
            circos.genomicInitialize(df_cuc, sector.names = df_cuc$name, factors = df_cuc$name, axis.labels.cex = 0.4, plotType = "axis")
            
            
            circos.trackHist( factors = csv_f1[["Found_chromosome"]], x=csv_f1[["Found_start_position_on_chromosome"]], bin.size = 200000, col = input$color1, border = input$color1)
            circos.yaxis(side = "right", sector.index =CELL_META$sector.index,
                         track.index = get.cell.meta.data("track.index"), labels.cex = 0.4)
            circos.trackHist( factors =csv_f2[["Found_chromosome"]], x=csv_f2[["Found_start_position_on_chromosome"]], bin.size = 200000, col = input$color2, border = input$color2)
            circos.yaxis(side = "right", sector.index =CELL_META$sector.index,
                         track.index = get.cell.meta.data("track.index"), labels.cex = 0.4)
            circos.trackHist( factors =csv_f3[["Found_chromosome"]], x=csv_f3[["Found_start_position_on_chromosome"]], bin.size = 200000, col = input$color3, border = input$color3)
            circos.yaxis(side = "right", sector.index =CELL_META$sector.index,
                         track.index = get.cell.meta.data("track.index"), labels.cex = 0.4)
            circos.trackHist( factors =csv_f4[["Found_chromosome"]], x=csv_f4[["Found_start_position_on_chromosome"]], bin.size = 200000, col = input$color4, border = input$color4)
            circos.yaxis(side = "right", sector.index =CELL_META$sector.index,
                         track.index = get.cell.meta.data("track.index"), labels.cex = 0.4)
            
            
            circos.track(ylim = c(0, 1), 
                         bg.col = "#edece8", 
                         bg.border = "black", track.height = 0.05, panel.fun = function(x, y) {
                           circos.text(CELL_META$xcenter, CELL_META$ycenter, CELL_META$sector.index, facing = "inside", niceFacing = TRUE, cex = 0.5)})
            
            dev.off()
          }
        )
        
        
      }
      if (input$num_files == 5) {
        req(get("csv_file_2", envir = .GlobalEnv))
        req(get("csv_file_1", envir = .GlobalEnv))
        req(get("csv_file_3", envir = .GlobalEnv))
        req(get("csv_file_4", envir = .GlobalEnv))
        req(get("csv_file_5", envir = .GlobalEnv))
        
        csv_f1 <- get("csv_file_1", envir = .GlobalEnv)
        csv_f2 <- get("csv_file_2", envir = .GlobalEnv)
        csv_f3 <- get("csv_file_3", envir = .GlobalEnv)
        csv_f4 <- get("csv_file_4", envir = .GlobalEnv)
        csv_f5 <- get("csv_file_5", envir = .GlobalEnv)
        
        
        circos.clear()
        circos.par("start.degree" = 83.5, "gap.degree" = c(2,2,2,2,2,2,12), "track.height" = 0.1, "track.margin" = c(0, 0.01))
        circos.genomicInitialize(df_cuc, sector.names = df_cuc$name, factors = df_cuc$name, axis.labels.cex = 0.4, plotType = "axis")
        
        
        circos.trackHist( factors = csv_f1[["Found_chromosome"]], x=csv_f1[["Found_start_position_on_chromosome"]], bin.size = 200000, col = input$color1, border = input$color1)
        circos.yaxis(side = "right", sector.index =CELL_META$sector.index,
                     track.index = get.cell.meta.data("track.index"), labels.cex = 0.4)
        circos.trackHist( factors =csv_f2[["Found_chromosome"]], x=csv_f2[["Found_start_position_on_chromosome"]], bin.size = 200000, col = input$color2, border = input$color2)
        circos.yaxis(side = "right", sector.index =CELL_META$sector.index,
                     track.index = get.cell.meta.data("track.index"), labels.cex = 0.4)
        circos.trackHist( factors =csv_f3[["Found_chromosome"]], x=csv_f3[["Found_start_position_on_chromosome"]], bin.size = 200000, col = input$color3, border = input$color3)
        circos.yaxis(side = "right", sector.index =CELL_META$sector.index,
                     track.index = get.cell.meta.data("track.index"), labels.cex = 0.4)
        circos.trackHist( factors =csv_f4[["Found_chromosome"]], x=csv_f4[["Found_start_position_on_chromosome"]], bin.size = 200000, col = input$color4, border = input$color4)
        circos.yaxis(side = "right", sector.index =CELL_META$sector.index,
                     track.index = get.cell.meta.data("track.index"), labels.cex = 0.4)
        circos.trackHist( factors =csv_f5[["Found_chromosome"]], x=csv_f5[["Found_start_position_on_chromosome"]], bin.size = 200000, col = input$color5, border = input$color5)
        circos.yaxis(side = "right", sector.index =CELL_META$sector.index,
                     track.index = get.cell.meta.data("track.index"), labels.cex = 0.4)
        
        
        circos.track(ylim = c(0, 1), 
                     bg.col = "#edece8", 
                     bg.border = "black", track.height = 0.05, panel.fun = function(x, y) {
                       circos.text(CELL_META$xcenter, CELL_META$ycenter, CELL_META$sector.index, facing = "inside", niceFacing = TRUE, cex = 0.5)})
        
        
        output$downloadPlot <- downloadHandler(
          filename = function() {
            paste("plot-", Sys.Date(), ".pdf", sep = "")
          },
          content = function(file) {
            pdf(file)
            circos.clear()
            circos.par("start.degree" = 83.5, "gap.degree" = c(2,2,2,2,2,2,12), "track.height" = 0.1, "track.margin" = c(0, 0.01))
            circos.genomicInitialize(df_cuc, sector.names = df_cuc$name, factors = df_cuc$name, axis.labels.cex = 0.4, plotType = "axis")
            
            
            circos.trackHist( factors = csv_f1[["Found_chromosome"]], x=csv_f1[["Found_start_position_on_chromosome"]], bin.size = 200000, col = input$color1, border = input$color1)
            circos.yaxis(side = "right", sector.index =CELL_META$sector.index,
                         track.index = get.cell.meta.data("track.index"), labels.cex = 0.4)
            circos.trackHist( factors =csv_f2[["Found_chromosome"]], x=csv_f2[["Found_start_position_on_chromosome"]], bin.size = 200000, col = input$color2, border = input$color2)
            circos.yaxis(side = "right", sector.index =CELL_META$sector.index,
                         track.index = get.cell.meta.data("track.index"), labels.cex = 0.4)
            circos.trackHist( factors =csv_f3[["Found_chromosome"]], x=csv_f3[["Found_start_position_on_chromosome"]], bin.size = 200000, col = input$color3, border = input$color3)
            circos.yaxis(side = "right", sector.index =CELL_META$sector.index,
                         track.index = get.cell.meta.data("track.index"), labels.cex = 0.4)
            circos.trackHist( factors =csv_f4[["Found_chromosome"]], x=csv_f4[["Found_start_position_on_chromosome"]], bin.size = 200000, col = input$color4, border = input$color4)
            circos.yaxis(side = "right", sector.index =CELL_META$sector.index,
                         track.index = get.cell.meta.data("track.index"), labels.cex = 0.4)
            circos.trackHist( factors =csv_f5[["Found_chromosome"]], x=csv_f5[["Found_start_position_on_chromosome"]], bin.size = 200000, col = input$color5, border = input$color5)
            circos.yaxis(side = "right", sector.index =CELL_META$sector.index,
                         track.index = get.cell.meta.data("track.index"), labels.cex = 0.4)
            
            
            circos.track(ylim = c(0, 1), 
                         bg.col = "#edece8", 
                         bg.border = "black", track.height = 0.05, panel.fun = function(x, y) {
                           circos.text(CELL_META$xcenter, CELL_META$ycenter, CELL_META$sector.index, facing = "inside", niceFacing = TRUE, cex = 0.5)})
            
            dev.off()
          }
        )
        
        
      }
      
      if (input$num_files == 6) {
        req(get("csv_file_2", envir = .GlobalEnv))
        req(get("csv_file_1", envir = .GlobalEnv))
        req(get("csv_file_3", envir = .GlobalEnv))
        req(get("csv_file_4", envir = .GlobalEnv))
        req(get("csv_file_5", envir = .GlobalEnv))
        req(get("csv_file_6", envir = .GlobalEnv))

        csv_f1 <- get("csv_file_1", envir = .GlobalEnv)
        csv_f2 <- get("csv_file_2", envir = .GlobalEnv)
        csv_f3 <- get("csv_file_3", envir = .GlobalEnv)
        csv_f4 <- get("csv_file_4", envir = .GlobalEnv)
        csv_f5 <- get("csv_file_5", envir = .GlobalEnv)
        csv_f6 <- get("csv_file_6", envir = .GlobalEnv)


        circos.clear()
        circos.par("start.degree" = 83.5, "gap.degree" = c(2,2,2,2,2,2,12), "track.height" = 0.1, "track.margin" = c(0, 0.01))
        circos.genomicInitialize(df_cuc, sector.names = df_cuc$name, factors = df_cuc$name, axis.labels.cex = 0.4, plotType = "axis")


        circos.trackHist( factors = csv_f1[["Found_chromosome"]], x=csv_f1[["Found_start_position_on_chromosome"]], bin.size = 200000, col = input$color1, border = input$color1)
        circos.yaxis(side = "right", sector.index =CELL_META$sector.index,
                     track.index = get.cell.meta.data("track.index"), labels.cex = 0.4)
        circos.trackHist( factors =csv_f2[["Found_chromosome"]], x=csv_f2[["Found_start_position_on_chromosome"]], bin.size = 200000, col = input$color2, border = input$color2)
        circos.yaxis(side = "right", sector.index =CELL_META$sector.index,
                     track.index = get.cell.meta.data("track.index"), labels.cex = 0.4)
        circos.trackHist( factors =csv_f3[["Found_chromosome"]], x=csv_f3[["Found_start_position_on_chromosome"]], bin.size = 200000, col = input$color3, border = input$color3)
        circos.yaxis(side = "right", sector.index =CELL_META$sector.index,
                     track.index = get.cell.meta.data("track.index"), labels.cex = 0.4)
        circos.trackHist( factors =csv_f4[["Found_chromosome"]], x=csv_f4[["Found_start_position_on_chromosome"]], bin.size = 200000, col = input$color4, border = input$color4)
        circos.yaxis(side = "right", sector.index =CELL_META$sector.index,
                     track.index = get.cell.meta.data("track.index"), labels.cex = 0.4)
        circos.trackHist( factors =csv_f5[["Found_chromosome"]], x=csv_f5[["Found_start_position_on_chromosome"]], bin.size = 200000, col = input$color5, border = input$color5)
        circos.yaxis(side = "right", sector.index =CELL_META$sector.index,
                     track.index = get.cell.meta.data("track.index"), labels.cex = 0.4)
        circos.trackHist( factors =csv_f6[["Found_chromosome"]], x=csv_f6[["Found_start_position_on_chromosome"]], bin.size = 200000, col = input$color6, border = input$color6)
        circos.yaxis(side = "right", sector.index =CELL_META$sector.index,
                     track.index = get.cell.meta.data("track.index"), labels.cex = 0.4)


        circos.track(ylim = c(0, 1),
                     bg.col = "#edece8",
                     bg.border = "black", track.height = 0.05, panel.fun = function(x, y) {
                       circos.text(CELL_META$xcenter, CELL_META$ycenter, CELL_META$sector.index, facing = "inside", niceFacing = TRUE, cex = 0.5)})


        output$downloadPlot <- downloadHandler(
          filename = function() {
            paste("plot-", Sys.Date(), ".pdf", sep = "")
          },
          content = function(file) {
            pdf(file)
            circos.clear()
            circos.par("start.degree" = 83.5, "gap.degree" = c(2,2,2,2,2,2,12), "track.height" = 0.1, "track.margin" = c(0, 0.01))
            circos.genomicInitialize(df_cuc, sector.names = df_cuc$name, factors = df_cuc$name, axis.labels.cex = 0.4, plotType = "axis")


            circos.trackHist( factors = csv_f1[["Found_chromosome"]], x=csv_f1[["Found_start_position_on_chromosome"]], bin.size = 200000, col = input$color1, border = input$color1)
            circos.yaxis(side = "right", sector.index =CELL_META$sector.index,
                         track.index = get.cell.meta.data("track.index"), labels.cex = 0.4)
            circos.trackHist( factors =csv_f2[["Found_chromosome"]], x=csv_f2[["Found_start_position_on_chromosome"]], bin.size = 200000, col = input$color2, border = input$color2)
            circos.yaxis(side = "right", sector.index =CELL_META$sector.index,
                         track.index = get.cell.meta.data("track.index"), labels.cex = 0.4)
            circos.trackHist( factors =csv_f3[["Found_chromosome"]], x=csv_f3[["Found_start_position_on_chromosome"]], bin.size = 200000, col = input$color3, border = input$color3)
            circos.yaxis(side = "right", sector.index =CELL_META$sector.index,
                         track.index = get.cell.meta.data("track.index"), labels.cex = 0.4)
            circos.trackHist( factors =csv_f4[["Found_chromosome"]], x=csv_f4[["Found_start_position_on_chromosome"]], bin.size = 200000, col = input$color4, border = input$color4)
            circos.yaxis(side = "right", sector.index =CELL_META$sector.index,
                         track.index = get.cell.meta.data("track.index"), labels.cex = 0.4)
            circos.trackHist( factors =csv_f5[["Found_chromosome"]], x=csv_f5[["Found_start_position_on_chromosome"]], bin.size = 200000, col = input$color5, border = input$color5)
            circos.yaxis(side = "right", sector.index =CELL_META$sector.index,
                         track.index = get.cell.meta.data("track.index"), labels.cex = 0.4)
            circos.trackHist( factors =csv_f6[["Found_chromosome"]], x=csv_f6[["Found_start_position_on_chromosome"]], bin.size = 200000, col = input$color6, border = input$color6)
            circos.yaxis(side = "right", sector.index =CELL_META$sector.index,
                         track.index = get.cell.meta.data("track.index"), labels.cex = 0.4)


            circos.track(ylim = c(0, 1),
                         bg.col = "#edece8",
                         bg.border = "black", track.height = 0.05, panel.fun = function(x, y) {
                           circos.text(CELL_META$xcenter, CELL_META$ycenter, CELL_META$sector.index, facing = "inside", niceFacing = TRUE, cex = 0.5)})

            dev.off()
          }
        )


      }
      
      if (input$num_files == 7) {
        req(get("csv_file_2", envir = .GlobalEnv))
        req(get("csv_file_1", envir = .GlobalEnv))
        req(get("csv_file_3", envir = .GlobalEnv))
        req(get("csv_file_4", envir = .GlobalEnv))
        req(get("csv_file_5", envir = .GlobalEnv))
        req(get("csv_file_6", envir = .GlobalEnv))
        req(get("csv_file_7", envir = .GlobalEnv))
        
        
        csv_f1 <- get("csv_file_1", envir = .GlobalEnv)
        csv_f2 <- get("csv_file_2", envir = .GlobalEnv)
        csv_f3 <- get("csv_file_3", envir = .GlobalEnv)
        csv_f4 <- get("csv_file_4", envir = .GlobalEnv)
        csv_f5 <- get("csv_file_5", envir = .GlobalEnv)
        csv_f6 <- get("csv_file_6", envir = .GlobalEnv)
        csv_f7 <- get("csv_file_7", envir = .GlobalEnv)
        
        circos.clear()
        circos.par("start.degree" = 83.5, "gap.degree" = c(2,2,2,2,2,2,12), "track.height" = 0.1, "track.margin" = c(0, 0.01))
        circos.genomicInitialize(df_cuc, sector.names = df_cuc$name, factors = df_cuc$name, axis.labels.cex = 0.4, plotType = "axis")
        
        
        circos.trackHist( factors = csv_f1[["Found_chromosome"]], x=csv_f1[["Found_start_position_on_chromosome"]], bin.size = 200000, col = input$color1, border = input$color1)
        circos.yaxis(side = "right", sector.index =CELL_META$sector.index,
                     track.index = get.cell.meta.data("track.index"), labels.cex = 0.4)
        circos.trackHist( factors =csv_f2[["Found_chromosome"]], x=csv_f2[["Found_start_position_on_chromosome"]], bin.size = 200000, col = input$color2, border = input$color2)
        circos.yaxis(side = "right", sector.index =CELL_META$sector.index,
                     track.index = get.cell.meta.data("track.index"), labels.cex = 0.4)
        circos.trackHist( factors =csv_f3[["Found_chromosome"]], x=csv_f3[["Found_start_position_on_chromosome"]], bin.size = 200000, col = input$color3, border = input$color3)
        circos.yaxis(side = "right", sector.index =CELL_META$sector.index,
                     track.index = get.cell.meta.data("track.index"), labels.cex = 0.4)
        circos.trackHist( factors =csv_f4[["Found_chromosome"]], x=csv_f4[["Found_start_position_on_chromosome"]], bin.size = 200000, col = input$color4, border = input$color4)
        circos.yaxis(side = "right", sector.index =CELL_META$sector.index,
                     track.index = get.cell.meta.data("track.index"), labels.cex = 0.4)
        circos.trackHist( factors =csv_f5[["Found_chromosome"]], x=csv_f5[["Found_start_position_on_chromosome"]], bin.size = 200000, col = input$color5, border = input$color5)
        circos.yaxis(side = "right", sector.index =CELL_META$sector.index,
                     track.index = get.cell.meta.data("track.index"), labels.cex = 0.4)
        circos.trackHist( factors =csv_f6[["Found_chromosome"]], x=csv_f6[["Found_start_position_on_chromosome"]], bin.size = 200000, col = input$color6, border = input$color6)
        circos.yaxis(side = "right", sector.index =CELL_META$sector.index,
                     track.index = get.cell.meta.data("track.index"), labels.cex = 0.4)
        circos.trackHist( factors =csv_f7[["Found_chromosome"]], x=csv_f7[["Found_start_position_on_chromosome"]], bin.size = 200000, col = input$color7, border = input$color7)
        circos.yaxis(side = "right", sector.index =CELL_META$sector.index,
                     track.index = get.cell.meta.data("track.index"), labels.cex = 0.4)
        
        
        circos.track(ylim = c(0, 1),
                     bg.col = "#edece8",
                     bg.border = "black", track.height = 0.05, panel.fun = function(x, y) {
                       circos.text(CELL_META$xcenter, CELL_META$ycenter, CELL_META$sector.index, facing = "inside", niceFacing = TRUE, cex = 0.5)})
        
        
        output$downloadPlot <- downloadHandler(
          filename = function() {
            paste("plot-", Sys.Date(), ".pdf", sep = "")
          },
          content = function(file) {
            pdf(file)
            circos.clear()
            circos.par("start.degree" = 83.5, "gap.degree" = c(2,2,2,2,2,2,12), "track.height" = 0.1, "track.margin" = c(0, 0.01))
            circos.genomicInitialize(df_cuc, sector.names = df_cuc$name, factors = df_cuc$name, axis.labels.cex = 0.4, plotType = "axis")
            
            
            circos.trackHist( factors = csv_f1[["Found_chromosome"]], x=csv_f1[["Found_start_position_on_chromosome"]], bin.size = 200000, col = input$color1, border = input$color1)
            circos.yaxis(side = "right", sector.index =CELL_META$sector.index,
                         track.index = get.cell.meta.data("track.index"), labels.cex = 0.4)
            circos.trackHist( factors =csv_f2[["Found_chromosome"]], x=csv_f2[["Found_start_position_on_chromosome"]], bin.size = 200000, col = input$color2, border = input$color2)
            circos.yaxis(side = "right", sector.index =CELL_META$sector.index,
                         track.index = get.cell.meta.data("track.index"), labels.cex = 0.4)
            circos.trackHist( factors =csv_f3[["Found_chromosome"]], x=csv_f3[["Found_start_position_on_chromosome"]], bin.size = 200000, col = input$color3, border = input$color3)
            circos.yaxis(side = "right", sector.index =CELL_META$sector.index,
                         track.index = get.cell.meta.data("track.index"), labels.cex = 0.4)
            circos.trackHist( factors =csv_f4[["Found_chromosome"]], x=csv_f4[["Found_start_position_on_chromosome"]], bin.size = 200000, col = input$color4, border = input$color4)
            circos.yaxis(side = "right", sector.index =CELL_META$sector.index,
                         track.index = get.cell.meta.data("track.index"), labels.cex = 0.4)
            circos.trackHist( factors =csv_f5[["Found_chromosome"]], x=csv_f5[["Found_start_position_on_chromosome"]], bin.size = 200000, col = input$color5, border = input$color5)
            circos.yaxis(side = "right", sector.index =CELL_META$sector.index,
                         track.index = get.cell.meta.data("track.index"), labels.cex = 0.4)
            circos.trackHist( factors =csv_f6[["Found_chromosome"]], x=csv_f6[["Found_start_position_on_chromosome"]], bin.size = 200000, col = input$color6, border = input$color6)
            circos.yaxis(side = "right", sector.index =CELL_META$sector.index,
                         track.index = get.cell.meta.data("track.index"), labels.cex = 0.4)
            circos.trackHist( factors =csv_f7[["Found_chromosome"]], x=csv_f7[["Found_start_position_on_chromosome"]], bin.size = 200000, col = input$color7, border = input$color7)
            circos.yaxis(side = "right", sector.index =CELL_META$sector.index,
                         track.index = get.cell.meta.data("track.index"), labels.cex = 0.4)
            
            
            circos.track(ylim = c(0, 1),
                         bg.col = "#edece8",
                         bg.border = "black", track.height = 0.05, panel.fun = function(x, y) {
                           circos.text(CELL_META$xcenter, CELL_META$ycenter, CELL_META$sector.index, facing = "inside", niceFacing = TRUE, cex = 0.5)})
            
            dev.off()
          }
        )
        
        
      }
      

      
    }) 
  })
}



# Combined UI
ui <- navbarPage("Create Circle Plot for B10 Cucumber ",
                 tabPanel("Convert ctgs to chr", ui1),
                 tabPanel("Create circle plot", ui2)
)

# Combined Server
server <- function(input, output, session) {
  #callModule(server1, "app1")  # optional: if modularized
  server1(input, output, session)
  server2(input, output, session)
}

shinyApp(ui, server)
