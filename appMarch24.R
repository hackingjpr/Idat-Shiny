source("./source_functionsNov.R")


ui <- shiny::fluidPage(
  
  #Loading bar
  useAttendant(),
  #Dashboard looks etc
  dashboardPage(dashboardHeader(title= "", dropdownMenuOutput("messageMenu")), 
                dashboardSidebar(
                  #logo
                  # (img(src='Free_Sample_By_Wix (5).jpg', align = "center")),
                  (img(src='Free_Sample_By_Wix%20(5).jpg', align = "center")),
                  #uploading idat
                  h2("Step 1"),
                  fileInput(
                    "idatFile",
                    "Upload idat file:",
                    multiple = TRUE,
                    accept = ".idat"),
                  #Step 2
                  h2("Step 2"),
                  selectInput("metagenes", "Select Metagene Set", c("MRT (ATRT & ECRT)", "ATRT", "ECRT")),
                  #set it away
                  h2("Step 3"),

                  
                  actionBttn(
                    inputId = "bttn1",
                    label = "Generate Risk Values",
                    color = "primary",
                    style = "jelly"
                  ),
                  
                  
                  #Loading Bar
                  attendantBar("progress-bar"),
                  #Reset button
                  column(width = 12,
                         "Want to upload a different data set?"),
                  column(
                    width = 12,
                    actionBttn(
                      inputId = "bttn2",
                      label = "Reset",
                      color = "danger",
                      style = "fill",
                      size = "sm"
                    ),
                    align = "center"
                  ),
                  
                  hr(),
                  hr(),
                  hr(),
                  hr(),

                  tags$footer(
                    hr(),
                    HTML(
                    paste(
                      "Author: James Hacking,", "<br/>",  
                      "Date Created: 28-04-2022,", "<br/>", 
                      "Copyright (c) James Hacking, 2022,", "<br/>",
                      h5(
                        a("Email: james.hacking@ncl.ac.uk", href = "mailto:james.hacking@ncl.ac.uk")
                      ),
                      h4(
                        "Disclaimer : This app is designed exclusively for research purposes and is strictly not for diagnostic or clinical use."
                      )
                      
                    )))
                ), dashboardBody((
                  tabsetPanel( id = "tabs",
                               #Information Tab
                               tabPanel(
                                 "Info",
                                 (fluidRow(
                                   includeMarkdown("./introduction.md")
                                  )
                                 )
                               ),
                               #Tutorial Tab
                               tabPanel("Tutorial",
                                        includeMarkdown("./Tutorial/tutorial.md")
                                          ),
                               #Results Tab
                               tabPanel("Results Table",
                                        (fluidRow(
                                          box(
                                            width = 12,
                                            title = "Risk Values",
                                            status = "info",
                                            solidHeader = TRUE,
                                            DTOutput('Mval')
                                          ))),
                                        (fluidRow(
                                          box(
                                            title = "Risk plot", 
                                            status = "warning", 
                                            solidHeader = TRUE,
                                            collapsible = FALSE,
                                            plotOutput("figure")),
                                          box(
                                            title = "Selections",
                                            status = "success",
                                            solidHeader = TRUE,
                                            h3("Metagene Set:"),
                                            textOutput("metagenechoice"),
                                            h3("Sample Selected:"),
                                            textOutput("sample"),
                                            h3("Patient's Risk Percentile:"),
                                            textOutput("percentages")
                                          ))
                                        )),
                               
                               #Download Tab
                               tabPanel("Download",
                                        textInput("filename", "Please insert desired filename", "M-values"),
                                        radioButtons(inputId = "download", label = "Select file type", choices = c("csv", "pdf")),
                                        downloadButton("down", "Download the results"),
                               )
                  )
                )),
                #Font Selection            
                tags$head(tags$style(HTML('* {font-family: "Courier New"};')))),
  
)

server <- function(session, input, output) {
  # Allow larger files to be uploaded
  options(shiny.maxRequestSize = 30 * 1024 ^ 2)
  #Hide tabs until graph made
  hideTab(inputId = "tabs", target = "Results Table")
  hideTab(inputId = "tabs", target = "Download")
  
  #Loading bar
  att <- Attendant$new("progress-bar")
  unlink("./temp/", recursive = T)
  if(!dir.exists("./temp")){
    dir.create("./temp")
  }
  #fix the overwriting issue
  fs::path("./temp/",createRandString()) -> tempDIR
  dir.create(tempDIR)
  
  Mvals <- reactiveValues(data = NULL)
    observeEvent(input$bttn1, {
      message("Loading Start")
      att$set(10, text = "Loading") #Start at 10% 
      message("Loading Move")
      att$auto(ms = 1600, value = 0.01) # automatically increment
      
      #need idat file to run
      req(input$idatFile)
      require(R.utils)
      
      message("input")
      input$idatFile$datapath -> in.files
      message(in.files)
      paste0(tempDIR,"/", input$idatFile$name) -> out.files
      message(out.files)
    
      copied <- file.copy(in.files, out.files)
      message("Copied")
      
      temp.base <- get_basenames(tempDIR)
      message("get_basenames done")
      cat("Timing start\n")
      ptm <- proc.time()
      
      temp.processed <- process_idats(temp.base)
      on.exit({
        att$done()
      })
      message("process_idats done")
      # test.res <-(  
      #   
      # if (input$metagenes == "MRT (ATRT & ECRT)") {
      #     ALL -> meta
      #   } %>%
      # if (input$metagenes == "ATRT") {
      #   ATRT -> meta
      # } %>%
      # if (input$metagenes == "ECRT") {
      #   ECRT -> meta
      # } %>%
      # #this needs to be reactive
      #  #output$test.res <- reactive({extract.metagene(
      #  extract.metagene(
      #   as.character(meta[[1]]$genes),
      #   as.numeric(meta[[1]]$weights),
      #   beta2m(temp.processed$betas),
      #   as.numeric(meta[[2]])
      # ))
     #)}))
      meta <- reactive({
        if (input$metagenes == "MRT (ATRT & ECRT)") {
          return(ALL)
        } else if (input$metagenes == "ATRT") {
          return(ATRT)
        } else if (input$metagenes == "ECRT") {
          return(ECRT)
        } else {
          # Handle other cases if needed
          return(NULL)
        }
      })
      
      # test.res <- reactive({
      test.res <- if (!is.null(meta())) {
          extract.metagene(
            as.character(meta()[[1]]$genes),
            as.numeric(meta()[[1]]$weights),
            beta2m(temp.processed$betas),
            as.numeric(meta()[[2]])
          )
        }
      # })
      message("test.res done")
      message(test.res)
      # as.numeric(test.res)
      round(test.res, digits = 3) -> test.res
      output$resultText <- renderText({
        # Display the result in the UI
        req(test.res())  # Ensure that test.res() is not NULL
        result <- round(test.res(), digits = 3)
        paste("Result: ", result)
      })
    
      message("round test.res done")
      output$Mval <- renderDT (({test.res
      }),
      options = list(
        pageLength = 10, 
        processing=FALSE),
      selection = "single"
      )
      message("render dt done")
      print(test.res)
      figure.input <- test.res$Risk_Value
      names(figure.input) <- rownames(test.res)
      print(figure.input)
      
      att$done(text = "Complete")


      output$time <- renderText({proc.time() - ptm})
      rowSelect <- reactive({input$Mval_rows_selected})
      
      output$figure <-
        renderPlot(
          
          # generate_figure_highlight_mrt(
          #   figure.input
          #   ,1)

# generate_figure_highlight_mrt(
#   figure.input
#   ,input$Mval_row_last_clicked)

          if (input$metagenes == "MRT (ATRT & ECRT)") {
            return(
            generate_figure_highlight_mrt(
              figure.input
              ,input$Mval_row_last_clicked))
          }
          else if (input$metagenes == "ATRT") {
            return(
            generate_figure_highlight_atrt(
              figure.input
              ,input$Mval_row_last_clicked))
          }
          else if (input$metagenes == "ECRT") {
            return(
            generate_figure_highlight_ecrt(
              figure.input
              ,input$Mval_row_last_clicked))
          }

          # https://rstudio.github.io/DT/shiny.html


        )
      output$percentages <- renderText (
        
        # input <- "./mrt54.dist.rds",
        # generate_figure_percentage_mrt(
        #   figure.input,
        #   input$Mval_row_last_clicked)
          # ,input$Mval_row_last_clicked)
        # "hello",
        
        if (input$metagenes == "MRT (ATRT & ECRT)") {
          figureFile <- "./mrt54.dist.rds"
          return(

          generate_figure_percentage_mrt(
            figure.input
            ,input$Mval_row_last_clicked))
        }
        else if (input$metagenes == "ATRT") {
          figureFile <- "./atrt8.dist.rds"
          return(
          
          generate_figure_percentage_atrt(
            figure.input
            ,input$Mval_row_last_clicked))
        }
        else if (input$metagenes == "ECRT") {
          figureFile <- "./ecrt20.dist.rds"
          return(

          generate_figure_percentage_ecrt(
            figure.input
            ,input$Mval_row_last_clicked))
        }
      )
        # {generate_figure_percentage(figure.input, input$Mval_row_last_clicked)})
      
      output$metagenechoice <- renderText({input$metagenes})
      
      showTab(inputId = "tabs", target = "Results Table", select = TRUE)
      showTab(inputId = "tabs", target = "Download")
      
      #For displaying currently selected sample

      output$sample <- renderText(input$Mval_cell_clicked$value) 
      Mvals$data <- if (!is.null(meta())) {
        extract.metagene(
          as.character(meta()[[1]]$genes),
          as.numeric(meta()[[1]]$weights),
          beta2m(temp.processed$betas),
          as.numeric(meta()[[2]])
        )
      }
    })
  #Reset session and delete tempDIR
  observeEvent(input$bttn2, {
    unlink(tempDIR, recursive = T)
    session$reload()
    return()
    print("session reload not working")
  })
  renderText(output$metagenes <- input$metagenes)
  # output$down <- downloadHandler(
  #   filename = function() {
  #     paste(input$filename,Sys.time(), ".csv", sep="_")
  #   },
  #   content ={ 
  #     function(file) {
  #       write.csv(Mvals(), file, row.names = TRUE)
  #     }
  #     # if (input$download == "csv")
  #     #   write.csv(output$Mval)
  #     # else
  #     #   pdf(output$Mval,
  #     #       width = 14
  #     #   )
  #     # dev.off()
  #   })
  output$down <- downloadHandler(
    filename = function() {
      paste(input$filename, Sys.time(), ifelse(input$download == "csv", ".csv", ".pdf"), sep = "_")
    },
    content = function(file) {
      if (input$download == "csv") {
        write.csv(Mvals$data, file, row.names = TRUE)
      } else if (input$download == "pdf") {
        pdf(file, width = 14, height = 8)
        
        Mvals$data <- round(Mvals$data, digits =3)
        grid.table(Mvals$data)
        # 
        # message("input")
        # input$idatFile$datapath -> in.files
        # message(in.files)
        # paste0(tempDIR,"/", input$idatFile$name) -> out.files
        # message(out.files)
        # 
        # copied <- file.copy(in.files, out.files)
        # message("Copied")
        # 
        # temp.base <- get_basenames(tempDIR)
        # message("get_basenames done")
        # cat("Timing start\n")
        # ptm <- proc.time()
        # 
        # temp.processed <- process_idats(temp.base)
        # on.exit({
        #   att$done()
        # })
        # message("process_idats done")
        # # test.res <-(  
        # #   
        # # if (input$metagenes == "MRT (ATRT & ECRT)") {
        # #     ALL -> meta
        # #   } %>%
        # # if (input$metagenes == "ATRT") {
        # #   ATRT -> meta
        # # } %>%
        # # if (input$metagenes == "ECRT") {
        # #   ECRT -> meta
        # # } %>%
        # # #this needs to be reactive
        # #  #output$test.res <- reactive({extract.metagene(
        # #  extract.metagene(
        # #   as.character(meta[[1]]$genes),
        # #   as.numeric(meta[[1]]$weights),
        # #   beta2m(temp.processed$betas),
        # #   as.numeric(meta[[2]])
        # # ))
        # #)}))
        # meta <- reactive({
        #   if (input$metagenes == "MRT (ATRT & ECRT)") {
        #     return(ALL)
        #   } else if (input$metagenes == "ATRT") {
        #     return(ATRT)
        #   } else if (input$metagenes == "ECRT") {
        #     return(ECRT)
        #   } else {
        #     # Handle other cases if needed
        #     return(NULL)
        #   }
        # })
        # 
        # # test.res <- reactive({
        # test.res <- if (!is.null(meta())) {
        #   extract.metagene(
        #     as.character(meta()[[1]]$genes),
        #     as.numeric(meta()[[1]]$weights),
        #     beta2m(temp.processed$betas),
        #     as.numeric(meta()[[2]])
        #   )
        # }
        # # })
        # message("test.res done")
        # message(test.res)
        # # as.numeric(test.res)
        # round(test.res, digits = 3) -> test.res
        # output$resultText <- renderText({
        #   # Display the result in the UI
        #   req(test.res())  # Ensure that test.res() is not NULL
        #   result <- round(test.res(), digits = 3)
        #   paste("Result: ", result)
        # })
        # 
        # message("round test.res done")
        # output$Mval <- renderDT (({test.res
        # }),
        # options = list(
        #   pageLength = 10, 
        #   processing=FALSE),
        # selection = "single"
        # )
        # message("render dt done")
        # print(test.res)
        # figure.input <- test.res$Risk_Value
        # names(figure.input) <- rownames(test.res)
        # print(figure.input)
        # 
        # att$done(text = "Complete")
        # 
        # 
        # output$time <- renderText({proc.time() - ptm})
        # rowSelect <- reactive({input$Mval_rows_selected})
        # 
        #  x <-
        #   grid.draw(
        #     
        #     # generate_figure_highlight_mrt(
        #     #   figure.input
        #     #   ,1)
        #     
        #     generate_figure_highlight_mrt(
        #       figure.input
        #       ,input$Mval_row_last_clicked)
        #     # 
        #     # if (input$metagenes == "MRT (ATRT & ECRT)") {
        #     #   generate_figure_highlight_mrt(
        #     #     figure.input
        #     #     ,input$Mval_row_last_clicked)
        #     # },
        #     # elseif (input$metagenes == "ATRT") {
        #     #   generate_figure_highlight_atrt(
        #     #     figure.input
        #     #     ,input$Mval_row_last_clicked)
        #     # },
        #     # elseif (input$metagenes == "ECRT") {
        #     #   generate_figure_highlight_ecrt(
        #     #     figure.input
        #     #     ,input$Mval_row_last_clicked)
        #     # }
        #     
        #     # https://rstudio.github.io/DT/shiny.html
        #     
        #     
        #   )
        #  x
        # 
        # # Print DataTable
        # # print(DT::datatable(Mvals$data))
        # # 
        # # # Print ggplot figure
        # # print(generate_figure_highlight_mrt(Mvals$data, 1))
        # 
        
        dev.off()
      }
    }
  )
  
  
  
  
  
}

shinyApp(server = server, ui = ui)
