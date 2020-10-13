#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinycssloaders)
library(shinyjs)
library(DT)
library(plotly)
library(jsonlite)

source("setup_python_env.R")

# Define any Python packages needed for the app here:
PYTHON_DEPENDENCIES <- c('numpy', 'gudhi', 'plotly')

# Define the js method that resets the page
jsResetCode <- "shinyjs.reset = function() {history.go(0)}"

# Define UI for application that draws a histogram
ui <- fluidPage(
    # Include shinyjs in the UI
    useShinyjs(),
    # Add the js code to the page
    extendShinyjs(text = jsResetCode, functions = "reset"),

    titlePanel("A visualisation tool for TDA using GUDHI"),

    sidebarLayout(

        # ---------------- Sidebar panel with changeable inputs ----------------- #
        sidebarPanel(

            fileInput(
                inputId = "file",
                label = "Data"
            ),
            sliderInput(
                inputId = 'alpha',
                label = "Filtration value",
                value = 0.005,
                min = 0.0,
                max = 0.01,
                step = 0.0001
            )

        ),

        # ---------------- Sidebar panel with changeable inputs ----------------- #
        mainPanel(

            tabsetPanel(
                type = 'tabs',
                tabPanel(
                    'Interactive 3D view',
                    withSpinner(plotly::plotlyOutput(
                        outputId = 'filtrationPlot',
                        height = "800px"
                    ))
                ),
                tabPanel(
                    'Architecture Info',
                    h3('Current architecture info'),
                    hr(),
                    withSpinner(DT::dataTableOutput('sysinfo')),
                    actionButton("reset_button", "Reset App")
                )
            )

        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {

    # Call the method from
    # somewhere within the server
    observeEvent(input$reset_button, {js$reset()})

    # ------------------ App virtualenv setup (Do not edit) ------------------- #

    virtualenv_dir <- Sys.getenv('VIRTUALENV_NAME')
    python_path <- Sys.getenv('PYTHON_PATH')

    # Create virtual env and install dependencies
    reticulate::virtualenv_create(
        envname = virtualenv_dir,
        python = python_path
    )

    reticulate::virtualenv_install(
        envname = virtualenv_dir,
        packages = PYTHON_DEPENDENCIES,
        ignore_installed = TRUE
    )

    reticulate::use_virtualenv(
        virtualenv = virtualenv_dir,
        required = TRUE
    )

    # ------------------ App server logic (Edit anything below) --------------- #

    # Import python functions to R
    reticulate::source_python('gudhi-viz.py')

    data <- reactive({
        compute_filtration(req(input$file$datapath), input$alpha)
    })

    output$filtrationPlot <- plotly::renderPlotly({
        plotly::as_widget(
            jsonlite::fromJSON(
                txt = compute_figure(data()),
                simplifyVector = FALSE
            )
        )
    })

    # Display info about the system running the code
    output$sysinfo <- DT::renderDataTable({
        s <- Sys.info()
        df <- data.frame(
            Info_Field = names(s),
            Current_System_Setting = as.character(s)
        )
        DT::datatable(
            data = df,
            rownames = FALSE,
            selection = 'none',
            style = 'bootstrap',
            filter = 'none',
            options = list(dom = 't')
        )
    })

}

# Run the application
shinyApp(ui = ui, server = server)
