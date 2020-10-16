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

source("setup_python_env.R")

# Define any Python packages needed for the app here:
PYTHON_DEPENDENCIES <- c('numpy', 'gudhi', 'plotly')

# Define the js method that resets the page
jsResetCode <- "shinyjs.reset = function() {history.go(0)}"

# Define UI for application that draws a histogram
ui <- fluidPage(
    # Include shinyjs in the UI
    shinyjs::useShinyjs(),
    # Add the js code to the page
    shinyjs::extendShinyjs(
        text = jsResetCode,
        functions = "reset"
    ),

    titlePanel("A visualisation tool for TDA using GUDHI"),

    sidebarLayout(

        # ---------------- Sidebar panel with changeable inputs ----------------- #
        sidebarPanel(

            fileInput(
                inputId = "file",
                label = "Data"
            ),
            sliderInput(
                inputId = 'fps',
                label = "Frames per second",
                value = 25,
                min = 0,
                max = 50,
                step = 1
            ),
            sliderInput(
                inputId = 'alpha',
                label = "Maximal simplex diameter",
                value = 0.5,
                min = 0.0,
                max = 1.0,
                step = 0.1
            )

        ),

        # ---------------- Sidebar panel with changeable inputs ----------------- #
        mainPanel(

            tabsetPanel(
                type = 'tabs',
                tabPanel(
                    'Simplicial Complexes',
                    withSpinner(htmlOutput(
                        outputId = 'complexPlot',
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

    # Import C++ functions to R
    Rcpp::sourceCpp("gudhi-viz.cpp")

    point_cloud <- reactive({
        load_off_file(req(input$file$datapath))
    })

    distance_matrix <- reactive({
        dist(point_cloud())
    })

    diameter_ub <- reactive({
        dm <- distance_matrix()
        n <- attr(dm, "Size")
        round(
            x = get_diameter_upper_bound(dm, n, 2),
            digits = 3
        )
    })

    simplicial_complexes <- reactive({
        build_complexes(
            pts = point_cloud(),
            max_diameter = diameter_ub(),
            dimension = 2
        )
    })

    diameter_lb <- reactive({
        round(
            x = get_diameter_lower_bound(
                complexes = simplicial_complexes(),
                dimension = 2
            ),
            digits = 3
        )
    })

    observe({
        lb <- diameter_lb()
        ub <- diameter_ub()
        updateSliderInput(
            session = session,
            inputId = "alpha",
            value = (lb + ub) / 2,
            min = lb,
            max = ub,
            step = round((ub - lb) / 99, digits = 3),
        )
    })

    output$complexPlot <- renderUI({
        sync_figures(simplicial_complexes(), input$alpha, input$fps)
        includeHTML("temp-plot.html")
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
