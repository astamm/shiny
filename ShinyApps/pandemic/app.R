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
library(ggplot2)

source("setup_python_env.R")

# Define any Python packages needed for the app here:
PYTHON_DEPENDENCIES <- c('numpy', 'scipy', 'matplotlib')

# Define the js method that resets the page
jsResetCode <- "shinyjs.reset = function() {history.go(0)}"

# Define UI for application that draws a histogram
ui <- fluidPage(
    # Include shinyjs in the UI
    useShinyjs(),
    # Add the js code to the page
    extendShinyjs(text = jsResetCode, functions = "reset"),

    titlePanel("Modèles de propagation d'épidémie"),

    sidebarLayout(

        # ---------------- Sidebar panel with changeable inputs ----------------- #
        sidebarPanel(

            h4('Choix du modèle'),
            br(),

            radioButtons(
                inputId = 'dist',
                label = "Modèle pour le nombre d'infectés:",
                choices = c(
                    'Loi binomiale négative' = 'negbin',
                    'Loi de Poisson' = 'pois',
                    'Constante' = 'const'
                )
            ),
            sliderInput(
                inputId = 'mean',
                label = "Nombre moyen d'infectés sur le long terme",
                value = 2.5,
                min = 0,
                max = 10,
                step = 0.5
            ),
            sliderInput(
                inputId = 'variance',
                label = "Variabilité autour du nombre moyen d'infectés",
                value = 65,
                min = 0, max = 100,
                step = 0.5
            ),
            sliderInput(
                inputId = 'support',
                label = "Taille du support",
                value = 10,
                min = 0, max = 100,
                step = 1
            ),
            actionButton(
                inputId = "run",
                label = "Lancer la simulation"
            ),
            hr(),

            h4("Simulation de l'effet des mesures sanitaires"),
            br(),

            sliderInput(
                inputId = 'trunc',
                label = "Nombre maximal de personnes qu'un infecté peut infecter à son tour",
                value = 20,
                min = 0,
                max = 20,
                step = 1
            ),
            checkboxInput(
                inputId = "gotrunc",
                label = "Mise en vigueur des mesures",
                value = FALSE
            )
        ),

        # ---------------- Sidebar panel with changeable inputs ----------------- #
        mainPanel(

            # Output: Tabset w/ plot, summary, and table ----
            tabsetPanel(
                type = 'tabs',
                tabPanel(
                    'Modèle',
                    h3('Visualisation du modèle choisi'),
                    hr(),
                    withSpinner(plotOutput('distPlot')),
                    h3('Chaine de propagation'),
                    hr(),
                    htmlOutput("graphSumm"),
                    withSpinner(imageOutput('graphPlot'))
                ),
                tabPanel(
                    'Architecture Info',
                    h3('Current architecture info'),
                    '(These values will change when app is run locally vs on Shinyapps.io)',
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
    reticulate::source_python('pandemic.py')

    observeEvent(input$dist, {
        if (input$dist != "negbin")
            updateSliderInput(
                session = session,
                inputId = "variance",
                value = if (input$dist == "pois")
                    input$mean
                else if (input$dist == "const")
                    0
            )
        shinyjs::toggleState(id = "variance", condition = input$dist == "negbin")
    })

    # Generate the requested distribution
    rfunc <- reactive({
        if (input$dist == "negbin")
            stopifnot(input$mean < input$variance)
        switch(
            EXPR   = input$dist,
            negbin = function(n) as.integer(rnbinom(
                n = n,
                size = input$mean^2 / (input$variance - input$mean),
                mu = input$mean
            )),
            pois   = function(n) as.integer(rpois(
                n = n,
                lambda = input$mean
            )),
            const  = function(n) as.integer(rep(round(input$mean), n))
        )
    })

    dfunc <- reactive({
        if (input$dist == "negbin")
            stopifnot(input$mean < input$variance)
        switch(
            EXPR   = input$dist,
            negbin = function(x) dnbinom(
                x = x,
                size = input$mean^2 / (input$variance - input$mean),
                mu = input$mean
            ),
            pois   = function(x) dpois(
                x = x,
                lambda = input$mean
            ),
            const  = function(x) ifelse(x == input$mean, 1, 0)
        )
    })

    output$distPlot <- renderPlot({
        Value <- 0:input$support
        Probability <- dfunc()(Value)
        ggplot(
            data = data.frame(Value, Probability),
            mapping = aes(as.factor(Value), Probability)
        ) +
            geom_col() +
            theme_bw() +
            labs(
                x = "Nombre de personnes qu'un infecté peut infecter à son tour",
                y = "Probabilité"
            )
    })

    tree <- eventReactive(input$run, {
        generate_tree(
            maxgen = 6,
            minpatients = 3,
            law = function() rfunc()(1)
        )
    })

    output$graphPlot <- renderImage(
        expr = {
            file <- tempfile(fileext = ".svg")
            if (!input$gotrunc)
                plot_truncated_tree(tree(), maxgen = 6, file = file)
            else
                plot_truncated_tree(tree()$copy(input$trunc), maxgen = 6, file = file)
            list(src = file)
        },
        deleteFile = TRUE
    )

    output$graphSumm <- renderUI({
        if (!input$gotrunc)
            info <- summarize_tree(tree())
        else
            info <- summarize_tree(tree()$copy(input$trunc))
        info1 <- paste0("Nombre d'infectés : ", info[1])
        info2 <- paste0("Nombre de générations : ", info[2], " / 6")
        HTML(paste(info1, info2, sep = "<br/>"))
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
