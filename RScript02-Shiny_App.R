################################################################################
#        1         2         3         4         5          6        7         8                   
#2345678901234567890123456789012345678901234567890123456789012345678901234567890
################################################################################
# Shiny App to calculate multipliers applied to offset projects to account for
# uncertainty in the efficacy of the proposed projects
# the App is based on Brandford (2017) and vanderlee et al. (in review)
# the method used Monte Carlo simulations of project values with associated
# uncertainty to provide a ratio of relative value for impacts and offset
# by taking some high percentile (> 80%) of the ratio distribution give and high
# likelihood of achieving no Net Loss (NNL)
# The app can take multiple impact and offset projects with different values
# and uncertainty as well as different weights (ability to apply multiplers)
# for offsets - resulting in multiple offset value to be calculated and applied
# to the different offset projects.
# User controls:
# - the number of simulations draw and risk tolerance
# - the number of impact and offset projects
# - the method of uncertainty: CV, SD, Categorical, or Pert dist'n
# - the parameterization of each project conditional of uncertainty method
# User has control over the number of impacts and offsets
#-------------------------------------------------------------------------------

# Clear workspace
rm(list=ls())

# load libraries
library(shiny)
library(bslib)
library(bsicons)
library(shinyWidgets)
library(shinyvalidate)
library(data.table)
library(ggplot2)

#-------------------------------------------------------------------------------
# Functions and settings
#-------------------------------------------------------------------------------

# ---- SETTINGS ----------------------------------------------------------------

# Initial values - default values for model inputs
init_val <- list(
  proj_label = "1",# project name
  err_type = "CV", # error type: CV, CS, pert, or u_cat
  mean = 1,        # mean
  cv = 0.5,        # coefficient of variation
  sd = 0.5,        # standard deviation
  mode = 0.5,      # mode
  min = 0,         # minimum
  max = 3.5,       # maximum
  u_cat = "Medium",# uncertainty category: Low, Medium, High = 0.25, 0.5, 1.0 CV 
  CR_w = 1         # offset weight 0-1 
)

# tool tip icon
tool_tip_icon <- bs_icon("info-circle")

# ggplot theme
theme_me <- theme_bw() +
  theme(axis.title = element_text(size = 14, family = "sans", face = "bold"),
        axis.text.x = element_text(size = 12, family = "sans", colour = "black"),
        axis.text.y = element_text(size = 12, family = "sans", hjust = 0.6,
                                   angle = 90, colour = "black"),
        legend.title = element_text(size = 14, family = "sans"),
        legend.text = element_text(size = 14, family = "sans"),
        strip.text = element_text(size = 14, family = "sans", face = "bold"))

# ---- SIMULATION FUNCTIONS ----------------------------------------------------

# function to generate distribution data for impacts or offset 
# function takes user input and simulation a distribution of potential 
# values conditional on inputs
# user can choose: CV, SD, qualitative or Pert option for generation and must
# supply the necessary parameters for each
# teh function with take the parameters and apply a log-normal or pert distn
# to generate simulated value after applying necessary conversions
dist_data.f <- function(data,  # parameterization data
                        n.sim) # number of simulations 
  {
  # number of projects
  n_val <- length(data)
  
  # loop through projects 
  # - extract parameters and create n_val random values conditional on err_type
  # - organize in to matrix - dimension n.sim, n_val
  do.call(cbind, lapply(seq_len(n_val), function(i) {
    
    # extract err_type : "Coefficient of Variation", "Standard Deviation", 
    #                    "Qualitative", or "Pert Distribution"
    # - "Coefficient of Variation" is the default
    err_type <- data[[i]]$par_data$err_type %||% init_val$err_type
    
    # if Standard Deviation
    if (err_type == "SD") {
      
      mean <- data[[i]]$par_data$mean # extract mean
      SD <- data[[i]]$par_data$sd     # extract SD
      
      # Convert parameters for input into rlnorm
      CV <- SD / mean             # calc CV
      sigma <- sqrt(log(CV^2+1))  # convert CV to log-sd of distributions
      mu <- log(mean) - sigma^2/2 # convert mean to log-mean of distributions
      
      # generate random values
      d = rlnorm(n.sim, meanlog = mu,  sdlog = sigma)
   
    # if Coefficient of variation     
    } else if (err_type == "CV") {
      
      mean <- data[[i]]$par_data$mean # extract mean
      CV <- data[[i]]$par_data$cv     # extract CV
      
      # Convert parameters for input into rlnorm
      sigma <- sqrt(log(CV^2+1))   # convert CV to log-sd of distributions
      mu <- log(mean) - sigma^2/2  # convert mean to log-mean of distributions
      
      # generate random values
      d = rlnorm(n.sim, meanlog = mu,  sdlog = sigma)
      
    # if Pert distribution  
    } else if (err_type == "pert") {
      
      mode <- data[[i]]$par_data$mode # extract mode
      min <- data[[i]]$par_data$min   # extract min
      max <- data[[i]]$par_data$max   # extract max
     
      # generate random values - mc2d package
      d = mc2d::rpert(n.sim,       # number of random values
                      mode = mode, # mode
                      min = min,   # min
                      max = max,   # max
                      shape = 4)   # shape parameter - 4 default
                                   # larger -> more skewed distn
    
    # if categorical definition of uncertainty    
    } else if (err_type == "u_cat") {
      
      mean <- data[[i]]$par_data$mean   # extract mean
      u_cat <- data[[i]]$par_data$u_cat # extract uncertainty category
      
      # Set CV based on qualitative value - *TO DO* - improve categories
      if(u_cat == "High") {
        CV = 1.0
      } else if(u_cat == "Medium") { 
        CV = 0.5
      }  else if(u_cat == "Low") { 
        CV = 0.25
      }
      
      # Convert parameters for input into rlnorm
      sigma <- sqrt(log(CV^2+1))  # convert CV to log-sd of distributions
      mu <- log(mean) - sigma^2/2 # convert mean to log-mean of distributions
      
      # generate random values
      d = rlnorm(n.sim, meanlog = mu,  sdlog = sigma)
      
    }
    
    # OUTPUT
    d
  }))
}

# Function to execute dist_data.f for impact(s) and offset(s)
dist.f <- function(
    impact_data,        # Impact data - distribution data and type
    offset_data,         # Offset data - distribution data and type
    n.sim = 10000  # number of draws
) {
  
  # Impact distributions
  d.i <- dist_data.f(data = impact_data, n.sim = n.sim)
  
  # offset distributions
  d.o <- dist_data.f(data = offset_data, n.sim = n.sim)
  
  # output
  list(impacts = d.i, # impact distributions
       offsets = d.o) # offset distributions
}

# Function to simulated and calculated Compensation Ratios
CRu.f <- function(
    d.i,     # impact distributions
    d.o,     # offset distributions
    p = 0.8, # risk tolerance threshold,
    CR_weight
) {
  
  # sum offsets across 
  d_i <- apply(d.i, 1, sum)
  
  # sum offsets across replicates
  d_o <- apply(d.o, 1, sum)
  
  # CR distributions 
  M = d_i/pmax(d_o, 1e-8)  # CR assuming applied across all offsets 
  # CR dist for weighted offsets
  M.weight = 1 + (d_i - d_o) / apply(d.o, 1, function(x) sum(CR_weight * x)) 
  
  # Compensation Ratios
  # Assuming applied across all offsets   
  CR <- round(quantile(M, p) ,2)[[1]]
  
  # Weighted CR by available offsets
  CR.adjusted <- round(quantile(M.weight, p) , 2)[[1]] # values
  CR.adjusted <- (CR.adjusted - 1) * CR_weight + 1     # turn into vector
  
  # Output
  list("EQ" = p,                      # risk tolerance   
       "CR - All" = CR,               # offset multiplier
       "CR - Adjusted" = CR.adjusted, # offset multiplier adjusted by weights
       "M" = M.weight                 # multiplier distribution
  )
}

# ---- VALIDATION --------------------------------------------------------------

# function to generate validation rules for dynamic UI elements - e.g. mean_i.1

add_validation_rules <- function(iv, input, prefix, ns_fun, n) {
  
  # loop through number of impacts/offsets
  for (id in seq_len(n)) {
    
    # extract input elements
    mean <- ns_fun(paste0("mean_", prefix, ".", id))
    cv   <- ns_fun(paste0("cv_", prefix, ".", id))
    sd   <- ns_fun(paste0("sd_", prefix, ".", id))
    mode <- ns_fun(paste0("mode_", prefix, ".", id))
    min  <- ns_fun(paste0("min_", prefix, ".", id))
    max  <- ns_fun(paste0("max_", prefix, ".", id))
    
    # extract error type
    err_id  <- ns_fun(paste0("err_", prefix, ".", id))
    err_type <- input[[err_id]] %||% "CV"
    
    # ---- MEAN RULES ----
    iv$add_rule(mean, sv_required()) # required
    iv$add_rule(mean, sv_numeric())  # numeric
    iv$add_rule(mean, sv_gt(0))      # > 0
    
    # ---- CV RULES ----
    if (err_type == "CV") {
      iv$add_rule(cv, sv_required()) # required
      iv$add_rule(cv, sv_numeric())  # numeric
      iv$add_rule(cv, sv_gte(0))     # >= 0
    }
    
    # ---- SD RULES ----
    if (err_type == "SD") {
      iv$add_rule(sd, sv_required()) # required
      iv$add_rule(sd, sv_numeric())  # numeric
      iv$add_rule(sd, sv_gte(0))     # >= 0
    }
    
    # ---- PERT RULES ----
    if (err_type == "pert") {
      
      ## - MODE -
      iv$add_rule(mode, sv_required()) # required
      iv$add_rule(mode, sv_numeric())  # numeric
      
      ## - Min -
      iv$add_rule(min, sv_required()) # required 
      iv$add_rule(min, sv_numeric())  # numeric
      iv$add_rule(min, sv_gte(0))     # >= 0
      
      ## - Max -
      iv$add_rule(max, sv_required()) # required
      iv$add_rule(max, sv_numeric())  # numeric
      
      ## - Custom Rules -
      # min < mean < max
      # - rules for mode
      iv$add_rule(mode, function(value) {
        
        # extract values - ensure no warnings
        mode_val <- suppressWarnings(as.numeric(value))
        min_val  <- suppressWarnings(as.numeric(input[[min]]))
        max_val  <- suppressWarnings(as.numeric(input[[max]]))
        
        # guards - all scalar values & all not NAs
        if (length(mode_val) != 1 || length(min_val) != 1 ||
            length(max_val) != 1) return(NULL)
        if (any(is.na(c(mode_val, min_val, max_val)))) return(NULL)
        
        # constraints
        if (min_val > mode_val) return("Mode ≥ Min") # mode > min
        if (mode_val > max_val) return("Mode ≤ Max") # mode < max
        
        NULL
      })
      
      # rules for min
      iv$add_rule(min, function(value) {
        
        # extract values - ensure no warnings
        min_val  <- suppressWarnings(as.numeric(value))
        mode_val <- suppressWarnings(as.numeric(input[[mode]]))
        
        # guards - all scalar values & all not NAs
        if (length(min_val) != 1 || length(mode_val) != 1) return(NULL)
        if (any(is.na(c(min_val, mode_val)))) return(NULL)
        
        # constraints
        if (min_val > mode_val) return("Min ≤ Mode") # min < mode
        
        NULL
      })
      
      # Rules for max
      iv$add_rule(max, function(value) {
        
        # extract values - ensure no warnings
        max_val  <- suppressWarnings(as.numeric(value))
        mode_val <- suppressWarnings(as.numeric(input[[mode]]))
        
        # guards - all scalar values & all not NAs
        if (length(max_val) != 1 || length(mode_val) != 1) return(NULL)
        if (any(is.na(c(max_val, mode_val)))) return(NULL)
        
        # constraints
        if (max_val < mode_val) return("Max ≥ Mode") # max > mode
        
        NULL
      })
    }
  }
}

#-------------------------------------------------------------------------------
# MODULES
#-------------------------------------------------------------------------------

# UI Module function for dynamic UI - number of Impact/offsets
projectModule_UI <- function(id) {
  ns <- NS(id)
  uiOutput(ns("ui"))
}

# Server module function for projects
# controls dynamic UI - user inputs for number of impact/offsets
# Initialized parameter data & updates with user input
# organize parameter data for use in simulation
projectModule_server <- function(
    id,                # Module namespace ID
    n,                 # number of elements - reactive val. n.i|n.o
    label,             # "impact" or "offset
    prefix,            # "i" or "o"
    has_weight = FALSE # include CR-weight for offsets
) {
  moduleServer(id, function(input,  # Read user input values (module-scoped)
                            output, # Create reactive outputs (module-scoped)
                            session
  ) {
    
    # server session
    ns <- session$ns
    
    # ---- PARAMETER VALUES ----------------------------------------------------
    # initialized and assign parameter values 
    # - takes init_val if user input isn't availbel yet.
    # Parameter values
    par_vals <- reactive({
      
      # require number or projects
      req(n())
      
      # loop over projects
      lapply(seq_len(n()), function(id) {
        
        # list of parameter values
        list(
          
          err_type = input[[paste0("err_", prefix, ".", id)]] %||% # error type
            init_val$err_type,
          mean     = input[[paste0("mean_", prefix, ".", id)]] %||% # mean
            init_val$mean,
          cv       = input[[paste0("cv_", prefix, ".", id)]] %||% # CV
            init_val$cv,
          sd       = input[[paste0("sd_", prefix, ".", id)]] %||% # SD
            init_val$sd,
          mode     = input[[paste0("mode_", prefix, ".", id)]] %||% # mode
            init_val$mode, 
          min      = input[[paste0("min_", prefix, ".", id)]] %||% # min
            init_val$min,
          max      = input[[paste0("max_", prefix, ".", id)]] %||% # max
            init_val$max,
          u_cat    = input[[paste0("u_cat_", prefix, ".", id)]] %||% # category
            init_val$u_cat,
          name     = input[[paste0("name_", prefix, ".", id)]] %||% # name
            paste(label, id),
          
          # for offset only
          CR_w = if (has_weight)
            input[[paste0("CR_w.", id)]] %||% init_val$CR_w # capacity
          else NULL
          
        ) # close list
      }) # close lapply
    }) # close reactive
    
    # ---- UI ------------------------------------------------------------------
    # UI displayed for n impacts/offsets and uncertainty type
    # changes the number of entry fields with changes to n.i or n.o 
    # changes the error field with selection of uncertainty type
    
    # - FUNCTION - 
    # function factory (function that return function) that produce
    # numeric input UI for different variables
    # tip - text to enter as a tooltip for input
    makeNumericInput <- function(name,      # input name
                                 par_label, # input label in UI
                                 tip = NULL # tool tip text
                                 ) {
      function(id, # project id - i.e. "1" (part of loop)
               pv  # project data par_vals[[id]]
               ) {
        numericInput( # create numeric input
          ns(paste0(name, "_", prefix, ".", id)),  # input id
          label = if (!is.null(tip)) {             # label with tooltip
              tooltip(
                trigger = list(par_label, tool_tip_icon),
                tip,
                placement = "right"
              )
            } else {
              par_label                            # label without tool tip
            },
          # value with gaurds for NULL entries
          value =  pv[[name]],
          min = 0,   # minimum value
          step = 0.1 # step interval
        )
      }
    }
    
    # helper function for tool tip text - makes HTML
    makeTip <- function(...) {
      HTML(paste(..., sep = "<br>"))
    }
    
    # Generate numeric input functions  
    # mean function
    meanInput <- makeNumericInput("mean", "Mean", 
                                  tip = makeTip("Arithmetic mean value")
    )
    # CV function
    cvInput   <- makeNumericInput("cv", "CV", 
                                  tip = makeTip(
                                    "Variability as a proportion of the mean.",
                                    "CV = SD / mean" )
    )
    # SD function                         
    sdInput   <- makeNumericInput("sd", "SD", 
                                  tip = 
                                    HTML("Variability around the mean.
                                         <br> SD = square root of variance")
    )
    # mode function
    modeInput <- makeNumericInput("mode", "Mode", 
                                  tip = makeTip("Most likely value")
    ) 
    # min function
    minInput  <- makeNumericInput("min", "Min",  
                                  tip = makeTip("Smallest possible value") 
                                              
    ) 
    # max function
    maxInput  <- makeNumericInput("max", "Max",  
                                  tip = makeTip("Largest possible value")
    ) 
    
    # - RENDER UI - 
    output$ui <- renderUI({
      
      # require number of projects
      req(n())
      tagList(
        # loop through number of projects 
        # - generate card containing required inputs
        lapply(seq_len(n()), function(id){
          
          # subset par_vals or if doesn't exist empty list
          pv <- par_vals()[[id]]
          
          # Extract error type from user input or pv - whichever is newer
          # or set to default: "Coefficient of Variation"
          err_type <- pv$err_type
          
          # Generate UI
          card(
            # - HEADER -
            # project name and access to select error type in popover
            card_header(
              
              # - PROJECT NAME INPUT -
              popover(
                
                # ---- DISPLAYED TITLE ----
                tags$span(pv$name,
                          bs_icon("pencil-square", class = "ms-1")
                ),
                
                # ---- POPOVER CONTENT ----
                title = "Project name",
                placement = "right",
                layout_column_wrap(
                  textInput(
                    ns(paste0("name_", prefix, ".", id)),
                    label = NULL,
                    value = pv$name
                  ),
                  
                  actionButton(
                    ns(paste0("name_apply_", prefix, ".", id)),
                    "Apply"
                )
                )
              ),
              
              popover(          # create popover
                bs_icon("list"),# icon 
                radioButtons(   # option to pick error type
                  ns(paste0("err_", prefix, ".", id)), # input id err_i/o.*
                  label = NULL,                        # no title
                  # displayed options
                  choiceNames = c("Standard Deviation (SD)",
                                  "Coefficient of Variation (CV)",
                                  "Categorical", 
                                  "Pert Distribution"
                                  ),
                  # value of options 
                  choiceValue = c( "SD","CV", "u_cat", "pert"),
                  selected = err_type # initial choice
                  ), # close radio button
                
                title = "Uncertainty Input", # popover title
                placement = "right",         # pop0ver placed on right
                options = list(trigger = "click"), # close when click off
                # popover text - helpText() controls formatting
                helpText(HTML(paste0("
                   Choose how uncertainty for this <b>", tolower(label), 
                   "</b> is defined:<br><br>

                   <b>Log-normal distribution (default):</b><br>
                   The following options define a log-normal distribution using 
                   the mean and a measure of variability:<br>
                   • <b>Standard Deviation (SD)</b>: absolute variability 
                     around the mean<br>
                   • <b>Coefficient of Variation (CV)</b>: variability relative 
                     to the mean (CV = SD / mean)<br>
                   • <b>Qualitative</b>: select a category that corresponds to
                     a CV value:<br>
                    &nbsp;&nbsp;&nbsp;– High: CV = 1.0<br>
                    &nbsp;&nbsp;&nbsp;– Medium: CV = 0.5<br>
                    &nbsp;&nbsp;&nbsp;– Low: CV = 0.25<br><br>

                    <b>Alternative distribution:</b><br>
                   • <b>Pert Distribution</b>: defines uncertainty using three 
                     points:<br>
                    &nbsp;&nbsp;&nbsp;– Minimum (Min)<br>
                    &nbsp;&nbsp;&nbsp;– Most likely value (Mode)<br>
                    &nbsp;&nbsp;&nbsp;– Maximum (Max)<br><br>

                    The PERT distribution provides a flexible alternative when 
                    uncertainty is better described by bounds and a most likely 
                    value rather than log-normal variability.
                    ")))
              ), # close popover
              class = "d-flex justify-content-between" # CSS control for popover
            ),
            # - USER PARAMETER INPUTS -
            card_body(     # create card to display user inputs             
              fill = TRUE, # grows with container
              if (err_type == "CV") {         # If coefficient of variation
                layout_columns(      # column layout - side-by-side
                  meanInput(id, pv), # input mean
                  cvInput(id, pv)    # input CV
                )
              } else if (err_type == "SD") { # fi standard deviation
                layout_columns(      # column layout - side-by-side
                  meanInput(id, pv), # input mean
                  sdInput(id, pv)    # input SD
                )
              } else if (err_type == "pert"){ # if pert distribution
                tagList(             # row layout 
                  modeInput(id, pv), # input mode
                  layout_columns(    # column layout - side-by-side
                    minInput(id, pv),# input min
                    maxInput(id, pv))# input max
                )
              } else {                         # if qualitative
                layout_columns(      # column layout  - side-by-side
                  meanInput(id, pv), # input mean
                  radioButtons(      # radio button selection for category
                    ns(paste0("u_cat_",prefix,".", id)), # input name
                    label =  tooltip( # label with tooltip
                      trigger = 
                        list(
                          tags$span( 
                            "Uncertainty", # label
                            tool_tip_icon, # tooltip icon
                            # control span to get icon inline iwth label
                            style = "display: inline-flex; align-items: center; 
                                     gap: 4px;")
                        ),
                      # tooltip text
                      makeTip(
                        "Qualitative measure of uncertainty.",
                        "corresponds to categorical choices to a coefficient of 
                        variation (CV):",
                        "High → CV = 1.0",
                        "Medium → CV = 0.5",
                        "Low → CV = 0.25"
                      ),
                       placement = "right"), 
                    choices = c("High", "Medium", "Low"), # categories
                    # selected value with gaurds for NULL
                    selected = pv$u_cat)
                )
              },
              # for offsets include a slider for offset weights for multipliers
            if(has_weight) {
              sliderInput(
                ns(paste0("CR_w.", id)), # input nname
                label = tooltip(
                  trigger = list("Offset Capacity",
                                 tool_tip_icon),
                  makeTip(
                    "Available capacity for this offset to provide additional
                    compensation.",
                    "Represents how much more can be used if higher compensation
                    is needed.",
                    "",
                    "Weight = 1 → full capacity available (can fully scale up)",
                    "Weight = 0 → no remaining capacity (cannot provide more)",
                    "",
                    "Intermediate values reflect partial availability"
                  ),
                  placement = "right"),
                min = 0, 
                max = 1,
                # initial value with guards for NULLs
                value = pv$CR_w, 
                step = 0.05)
            }
          ) # cardbody
          ) # card
        }) # lapply
      ) # tagList
    # close renderUI
    }) |> bindEvent( # only react to change in n() of err_type
      n(),           # change in number of projects
      lapply(        # change in error type
        grep(paste0("^err_", prefix), names(input), value = TRUE),
        function(id) input[[id]]
      ),
      lapply(        # change in project name
        grep(paste0("^name_apply_", prefix), names(input), value = TRUE),
        function(id) input[[id]]
      ),
      ignoreInit = FALSE
    )
    
    # ---- RETURN --------------------------------------------------------------
    # Output parameters 
    # outputs a parameter list updated with user inputs
    params <- reactive({
      
      # loop over projects
      lapply(par_vals(), function(pv) {
        
        # output list of parameters
        out <- list(
          par_data = pv
        )
        
        # if offset - output capacity weight
        if (has_weight) {
          out$CR_weight <- pv$CR_w
        }
        
        out
      })
      
    })
    
    # output parameters
    return(params)
    
  }) # close module
} # close function

#-------------------------------------------------------------------------------
# Module to create distribution plot output and UI Display
# UI plots distributions of impact and offset projects with plot options:
# - plot type - frequency polygon or histogram
# - display - combine or by project
# - slow legend
# 
# UI Function
distributionPlotUI <- function(id, title, color) {
  ns <- NS(id) # server
  
  # creat card
  card(
    # border colour of card
    class = paste0("border-", color),
    
    # card header - included popover for plot options
    card_header(
      
      # card title
      title, 
      
      # colour of card header
      class = paste0("bg-", color, " text-white d-flex justify-content-between"),
      
      # create popover
      popover(
        bs_icon("gear-fill"),   # toggle icon - gear
        title = "Plot Options", # title in popover
        placement = "auto",    # location to icon
        options = list(trigger = "click",          # trigger option - click icon
                       customClass = "narrow-pop"),# name narrow
        class = "narrow-pop",
        
        # set layout as grid with spacing - as control panel
        div(class = "d-grid gap-3",
            
            # ---- VIEW ----
            # User choice of display - total or individual
            div(
              tags$div(class = "fw-bold mb-1", "View"), # section name
              radioButtons(                             # radio button 
                ns("view"),                             # input ID
                NULL,                                   # no label
                # displayed names - include icon
                choiceNames = list(
                  tagList(icon("layer-group"), " Combined"),
                  tagList(icon("object-ungroup"), " By Project")
                ),
                # choice values - for server
                choiceValues = c("Combined", "By Project")
              )
            ),
            
            # ---- TYPE ----
            # user choice of plot type - freq. polygon or histogram
            div(
              tags$div(class = "fw-bold mb-1", "Plot Type"), # section name
              radioButtons(                                  # radio button
                ns("plot"),                                  # input ID
                NULL,                                        # no label
                # displayed names - include icon
                choiceNames = list(
                  tagList(icon("chart-line"), " Line"),
                  tagList(icon("chart-bar"), " Histogram")
                ),
                # choice values - for server
                choiceValues = c("Line", "Histogram")
              )
            ),
            
            # ---- LEGEND ----
            # user choice to display legend - shinyWidget::prettySwitch
            div(
              tags$div(class = "fw-bold mb-1", "Show Legend"),
              prettySwitch(ns("legend"), label = NULL, value = TRUE)
            )
        )
      ) # close popover
    ),
    
    # display plot
    plotOutput(ns("plot"), height = "300px")
  )
}

# Server function
distributionPlotServer <- function(
    id,           # modUle ID
    dist_data,    # distriubtion data
    project_names,# vector of project names
    n_sim,        # number of simulations
    p_val = NULL  # equivalency threshold (1-risk tolerance
) {
  
  moduleServer(id, function(input, output, session) {
    
    # Plot output
    output$plot <- renderPlot({
      
      # ---- DATA ----
      # distribution data -  from dist_f()
      d <- dist_data()
      
      # names of projects - User input
      names <- project_names()
      
      # if plotting combine totals
      if (input$view == "Combined") {
        df <- data.table(
           x = apply(d, 1, sum),
           type = "Combined"
        )
      # if plotting by project    
      } else {
        df <- data.table(
          x = as.vector(d),
            type = rep(names, each = nrow(d))
        )
      }
       
      # ---- PLOT TYPE ----
      # frequency polygon or histogram
      # if frequency polygon
      if (input$plot == "Line") {
        p <- ggplot(df) +
          geom_freqpoly(aes(df$x, colour = type), bins = 40) 
        # set colour to black if combined
        if (input$view == "Combined") {
          p <- p + scale_colour_manual(values = c("Combined" = "black"))
        }
      # if histogram
      } else {
        p <- ggplot(df) +
          geom_histogram(
            aes(df$x, fill = type),
            bins = 40,
            alpha = 0.3,
            position = "identity"
          )        
        # set colour to black if combined
        if (input$view == "Combined") { 
            p <- p + scale_fill_manual(values = c("Combined" = "black"))
        }
      }
      
      # ---- LEGEND ----
      # Legend if disired
      p <- p + theme_me +
        theme(
          legend.position =
            if (isTRUE(input$legend)) "inside" else "none",
          legend.position.inside = c(0.99, 0.99),
          legend.justification = c(1, 1)
        )
      
      # ---- LABELS ----
      p + labs(
        x = "Value",
        y = "Frequency",
        fill = NULL,
        colour = NULL
      )
    })
  })
}

#-------------------------------------------------------------------------------
# USER INTERFACE
#-------------------------------------------------------------------------------

# initialize sidebar layout
ui <- page_sidebar( 
  # Bootstrap version
  theme = bs_theme(version = 5,
                   bootswatch = "cosmo") |> 
    
    bs_add_rules("
      .popover.narrow-pop {
        max-width: 150px !important;
      }
    
      .popover.narrow-pop .popover-body {
        max-width: 150x !important;
        width: 150px !important;
      }
    "),
  
  # App title
  title = "Offset Uncertainty Multiplier Calculator",
  
  # ---- SIDEBAR ---------------------------------------------------------------
  sidebar = sidebar(
    
    width = 300, # set default sibdebar width
    
    # Go Button - execute CR calculation
    input_task_button("go", "Run", icon("play"), type = "default"),
    
    # User input - accordion panel so can be minimized
    accordion(
      open = TRUE,
      # - Simulation settings - 
      accordion_panel("Simulation Settings", # title
                      # replicates
                      numericInput("n", 
                                   label = tooltip(
                                     trigger = list("N",
                                                    tool_tip_icon),
                                     "Number of replicates for Monte Carlo 
                                     simulations",
                                     placement = "right"),
                                   value = 10000, step = 1000),
                      # Risk tolerance 
                      sliderInput("p", 
                                   label = tooltip(
                                     trigger = list("Risk Tolerance Threshold",
                                                    tool_tip_icon),
                                     "Probability that the offset does not fully 
                                     compensate for the impact. Lower values 
                                     indicate greater certainty. Recommended ≤ 
                                     0.2",
                                     placement = "right"),
                                   value = 0.2, min = 0, max = 0.5, step = 0.05)
      ),
      # - Impact projects - 
      accordion_panel("Impacts",
                      numericInput("n.i", 
                                   label = tooltip(
                                     trigger = list("Number of Impacts",
                                                    tool_tip_icon),
                                     "Specifies how many separate impact 
                                     projects are included. Each impact is 
                                     assigned its own value and uncertainty, 
                                     and the total impact is calculated by 
                                     combining all projects.",
                                     placement = "right"),
                                   value = 1,
                                   min = 1,
                                   step = 1),
                      projectModule_UI("impact") # dynamic UI - see module
      ),
      # - Offset projects - 
      accordion_panel("Offsets",
                      numericInput("n.o", 
                                   label = tooltip(
                                     trigger = list("Number of Offsets",
                                                    tool_tip_icon),
                                     "Specifies how many separate offset 
                                     projects are included. Each offset is 
                                     assigned its own value and uncertainty, 
                                     and the total offset is calculated by 
                                     combining all projects.",
                                     placement = "right"),
                                   value = 1,
                                   min = 1,
                                   step = 1),
                      projectModule_UI("offset") # dynamic UI - see module
                      
      ) # close panels
    )  # close accordion
  ), # close sidebar
  
  # ---- MAIN PANEL ------------------------------------------------------------
  
  # Offset Multiplier VALUE
  # Value box output for CR - one for each offset project
  uiOutput("CR.text"), # dynamic UI - see server
  
  # DISTRIBUTION PLOTS
  # - call distribution plot Module - see above
  
  # Column layout
  layout_column_wrap(
    
    ## - IMPACT PLOT - 
    distributionPlotUI("impact_plot", "Impact Distribution", "danger"),
    
    ## - OFFSET PLOT -
    distributionPlotUI("offset_plot", "Offset Distribution", "success")
    
  ), # Close layout_column
  
  ## - MULTIPLIER PLOT - 
  # Create similar card to distributionPlotUI Module - less plot options
  card( # create card
    
    # border colour
    class = "border-primary", 
    
    # card header
    card_header(
      "Offset Distribution", # title
      
      # header colour - text colour - flexbox - justify
      class = "bg-primary text-white d-flex justify-content-between",
      
      # create popover with plot control
      popover(
        bs_icon("gear-fill"),   # icon
        title = "Plot Options", # popover title
        placement = "left",     # place to right of icon
        options = list(trigger = "click",          # trigger option - click icon
                       customClass = "narrow-pop"),# name narrow
        
        # grid layout with spacing
        div(class = "d-grid gap-3",
            
          # ---- TYPE ----
          # user choice of plot type - freq. polygon or histogram
          div(
            tags$div(class = "fw-bold mb-1", "View"), # secont title
            radioButtons(                             # radio button
              "M_plot",                               # input ID
              label = NULL,                           # no label
              # choice display name with icons
              choiceNames = list(
                tagList(icon("chart-line"), " Line"),
                tagList(icon("chart-bar"), " Histogram")
              ),
              # coice values
              choiceValues = c("Line", "Histogram"),
              inline = FALSE
            ) # close radio button
          ) # close div - section
        ) # close div - layout
      ) # close popover
    ),# close card header
    
    # ---- PLOT ----
    plotOutput("plot.M", height = "300px")
  ) # close card
) # close page_sider - UI

#-------------------------------------------------------------------------------
# Server
#-------------------------------------------------------------------------------

## server.R
server <- function(input, output, session){
  
  #-----------------------------------------------------------------------------
  # VALIDATION
  #-----------------------------------------------------------------------------
  # validation rules using shingvalidate package
  
  iv <- InputValidator$new()      # set up validater
  
  # add rules
  iv$add_rule("n", sv_required())   # n sim is required
  iv$add_rule("n", sv_numeric())    # n sim must be numeric
  iv$add_rule("n", sv_gte(1))       # n sim must be >=1
  
  iv$add_rule("n.i", sv_required()) # n impacts is required
  iv$add_rule("n.i", sv_numeric())  # n impacts must be numeric
  iv$add_rule("n.i", sv_gte(1))     # n impacts must be >=1
  
  iv$add_rule("n.o", sv_required()) # n offset is required
  iv$add_rule("n.o", sv_numeric())  # n offset must be numeric
  iv$add_rule("n.o", sv_gte(1))     # n offset must be >=1
  
  # rules for dynamic impact UI inputs - see functions
  observe({
    
    # require 
    req(input$n.i, input$n.o)
    
    # call name space of module
    ns_impact <- NS("impact")
    
    # run validation rule function
    add_validation_rules(
      iv = iv,            # validater
      input = input,      # UI input
      prefix = "i",       # indicator for impacts
      ns_fun = ns_impact, # name space
      n = input$n.i       # number of impacts
    )
    
  })
  
  # rules for dynamic offset UI inputs - see functions
  observe({
    
    # require
    req(input$n.i, input$n.o)
    
    # call name space for module
    ns_offset <- NS("offset")
    
    # run validation rule function
    add_validation_rules(
      iv = iv,            # validater
      input = input,      # UI input
      prefix = "o",       # indicator for impacts
      ns_fun = ns_offset, # name space
      n = input$n.o       # number of impacts
    )
    
  })
  
  # enable rules
  iv$enable()
  
  #-----------------------------------------------------------------------------
  # DYNAMIC UI
  #-----------------------------------------------------------------------------
  # Call modules to produce dynamic UI as the number of Impact for offset
  # projects changes
  
  # IMPACTS 
  impact_params <- projectModule_server(
    "impact",
    n = reactive(input$n.i),
    label = "Impact",
    prefix = "i",
    has_weight = FALSE
  )
  
  offset_params <- projectModule_server(
    "offset",
    n = reactive(input$n.o),
    label = "Offset",
    prefix = "o",
    has_weight = TRUE
  )
  
  #-----------------------------------------------------------------------------
  # DATA
  #-----------------------------------------------------------------------------
  
  # Create simulation data for impacts
  
  impact_projects <- eventReactive(input$go, {
    req(iv$is_valid())
    impact_params()
  })
  
  # Create simulation data
  
  offset_projects <- eventReactive(input$go, {
    req(iv$is_valid())
    offset_params()
  })
  
  #-----------------------------------------------------------------------------
  # SIMULATIONS
  #-----------------------------------------------------------------------------
  
  # generate distributions - call dist.f() function 
  # outputs dataframes of simulated impact and offset values in a list
  dist.data <- eventReactive(input$go, {
    req(iv$is_valid(), impact_projects(), offset_projects())

    # generate simulated data
    dist.f(
      impact_data = impact_projects(), # mean impact value
      offset_data = offset_projects(), # mean offset value
      n.sim = input$n                  # number of draws
    )
  })
  
  # Calc offset multipler
  CR_sim <- eventReactive(input$go, {
    req(iv$is_valid(), dist.data())
    CRu.f( 
      d.i = dist.data()$impacts,  # distributions of n.i impacts
      d.o = dist.data()$offsets,  # distributions of n.o offsets
      p = 1-input$p,               # risk tolerance level
      CR_weight = sapply(offset_projects(), function(x) x$CR_w) # offset weights
    ) 
  })
  
  #-----------------------------------------------------------------------------
  # OUTPUTS
  #-----------------------------------------------------------------------------
  
  ## - IMPACTS -
  distributionPlotServer(
    "impact_plot",
    dist_data = reactive(dist.data()$impacts),
    project_names = reactive(
      sapply(impact_projects(), function(x) x$par_data$name)
    ),
    n_sim = reactive(input$n)
  )
  
  ## - OFFSET-
  distributionPlotServer(
    "offset_plot",
    dist_data = reactive(dist.data()$offsets),
    project_names = reactive(
      sapply(offset_projects(), function(x) x$par_data$name)
    ),
    n_sim = reactive(input$n)
  )
  
  ## ---- MULTIPLIERS ----
  # Multiple data to display in UI
  CR_display <- eventReactive(input$go, {
    list(
      cr = CR_sim(),        # M distribution
      n_offsets = input$n.o,# number of offset projects
      p = 1 - input$p       # risk tolerance convert to ET
    )
  })
  
  # output plot distributions of M - Compensation Ratio
  output$plot.M <- renderPlot({ 
    
    # reactive data
    res <- CR_display()
    
    # plot data
    df <- data.table(
      M = res$cr$M
    )
    
    # plot limit - cut off ends of large dist'n
    limits <- quantile(df$M, c(0,0.99))
    if(limits[1] < 0) limits[1] <- 0
    
    # - PLOT TYPE -
    # user input option for plotting
    # if frequency polygon
    if (input$M_plot == "Line") {
      p <- ggplot(df) +
        geom_freqpoly(aes(M),bins = 40)
    # if histogram
    } else {
      p <- ggplot(df) +
        geom_histogram(aes(M),
                       position = "identity",
                       bins = 40, alpha = 0.3)
    } 
    
    
    p + 
      scale_x_continuous(limits = limits)+                   # apply limits
      labs(x = "Offset Multiplier",y="Frequency")+           # labels
      geom_vline(xintercept=round(quantile(df$M, res$p),2),  # risk tolerance
                 linetype = 2) +                             # reference line
      theme_me                                               # my theme
    
  })
  
  # Multipler value output as text
  output$CR.text <- renderUI({
    
    # create emtry value box to display when App opens
    if (input$go == 0) {
      return(        
        value_box(
          title = "Offset Multiplier",
          value = "—",
          theme = value_box_theme(bg = "#f79e25")
      ))
    }
    
    # data
    res <- CR_display()
    
    # extract values
    cr_vals <- res$cr[["CR - Adjusted"]] # multiplier dist'n
    n_offsets <- res$n_offsets           # number of offsets
    
    # create value box for each offset projects & display M[r] value
    layout_columns(
      !!!lapply(seq_len(n_offsets), function(i) {
        value_box(
          title = HTML(paste("Offset Multiplier <br> Offset", i)),
          value = round(cr_vals[i], 3),
          theme = value_box_theme(bg = "#f79e25")
        )
      })
    )
  })
  
} # close server

#-------------------------------------------------------------------------------
# RUN
#-------------------------------------------------------------------------------

### Run Application
shinyApp(ui, server)

################################################################################