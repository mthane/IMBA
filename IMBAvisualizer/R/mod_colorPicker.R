
mod_colorPicker_ui <- function(id) {
  ns <- NS(id)
  fluidRow(
    column(2,
           
           uiOutput(ns("groupNamePickers"))
    ),
    column(2,
           
           uiOutput(ns("colorPickers"))
    ),
    column(2,
           uiOutput(ns("groupSelectionsCheckboxes"))
    )
  )
}

mod_colorPicker_server <- function(id,summaryData){
  moduleServer(
    id,
    function(input, output, session) {
      
      ### colorPicker ###
      groups <- reactive(unique(summaryData()$group_condition))
      
      pal <- reactive(palette.colors(n=length(groups()),recycle = T))
      
      cols <- reactive({
        ns <- session$ns
        lapply(1:length(groups()), function(i) {
          spectrumInput(
            inputId = ns(paste("col", i, sep="_")),
            label = groups()[i],
            selected = pal()[i]
          )
        })
      })
      
      groupNames <- reactive({
        ns <- session$ns
        lapply(groups(), function(gname) {
          textInput(
            inputId = ns(paste("gname", gname, sep="_")),
            label = gname,
            value = gname
          )
        })
      })
      
      groupSelections <- reactive({
        ns <- session$ns
        lapply(groups(), function(gname) {
          checkboxInput(
            inputId = ns(paste("sel", gname, sep="_")),
            label = gname,
            value = T
          )
        })
      })
      
      groupColorValues <- reactive({
        lapply(1:length(groups()), function(i) {
          input[[paste("col", i, sep="_")]]
        })
      })
      
      groupNameValues <- reactive({
        lapply(groups(), function(gname) {
          input[[paste("gname", gname, sep="_")]]
        })
      })
      
      groupSelectionValues <- reactive({
        lapply(groups(), function(gname) {
          input[[paste("sel", gname, sep="_")]]
        })
        
      })
      
      
      output$colorPickers <- renderUI({cols()})
      
      output$groupNamePickers <- renderUI({groupNames()})
      
      output$groupSelectionsCheckboxes <- renderUI({groupSelections()})
      
      
      
      output$res <- renderPrint({groupNameValues()})
      
      result <- reactive({
        idx = unlist(groupSelectionValues())
        list(
          groupPrevValues = unlist(groups()[idx]),
          groupNameValues = unlist(groupNameValues()[idx]),
          groupColorValues = unlist(groupColorValues()[idx])
        )
      })
      result
      
      
    }
  )
}