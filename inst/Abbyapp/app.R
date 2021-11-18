
##### AbbyApp  by slphyx

app_version <- "V0.1"


library(shiny)
library(shinythemes)
library(stringr)
library(ggplot2)
library(DT)

# min/max columns
col.order <- c(3,4)
col.family <- c(7,8)
col.genus <- c(11,12)
col.species <- c(15,16)


## get the list of markers
Markers.list <- function(data){
    len <- nrow(data)
    data[2:len,1]
}

## get row data from input maker
Marker.data <- function(marker, data, OutTable = F){
    tmp <- data[data$Genetic.marker==marker,]
    if(OutTable == F)
        return(as.numeric(tmp[1,2:ncol(tmp)]))
    else
        return(tmp)
    
}

## get the min/max range
MinMax <- function(level, marker, data){
    tmp <- switch(level, 
                  Order = Marker.data(marker, data)[col.order],
                  Family = Marker.data(marker, data)[col.family],
                  Genus = Marker.data(marker, data)[col.genus],
                  Species = Marker.data(marker, data)[col.species])
    
    return(c(Min = tmp[1], Max = tmp[2]))
}

## check if x is between left and right
Between <- function(x, left, right){
    x >= left & x <= right
} 


## show available levels ,i.e., MinMax not return NA
Level.Available <- function(data, marker, level.out = FALSE){
    df <- data.frame(level=NULL,Min=NULL,Max=NULL)
    for(or in c("Order","Family","Genus","Species")){
        mm <- MinMax(data = data,level = or,marker = marker)
        df <- rbind(df,data.frame(level=or,Min=mm[1],Max=mm[2]))
    }
    row.names(df) <- NULL
    
    if(level.out == TRUE)
        return(df[!is.na(df$Min),])
    
    return(df[!is.na(df$Min),]$level)
}


Convert.group <- function(group){
    if(group == "Nematode (Ascaridida and Spirurida)")
        return("NAS")
    if(group == "Nematode (Strongylida)")
        return("NS")
    if(group == "Nematode (Trichocephalida)")
        return("NT")
    if(group == "Trematode")
        return("TR")
    if(group == "Cestode")
        return("CE")
}

### load data
Load.data <- function(group){
    data <- switch(group,
                   NAS = read.csv('./data/nematode1.csv'),
                   NS = read.csv('./data/nematode2.csv'),
                   NT = read.csv('./data/nematode3.csv'),
                   TR = read.csv('./data/trematode.csv'),
                   CE = read.csv('./data/cestode.csv'))
    return(data)
}


## check if distance between A and B val from helminth group using marker are the same at level
# group NAS = Nematode (Ascaridida and Spirurida)
# group NS = Nematode (Strongylida)
# group NT = Nematode (Trichocephalida)
# group TR = Trematode
# group CE = Cestode 
# return TRUE if val is between min/max
Abby <- function(group, level, marker, val){
    
    data <- Load.data(group)
    
    minmax <- MinMax(data = data, level = level, marker = marker)
    q<-Between(val, minmax[1], minmax[2])
    names(q) <- NULL
    return(q)
}

ShowRanges <- function(ranges, distance=NULL, group = NULL,marker=NULL){
    len <- nrow(ranges)
    # print(ranges)
    grp <- ggplot(data = ranges, aes(x=Min, xend=Max, y=level, colour=level)) +
        geom_segment(aes(xend=Max, yend=level), linetype=1, size=2) +
        geom_label(aes(label=str_wrap(level,12), x=(Min + Max)/2), size=5) +
        scale_colour_brewer(palette = "Set1") +
        xlab("Genetic Distance")+ ylab("Level")+
        theme_bw() + 
        theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), 
              aspect.ratio = .2,
              legend.position="none",
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank()) +
        ggtitle(paste(group,":" ,marker))
    
    if(!is.null(distance)){
        grp <- grp + geom_vline(xintercept = distance,linetype="longdash", size=1, col="gray")
    }
    
    return(grp)
    
}



ui.bak <- fluidPage(
    titlePanel("Abby App"),
    sidebarLayout(
        sidebarPanel(
            selectInput(inputId = 'group', label = "Please choose your helminth group",
                    choices = c("Nematode (Ascaridida and Spirurida)",
                          "Nematode (Strongylida)",
                          "Nematode (Trichocephalida)",
                          "Trematode", "Cestode")
                        ),
            
            uiOutput('geninput.marker'),
            uiOutput('geninput.level'),
            actionButton('submit',"Submit")
            
        ),
        mainPanel(
            htmlOutput('Abby.ans'),
            plotOutput('distplot')
        )
    )
)

ui <- fluidPage(
    title = "AbbyApp",
    theme = shinytheme("readable"),
    navbarPage("AbbyApp",
        tabPanel("About",
            includeMarkdown("./www/text/about.md") 
        ),
        tabPanel("Application",
            fluidRow(
                column(4,
                    numericInput(inputId = 'distance', label = "Genetic distance between taxa of interest (p-distance:0-1):
", min = 0,max = 1, step = 0.001, value = 0)       
                ),
                column(4,
                    selectInput(inputId = 'group', label = "Helminth group :",
                       choices = c("Nematode (Ascaridida and Spirurida)",
                                   "Nematode (Strongylida)",
                                   "Nematode (Trichocephalida)",
                                   "Trematode", 
                                   "Cestode")
                    )    
                ),
                column(4,
                    uiOutput('geninput.marker')
                ),
            ),
            actionButton('submit',"Go!"),
            hr(),
            fluidRow(
                tabsetPanel(
                    tabPanel("Graph",
                        column(8,
                            plotOutput('distplot')
                        ),
                        column(4,
                            htmlOutput('Abby.ans')
                        )
                        
                    ),
                    tabPanel("Data",
                        DTOutput('tbl')
                    )
                )
            )
            
        )
        
    )
    
)


server <- function(input, output) {
    
    values <- reactiveValues()

    # values$group <- Convert.group(input$group)
    # values$data <- Load.data(group = values$group)
    
    output$geninput.marker <- renderUI({
        
        values$group <- Convert.group(input$group)
        values$data <- Load.data(group = values$group)
        values$marker.list <- Markers.list(data = values$data)
        
        
        div(
            selectInput('marker',label = "Genetic marker :", choices = values$marker.list, selected = NULL)
        )
    })
    
    
    # output$geninput.level <- renderUI({
    #     
    #     values$marker <- input$marker
    #     values$level.list <- Level.Available(data =values$data,marker = values$marker)
    #     
    #     #print(values$level.list)
    #     
    #     div(
    #         radioButtons(inputId = 'level', label = 'Level :', choices = values$level.list),
    #         
    #         numericInput(inputId = 'distance', label = "Genetic Distance (0 - 1)", min = 0,max = 1,
    #                      step = 0.001, value = 0)
    #         
    #     )
    # })
    # 
    
    observeEvent(input$submit,{
        
        # values$level <- input$level
        
        values$distance <- input$distance
        
        # values$abbyQ <- Abby(group = values$group, level = values$level,marker = values$marker, val = values$distance)
        values$ranges <- Level.Available(data = values$data, marker = input$marker, level.out = TRUE)

        
        # output$Abby.ans <- renderText({
        #     if(values$abbyQ == TRUE){
        #         paste(h2("They are different in", values$level,"level."))
        #     }else{
        #         paste(h2("They are the same in", values$level,"level."))
        #     }
        #     
        # })
        
        output$distplot <- renderPlot({
            ShowRanges(ranges = values$ranges,distance = values$distance, group = values$group, marker = values$marker)
            
        })
        
        
        
    })
    
    output$tbl <- renderDT({ Marker.data(data = values$data, marker = input$marker, OutTable = TRUE)})
    
    
}

shinyApp(ui = ui, server = server)

