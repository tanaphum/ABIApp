
##### ABIyApp  by slphyx

app_version <- "V0.4"

if(!require("ggtree", quietly = TRUE)){
  if (!require("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
  }
  BiocManager::install("ggtree")
}

library(shiny)
library(shinythemes)
library(shinyjs)
library(shinydashboard)
library(stringr)
library(ggplot2)
library(htmltools)
library(markdown)
library(shinyWidgets)
library(bsplus)
library(seqinr)
library(ape)
library(ggtree)

flag.choose <- function(language) {
  switch(language,
         EN = "./img/gb.png",
         TH = "./img/th.png"
  )
}

#changeLanguage EN to TH include Document path
changeLanguage <-function(language) {
  switch(language,
         EN = "TH",
         TH = "EN")
}

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
  colnames(tmp) <- c( "Genetic.marker", "Order Suborder \n Mean", "SD" ,"Min" ,"Max",
                      "Family \n Mean", "SD" ,"Min" ,"Max",
                      "Genus \n Mean" , "SD" ,"Min" ,"Max",
                      "Species \n Mean" ,"SD" ,"Min" ,"Max")
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
  if(x < 0) return(F)
  if(is.na(x >= left & x <= right)) return(F)
  else if(length(x >= left & x <= right) ==0) return(F)
  else return(x >= left & x <= right)
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
  if(group == "Trematode (Plagiorchiida)")
    return("TR")
  if(group == "Trematode (Diplostomida)")
    return("TRD")
  if(group == "Cestode")
    return("CE")
}

### load data
Load.data <- function(group){
  data <- switch(group,
                 NT = read.csv('./data/nematode1.csv'),
                 NAS = read.csv('./data/nematode2.csv'),
                 NS = read.csv('./data/nematode3.csv'),
                 TR = read.csv('./data/trematode1.csv'),
                 TRD = read.csv('./data/trematode2.csv'),
                 CE = read.csv('./data/cestode.csv'))
  return(data)
}


## check if distance between A and B val from helminth group using marker are the same at level
# Between(x, left, right)
Check.distance.between <- function(group,distance,marker){

  data <- Load.data(group)
  markerNumber <-as.numeric(rownames(data[data$Genetic.marker == marker,]))
  #min max Species
  if(Between(distance,as.numeric(data$X.10[markerNumber]),as.numeric(data$X.11[markerNumber]))){
    #Species
    return("Species")

    #max Species & min Genus
  }else if(Between(distance,as.numeric(data$X.11[markerNumber]),as.numeric(data$X.7[markerNumber]))){
    #Between Genus-Species
    return("Genus-Species")

    #min max Genus
  }else if(Between(distance,as.numeric(data$X.7[markerNumber]),as.numeric(data$X.8[markerNumber]))){
    #Genus
    return("Genus")

    #max Genus & min Family
  }else if(Between(distance,as.numeric(data$X.8[markerNumber]),as.numeric(data$X.4[markerNumber]))){
    #Between Family-Genus
    return("Family-Genus")

    #min max Family
  }else if(Between(distance,as.numeric(data$X.4[markerNumber]),as.numeric(data$X.5[markerNumber]))){
    #
    return("Family")

    #max Family & min Order
  }else if(Between(distance,as.numeric(data$X.5[markerNumber]),as.numeric(data$X.1[markerNumber]))){
    #Between Order-Family
    return("Order-Family")
  }else if(Between(distance,as.numeric(data$X.1[markerNumber]),as.numeric(data$X.2[markerNumber]))){
    #Order
    return("Order")
  }else{
    #out of bounds
    return("Out")
  }

}

ShowRanges <- function(ranges, distance=NULL, group = NULL,marker=NULL){
  len <- nrow(ranges)

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

ui <- dashboardPage(
  # Header
  dashboardHeader(
    title = "ABIApp",
    tags$li(
      actionLink(
        inputId = "goto_Contact",
        label = HTML("<i class='fa fa-envelope'></i> Contact Us"), # Add icon using HTML
        class = "btn btn-default"
      ),
      class = "dropdown"
    )
  ),

  # Sidebar
  dashboardSidebar(
    sidebarMenu(
      id = "sidebarMenu",

      # About Section with Dropdown
      menuItem("About", icon = icon("info-circle"), startExpanded = F,
               menuSubItem("General Information", tabName = "general_information"),
               menuSubItem("How ABI App Works", tabName = "how_abi_app_works"),
               menuSubItem("Benefits", tabName = "benefits"),
               menuSubItem("Assumptions", tabName = "assumptions"),
               menuSubItem("References", tabName = "references")
      ),
      menuItem("Quick guidelines", tabName = "quick_guidelines_tab", icon = icon("book")),
      menuItem("Application",tabName = "application_tab", icon = icon("flask"), startExpanded = F,
               menuSubItem("Fasta file", tabName = "fasta_file"),
               menuSubItem("Genetic distance", tabName = "genetic_distance")
      )
    )
  ),

  # Body
  dashboardBody(
    useShinyjs(),
    includeCSS("./www/css/style.css"),
    tags$link(
      href = "https://fonts.googleapis.com",
      rel = "preconnect"
    ),
    tags$link(
      href = "https://fonts.gstatic.com",
      rel = "preconnect",
      crossorigin = ""
    ),
    tags$link(
      href = "https://fonts.googleapis.com/css2?family=Noto+Sans+JP:wght@100..900&display=swap",
      rel = "stylesheet"
    ),
    fluidRow(id = "tabItems_fluidRow",
             tabItems(
               # Intro Tab
               tabItem(tabName = "general_information",
                       h2("General Information"),
                       div(
                         tags$img(id = "ABIApp", src = "./img/ABI_app.png", alt = "ABI App"),
                         tags$hr(),
                         includeMarkdown("./www/text/General Information.md")
                       )
               ),

               # How ABI App Works Tab
               tabItem(tabName = "how_abi_app_works",
                       h2("How ABI App Works"),
                       includeMarkdown("./www/text/app work.md")
               ),

               # Benefits Tab
               tabItem(tabName = "benefits",
                       h2("Benefits"),
                       includeMarkdown("./www/text/Benefits.md")
               ),

               # Assumptions Tab
               tabItem(tabName = "assumptions",
                       h2("Assumptions"),
                       includeMarkdown("./www/text/Assumptions.md")
               ),

               # References Tab
               tabItem(tabName = "references",
                       h2("References"),
                       includeMarkdown("./www/text/References.md")
               ),

               # Quick Guidelines Tab
               tabItem(
                 tabName = "quick_guidelines_tab",
                 fluidRow(
                   column(10, tags$h2("Quick guidelines")),
                   column(2, uiOutput("button.flag"))
                 ),
                 tabsetPanel(
                   tabPanel("Before you start", uiOutput("beforeStart")),
                   tabPanel("Input data", uiOutput("inputData")),
                   tabPanel("Output visualization",
                            tabsetPanel(
                              tabPanel("Fasta file", uiOutput("outputVisualization_fasta")),
                              tabPanel("Genetic distance", uiOutput("outputVisualization"))
                            )
                   ),
                   tabPanel("Suggested primers",
                            uiOutput("suggestedPrimers"),
                            tags$a("Suggested PCR primers for ABIApp", target = "_blank", href = "pdf/Suggested_PCR_primers_for_ABIapp.pdf")
                   )
                 )
               ),
               # Fasta File Tab
               tabItem(
                 tabName = "fasta_file",
                 h2("Fasta file"),
                 tags$br(),
                 fluidRow(
                 column(9, fileInput("fastaFile", "Upload File:", width = "500px", accept = ".fasta")),
                 column(3, uiOutput("distanceText")),
                 column(6, selectizeInput("fastaSelect1", "Sequence ID of queried sequence (ID name should match name in fasta file):", choices = c(), width = "100%")),
                 column(6, selectizeInput("fastaSelect2", "Sequence ID of sequence to be compared with (ID name should match name in fasta file):", choices = c(), width = "100%")),

                 )
               ),

               # Genetic Distance Tab
               tabItem(
                 tabName = "genetic_distance",
                 h2("Genetic distance")
               )
             ),
             # Common Elements for Both Tabs
             div(
               id = "application_tab",
               fluidRow(
                 column(6,
                        div(id="genetic_distance_input",
                            column(12, numericInput("distance", "Genetic distance between taxa of interest (p-distance: 0-1):", min = 0, max = 1, step = 0.001, width = "100%", value = 0))
                            ),
                   column(12, selectInput("group", "Please choose your helminth group",
                                         choices = c("Nematode (Trichocephalida)",
                                                     "Nematode (Ascaridida and Spirurida)",
                                                     "Nematode (Strongylida)",
                                                     "Trematode (Plagiorchiida)",
                                                     "Trematode (Diplostomida)",
                                                     "Cestode"),
                                         width = "100%")),
                   column(12, uiOutput("geninput.marker")),
                   column(12, actionButton("submit", "Submit!")),
                   ),
                 column(6,
                        div(
                   # Output Table and Tabs Output
                   column(12,tableOutput("tbl")),
                    )

                 ),
               ),
             ),
             div(id="output_tabs",
                 uiOutput("tabsOutput")
             )
    ),
    fluidRow(id="Contact_us",
      column(12,
             # Title
             h1("Contact us"),
             br(),
             h4(
               # Contact Information
               tags$p(
                 h4("Dr. Abigail Hui En Chan"),
                 "2nd floor Santasiri Sornmani Building,Department of Helminthology,",
                 br(),
                 "Faculty of Tropical Medicine,  Mahidol University,",
                 br(),
                 "420/6  Ratchawithi Rd., Ratchathewi, Bangkok 10400. Thailand",
                 br(),
                 br(),
                 "For any inquiry about the dashboard, please contact",
                 br(),
                 "Email: ", tags$a(href="mailto:abigail.cha@mahidol.edu", "abigail.cha@mahidol.edu")
               )
             ),
             br(),
             column(6,
                    tags$img(
                      src = "img/mahidol_logo.png", width="70%" ,height="70%"
                    )
             ),
             column(6,
                    tags$img(
                       src = "img/MORU_FS_RGB.png", width="70%" ,height="70%"
                    )
             ),
             column(12,
             hr(),
             ),
             h4("Team"),
             tags$table(
               style = "width: 100%; border-collapse: collapse; border-style: none;",
               tags$tbody(
                 tags$tr(
                   tags$td(tags$b("Abigail Hui En Chan"), style = "text-align: left; padding: 8px;"),
                   tags$td("Department of Helminthology, Faculty of Tropical Medicine", style = "padding: 8px;")
                 ),
                 tags$tr(
                   tags$td(tags$b("Urusa Thaenkham"), style = "padding: 8px;"),
                   tags$td("Department of Helminthology, Faculty of Tropical Medicine", style = "padding: 8px;")
                 ),
                 tags$tr(
                   tags$td(tags$b("Tanaphum Wichaita"), style = "padding: 8px;"),
                   tags$td("Mahidol Oxford Tropical Medicine Research Unit", style = "padding: 8px;")
                 ),
                 tags$tr(
                   tags$td(tags$b("Sompob Saralamba"), style = "padding: 8px;"),
                   tags$td("Mahidol Oxford Tropical Medicine Research Unit", style = "padding: 8px;")
                 )
               )
             )

      )
    )
  )
)

server <- function(input, output,session) {
  useShinyjs()
  values <- reactiveValues()
  values$distance <- NULL
  values$heighttree <- "400px"
  values$displayTable <- F
  values$language <- "EN"

  hide("application_tab")
  hide("genetic_distance_input")
  hide("Contact_us")


  output$geninput.marker <- renderUI({

    values$group <- Convert.group(input$group)
    values$data <- Load.data(group = values$group)
    values$marker.list <- Markers.list(data = values$data)
    output$Abby.ans <-NULL
    values$displayTable <- F

    div(
      selectInput('marker',label = "Genetic marker :", choices = values$marker.list, selected = NULL,width = "100%")
    )
  })


  values$MD_Before <- "./www/text/Before you start.md"
  values$MD_InputData <- "./www/text/Input data.md"
  values$MD_Visualization_fasta <-"./www/text/Output visualization_fasta.md"
  values$MD_Visualization <-"./www/text/Output visualization.md"
  values$MD_SuggestedPrimers <- "./www/text/Suggested primers.md"


  output$beforeStart<- renderUI({
    includeMarkdown(values$MD_Before)
  })
  output$inputData <- renderUI({
    includeMarkdown(values$MD_InputData)
  })
  output$outputVisualization_fasta <- renderUI({
    includeMarkdown(values$MD_Visualization_fasta)
  })
  output$outputVisualization <- renderUI({
    includeMarkdown(values$MD_Visualization)
  })
  output$suggestedPrimers<- renderUI({
    includeMarkdown(values$MD_SuggestedPrimers)
  })

  output$button.flag <- renderUI({
    tags$button(
      id = "language_button",
      class = "btn action-button",
      style='background-color:transparent',
      tags$img(id = "flag",src = flag.choose(values$language),height = "25px")
    )
  })

  ### Change Language
  observeEvent(input$language_button, {
    values$MD_Before <-   switch(values$language,
                                 "EN" = "./www/text/Before you start_TH.md",
                                 "TH" = "./www/text/Before you start.md")
    values$MD_InputData <- switch(values$language,
                                  "EN" = "./www/text/Input data_TH.md",
                                  "TH" = "./www/text/Input data.md")
    values$MD_Visualization <- switch(values$language,
                                      "EN" = "./www/text/Output visualization_TH.md",
                                      "TH" = "./www/text/Output visualization.md")
    values$MD_Visualization_fasta <- switch(values$language,
                                            "EN" = "./www/text/Output visualization_fasta_TH.md",
                                            "TH" = "./www/text/Output visualization_fasta.md")
    values$MD_SuggestedPrimers <- switch(values$language,
                                         "EN" = "./www/text/Suggested primers_TH.md",
                                         "TH" = "./www/text/Suggested primers.md")
    values$language <- changeLanguage(values$language)

  })

  observeEvent(input$goto_Contact,{
    hide("tabItems_fluidRow")
    show("Contact_us")
  })


  observeEvent(input$sidebarMenu, {
    show("tabItems_fluidRow")
    hide("Contact_us")
    if(input$sidebarMenu=="fasta_file" | input$sidebarMenu == "genetic_distance"){
      show("application_tab")
      show("output_tabs")
      if(input$sidebarMenu == "genetic_distance"){
        show("genetic_distance_input")
      }else
      {
        hide("genetic_distance_input")
      }

    }else{
      hide("application_tab")
      hide("output_tabs")
      hide("genetic_distance_input")
    }
  })
  # read Fasta File
  observeEvent(input$fastaFile,{
    fastaFile <- input$fastaFile
    sequences <- read.dna(fastaFile$datapath,
                          format = "fasta")


    # Convert sequences to distance matrix
    dist_matrix <- dist.dna(sequences, model = "raw")

    # Construct neighbor-joining tree
    values$dna_tree <- nj(dist_matrix)

    # Calculate genetic distance between each pair of sequences
    values$genetic_distance <- cophenetic(values$dna_tree)
    if(length(values$genetic_distance[1,]) > 16){
      h <- length(values$genetic_distance[1,])*25
      values$heighttree <- paste0(h,"px")
    }

    colname_genetic <- colnames(values$genetic_distance)
    rowname_genetic <- rownames(values$genetic_distance)
    # update Selectize
    updateSelectizeInput(session,
                         "fastaSelect1",
                         "Sequence ID of queried sequence(ID name should match name in fasta file):",
                         choices=c(colname_genetic),

    )

    # update Selectize
    updateSelectizeInput(session,
                         "fastaSelect2",
                         "Sequence ID of sequence to be compared with(ID name should match name in fasta file):",
                         choices=c(rowname_genetic),
                         selected = tail(rowname_genetic, 1)
    )

  })


  # plot tree
  output$treeplot <- renderPlot({
    # plot(nj_tree)
    taxon1 <- match(input$fastaSelect1,values$dna_tree$tip.label)
    taxon2 <- match(input$fastaSelect2,values$dna_tree$tip.label)
    num_char <- nchar(input$fastaSelect1)
    branch_length <- "none"
    if(input$branch_length) branch_length <- "branch.length"
    ggtree(values$dna_tree,layout=input$treelayout,branch.length=branch_length)+
      geom_tiplab(size=5)+
      geom_highlight(node=taxon1,alpha=0.3,fill="forestgreen",extend = num_char) +
      geom_hilight(node=taxon2,alpha=0.3,fill="forestgreen", extend = num_char) +
      hexpand(.35)+
      geom_taxalink(input$fastaSelect1, input$fastaSelect2, color="orange2", curvature=-.9)
  })

  output$treeOutput <- renderUI({
    req(!is.null(input$fastaFile))
    tagList(
      fluidRow(

        column(3,offset = 1,
               checkboxInput("branch_length","branch length scaling")
        ),
        column(8,
               selectInput("treelayout", "Layout of tree :",
                           c("rectangular" = "rectangular",
                             "slanted" = "slanted",
                             "circular" = "circular"
                           ),
                           selected = "rectangular"
               ),
        ),

      ),
      plotOutput('treeplot',width = "100%",height =values$heighttree)
    )
  })

  observeEvent(input$sidebarMenu,{
    values$displayTable <-F
  })
  output$taxa <- renderText({
    req(input$sidebarMenu == "fasta_file" && values$displayTable)
    paste0("Based on a genetic distance of ", round(values$distance,3),
           " between ",input$fastaSelect1,
           " and " , input$fastaSelect2)
  })

  observeEvent(c(input$submit,input$fastaSelect1,
                 input$fastaSelect2,input$group,input$marker),{
                   if(input$sidebarMenu == "fasta_file" && !is.null(input$fastaFile)){
                     HTML(paste("<h4>Error
                         </h4>"))
                   }
                   req((input$sidebarMenu == "fasta_file" && !is.null(input$fastaFile)) || input$sidebarMenu == "genetic_distance")
                   # values$level <- input$level
                   values$displayTable <-T
                   values$distance <- 0
                   values$distance_fasta_max <- 0
                   #Fasta file
                   if(input$sidebarMenu == "fasta_file"){

                     col_Select <- input$fastaSelect1
                     row_Select <- input$fastaSelect2
                     values$distance <- values$genetic_distance[row_Select,col_Select]
                   }else{
                     # input distance
                     values$distance <- input$distance
                   }



                   # values$abbyQ <- Check.distance.ranges(group = values$group, level = values$level,marker = values$marker, val = values$distance)
                   values$ranges <- Level.Available(data = values$data, marker = input$marker, level.out = TRUE)
                   values$distanceBetween <- Check.distance.between(group=values$group,distance=values$distance,marker=input$marker)

                   #Species Genus Family Order
                   #NAS NS NT TR CE
                   output$AbbyText <- renderText({
                     req(values$displayTable)
                     if(values$distance==0){
                       if(values$group == "NAS" && (input$marker =="18S rRNA" || input$marker =="28S rRNA" || input$marker =="ITS1" || input$marker =="ITS2" || input$marker =="COII" ||input$marker =="12S rRNA")){
                         HTML(paste("<h4><ul><li>Suggest to use mt 16S rRNA gene or mt COII as an alternative genetic marker</li>
                         <li>Although 18S has the smallest gap, but the low sequence variation at the genus-species level is challenging for species delimitation</li>
                         <li>Suggest nematode 16S primer from nematode systematics paper</li><ul>
                         </h4>"))
                       }else if (values$group == "NT" && (input$marker =="18S rRNA" || input$marker =="ITS2")){
                         paste(h3("Suggest to use mt 12S"))
                       }else if (values$group == "TR"&& (input$marker =="18S rRNA" || input$marker =="ITS1" || input$marker =="ITS2")){
                         HTML(paste("<h4><ul><li>Suggest to use mt 16S </li>
                         <li>Although 18S has small gap (same as 16S), but the low sequence variation at the genus-species level is challenging for species delimitation</li></ul>
                         </h4>"))
                       }else if (values$group == "CE"&& input$marker =="12S rRNA"){
                         paste(h3("Recommendation to use another mt genetic marker"))
                       }
                     }else{

                       if(values$distanceBetween == "Genus-Species"){
                         if(values$group == "NAS" || values$group == "NS"){
                           HTML(paste("<h4><ul> <li>Suggest using the mt 16S rRNA or COII gene as an alternative genetic marker</li>
                                <li>	Refer to the suggested PCR primer list for the 16S primer for nematodes</li></ul>
                           </h4>"))
                         }else if (values$group == "NT"){
                           HTML(paste("<h4><ul><li>Suggest using the mt 12S rRNA gene as an alternative genetic marker</li>
                          <li>Refer to the suggested PCR primer list for the 12S primer for nematodes</li></ul>
                          </h4>"))
                         }else if (values$group == "TR"){
                           HTML(paste("<h4><ul><li>Suggest using the mt 16S rRNA gene as an alternative genetic marker</li>
                                <li>Refer to the suggested PCR primer list for the 16S primer for platyhelminths</li>
                            </h4>"))
                         }else if (values$group == "CE"){
                           HTML(paste("<h4><ul><li>Suggest using the mt 12S rRNA gene as an alternative genetic marker</li>
                          <li>Refer to the suggested PCR primer list for the 12S primer for platyhelminths</li><ul></h4>
                        "))
                         }else{
                           paste(h2("They are different in", values$distanceBetween,"level."))
                         }
                       }else if(values$distanceBetween == "Family-Genus"){
                         if(values$group == "NAS" ){
                           HTML(paste("<h4><ul><li>Suggest using the mt 12S rRNA gene as an alternative genetic marker</li>
                          <li>Refer to the suggested PCR primer list for the 12S primer for nematodes</li></ul></h4>
                         "))
                         }else if (values$group == "NS"){
                           HTML(paste("<h4><ul><li>Suggest using the mt COI, 12S, or 16S rRNA gene as an alternative genetic marker</li>
                                <li>Refer to the suggested PCR primer list</li></ul>
                            </h4>"))
                         }else if (values$group == "NT"){
                           paste(h2("They are between in", values$distanceBetween))
                         }else if (values$group == "TR"){
                           HTML(paste("<h4><ul><li>Suggest using the nuclear 28S rRNA gene as an alternative genetic marker</li>
                                <li>Refer to the suggested PCR primer list for the 28S primer for platyhelminths</li></ul>
                            </h4>"))
                         }else if (values$group == "CE"){
                           HTML(paste("<h4><ul><li>Suggest using the mt 16S rRNA or COI gene as an alternative genetic marker</li>
                          <li>Refer to the suggested PCR primer list</li></ul>
                          </h4>"))
                         }else{
                           paste(h2("They are different in", values$distanceBetween,"level."))
                         }
                       }else if(values$distanceBetween == "Order-Family"){
                         paste(h2("They are between in", values$distanceBetween))
                       }else if(values$distanceBetween == "Out"){
                         #out of bounds
                         paste(h1("Out Of Bounds"))
                       }else{
                         paste(h2("They are different in", values$distanceBetween,"level."))
                       }

                     }



                   })


                   output$distplot <- renderPlot({
                     req(values$displayTable)
                     ShowRanges(group = values$group,ranges = values$ranges,distance = values$distance, marker = values$marker)

                   })



                 })

  output$distanceText <- renderUI({
    req(!is.null(input$fastaFile))
    distance <- NULL
    if(!is.null(values$distance)) distance <- round(values$distance,3)
    paste0()
    tagList(
      div(id="dtext",
          p(h5("Distance between taxons:"), distance)
      )
    )
  })

  output$tabsOutput <- renderUI({
    if(input$sidebarMenu == "fasta_file" && values$displayTable){
      req(!is.null(input$fastaFile))
      tagList(

        tabsetPanel(
          tabPanel("Distance Plot",
                   fluidRow(
                   column(12,
                          h4(textOutput("taxa")),



                            column(8,
                                   plotOutput('distplot')
                            ),
                            column(4,
                                   htmlOutput('AbbyText')
                            ),


                          )
                   ),
          ),
          tabPanel("Tree Plot",
                   uiOutput('treeOutput'),
          )

        ),

      )
    }else{
      tagList(

        fluidRow(
        column(12,
               h4(textOutput("taxa")),

                 column(8,
                        plotOutput('distplot')
                 ),
                 column(4,
                        htmlOutput('AbbyText')
                 ),
               )
        ),
      )
    }
  })

  output$tbl <- renderTable({

    tmp <- Marker.data(data = values$data, marker = input$marker, OutTable = TRUE) |> unlist()
    df <- data.frame(x=c("Order","Family","Genus","Species"),
                     Mean= c(tmp[2],tmp[6],tmp[10],tmp[14]),
                     SD= c(tmp[3],tmp[7],tmp[11],tmp[15]),
                     Min= c(tmp[4],tmp[8],tmp[12],tmp[16]),
                     Max= c(tmp[5],tmp[9],tmp[13],tmp[17])
    )
    names(df)[1] <- c(tmp[1])
    df
  })


}

shinyApp(ui = ui, server = server)

