
##### ABIyApp  by slphyx

app_version <- "V0.3"


library(shiny)
library(shinythemes)
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
                   NAS = read.csv('./data/nematode1.csv'),
                   NS = read.csv('./data/nematode2.csv'),
                   NT = read.csv('./data/nematode3.csv'),
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



# ui.bak <- fluidPage(
#     titlePanel("ABI App"),
#     sidebarLayout(
#         sidebarPanel(
#             selectInput(inputId = 'group', label = "Please choose your helminth group",
#                     choices = c("Nematode (Ascaridida and Spirurida)",
#                           "Nematode (Strongylida)",
#                           "Nematode (Trichocephalida)",
#                           "Trematode", "Cestode")
#                         ),
#
#             uiOutput('geninput.marker'),
#             uiOutput('geninput.level'),
#             actionButton('submit',"Submit")
#
#         ),
#         mainPanel(
#             htmlOutput('Abby.ans'),
#             plotOutput('distplot')
#         )
#     )
# )


ui <- fluidPage(
  
    title = "ABI App",
    theme = shinytheme("readable"),

    navbarPage("ABIApp",
        tabPanel("About",
            # includeMarkdown("./www/text/about.md"),
            # tags$img(id="Diagram",src="./img/Abi_Diagram.png" ,alt="Abi Diagram")
            includeCSS("./www/css/style.css"),
            div(
              tags$img(id="ABIApp",src="./img/ABI_app.png" ,alt="ABI App"),
              tags$hr()
            ),
            bs_accordion(id = "about") %>%
              bs_set_opts(panel_type = "default", use_heading_link = TRUE) %>%
              bs_append(title = "General Information", content = includeMarkdown("./www/text/General Information.md")) %>%
              bs_append(title = "How ABI app works", content = includeMarkdown("./www/text/app work.md")) %>%
              bs_append(title = "Benefits", content = includeMarkdown("./www/text/Benefits.md")) %>%
              bs_append(title = "Assumptions", content = includeMarkdown("./www/text/Assumptions.md")) %>%
              bs_append(title = "References", content = includeMarkdown("./www/text/References.md"))


        ),
        tabPanel("Quick guidelines",
                 fluidRow(
                 column(10,tags$h4("Quick guidelines")),
                 column(2,uiOutput("button.flag"))
                 ),
                 navbarPage("",
                   tabPanel("Before you start",
                            uiOutput("beforeStart")
                            ),
                   tabPanel("Input data",
                            uiOutput("inputData")
                            ),
                   tabPanel("Output visualization",
                            tabsetPanel(
                              tabPanel("Fasta file",
                                uiOutput("outputVisualization_fasta")
                              ),
                              tabPanel("Genetic distance",
                                uiOutput("outputVisualization")
                              )
                            )

                            ),
                   tabPanel("Suggested primers",
                            uiOutput("suggestedPrimers"),
                              tags$a("Suggested PCR primers for ABIapp",target="_blank",href="pdf/Suggested_PCR_primers_for_ABIapp.pdf")
                            )
                 )
        ),
        tabPanel("Application",
            fluidRow(  
              tabsetPanel(id = "tabs",
              tabPanel("Fasta file", 
                       value = 1,
                       tags$br(),
                       column(12,
                       fileInput("fastaFile", 
                       "Upload File :"
                       ,width = "500px"
                       , accept = ".fasta"),
                       ),
                       column(5,
                              selectizeInput("fastaSelect1", 
                                             "Sequence ID of queried sequence(ID name should match name in fasta file):", 
                                             choices=c(),
                                             width = "100%"),
                       ),
                       column(5,
                              selectizeInput("fastaSelect2", 
                                             "Sequence ID of sequence to be compared with(ID name should match name in fasta file):", 
                                             choices=c(),
                                             width = "100%")
                       ),
                       column(2,
                            uiOutput('distanceText')
                       ),
                       ),
              tabPanel("Genetic distance",
                       value = 2,
                       tags$br(),
                       column(4,
                              numericInput(inputId = 'distance', label = "Genetic distance between taxa of interest (p-distance:0-1):",
                                           min = 0,max = 1,width = "100%", step = 0.001, value = 0),

                       ),
              ),
              ),#end tab
                       column(4,
                              selectInput(inputId = 'group', label = "Please choose your helminth group",
                                          choices = c("Nematode (Ascaridida and Spirurida)",
                                                      "Nematode (Strongylida)",
                                                      "Nematode (Trichocephalida)",
                                                      "Trematode (Plagiorchiida)",
                                                      "Trematode (Diplostomida)",
                                                      "Cestode"),
                                          width = "100%"
                              ),
                       ),
                       column(4,
                              uiOutput('geninput.marker')
                       ),
                       column(12,
                              actionButton('submit',"Go!"),
                              ),
            ),
            tableOutput('tbl'),
            hr(),
            uiOutput("tabsOutput"),


        )

    )

)


server <- function(input, output,session) {

    values <- reactiveValues()
    values$distance <- NULL
    values$heighttree <- "400px"
    values$displayTable <- F
    values$language <- "EN"
    

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
      plotOutput('treeplot',width = "80%",height =values$heighttree)
      )
    })
    
    observeEvent(input$tabs,{
      values$displayTable <-F
    })
    output$taxa <- renderText({
      req(input$tabs == 1 && values$displayTable)
      paste0("Based on a genetic distance of ", round(values$distance,3),
             " between ",input$fastaSelect1, 
             " and " , input$fastaSelect2)
    })
    
    observeEvent(c(input$submit,input$fastaSelect1,
                   input$fastaSelect2),{
      if(input$tabs == 1 && !is.null(input$fastaFile)){
        HTML(paste("<h4>Error
                         </h4>"))
      }
        req((input$tabs == 1 && !is.null(input$fastaFile)) || input$tabs == 2)
        # values$level <- input$level
        values$displayTable <-T
        values$distance <- 0
        values$distance_fasta_max <- 0
        #Fasta file
        if(input$tabs == 1){

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

      if(input$tabs == 1){
        req(!is.null(input$fastaFile))
      tagList(

        tabsetPanel(
          tabPanel("Tree Plot",
            uiOutput('treeOutput'),
          ),
          tabPanel("Distance Plot",
            column(12,
                   h4(textOutput("taxa")),
                   fluidRow(
                     
                     
                     column(8,
                            plotOutput('distplot')
                     ),
                     column(4,
                            htmlOutput('AbbyText')
                     ),
                     
                     
                   )
            ),
          ),

        ),

      )
      }else{
        tagList(
          column(12,
                 h4(textOutput("taxa")),
                 fluidRow(
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
      Marker.data(data = values$data, marker = input$marker, OutTable = TRUE)
      })


}

shinyApp(ui = ui, server = server)

