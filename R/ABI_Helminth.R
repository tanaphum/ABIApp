#' A tool for helminth species delimitation at various taxonomic levels
#'
#' @param distance Number >= 0 (From our database the number sh ould be 0.000 - 0.909)
#' @param group group of Helminth ("NAS","NS","NT","TR","TRD","CE") /n
#' "NAS" is "Nematode (Ascaridida and Spirurida)" /n
#' "NS"  is "Nematode (Strongylida)") /n
#' "NT"  is "Nematode (Trichocephalida)" /n
#' "TR"  is "Trematode (Plagiorchiida)" /n
#' "TRD" is "Trematode (Diplostomida)" /n
#' "CE"  is "Cestode"
#' @param marker Helminth Genetic Markers
#' @param Fastafile directory of Fasta file
#' @param fastaSelect1 Selection of 1st taxon from Fasta file (number or text(label))
#' @param fastaSelect2 Selection of 2nd taxon from Fasta file (number or text(label))
#' ("18S rRNA","28S rRNA","ITS1","ITS2","COI","COII","cytB","NAD1","12S rRNA","16S rRNA")
#' @return
#' Plot of ggplot
#' @export ABI_Helminth
#'
#' @examples
#' ABI_Helminth()
#' ABI_Helminth(0.06)
#' ABI_Helminth(0.02,"NS","18S rRNA")
#' ABI_Helminth(distance = 0.5,group = "CE",marker = "ITS2")
#'
#' ##### Fasta file #####
#' # Select with numbers
#' ABI_Helminth(Fastafile="dir/fastaFile.fasta","NT","18S rRNA"
#' fastaSelect1 = 1, fastaSelect2 = 2
#' )
#' # Select with number and text
#' ABI_Helminth(Fastafile="dir/fastaFile.fasta","NT","18S rRNA"
#' fastaSelect1 = "Label 1", fastaSelect2 = 2
#' )
#'# Select with texts
#' ABI_Helminth(Fastafile="dir/fastaFile.fasta","NT","18S rRNA"
#' fastaSelect1 = "Label 1", fastaSelect2 = "Label 2"
#' )
#'
#' @import ggplot2 stringr utils ape
#
ABI_Helminth <- function(distance=0,Fastafile = NULL,fastaSelect1=1,fastaSelect2=2, group = "NAS",marker="18S rRNA"){
  if(distance < 0 && is.null(Fastafile)){
    warning(("• Number must more than 0"))
  }else{
    data <- Load.data(group = group)
    ranges <- Level.Available(data = data, marker = marker, level.out = TRUE)
    if(is.null(Fastafile)){
      distanceBetween <- Check.distance.between(group=group,distance=distance,marker=marker)
      ABI_warning(distance = distance,group = group, marker = marker,distanceBetween=distanceBetween)
      ShowRanges(ranges = ranges,group = group,distance = distance, marker = marker)
    }else{
        if(is.null(fastaSelect1)){
          warning(("• fastaSelect1 is NULL. System set fastaSelect1=1"))
          fastaSelect1 <- 1
        }
          if(is.null(fastaSelect2)){
          warning(("• fastaSelect2 is NULL. System set fastaSelect2=2"))
          fastaSelect2 <- 2
        }
        #read Fasta file
        sequences <- read.dna(Fastafile,format = "fasta")
        # Convert sequences to distance matrix
        dist_matrix <- dist.dna(sequences, model = "raw")
        # Construct neighbor-joining tree
        dna_tree <- nj(dist_matrix)
        # Calculate genetic distance between each pair of sequences
        genetic_distance <- cophenetic(dna_tree)

        # 1st Select
        col_Select <- fastaSelect1
        # 2nd Select
        row_Select <- fastaSelect2
        distance <-   genetic_distance[row_Select,col_Select]

        distanceBetween <- Check.distance.between(group=group,distance=distance,marker=marker)
        ABI_warning(distance = distance,group = group, marker = marker,distanceBetween=distanceBetween)
        ShowRanges(ranges = ranges,group = group,distance = distance, marker = marker)
      }
    }
}

# min/max columns
col.order <- c(3,4)
col.family <- c(7,8)
col.genus <- c(11,12)
col.species <- c(15,16)


Between <- function(x, left, right){
  if(x < 0) return(F)
  if(is.na(x >= left & x <= right)) return(F)
  else if(length(x >= left & x <= right) ==0) return(F)
  else return(x >= left & x <= right)
}

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

MinMax <- function(level, marker, data){
  tmp <- switch(level,
                Order = Marker.data(marker, data)[col.order],
                Family = Marker.data(marker, data)[col.family],
                Genus = Marker.data(marker, data)[col.genus],
                Species = Marker.data(marker, data)[col.species])

  return(c(Min = tmp[1], Max = tmp[2]))
}

Load.data <- function(group){
  dir <- system.file(package='ABI')
  data <- switch(group,
                 NT = read.csv(paste0(dir,'/ABI/data/nematode1.csv')),
                 NAS = read.csv(paste0(dir,'/ABI/data/nematode2.csv')),
                 NS = read.csv(paste0(dir,'/ABI/data/nematode3.csv')),
                 TR = read.csv(paste0(dir,'/ABI/data/trematode1.csv')),
                 TRD = read.csv(paste0(dir,'/ABI/data/trematode2.csv')),
                 CE = read.csv(paste0(dir,'/ABI/data/cestode.csv')))
  return(data)
}

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
ABI_warning <-function(distance,group,marker,distanceBetween){
    if(distance==0){
      if(group == "NAS" && (marker =="18S rRNA" || marker =="28S rRNA" || marker =="ITS1" || marker =="ITS2" || marker =="COII" ||marker =="12S rRNA")){
        warning(("• Suggest to use mt 16S rRNA gene or mt COII as an alternative genetic marker /n
                  • Although 18S has the smallest gap, but the low sequence variation at the genus-species level is challenging for species delimitation /n
                  • Suggest nematode 16S primer from nematode systematics paper
                   "))
      }else if (group == "NT" && (marker =="18S rRNA" || marker =="ITS2")){
        warning("• Suggest to use mt 12S")
      }else if (group == "TR"&& (marker =="18S rRNA" || marker =="ITS1" || marker =="ITS2")){
        warning(("• Suggest to use mt 16S /n
                       • Although 18S has small gap (same as 16S), but the low sequence variation at the genus-species level is challenging for species delimitation
                         "))
      }else if (group == "CE"&& marker =="12S rRNA"){
        warning(("• Recommendation to use another mt genetic marker"))
      }
    }else{

      if(distanceBetween == "Genus-Species"){
        if(group == "NAS" || group == "NS"){
          warning(("• Suggest using the mt 16S rRNA or COII gene as an alternative genetic marker /n
                         • Refer to the suggested PCR primer list for the 16S primer for nematodes
                          "))
        }else if (group == "NT"){
          warning(("• Suggest using the mt 12S rRNA gene as an alternative genetic marker /n
                    • Refer to the suggested PCR primer list for the 12S primer for nematodes
                    "))
        }else if (group == "TR"){
          warning(paste(" • Suggest using the mt 16S rRNA gene as an alternative genetic marker /n
                       • Refer to the suggested PCR primer list for the 16S primer for platyhelminths /n
                       "))
        }else if (group == "CE"){
          warning(paste("• Suggest using the mt 12S rRNA gene as an alternative genetic marker /n
                         • Refer to the suggested PCR primer list for the 12S primer for platyhelminths
                        "))
        }else{
          cat(paste("• They are different in", distanceBetween,"level."))
        }
      }else if(distanceBetween == "Family-Genus"){
        if(group == "NAS" ){
          warning(paste("• Suggest using the mt 12S rRNA gene as an alternative genetic marker/n
                         • Refer to the suggested PCR primer list for the 12S primer for nematodes
                         "))
        }else if (group == "NS"){
          warning(paste("• Suggest using the mt COI, 12S, or 16S rRNA gene as an alternative genetic marker/n
                         • Refer to the suggested PCR primer list
                          "))
        }else if (group == "NT"){
          warning(h2("They are between in", distanceBetween))
        }else if (group == "TR"){
          warning(("• Suggest using the nuclear 28S rRNA gene as an alternative genetic marker/n
                 • Refer to the suggested PCR primer list for the 28S primer for platyhelminths
                      "))
        }else if (group == "CE"){
          warning(("• Suggest using the mt 16S rRNA or COI gene as an alternative genetic marker/n
                    • Refer to the suggested PCR primer list
                    • "))
        }else{
          cat(paste("• They are different in", distanceBetween,"level."))
        }
      }else if(distanceBetween == "Order-Family"){
        warning(paste("• They are between in", distanceBetween))
      }else if(distanceBetween == "Out"){
        #out of bounds
        warning(("• Out Of Bounds"))
      }else{
        cat(paste("• They are different in", distanceBetween,"level."))
      }

    }
}
