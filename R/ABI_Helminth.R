#' A tool for helminth species delimitation at various taxonomic levels
#'
#' @param distance Number >= 0 (Number should be 0.000 - 0.735)
#' @param group group of Helminth ("NAS","NS","NT","TR","TRD","CE") \cr
#' "NAS" is "Nematode (Ascaridida and Spirurida)" \cr
#' "NS"  is "Nematode (Strongylida)" \cr
#' "NT"  is "Nematode (Trichocephalida)" \cr
#' "TR"  is "Trematode (Plagiorchiida)" \cr
#' "TRD"  is "Trematode (Diplostomida)" \cr
#' "CE"  is "Cestode" \cr
#' @param marker Helminth Genetic Markers
#' ("18S rRNA","28S rRNA","ITS1","ITS2","COI","COII","cytB","NAD1","12S rRNA","16S rRNA")
#' @return
#' Plot of ggplot
#' @export ABI_Helminth
#' @examples
#' ABI_Helminth()
#' ABI_Helminth(0.06)
#' ABI_Helminth(0.02,"NS","18S rRNA")
#' ABI_Helminth(distance = 0.5,group = "CE",marker = "ITS2")
#' 
#' ##### Fasta file #####
#' library(ape)
#' sequences <- read.dna("file.fasta",format = "fasta")
#' # Convert sequences to distance matrix
#' dist_matrix <- dist.dna(sequences, model = "raw")
#'
#' # Construct neighbor-joining tree
#' nj_tree <- nj(dist_matrix)
#' genetic_distance <- cophenetic(nj_tree)
#' 
#' # set column and row name
#' colM <- colnames(genetic_distance)
#' rowM <- rownames(genetic_distance)
#'
#' # select taxa
#' distance_selected  <- genetic_distance[rowM[1],colM[2]]
#' # Use ABI_Helminth()
#' ABI_Helminth(distance = distance_selected,
#' group = "CE",marker = "ITS2")
#' 
#' #############################################
#' 
#' Warning!Some groups don\'t some markers
#' plot will show Noting
#' #ABI_Helminth(0.1,"NAS","28S rRNA")
#'
#' 
#' @import ggplot2 stringr utils
# 
ABI_Helminth <- function(distance=0, group = "NAS",marker="18S rRNA"){
  data <- Load.data(group = group)
  ranges <- Level.Available(data = data, marker = marker, level.out = TRUE)
  distanceBetween <- Check.distance.between(group=group,distance=distance,marker=marker)
  ABI_warning(distance = distance,group = group, marker = marker,distanceBetween=distanceBetween)
  ShowRanges(ranges = ranges,group = group,distance = distance, marker = marker)
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
      warning("• Suggest to use mt 16S rRNA gene or mt COII as an alternative genetic marker\n• Although 18S has the smallest gap, but the low sequence variation at the genus-species level is challenging for species delimitation\n• Suggest nematode 16S primer from nematode systematics paper")
    }else if (group == "NT" && (marker =="18S rRNA" || marker =="ITS2")){
      warning("• Suggest to use mt 12S")
    }else if (group == "TR"&& (marker =="18S rRNA" || marker =="ITS1" || marker =="ITS2")){
      warning("• Suggest to use mt 16S \n•	Although 18S has small gap (same as 16S), but the low sequence variation at the genus-species level is challenging for species delimitation")
    }else if (group == "CE"&& marker =="12S rRNA"){
      warning("• Recommendation to use another mt genetic marker")
    }
  }else{

    if(distanceBetween == "Genus-Species"){
      if(group == "NAS" || group == "NS"){
        warning("• Suggest to use mt 16S rRNA gene or mt COII as an alternative genetic marker \n• Although 18S has the smallest gap, but the low sequence variation at the genus-species level is challenging for species delimitation \n• Suggest nematode 16S primer from nematode systematics paper")
      }else if (group == "NT"){
        warning("• Suggest to use mt 12S")
      }else if (group == "TR"){
        warning("• Suggest to use mt 16S \n •Although 18S has small gap (same as 16S), but the low sequence variation at the genus-species level is challenging for species delimitation")
      }else if (group == "CE"){
        warning("• Suggest to use cytB or 12S (but rarely use cytB for cestodes)")
      }else{
        cat("They are between in", distanceBetween)
      }
    }else if(distanceBetween == "Family-Genus"){
      if(group == "NAS" ){
        warning("• Suggest to use mt 12S rRNA gene as an alternative genetic marker, with primer from nematode systematics paper")
      }else if (group == "NS"){
        warning("• Suggest to use mt COI or mt 12S or mt 16S\n•	But caution for COI primers (universal primers might not be able to amplify, should use specific primers)")
      }else if (group == "NT"){
        cat("They are between in", distanceBetween)
      }else if (group == "TR"){
        warning("• Suggest 28S, but need to caution (low sequence variation, alignment region) \n•	Other alternative is the mt genetic markers")
      }else if (group == "CE"){
        warning("• Suggest COI or 16S")
      }else{
        cat("They are between in", distanceBetween)
      }
    }else if(distanceBetween == "Order-Family"){
      cat("They are between in", distanceBetween)
    }else if(distanceBetween == "Out"){
      #out of bounds
      warning("• Out Of Bounds")
    }else{
      cat("They are different in", distanceBetween,"level.")
    }

  }
}
