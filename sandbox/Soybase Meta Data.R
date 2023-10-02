##################################################
## Project: Oil Mapping Manuscript
## Script purpose: Various functions to scrape and format data from soybase
## Into more easily readable tables/plots.
## Date: October 12, 2020
## Author: jay Gillenwater
##################################################


# Base search URL https://www.soybase.org/search/index.php?searchterm=%22Seed+Oil%22&list=bi_parental_qtl_listview
# A script for scraping QTL metadata from soybase.org
# Also processing of scraped data to publication ready tables and plots
source("./packages.R")

# The function takes some QTL trait as it's only argument. For a complete list of potential traits, see this page: https://www.soybase.org/search/qtllist.php
# This scrapes a substantial amount of data for all QTL reported for the given trait.
# This function unfortunately got quite large but I've tried to document the steps within the function. 
Soybase_QTL_Info <- function(TraitName = "Seed Oil"){
  
  TraitName <- gsub(" ", "+", TraitName) # Replace spaces with "+" symbols
  
  # The webpage which lists the bipaental QTL for a given trait
  QTL_List_URL <- paste("https://www.soybase.org/search/index.php?searchterm=", 
                        TraitName, 
                        "&list=bi_parental_qtl_listview", 
                        sep = "")
  
  # List of qtl names
  QTL_List <- tryCatch(
    read_html(QTL_List_URL) %>%
      html_node('#bi_parental_qtl_listview') %>%
      html_table(fill = TRUE),
    error = function(e) return(NA)
  )
  
  if(!is.data.frame(QTL_List)){
    return(NA)
  }
  
  
  # Cleaning the table
  EliminateRows      <- which(QTL_List[, 3] == "QTL Name") - 1 # Remove headers that are placed within the table body
  QTL_List           <- QTL_List[-EliminateRows, ]
  colnames(QTL_List) <- QTL_List[1, ]                          # Use the first row for column names and then remove the first row
  QTL_List           <- QTL_List[-1, ]
  EliminateRows      <- which(QTL_List[, 3] == "QTL Name")     # Remove more headers
  
  if(length(EliminateRows)){
    QTL_List         <- QTL_List[-EliminateRows, ]
  }
  
  EliminateCols      <- which(is.na(colnames(QTL_List)) | colnames(QTL_List) == "NA" | colnames(QTL_List) == "") # Remove empty columns
  QTL_List           <- QTL_List[, -EliminateCols]
  
  # Line for testing a range of pages (middle page has empty tables), keep for testing future weird tables
  # QTL_List <- QTL_List %>% filter(`QTL Name` %in% c("Seed oil 6-2", "Seed oil 6-3", "Seed oil 6-4"))
  
  # The entries in the 'QTL Name' column are (unsurprisingly) the names of the QTLs
  # Each QTL has it's own page of additional information
  QTL_URL_Base <- "https://www.soybase.org/sbt/search/search_results.php?category=QTLName&search_term="          # Base URL used to find information for each named QTL
  
  
  # Work down each QTL in the table and scrape additional information
  Additional_QTL_Data <- vector('list', length = nrow(QTL_List))
  for(QTL in 1:length(Additional_QTL_Data)){
    CurrentQTL <- QTL_List$`QTL Name`[[QTL]]
    CurrentURL <- paste(QTL_URL_Base, gsub(" ", "+", CurrentQTL), sep = "")
    
    CurrentPage <- read_html(CurrentURL)
    
    # Some pages had empty tables that were messing up the scraping, this is a temporary
    # solution where the scraper will just return a NA if the scraping failed and move on.
    # This is a rough fix that ignores the underlying causes of why a scrape would fail for some
    # page and should be fixed in the future
    Additional_QTL_Data[[QTL]] <- tryCatch(
      {
        CurrentPage %>% html_nodes('table') %>% html_table(fill = TRUE)
      },
      error=function(e){
        message(paste("Scraping failed for QTL: ",  QTL_List$`QTL Name`[[QTL]]))
        return(NA)
      }
    )
    
    if(any(is.na(Additional_QTL_Data[[QTL]]))){
      next
    }
    
    Additional_QTL_Data[[QTL]] <- Additional_QTL_Data[[QTL]][-length(Additional_QTL_Data[[QTL]])] # Remove last table
    
    # Get table names
    Data_headers <- CurrentPage %>% html_nodes(".head_link") %>% html_text(trim = TRUE)
    First_Name <- CurrentPage %>% html_nodes('.top_bar') %>% html_text(trim = TRUE)
    Data_headers <- c(First_Name, Data_headers)
    
    
    # This section often has several tables which will interfere with scraping. The upshot is that data from the following
    # sections will be merged into a single table which will be copied the same number of times as there are tables in this section.
    # To fix this, I selected only the first few rows of the first tabe and removed the copied tables. See the error correction section
    # for more details
    if("Loci positively associated with the QTL" %in% Data_headers & length(Data_headers) < length(Additional_QTL_Data[[QTL]])){
      # Which table holds the associated QTL
      Assoc_Table <- which(Data_headers == "Loci positively associated with the QTL")
      
      # How many repeats to trim
      TrimNum <- length(Additional_QTL_Data[[QTL]]) - length(Data_headers)
      
      # Indexes of tables to remove
      Rem_Tables <-(Assoc_Table+1):(Assoc_Table+TrimNum)
      
      # Trim association table to include only desired data
      Additional_QTL_Data[[QTL]][[Assoc_Table]] <- Additional_QTL_Data[[QTL]][[Assoc_Table]][1:TrimNum, ]
      Additional_QTL_Data[[QTL]] <- Additional_QTL_Data[[QTL]][-Rem_Tables]
    }
    
    
    names(Additional_QTL_Data[[QTL]]) <- Data_headers # Set table names
    Additional_QTL_Data[[QTL]] <- Additional_QTL_Data[[QTL]][!is.na(names(Additional_QTL_Data[[QTL]]))] # Remove unknown tables (also a rough fix)
    
    print(paste("Scraped:", QTL_List$`QTL Name`[[QTL]]))
    Sys.sleep(10) # Pause system to prevent spamming the website
  }
  
  names(Additional_QTL_Data) <- QTL_List$`QTL Name`
  
  Results <- list("QTL List" = QTL_List, "Additional Data" = Additional_QTL_Data)
  Results
}

# A function to format the "Additional Data" section of the output into a tidier format
tidy_AdditionalData <- function(AdditionalData){
  
  # A function to extract a numeric portion from a string
  numextract <- function(string){ 
    str_extract(string, "\\-*\\d+\\.*\\d*")
  } 
  
  # Available fields in the current QTL
  AvailableData <- names(AdditionalData)
  
  # The name of the SNP
  QTLName <- AvailableData[[1]]
  
  # Name of the paper in which the QTL was published
  if('References for the QTL' %in% AvailableData){
    LitName <- AdditionalData[['References for the QTL']][[2]]
    LitName <- data.frame(LiteratureTitle = LitName, QTL = QTLName)
    LitName$PubYear   <- as.numeric(numextract(AdditionalData[['References for the QTL']][[1]]))
    LitName$ShortName <- AdditionalData[['References for the QTL']][[1]]
  }else{
    LitName <- data.frame(LitName = character(0), QTL = character(0), PubYear = numeric(0), ShortName = character(0))
  }
  
  # Methods for QTL detection
  if('Methods used to identify the QTL' %in% AvailableData){
    DetectMethods <- AdditionalData[['Methods used to identify the QTL']]
    colnames(DetectMethods) <- "MethodOfDetection"
    DetectMethods$QTL <- QTLName
  }else{
    DetectMethods <- data.frame(MethodofDetection = character(0), QTL = character(0))
  }
  
  # Mapping population used
  if('Population types used in identification of the QTL' %in% AvailableData){
    Poptype <- AdditionalData[['Population types used in identification of the QTL']]
    colnames(Poptype) <- "MappingPopulationType"
    Poptype$QTL <- QTLName
  }else{
    Poptype <- data.frame(MappingPopulationType = character(0), QTL = character(0))
  }
  
  # Associated Loci
  if('Loci associated with the QTL' %in% AvailableData){
    AssocLoci <- AdditionalData[['Loci associated with the QTL']]
    colnames(AssocLoci) <- "AssociatedLoci"
    AssocLoci$QTL <- QTLName
  }else{
    AssocLoci <- data.frame(AssociatedLoci = character(0), QTL = character(0))
  }
  
  # Phenotypic R2
  if('Loci positively associated with the QTL' %in% AvailableData){
    AssocLoci_pos <- AdditionalData[['Loci positively associated with the QTL']]
    colnames(AssocLoci_pos) <- c("AssociatedLoci", "LociFactor", "FactorValue")
    AssocLoci_pos[, 1:ncol(AssocLoci_pos)] <- apply(AssocLoci_pos[, 1:ncol(AssocLoci_pos)], 2, as.character)
    AssocLoci_pos$QTL <- QTLName
  }else{
    AssocLoci_pos <- data.frame(AssociatedLoci = character(0), LociFactor = character(0), FactorValue = character(0), QTL = character(0))
  }
  
  # Where is the QTL found on the composite map?
  MapString   <- paste('Maps containing', QTLName)
  if(MapString %in% AvailableData){
    MapPosition <- AdditionalData[[MapString]]
    MapPosition$QTL <- QTLName
  }else{
    MapPosition <- NA
  }
  
  
  
  # Parent information, percent variation explained for the SNP
  ParentalInfo           <- as.data.frame(t(AdditionalData[[1]]))
  colnames(ParentalInfo) <- sub(":", "", as.character(as.matrix(ParentalInfo[1, ])))
  
  ParentalInfo     <- ParentalInfo[-1, ]
  ParentalInfo$QTL <- QTLName
  
  results <- list("ParentalInfo"  = ParentalInfo,
                  "MapPosition"   = MapPosition,
                  "LitName"       = LitName,
                  "AssocLoci"     = AssocLoci,
                  "Detectmethods" = DetectMethods,
                  "Poptype"       = Poptype,
                  "AssocLociMeta" = AssocLoci_pos)
  results
}



# A function to create some summary tables for a scrape
MetaSummaries <- function(Scrape){
  
  # Clean up the scrape using the tidying function from above
  Scrape_Clean <- lapply(Scrape[["Additional Data"]], tidy_AdditionalData)
  
  # Allocate some lists to hold parts of the cleaned data
  ParentInfo <- AssocLoc <- PubYears <- PublicationsByYear <- vector("list", length = length(Scrape_Clean))
  
  AssocLoc <- lapply(Scrape_Clean, function(x) x[["AssocLoci"]])
  AssocLoc <- AssocLoc[which(!is.na(AssocLoc))]
  
  PosAssocLoc <- lapply(Scrape_Clean, function(x) x[["AssocLociMeta"]])
  PosAssocLoc <- PosAssocLoc[which(!is.na(PosAssocLoc))]
  PosAssocLoc <- do.call(bind_rows, PosAssocLoc)
  
  ParentInfo <- lapply(Scrape_Clean, function(x) x[["ParentalInfo"]])
  ParentInfo <- ParentInfo[which(!is.na(ParentInfo))]
  
  PubYears <- lapply(Scrape_Clean, function(x) x[["LitName"]])
  AllPublications <- PubYears
  AllPublications <- do.call(bind_rows, AllPublications)
  AllPublications <- left_join(AllPublications, Scrape$`QTL List`, by = c("QTL" = "QTL Name"))
  
  PubYears <- PubYears[which(!is.na(PubYears))]
  
  PublicationsByYear <- PubYears
  
  PubYears <- lapply(PubYears, function(x) x[["PubYear"]])
  
  # How many QTL are detected by year?
  QTL_Tabulate <- as.data.frame(table(unlist(PubYears)))
  colnames(QTL_Tabulate) <- c("Year", "QTLCount")
  
  # How many QTL are detected by each paper
  PublicationsByYear <- do.call(bind_rows, PublicationsByYear)
  
  QTL_Per_Paper <- PublicationsByYear %>% group_by(LiteratureTitle) %>% summarise(QTLCount = n()) %>% arrange(desc(QTLCount))
  
  # How many papers are published each year?
  Papers_Per_Year <- PublicationsByYear %>% group_by(PubYear) %>% summarise(PaperCount = n()) %>% arrange(PubYear)
  
  # How many papers, QTL in total?
  nPapers <- nrow(QTL_Per_Paper)
  nQTL    <- length(Scrape_Clean)
  
  # Data on parents, other metadata on a QTL-level
  ParentInfo <- do.call(bind_rows, ParentInfo)
  
  # All loci associated with QTLs
  AssocLoc <- do.call(bind_rows, AssocLoc)
  
  # Tabulate associated Loci
  LociTabulate <- as.data.frame(table(AssocLoc$AssociatedLoci)) %>% arrange(desc(Freq))
  colnames(LociTabulate) <- c("Loci", "QTLCount")
  
  # A table to hold a general summary for printing data about where QTL may be found, in which study were they found
  # and what methods were used to detect them, created by joining multiple tables using the QTL name as an identifier
  
  # Need a table of QTL names and papers, methods and QTL name, and position and QTL name
  
  results <- list("QTLByYear"        = QTL_Tabulate,
                  "QTL_Per_Paper"    = QTL_Per_Paper,
                  "Papers_Per_Year"  = Papers_Per_Year,
                  "QTLMetadata"      = ParentInfo,
                  "AssociatedLoci"   = AssocLoc,
                  "PosAssociatedLoci"= PosAssocLoc,
                  "LociTabulation"   = LociTabulate,
                  "PublicationTable" = AllPublications)
  
  results
}


# A function to download images showing the positions of QTL for various traits
# on the composite linkage map from soybase

# An example URL, split up to show its components
# baseURL <- "https://www.soybase.org/cmap/cgi-bin/cmap/viewer?dotplot=0&
# eliminate_orphans=0&
# mapMenu=&
# featureMenu=1&
# corrMenu=0&
# displayMenu=1&
# advancedMenu=&
# map_start_0_GmComposite2003_A1=-36.1&
# map_stop_0_GmComposite2003_A1=104.58&
# map_mag_0_GmComposite2003_A1=1&map_in_menu_0_GmComposite2003_A1=1&highlight=&
# ft_Gene=0&
# ft_QTL_fungal=0&
# ft_QTL_inorganic=0&
# ft_QTL_insect=0&
# ft_QTL_leaf-stem=0&
# ft_QTL_misc=0&
# ft_QTL_nematode=0&
# ft_QTL_oil=2&
# ft_QTL_other-seed=0&
# ft_QTL_pod=0&
# ft_QTL_protein=2&
# ft_QTL_reprod-period=0&
# ft_QTL_root=0&
# ft_QTL_viral=0&
# ft_QTL_whole-plant=0&
# ft_QTL_yield=2&
# ft_RFLP=0&
# ft_SNP=0&
# ft_SSR=0&
# ft_DEFAULT=0&
# label_features=all&
# collapse_features=0&
# evidence_type_ANB=1&
# ets_ANB=0&
# evidence_type_map_based=0&
# ets_map_based=0&
# aggregate=0&
# corrs_to_map=0&
# show_intraslot_corr=0&
# split_agg_ev=0&
# sub=Redraw&
# phrb=800&
# pixel_height=1750&
# font_size=medium&
# image_type=png&
# clean_view=0&
# hide_legend=0&
# dotplot_ps=1&
# scale_maps=1&
# omit_area_boxes=0&
# comp_menu_order=display_order&
# ignore_image_map_sanity=0&
# session_id=a3d5f33843ced3f73935ca533b46dfce&step=6&
# session_mod=&
# ref_map_set_acc=GmComposite2003_&ref_map_accs=GmComposite2003_A1&
# ref_species_acc=Glycine_max&
# flip=&data_source=sbt_cmap"
# 
# baseURL <- str_replace_all(baseURL, "[\r\n]" , "")
# 
# # Read the page cotaining the map
# MapPage <- read_html(baseURL) 
# 
# # Extract the url for the image
# ImgAttrs <- MapPage %>% html_nodes("#image_div") %>% html_attrs() %>% unlist()
# ImgUrl <- paste("https://www.soybase.org", ImgAttrs["src"], sep = "")
# download.file(ImgUrl, destfile = "LG_A1.png", mode = "wb")

# Each linkage group can be found by giving the linkage group name, along with start and
# ending positions to display in units of cM. I have already created a table to store
# this data to make it easier
LGRanges <- read.csv("https://raw.githubusercontent.com/jhgille2/SoybaseData/master/LinkageGroupRanges.csv")

# A function to download genetic maps for each chromosome to some destination with 
# requested features
getGeneticMaps <- function(ImgHeight = 1750, destDirectory = getwd()){
  # For each linkage group, get the starting and ending positions of the LGs
  for(i in 1:nrow(LGRanges)){
    CurrentLG <- as.character(LGRanges[i, "LG"])
    LGStart   <- LGRanges[i, "Pos_Start"]
    LGEnd     <- LGRanges[i, "Pos_End"]
    
    CurrentURL <- paste("https://www.soybase.org/cmap/cgi-bin/cmap/viewer?dotplot=0&
                        eliminate_orphans=0&
                        mapMenu=&
                        featureMenu=1&
                        corrMenu=0&
                        displayMenu=1&
                        advancedMenu=&
                        map_start_0_GmComposite2003_", CurrentLG, "=", LGStart, "&
                        map_stop_0_GmComposite2003_", CurrentLG, "=", LGEnd, "&
                        map_mag_0_GmComposite2003_", CurrentLG, "=1&
                        map_in_menu_0_GmComposite2003_", CurrentLG, "=1&
                        highlight=&
                        ft_Gene=0&
                        ft_QTL_fungal=0&
                        ft_QTL_inorganic=0&
                        ft_QTL_insect=0&
                        ft_QTL_leaf-stem=0&
                        ft_QTL_misc=0&
                        ft_QTL_nematode=0&
                        ft_QTL_oil=2&
                        ft_QTL_other-seed=0&
                        ft_QTL_pod=0&
                        ft_QTL_protein=2&
                        ft_QTL_reprod-period=0&
                        ft_QTL_root=0&
                        ft_QTL_viral=0&
                        ft_QTL_whole-plant=0&
                        ft_QTL_yield=2&
                        ft_RFLP=0&
                        ft_SNP=0&
                        ft_SSR=0&
                        ft_DEFAULT=0&
                        label_features=all&
                        collapse_features=0&
                        evidence_type_ANB=1&
                        ets_ANB=0&
                        evidence_type_map_based=0&
                        ets_map_based=0&
                        aggregate=0&
                        corrs_to_map=0&
                        show_intraslot_corr=0&
                        split_agg_ev=0&
                        sub=Redraw&
                        phrb=800&
                        pixel_height=", ImgHeight, "&
                        font_size=medium&
                        image_type=png&
                        clean_view=0&
                        hide_legend=0&
                        dotplot_ps=1&
                        scale_maps=1&
                        omit_area_boxes=0&
                        comp_menu_order=display_order&
                        ignore_image_map_sanity=0&
                        ref_map_set_acc=GmComposite2003_&
                        ref_map_accs=GmComposite2003_", CurrentLG, "&
                        ref_species_acc=Glycine_max&
                        flip=&data_source=sbt_cmap", sep = "") 
    
    CurrentURL      <- str_replace_all(CurrentURL, "[\r\n]" , "") %>% str_replace_all(fixed(" "), "")
    CurrentImgAttrs <- read_html(CurrentURL) %>% html_nodes("#image_div") %>% html_attrs() %>% unlist()
    
    CurrentImgUrl <- paste("https://www.soybase.org", CurrentImgAttrs["src"], sep = "")
    
    # Need to fix this part
    destpath <- paste(destDirectory, "/", paste("LG_", CurrentLG, ".png", sep = ""), sep = "")
    download.file(CurrentImgUrl, destfile = destpath, mode = "wb")
    
    Sys.sleep(5)
  }
}
setwd("C:/Users/Jay/Desktop/Oil MP Paper/Soybase Genetic Maps")
getGeneticMaps()

# Example for a different trait (Seed Yield)
# YieldScrape <- Soybase_QTL_Info("Seed yield")
# YieldMeta   <- MetaSummaries(YieldScrape)


# Collecting data for multiple traits
# Download QTL position data from soybase
QTLPos <- read.csv(here::here("data", "SoybaseQTL_2023_09_08.csv"))
QTLPos$Object.Type <- as.character(QTLPos$Object.Type)
QTLPos$QTLCategory <- gsub("\\d+-\\d+", "", QTLPos$Object.Name) %>% str_trim()
QTLPos$QTLCategory <- gsub("-\\d+", "", QTLPos$QTLCategory)
QTLPos$QTLCategory <- gsub("\\.\\d+", "", QTLPos$QTLCategory) %>% str_trim()

# Shorter names for the QTL general categories
LongQTLNames <- as.character(unique(QTLPos$Object.Type))

ShortQTLNames <- c("Other Seed",
                   "Whole Plant",
                   "Inorganic",
                   "Fungal",
                   "Insect",
                   "Leaf-stem",
                   "Misc",
                   "Oil",
                   "Protein",
                   "Reproductive Period",
                   "Yield",
                   "Viral",
                   "Nematode",
                   "Root",
                   "Pod")

QTLNameConversion <- data.frame(LongName = LongQTLNames, ShortName = ShortQTLNames)

# How many QTL in each major category?
QTL_MajorCategoryCount <- QTLPos %>% group_by(Object.Type) %>% summarise(GroupCount = n()) %>% arrange(desc(GroupCount)) %>% ungroup()

# How many QTL in each sub-ctegory?
QTL_SubcategoryCount <- QTLPos %>% group_by(Object.Type, QTLCategory) %>% summarise(GroupCount = n()) %>% arrange(Object.Type, desc(GroupCount)) 

# Filter to just protein, oil, and yield
QTLCounts_ProtOilYield <- QTL_SubcategoryCount %>% filter(Object.Type %in% c("QTL_oil", "QTL_protein", "QTL_yield"))

# The traits which I want meta-data for
QTL_Types <- QTLPos %>% filter(Object.Type %in% c("QTL_oil", "QTL_protein", "QTL_yield"))
QTL_Types$LG <- as.character(QTL_Types$LG)
QTL_Types$LG <- sub("D1A", "D1a", QTL_Types$LG)
QTL_Types$LG <- sub("D1B", "D1b", QTL_Types$LG)


QTL_MajorCategoryCount <- left_join(QTL_MajorCategoryCount, QTLNameConversion, by = c("Object.Type" = "LongName"))
QTL_MajorCategoryCount <- QTL_MajorCategoryCount %>% select(ShortName, GroupCount)
colnames(QTL_MajorCategoryCount) <- c("QTL Category", "QTL Found")



# I've already done this part, skip to the next section to load in the file
AllTraits <- gsub(" ", "+", QTL_Types$QTLCategory) %>% unique()
AllScrapes <- vector("list", length = length(AllTraits))
for(i in 1:length(AllScrapes)){
  AllScrapes[[i]] <- Soybase_QTL_Info(AllTraits[[i]])
  print(noquote(paste("Finished Scrape for trait:", AllTraits[[i]])))
}


# Load in the completed scrape
load("C:/Users/Jay/Desktop/Oil MP Paper/AllProtOilYield.RData")
names(AllScrapes) <- AllTraits

# Clean all the additional data from every scrape

# Return only those traits which had additional data
AllScrapes <- AllScrapes[which(!unlist(lapply(AllScrapes, function(x) all(is.na(x)))))]


AllScrapes_CleanAdditional <- pblapply(AllScrapes, MetaSummaries)
names(AllScrapes_CleanAdditional) <- gsub("\\+", " ", names(AllScrapes_CleanAdditional))



# Pull out the literature tables
AllScrapes_Lit <- lapply(AllScrapes_CleanAdditional, function(x) x[["PublicationTable"]])

QTL_Categories <- QTL_Types$Object.Type[match(names(AllScrapes_CleanAdditional), QTL_Types$QTLCategory)]
for(i in 1:length(AllScrapes_Lit)){
  AllScrapes_Lit[[i]]$QTL_Type <- QTL_Categories[[i]]
}

# Bind all to a single dataframe
AllScrapes_Lit <- do.call(bind_rows, AllScrapes_Lit)
AllScrapes_Lit$QTL_Type <- as.factor(AllScrapes_Lit$QTL_Type)

# Pull all loci associated with QTL
AllScrapes_AssocLoci <- lapply(AllScrapes_CleanAdditional, function(x) x[["AssociatedLoci"]])
AllScrapes_AssocLoci <- do.call(bind_rows, AllScrapes_AssocLoci) %>% unique()

AllScrapes_Lit <- left_join(AllScrapes_Lit, AllScrapes_AssocLoci, by = "QTL")


# The phenotypic R2 for each marker
PhenoRSq <- lapply(AllScrapes_CleanAdditional, function(x) x[["PosAssociatedLoci"]])
PhenoRSq <- do.call(bind_rows, PhenoRSq)
PhenoRSq <- PhenoRSq %>% filter(LociFactor == "Phenotypic_R2") %>% unique() %>% select(AssociatedLoci, FactorValue, QTL)
colnames(PhenoRSq) <- c("AssociatedLoci", "R2", "QTL")
PhenoRSq$R2 <- round(as.numeric(PhenoRSq$R2), 2)

# Some R2 values are already converted to a % scale, convert all to be between 0 and 1
PhenoRSq$R2[which(PhenoRSq$R2 > 1)] <- PhenoRSq$R2[which(PhenoRSq$R2 > 1)]/100

# Add this to the literature table
AllScrapes_Lit <- left_join(AllScrapes_Lit, PhenoRSq, by = c("QTL", "AssociatedLoci"))

# Counts of QTL detected for each trait, by paper
QTLCounts_ByPaper <- AllScrapes_Lit %>% 
  group_by(LiteratureTitle, PubYear, QTL_Type) %>% 
  summarise(Count = n()) %>%
  ungroup() %>%
  spread(QTL_Type, Count, fill = 0) %>%
  mutate(TotalQTL = QTL_oil + QTL_protein + QTL_yield, QTL_GroupCount = as.numeric(QTL_oil > 0) + as.numeric(QTL_protein > 0) + as.numeric(QTL_yield > 0)) %>%
  arrange(desc(QTL_GroupCount), desc(PubYear), desc(TotalQTL))

# papers which identified qtl in all 3 categories
AllThreeLit <- QTLCounts_ByPaper$LiteratureTitle[QTLCounts_ByPaper$QTL_GroupCount == 3]
AllThreeLit <- AllScrapes_Lit %>% filter(LiteratureTitle %in% AllThreeLit)
AllThreeLit_reduced <- AllThreeLit %>% 
  select(ShortName, PubYear, QTL_Type, QTL, `LG(GmComposite2003)`, `Start(cM)`, `End(cM)`, AssociatedLoci, R2) %>% 
  arrange(PubYear, ShortName, QTL_Type, QTL, `LG(GmComposite2003)`) %>%
  select(ShortName, QTL_Type, QTL, `LG(GmComposite2003)`, `Start(cM)`, `End(cM)`, AssociatedLoci, R2)

AllThreeLit_reduced$QTL_Type <- capitalize(sub("QTL_", "", AllThreeLit_reduced$QTL_Type))
colnames(AllThreeLit_reduced) <- c("Literature Source", "QTL Category", 'QTL Name', "LG(GmComposite2003)", "Start(cM)", "End(cM)", "Associated Loci", "R-Squared")
AllThreeLit_reduced <- unique(AllThreeLit_reduced)
rownames(AllThreeLit_reduced) <- NULL

# Top QTL Categories
TopQTL_ByCat <- QTLCounts_ProtOilYield %>% group_by(Object.Type) %>% top_n(5, GroupCount) %>% ungroup()

# Create a table with counts of published QTL in 5cM intervals on the soybean consensus map
QTL_Plotting <- QTL_Types
QTL_Plotting$LG <- as.factor(QTL_Plotting$LG)


# Position of last detected QTL for each LG
LastQTLPos <- QTL_Plotting %>% group_by(LG) %>% summarise(MaxLen = max(Stop.cM)) %>% ungroup
LastQTLPos$EndRound <- ceiling(LastQTLPos$MaxLen / 5)*5

LastQTLPos <- LastQTLPos %>% select(LG, EndRound)
LastQTLPos$LG <- as.character(LastQTLPos$LG)

QTLPos_Plotting <- left_join(QTL_Plotting, LastQTLPos, by = "LG")

QTLPos_Plotting$Start.cM <- round(QTLPos_Plotting$Start.cM, 0)
QTLPos_Plotting$Stop.cM  <- round(QTLPos_Plotting$Stop.cM, 0)

# Create a list of numeric vectors in cm steps of 1 to hold QTL counts
OilPos <- vector("list", length = length(unique(QTLPos_Plotting$LG)))
names(OilPos) <- unique(QTLPos_Plotting$LG)

for(i in names(OilPos)){
  OilPos[[i]] <- rep(0, LastQTLPos$EndRound[match(i, LastQTLPos$LG)])
}

YieldPos <- ProtPos <- OilPos


for(i in 1:nrow(QTLPos_Plotting)){
  CurrentObject <- QTLPos_Plotting$Object.Type[[i]]
  
  CurrentStartPos <- QTLPos_Plotting$Start.cM[[i]]
  CurrentEndPos   <- QTLPos_Plotting$Stop.cM[[i]]
  CurrentLG       <- as.character(QTLPos_Plotting$LG[[i]])
  
  if(CurrentObject == "QTL_oil"){
    OilPos[[CurrentLG]][CurrentStartPos:CurrentEndPos] <- OilPos[[CurrentLG]][CurrentStartPos:CurrentEndPos] + 1
    next
  }
  
  if(CurrentObject == "QTL_protein"){
    ProtPos[[CurrentLG]][CurrentStartPos:CurrentEndPos] <- ProtPos[[CurrentLG]][CurrentStartPos:CurrentEndPos] + 1
    next
  }
  
  if(CurrentObject == "QTL_yield"){
    YieldPos[[CurrentLG]][CurrentStartPos:CurrentEndPos] <- YieldPos[[CurrentLG]][CurrentStartPos:CurrentEndPos] + 1
  }
  
  
}

# A function to process this data
ProcessCountLists <- function(CountList, qtlName = "Oil"){
  CountDf <- data.frame(LG = rep(names(CountList), sapply(CountList, length)), cM = unlist(sapply(sapply(CountList, length), function(x) 1:x)), count = unlist(CountList))
  rownames(CountDf) <- NULL
  colnames(CountDf) <- c("LG", "cM", qtlName)
  CountDf
}

OilDF   <- ProcessCountLists(OilPos, qtlName = "Oil")
ProtDF  <- ProcessCountLists(ProtPos, qtlName = "Protein")
YieldDF <- ProcessCountLists(YieldPos, qtlName = "Yield")

AllDf <- bind_cols(OilDF, bind_cols(ProtDF, YieldDF))
colnames(AllDf)[1:2] <- c("LG", "cM")

AllDf <- AllDf %>% select(LG, cM, Oil, Protein, Yield)

AllDf_melt    <- melt(AllDf, measure.vars = c("Oil", "Protein", "Yield"))
ExpectedSizes <- read.csv("https://raw.githubusercontent.com/jhgille2/SoybaseData/master/ExpectedLGSizes.csv")

Factlevels    <- as.character(ExpectedSizes$LG)
AllDf_melt$LG <- factor(AllDf_melt$LG, levels = Factlevels)

loadfonts(device = "win")

# Subset for testing plotting
AllDf_melt_LG1 <- AllDf_melt %>% filter(LG == "D1a")
Chr1Plot <- ggplot(AllDf_melt_LG1, aes(fill = variable, x = cM, y = value)) + 
  geom_bar(stat = "identity", color = NA) + 
 # ggtitle("Chromosome D1a") + 
  theme_hc() + 
  theme(text = element_text(family = "LM Roman 10"), 
        plot.title = element_text(hjust = 0.5, size = 20), 
        legend.title = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 10)) + 
  ylab("Count of QTL")

# Counts of QTL Category by Chromosome
QTL_Types_bychr <- QTL_Types %>% group_by(LG, Object.Type) %>% summarise(Count = n())
QTL_Types_bychr$Object.Type <- capitalize(sub("QTL_", "", QTL_Types_bychr$Object.Type))
QTL_Types_bychr <- QTL_Types_bychr %>% spread(Object.Type, Count, fill = 0)

# Counts of confirmed QTL by category
ConfirmedQTL <- QTL_Types %>% filter(QTLCategory %in% c("cqSeed oil", "cqSeed protein", "cqSeed yield", "cqSeed weight"))

# Confirmed Lit Sources
ConfirmedStrings <- c("cqSeed oil", "cqSeed protein", "cqSeed yield", "cqSeed weight")
ConfirmedStrings <- paste(ConfirmedStrings, collapse = "|")

ConfirmedLit <- AllScrapes_Lit %>% filter(str_detect(QTL, ConfirmedStrings)) %>% unique() %>% arrange(QTL, ShortName)
ConfirmedLit$LitName <- NULL

ConfirmedLit <- ConfirmedLit %>% 
  select(ShortName, PubYear, QTL_Type, QTL, `LG(GmComposite2003)`, `Start(cM)`, `End(cM)`, AssociatedLoci, R2) %>% 
  arrange(PubYear, ShortName, QTL_Type, QTL, `LG(GmComposite2003)`) %>%
  select(ShortName, QTL_Type, QTL, `LG(GmComposite2003)`, `Start(cM)`, `End(cM)`, AssociatedLoci, R2)

ConfirmedLit$QTL_Type <- capitalize(sub("QTL_", "", ConfirmedLit$QTL_Type))
colnames(ConfirmedLit) <- c("Literature Source", "QTL Category", 'QTL Name', "LG(GmComposite2003)", "Start(cM)", "End(cM)", "Associated Loci", "R-Squared")


# Shorter version of the all literature dataframe, for export
AllScraped_reduced <- AllScrapes_Lit %>%   select(ShortName, PubYear, QTL_Type, QTL, `LG(GmComposite2003)`, `Start(cM)`, `End(cM)`, AssociatedLoci, R2) %>% 
  arrange(PubYear, ShortName, QTL_Type, QTL, `LG(GmComposite2003)`) %>%
  select(ShortName, QTL_Type, QTL, `LG(GmComposite2003)`, `Start(cM)`, `End(cM)`, AssociatedLoci, R2)

AllScraped_reduced$QTL_Type <- capitalize(sub("QTL_", "", AllScraped_reduced$QTL_Type))



colnames(AllScraped_reduced) <- c("Literature Source", "QTL Category", 'QTL Name', "LG(GmComposite2003)", "Start(cM)", "End(cM)", "Associated Loci", "R-Squared")

AllScraped_reduced <- unique(AllScraped_reduced)


# Number of publications in each category by year
AllScrapes_Lit_Plotting <- unique(AllScrapes_Lit)
AllScrapes_Lit_Plotting <-  AllScrapes_Lit_Plotting %>% 
  group_by(QTL_Type, PubYear) %>% 
  summarise(Count = n()) %>% spread(QTL_Type, Count, fill = 0) %>%
  ungroup() %>%
  melt(measure.vars = c("QTL_oil", 'QTL_protein', "QTL_yield"))

AllScrapes_Lit_Plotting$variable <- capitalize(sub("QTL_", "", AllScrapes_Lit_Plotting$variable ))

AllScrapes_Lit_Plotting <- AllScrapes_Lit_Plotting[-which(is.na(AllScrapes_Lit_Plotting$PubYear)), ]

QTLCountPlot_ByYear <- ggplot(AllScrapes_Lit_Plotting, aes(fill = variable, x = PubYear, y = value)) + 
  geom_bar(stat = "identity", color = NA) + 
  # ggtitle("Chromosome D1a") + 
  theme_hc() + 
  theme(text = element_text(family = "LM Roman 10"), 
        plot.title = element_text(hjust = 0.5, size = 20), 
        legend.title = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) + 
  scale_x_discrete("Year", limits = 1992:2018, labels = as.character(1992:2018)) + 
  ylab("Count of QTL") + 
  xlab("Year")
QTLCountPlot_ByYear

# Data for a line plot
AllScrapes_LinePlt <- AllScrapes_Lit_Plotting %>% 
  group_by(variable) %>% 
  mutate(TotalCount = cumsum(value)) %>% 
  ungroup()

BaseTheme <-   theme(legend.position = 'bottom',
                     text = element_text(family = 'LM Roman 10', margin = margin(0.5,0.5,0.5,0.5, unit = 'lines')),
                     legend.text = element_text(size = 12, face = 'bold'),
                     legend.title = element_blank(),
                     legend.text.align = 0.5,
                     legend.justification = 'center',
                     plot.title = element_text(size = 25, hjust = 0.5, margin = margin(20,0,0,0)),
                     plot.subtitle = element_text(size = 20, hjust = 0.5, margin = margin(10,0,20,0)),
                     panel.border = element_rect(colour = "black", fill=NA, size=1),
                     #legend.margin = margin(20,0,20,0),
                     legend.key.size = unit(2, 'lines')) 

# Total count of QTL lineplot
QTL.Lineplot <- ggplot(AllScrapes_LinePlt, aes(x = PubYear, y = TotalCount, colour = variable)) + 
  geom_line(lty = 4, size = 1) + 
  theme_bw() + 
  BaseTheme+ 
  xlab("Publication Year") + 
  ylab("Total QTL Count") +
  guides(colour = guide_legend(override.aes = list(size = 10))) +
  theme(axis.text.x = element_text(size = 15, angle = 0.90),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 17, margin = margin(0.5,0,0,0, 'lines')),
        axis.title.y = element_text(size = 17, margin = margin(0,0.5,0,0, 'lines')),
        legend.key = element_rect(fill = NA)) +
  labs(title    = "Total QTL listed on SoyBase",
       subtitle = "As of July, 2020") + 
  scale_x_continuous(breaks = round(seq(min(AllScrapes_LinePlt$PubYear), max(AllScrapes_LinePlt$PubYear), by = 5),1))

ggsave(filename = "QTLLineplot.svg",
       plot = QTL.Lineplot,
       device = 'svg',
       path = "C:\\Users\\Jay\\Desktop\\Oil MP Paper\\Presentation Images\\",
       width = 7,
       height = 6,
       units = 'in')



# Total number of publications per year
AllScrapes_PubYearCount <-  AllScrapes_Lit %>%
  unique() %>%
  group_by(PubYear) %>% 
  summarise(Count = n()) %>%
  ungroup()

# The directory where publication summary tables are stored
PaperSummaryDir <- "C:\\Users\\Jay\\Desktop\\Oil MP Paper\\paper summary tables\\"

# Directory to export latex tables to
TexDir <- "C:\\Users\\Jay\\Desktop\\Oil MP Paper\\Tex Tables\\"

# Directory where plots are stored
PlotDir <- "C:\\Users\\Jay\\Desktop\\Oil MP Paper\\Plots\\"

# Write QTL Counts, overall literature table to this directory
write.csv(AllScrapes_Lit, 
          paste(PaperSummaryDir, "ProtOilYieldLit.csv", sep = ""),
          row.names = FALSE)

write.csv(QTLCounts_ByPaper,
          paste(PaperSummaryDir, "ProtOilYieldQTLCounts.csv", sep = ""),
          row.names = FALSE)

write.csv(TopQTL_ByCat,
          paste(PaperSummaryDir, "TopQTLCat.csv", sep = ""),
          row.names = FALSE)


# Formatting latex tables
QTLcategory_tex <- kable(QTL_MajorCategoryCount, 
                         "latex", 
                         caption = 'Counts of QTL Discovered by Category', 
                         booktabs = TRUE) 
save_kable(QTLcategory_tex, file = paste(TexDir, "QTLCategory.pdf", sep = ""), keep_tex = TRUE)


AllThreeLit_reduced[is.na(AllThreeLit_reduced)] <- ''
AllThree_tex <- kable(AllThreeLit_reduced, "latex", caption = "Literature Identifying Oil, Protein, and Yield QTL Together",longtable = TRUE, align = paste(rep('c', ncol(AllThreeLit_reduced)), collapse = ''), row.names = FALSE) %>%
  #collapse_rows(columns = c(1, 2)) %>%
  kable_styling(latex_options = "repeat_header", font_size = 7)
save_kable(AllThree_tex, file = paste(TexDir, "ThreecategoryQTL.pdf", sep = ""), keep_tex = TRUE)


ConfirmedLit[is.na(ConfirmedLit)] <- ''
ConfirmedTex <- kable(ConfirmedLit, "latex", caption = "Literature Identifying Confirmed QTL for Seed Oil, Protein, and Yield traits",longtable = TRUE, align = paste(rep('c', ncol(ConfirmedLit)), collapse = ''), row.names = FALSE) %>%
  #collapse_rows(columns = c(1, 2)) %>%
  kable_styling(latex_options = "repeat_header", font_size = 7)
save_kable(ConfirmedTex, file = paste(TexDir, "ConfirmedQTL.pdf", sep = ""), keep_tex = TRUE)


AllScraped_reduced[is.na(AllScraped_reduced)] <- ''
AllLitTex <- kable(AllScraped_reduced, "latex", caption = "QTL Mapping Studies for Seed Oil, Protein, and Yield Traits",longtable = TRUE, align = paste(rep('c', ncol(ConfirmedLit)), collapse = ''), row.names = FALSE) %>%
  #collapse_rows(columns = c(1, 2)) %>%
  kable_styling(latex_options = "repeat_header", font_size = 7)
save_kable(AllLitTex, file = paste(TexDir, "AllLitQTL.pdf", sep = ""), keep_tex = TRUE)

  


# Saving plots

# Plot of QTL COunts by year
ggsave(filename = paste("QTLCountByYear.pdf", sep = ""),
       plot = QTLCountPlot_ByYear,
       device = "pdf",
       width = 6,
       height = 4,
       units = "in",
       dpi = 2000,
       path = PlotDir)


# Plot of QTL counts for each chromosome
for(i in unique(AllDf_melt$LG)){
  
  AllDf_melt_LG <- AllDf_melt %>% filter(LG == i)
  ChrPlot <- ggplot(AllDf_melt_LG, aes(fill = variable, x = cM, y = value)) + 
    geom_bar(stat = "identity", color = NA) + 
    # ggtitle("Chromosome D1a") + 
    theme_hc() + 
    theme(text         = element_text(family = "LM Roman 10"), 
          plot.title   = element_text(hjust = 0.5, size = 20), 
          legend.title = element_blank(),
          axis.text    = element_text(size = 10),
          axis.title   = element_text(size = 15),
          legend.text  = element_text(size = 10),
          panel.border = element_rect(colour = "black", fill=NA, size=1)) + 
    ylab("Count of QTL")
  
  ggsave(filename = paste(i, "QTLCount.pdf", sep = ""),
         plot = ChrPlot,
         device = "pdf",
         width = 6,
         height = 4,
         units = "in",
         dpi = 2000,
         path = PlotDir)
  

}
PresentationFolder <- "C:\\Users\\Jay\\Desktop\\Oil MP Paper\\Presentation Images\\"


AllDf_melt_20 <- AllDf_melt %>% filter(LG %in% c("I", "E", "B2"))

StableChrCOnvert <- list("B2" = "Chromosome 14",
                         "E"  = 'Chromosome 15',
                         "I"  = 'Chromosome 20') 

AllDf_melt_20$LG <- factor(as.character(StableChrCOnvert[as.character(AllDf_melt_20$LG)]), levels = as.character(StableChrCOnvert))

StableQTLPlot <- ggplot(AllDf_melt_20, aes(fill = variable, x = cM, y = value)) + 
  geom_bar(stat = "identity", size = 0.01, colour = NA) + 
  facet_wrap(~LG, nrow = 3) + 
  theme_bw() + 
  theme(legend.position    = "bottom",
        axis.text          = element_text(size = 45, face = 'bold', margin = margin(0,1,0,0, 'lines')),
        axis.title         = element_text(size = 55),
        legend.text        = element_text(size = 45, margin = margin(0, 1, 0, 0, 'lines')),
        legend.title       = element_blank(),
        legend.key.size    = unit(5,"line"),
        #aspect.ratio      = 1,
        text               = element_text(family = "LM Roman 10"),
        plot.title         = element_text(size = 65, hjust = 0.5),
        plot.subtitle      = element_text(size = 50, hjust = 0.5, margin = margin(0,0,1,0, 'lines')),
        strip.background   = element_blank(),
        strip.text.x       = element_text(size = 50, face = "bold"),
        panel.spacing      = unit(2, "lines"),
        panel.grid.major   = element_line(colour = 'black', linetype = 1),
        panel.grid.minor   = element_line(colour = 'black', linetype = 1),
        axis.title.y = element_text(margin = margin(0, 1, 0, 0, 'lines'))) + 
  scale_x_continuous(breaks = (seq(0, (max(AllDf_melt_20$cM) + 10), 20)))+
  ylab("Count of QTL") + 
  labs(title = "Previously found QTL", subtitle = "As of July 2020") + 
  scale_y_continuous(expand = c(0,0))




ggsave(filename = paste(PresentationFolder, "SatbleQTL.svg", sep = ''), 
       plot     = StableQTLPlot,
       device   = 'svg',
       dpi      = 1000,
       width    = 25,
       height   = 20,
       units    = "in")


# A different plot that puts all LGs in two images, 10 per page

First10 <- Factlevels[1:10]
Last10  <- Factlevels[11:20]

AllDf_melt1 <- AllDf_melt %>% filter(LG %in% First10)
AllDf_melt2 <- AllDf_melt %>% filter(LG %in% Last10)

LGPlot_Facet1 <- ggplot(AllDf_melt1, aes(fill = variable, x = cM, y = value)) + 
  geom_bar(stat = "identity", colour = NA) + 
  theme_hc() + 
  ylab("Count of QTL") + 
  facet_wrap(~LG, ncol = 5) + 
  theme(text         = element_text(family = "LM Roman 10"), 
        plot.title   = element_text(hjust = 0.5, size = 20), 
        legend.title = element_blank(),
        axis.text    = element_text(size = 10),
        axis.title   = element_text(size = 15),
        legend.text  = element_text(size = 10),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

LGPlot_Facet2 <- ggplot(AllDf_melt2, aes(fill = variable, x = cM, y = value)) + 
  geom_bar(stat = "identity", colour = NA) + 
  theme_hc() + 
  ylab("Count of QTL") + 
  facet_wrap(~LG, ncol = 5) + 
  theme(text         = element_text(family = "LM Roman 10"), 
        plot.title   = element_text(hjust = 0.5, size = 20), 
        legend.title = element_blank(),
        axis.text    = element_text(size = 10),
        axis.title   = element_text(size = 15),
        legend.text  = element_text(size = 10),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

ggsave(filename = 'AllLGCounts1.pdf',
       plot     = LGPlot_Facet1,
       device   = "pdf",
       width    = 11,
       height   = 8,
       units    = "in",
       dpi      = 8000,
       path     = PlotDir)

ggsave(filename = 'AllLGCounts2.pdf',
       plot     = LGPlot_Facet2,
       device   = "pdf",
       width    = 11,
       height   = 8,
       units    = "in",
       dpi      = 8000,
       path     = PlotDir)

# Alternatively, plot five groups of 4
for(i in seq(1, 20, 4)){
  
  # Filter data to the current selection
  AllDf_CurrentLevels    <- Factlevels[i:(i+3)]
  AllDf_CurrentSelection <- AllDf_melt %>% filter(LG %in% AllDf_CurrentLevels)
  
  CurrentFacetPlot <- ggplot(AllDf_CurrentSelection, aes(fill = variable, x = cM, y = value)) + 
    geom_bar(stat = "identity", colour = NA) + 
    theme_hc() + 
    ylab("Count of QTL") + 
    facet_wrap(~LG, ncol = 1) + 
    theme(text         = element_text(family = "LM Roman 10"), 
          plot.title   = element_text(hjust = 0.5, size = 20), 
          legend.title = element_blank(),
          axis.text    = element_text(size = 10),
          axis.title   = element_text(size = 15),
          legend.text  = element_text(size = 10),
          panel.border = element_rect(colour = "black", fill=NA, size=1))
  
  FileNam <- paste("Chr_", i, "_", i+3, '_facetPlot.pdf', sep = "")
  
  ggsave(filename = FileNam,
         plot     = CurrentFacetPlot,
         device   = "pdf",
         width    = 6,
         height   = 8,
         units    = "in",
         dpi      = 1000,
         path     = PlotDir)
}






# ERROR CORRECTION SECTION

# Fixing errors with scraping tables from extended data
# The problem is that sometimes the loci under the "Loci positively associated wih the QTL" div
# are stored as individual tables instead of being stored in a single table, as would be expected.

# A couple assumptions I could make are that is a page is to have additional tables, it will always have them in the positively associated QTl section
# and another is that the last table on the page will always be a reference to funding. Using these two assumptions,
# it may be possible to correct tables by removing extra rows from the associated QTL table, and then extra tables from the scrape
# according to the expected number of tables as provided by the '.head_link' selectors

Good_URL <- "https://www.soybase.org/sbt/search/search_results.php?category=QTLName&search_term=mqSeed+oil-013" # QTL with single table
Bad_URL  <- "https://www.soybase.org/sbt/search/search_results.php?category=QTLName&search_term=Seed+oil+2-7" # QTL stored in multiple tables


Test_URL     <- New_Url
Test_Data    <- read_html(Test_URL) %>% html_nodes("table") %>% html_table(fill = TRUE)
Test_Data2   <- Test_Data
Test_headers <- read_html(Test_URL) %>% html_nodes('.head_link') %>% html_text()
Test_headers <- c("General Information", Test_headers)

# Test-Data has 17 elements while Test_headers only has 12. Presumably Test_headers is the name of the tables, so where are the 
# extra elements coming from?
# Actually, I will make the process more precise
Test_page <- read_html(Test_URL) # Read in the page

# The first table has a different selector for its title than the rest of the tables
First_TableName  <- Test_page %>% html_nodes('.top_bar') %>% html_text(trim = TRUE) # Name of the first table
Other_TableNames <- Test_page %>% html_nodes('.head_link') %>% html_text(trim = TRUE) # Names for the other tables
All_TableNames   <- c(First_TableName, Other_TableNames) # All the table names

# In total, there should only be as many tables in the final results as there are distinct table names, but there are more.
# The culprit in this case is two-fold. First, the funding information is stored as a table, but it isn't really a table, and is also
# not very interesting to read. This data is stored in the last element of the table list so can be removed right away
Test_Data <- Test_Data[-length(Test_Data)]

# There are still 16 elements of Test_Data though, when only 12 are expected. Where are the extra tables coming from?
# In this case, the extra tables appear because tlements under the "loci positively associated with this QTL" section
# are (sometimes) each stored as a seperate table. I have not found the exact reason why, but this causes the scraper to 
# combine the remaining table text into a single long table, the first few rows of which belong to the associated QTL section
# and the table seems to be repeated the same number of times as there are repeated tables. 

# The assumption I will make here is that if there are repeated tables under the same heaser, they will appear only in this sectio
# using this logic, the table can be trimmed to an appropriate size, and extra tables eliminated from the data list.... I'll illustrate

# Which table holds the associated QTL
Assoc_Table <- which(All_TableNames == "Loci positively associated with the QTL")

# How many repeats to trim
TrimNum <- length(Test_Data) - length(All_TableNames)

# Indexes of tables to remove
Rem_Tables <-(Assoc_Table+1):(Assoc_Table+TrimNum)

# Trim association table to include only desired data
Test_Data[[Assoc_Table]] <- Test_Data[[Assoc_Table]][1:TrimNum, ]
Test_Data                <- Test_Data[-Rem_Tables]

names(Test_Data) <- Test_headers
Test_Data        <- Test_Data[!is.na(names(Test_Data))]

# Add names
names(Test_Data) <- All_TableNames
Test_Data

# Works for this case....How about another?
# Another 'bad' url
New_Url <- "https://www.soybase.org/sbt/search/search_results.php?category=QTLName&search_term=Seed+oil+6-3"

# Seems like this strategy works (for now), I'll go ahead and implement it in the main function

# Another problem, Sometimes the function will try to scrape an empty table
# and throw an error, I'll try to counter this with the 'possibly' function
# from purrr

# See this link https://community.rstudio.com/t/how-to-skip-empty-table-while-scraping-with-rvest-html-table/27597/2

scrape_table <- function(ScrapedPage, xpath){
  ScrapedPage %>%
    html_nodes(xpath = xpath) %>%
    html_table()
}

# This url has an empty table in it
Empty_url  <- "https://www.soybase.org/sbt/search/search_results.php?category=QTLName&search_term=Seed+oil+6-3"
Empty_Page <- read_html(Empty_url)

# I will need to modify the code from the webpage a bit
nTables <- read_html(Empty_url) %>% html_nodes('table') %>% length()
xpaths  <- paste("/html/body/div[3]/div/table[", 1:nTables, "]", sep = "")

scrape_table_possibly <- possibly(scrape_table, otherwise = NULL)
scraped_tables        <- map(xpaths, ~ scrape_table_possibly(ScrapedPage = Empty_Page, xpath = .x))

Oats.lmer <- lme4::lmer(yield ~ 1 + (1|Variety) + factor(nitro) + (1|Block/Variety), data = nlme::Oats)
