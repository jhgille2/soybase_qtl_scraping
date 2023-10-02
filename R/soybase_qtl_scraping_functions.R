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
    Sys.sleep(3) # Pause system to prevent spamming the website
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