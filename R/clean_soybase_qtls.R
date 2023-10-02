#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param soybase_qtl_positions_file
#' @return
#' @author Jay Gillenwater
#' @export
clean_soybase_qtls <- function(soybase_qtl_positions_file) {

  # Read in the soybase qtl positions file
  QTLPos <- read.csv(soybase_qtl_positions_file)
  
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
  
  # Combine both dataframes in a list and return the list
  res <- list("QTL_positions"      = QTLPos,
              "QTL_NameConversion" = QTLNameConversion)
  
  return(res)
}
