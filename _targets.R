## Load your packages, e.g. library(targets).
source("./packages.R")

## Load your R files
lapply(list.files("./R", full.names = TRUE), source)

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

# The traits which I want meta-data for
QTL_Types <- QTLPos %>% dplyr::filter(Object.Type %in% c("QTL_oil", "QTL_protein", "QTL_yield"))
QTL_Types$LG <- as.character(QTL_Types$LG)
QTL_Types$LG <- sub("D1A", "D1a", QTL_Types$LG)
QTL_Types$LG <- sub("D1B", "D1b", QTL_Types$LG)

# I've already done this part, skip to the next section to load in the file
AllTraits <- gsub(" ", "+", QTL_Types$QTLCategory) %>% unique()
AllTraits_tbl <- tibble(TraitName = AllTraits)

# # What traits do I want QTL data for?
# qtl_traits <- c("Seed+protein",
#                 "Seed+oil",
#                 "Seed+weight",
#                 "Seed+yield")
# 
# # Make a table that has every unique trait on soybase
# AllTraits_tbl <- here::here("data", "SoybaseQTL_2023_09_08.csv") %>% 
#   clean_soybase_qtls() %>% 
#   get_qtl_types() %>%
#   dplyr::filter(TraitName %in% qtl_traits)

## tar_plan supports drake-style targets and also tar_target()
tar_plan(
  
  # Scrape QTL results for every unique trait
  trait_scrapes <- tar_map(
    values = AllTraits_tbl, 
    
    tar_target(TraitScrape,
               Soybase_QTL_Info(TraitName))
    
  ),
  
  # Combine all the scrapes into a list
  tar_combine(AllScrapes,
              trait_scrapes[[1]],
              command = list(!!!.x))

)
