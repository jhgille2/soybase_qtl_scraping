#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param soybase_qtl_df
#' @return
#' @author Jay Gillenwater
#' @export
get_qtl_types <- function(soybase_qtl_df) {

  AllTraits <- soybase_qtl_df %>% 
    pluck("QTL_positions") %>% 
    pluck("QTLCategory") %>% 
    gsub(" ", "+", .) %>% 
    gsub(",", "", .) %>%
    unique() 
  
  trait_tbl <- tibble(TraitName = AllTraits)
  
  return(trait_tbl)
}
