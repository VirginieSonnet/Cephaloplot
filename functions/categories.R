###################################################
# Function to assign a category to the predictors #
###################################################

predictor_categories <- function(predictors){
  # --- Function to extract the input or selected predictors ---
  # Inputs: 
  #       - predictors: tibble, dataset with column "variable" containing the predictors 
  # Outputs: 
  #       - predictors: tibble, added column "category"
  
  predictors <- predictors %>% 
    mutate(name=tolower(variable)) %>% 
    mutate(category=case_when(str_detect(name,"eke|fsle")~"Courantology",
                              str_detect(name,"ice")~"Ice",
                              str_detect(name,"chla|pp|tot_poc|diato|dino|pico|prochlo|prokar|green|micro|nano|hapto")~"Productivity",
                              str_detect(name,"k490|bbp|aph|adg")~"Optics",
                              str_detect(name,"sfco2|spco2|co2|co3|hco3|dic|talk|omega|revelle|s_|s_m")~"Carbonate",
                              str_detect(name,"m_|t_|sst|par|ph_|a_200_300|o_")~"Latitude",
                              str_detect(name,"n_|p_|i_|a_0_50|a_m|i_m|n_m|p_m")~"Chemistry",
                              TRUE ~ "Uncategorized"
                              )) %>% 
    dplyr::select(-name)
  return(predictors)
}
  
