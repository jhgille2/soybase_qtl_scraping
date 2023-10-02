# Install pacman if it is not already installed
if (!require("pacman")) install.packages("pacman")

# Load required packages
pacman::p_load(conflicted,
               dotenv,
               targets,
               tarchetypes,
               here,
               fnmate)

# Data wrangling
pacman::p_load(dplyr,
               stringr,
               purrr,
               tidyr,
               Hmisc,
               reshape2,
               rvest)

# Plotting
pacman::p_load(ggpubr,
               ggthemes)

# Reporting
pacman::p_load(kableExtra,
               extrafont)

# Misc
pacman::p_load(pbapply)
