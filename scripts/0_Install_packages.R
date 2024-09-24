# To install all of the packages we will need: select all of the lines of code and run them at the same time 

# install pacman - this is a library manager
install.packages('pacman')

# load all of the packages that are on CRAN
pacman::p_load('remotes', 
               'tidyverse', 'CoordinateCleaner', 
               'sf', 'terra',
               'rgbif', 'rnaturalearth', 'rnaturalearthdata')

# load all of the packages that are on GitHub
# pacman::p_load_gh('sjevelazco/flexsdm', 'rlesur/klippy')