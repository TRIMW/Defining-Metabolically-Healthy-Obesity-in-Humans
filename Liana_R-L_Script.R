library(tidyverse)
library(magrittr)
library(liana)

show_resources()
show_methods()

##Downsample larger files if not connected to the cluster
OM_Small <- subset(x = OM, downsample = 500)

##Set up liana
liana_path <- system.file(package = "liana")

##Run Liana
liana_OM_MUO <- liana_wrap(OM_MUO_All)

# Liana returns a list of results, each element of which corresponds to a method
liana_OM_MUO %>% dplyr::glimpse()

# We can aggregate these results into a tibble with consensus ranks
liana_OM_MUO <- liana_OM_MUO %>%
  liana_aggregate()

dplyr::glimpse(liana_OM_MUO)

liana_OM_MUO %>%
  liana_dotplot(source_groups = c('Myeloid cells'),
                target_groups = c('PDGFRA+', 'PI16+', 'ECs', 'SMCs', 'NK cells', 'T cells', 'B cells', 'Plasma cells', 'MSLN+'),
                ntop = 20)

liana_trunc <- liana_OM_MUO %>%
  # only keep interactions concordant between methods
  filter(aggregate_rank <= 0.05) # note that these pvals are already corrected


p <- chord_freq(liana_trunc,
                source_groups = c('PDGFRA+', 'PI16+', 'ECs', 'SMCs', 'NK cells', 'T cells', 'B cells', 'Plasma cells', 'Myeloid cells'),
                target_groups = c('MSLN+'))
