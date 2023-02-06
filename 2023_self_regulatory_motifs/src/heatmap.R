library(dplyr)
library(ggplot2)
library(plotly)

Arabidopsis_thaliana_results <- read.csv("~/2023_self_regulatory_motifs/Arabidopsis_thaliana/Arabidopsis_thaliana_results.csv")
Caenorhabditis_elegans_results <- read.csv("~/2023_self_regulatory_motifs/Caenorhabditis_elegans/Caenorhabditis_elegans_results.csv")
Candida_albicans_results <- read.csv("~/2023_self_regulatory_motifs/Candida_albicans/Candida_albicans_results.csv")
Danio_rerio_results <- read.csv("~/2023_self_regulatory_motifs/Danio_rerio/Danio_rerio_results.csv")
Dictyostelium_discoideum_results <- read.csv("~/2023_self_regulatory_motifs/Dictyostelium_discoideum/Dictyostelium_discoideum_results.csv")
Drosophila_melanogaster_results <- read.csv("~/2023_self_regulatory_motifs/Drosophila_melanogaster/Drosophila_melanogaster_results.csv")
Escherichia_coli_results <- read.csv("~/2023_self_regulatory_motifs/Escherichia_coli/Escherichia_coli_results.csv")
Glycine_max_results <- read.csv("~/2023_self_regulatory_motifs/Glycine_max/Glycine_max_results.csv")
Homo_sapiens_results <- read.csv("~/2023_self_regulatory_motifs/Homo_sapiens/Homo_sapiens_results.csv")
Methanocaldococcus_jannaschii_results <- read.csv("~/2023_self_regulatory_motifs/Methanocaldococcus_jannaschii/Methanocaldococcus_jannaschii_results.csv")
Mus_musculus_results <- read.csv("~/2023_self_regulatory_motifs/Mus_musculus/Mus_musculus_results.csv")
Oryza_sativa_results <- read.csv("~/2023_self_regulatory_motifs/Oryza_sativa/Oryza_sativa_results.csv")
Rattus_norvegicus_results <- read.csv("~/2023_self_regulatory_motifs/Rattus_norvegicus/Rattus_norvegicus_results.csv")
Saccharomyces_cerevisiae_results <- read.csv("~/2023_self_regulatory_motifs/Saccharomyces_cerevisiae/Saccharomyces_cerevisiae_results.csv")
Schizosaccharomyces_pombe_results <- read.csv("~/2023_self_regulatory_motifs/Schizosaccharomyces_pombe/Schizosaccharomyces_pombe_results.csv")
Zea_mays_results <- read.csv("~/2023_self_regulatory_motifs/Zea_mays/Zea_mays_results.csv")

'All_results <- rbind(Arabidopsis_thaliana_results, Caenorhabditis_elegans_results, Candida_albicans_results,
                     Danio_rerio_results, Dictyostelium_discoideum_results, Drosophila_melanogaster_results,
                     Escherichia_coli_results, Glycine_max_results, Homo_sapiens_results,
                     Methanocaldococcus_jannaschii_results, Mus_musculus_results, Oryza_sativa_results,
                     Rattus_norvegicus_results, Saccharomyces_cerevisiae_results, Schizosaccharomyces_pombe_results,
                     Zea_mays_results)'

#Arabidopsis_thaliana_results <- data.frame(Distance.Domain.Motif=c(500, -200, 3, 1000), Domain=c('dreta', 'esquerra', 'aprop', 'vfar'))

df <- Arabidopsis_thaliana_results %>%
  group_by(Domain, Distance.Domain.Motif) %>%
  summarize(ocurrences = n())

#df$text <- df$Domain
## Add a variable for the display text with domain name, number of ocurrences and position

p <- ggplot(df, aes(x=Distance.Domain.Motif, y=Domain, fill=ocurrences, text=ocurrences)) +
  geom_tile() +
  scale_fill_gradient(low = 'pink', high = 'darkred') +
  theme(legend.position = 'None',
        panel.background = element_rect(fill='white', color = 'black'))
p
ggplotly(p, tooltip='text')
# Check this:
#htmlwidgets::saveWidget(p, 'Arabidopsis_heatmap.html')

unique(Arabidopsis_thaliana_results$Domain)
max(Arabidopsis_thaliana_results$Distance.Domain.Motif)
min(Arabidopsis_thaliana_results$Distance.Domain.Motif)
