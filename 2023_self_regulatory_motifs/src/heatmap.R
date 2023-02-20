library(dplyr)
library(ggplot2)
library(plotly)

SwissProt_results <- read.csv("~/2023_self_regulatory_motifs/SwissProt_results.csv")

df <- SwissProt_results %>%
  group_by(Domain, Distance.Domain.Motif) %>%
  summarize(Ocurrences = n())

df <- df %>% mutate(text = paste('Peptide at', Distance.Domain.Motif, 'residues from', Domain, Ocurrences, 'times'))

df <- df[abs(df$Distance.Domain.Motif)>20,]

p <- ggplot(df, aes(x=Distance.Domain.Motif, y=Domain, fill=Ocurrences, text=text)) +
  geom_tile() +
  geom_vline(aes(xintercept=0), color='black', linetype='longdash') +
  geom_vline(aes(xintercept=-20), color='black', linetype='longdash', alpha=0.4) +
  geom_vline(aes(xintercept=20), color='black', linetype='longdash', alpha=0.4) +
  scale_fill_gradient(low = 'pink', high = 'darkred') +
  theme(panel.background = element_rect(fill='white', color = 'black'),
        text = element_text(size = 18))
p
#ggplotly(p, tooltip='text')

interesting.df <- SwissProt_results[(SwissProt_results$Domain %in% unique(df$Domain)) & abs(SwissProt_results$Distance.Domain.Motif)>20,]
write.csv(interesting.df, "~/2023_self_regulatory_motifs/interesting.all.data.csv", row.names=FALSE)

interesting.df <- interesting.df %>%
  group_by(Domain, Motif.seq) %>%
  summarise(ocurrences = n())
write.csv(interesting.df, "~/2023_self_regulatory_motifs/interesting.peptides.ocurrences.csv", row.names=FALSE)

for (dom in unique(interesting.df$Domain)){
  individual.df <- interesting.df[interesting.df$Domain==dom,]
  write.csv(individual.df, paste("~/2023_self_regulatory_motifs/interesting.", dom, ".csv", sep=''), row.names=FALSE)
  write.csv(individual.df, paste("~/2023_self_regulatory_motifs/presentation/interesting.", dom, ".csv", sep=''), row.names=FALSE)
}


###############################################3

Arabidopsis_thaliana_results <- read.csv("~/2023_self_regulatory_motifs/Arabidopsis_thaliana/Arabidopsis_thaliana_results.csv")

'All_results <- rbind(Arabidopsis_thaliana_results, Caenorhabditis_elegans_results, Candida_albicans_results,
                     Danio_rerio_results, Dictyostelium_discoideum_results, Drosophila_melanogaster_results,
                     Escherichia_coli_results, Glycine_max_results, Homo_sapiens_results,
                     Methanocaldococcus_jannaschii_results, Mus_musculus_results, Oryza_sativa_results,
                     Rattus_norvegicus_results, Saccharomyces_cerevisiae_results, Schizosaccharomyces_pombe_results,
                     Zea_mays_results)'

df <- Arabidopsis_thaliana_results %>%
  group_by(Domain, Distance.Domain.Motif) %>%
  summarize(ocurrences = n())

df <- df %>% mutate(text = paste('Peptide at', Distance.Domain.Motif, 'residues from', Domain, ocurrences, 'times'))

df <- df %>% mutate(distances2 =
                      case_when(Distance.Domain.Motif <= -100 ~ -100,
                                Distance.Domain.Motif >= 100 ~ 100,
                                Distance.Domain.Motif > -20 & ocurrences < 20 ~ 0,
                                .default = Distance.Domain.Motif))
df <- df %>%
  group_by(Domain, distances2) %>%
  summarize(ocurrences2 = n())

write.csv(df, "~/2023_self_regulatory_motifs/Arabidopsis_thaliana/Arabidopsis_thaliana_ocurrences.csv", row.names=FALSE)

p <- ggplot(df, aes(x=distances2, y=Domain, fill=ocurrences2, text=Domain)) +
  geom_tile() +
  scale_fill_gradient(low = 'pink', high = 'darkred') +
  theme(legend.position = 'None',
        panel.background = element_rect(fill='white', color = 'black'))
p
ggplotly(p, tooltip='text')
# Check this:
#htmlwidgets::saveWidget(p, 'Arabidopsis_heatmap.html')

Caenorhabditis_elegans_results <- read.csv("~/2023_self_regulatory_motifs/Caenorhabditis_elegans/Caenorhabditis_elegans_results.csv")

df <- Caenorhabditis_elegans_results %>%
  group_by(Domain, Distance.Domain.Motif) %>%
  summarize(ocurrences = n())

#df$text <- df$Domain
## Add a variable for the display text with domain name, number of ocurrences and position

write.csv(df, "~/2023_self_regulatory_motifs/Caenorhabditis_elegans/Caenorhabditis_elegans_ocurrences.csv", row.names=FALSE)

p <- ggplot(df, aes(x=Distance.Domain.Motif, y=Domain, fill=ocurrences, text=ocurrences)) +
  geom_tile() +
  scale_fill_gradient(low = 'pink', high = 'darkred') +
  theme(legend.position = 'None',
        panel.background = element_rect(fill='white', color = 'black'))
p
ggplotly(p, tooltip='text')
# Check this:
#htmlwidgets::saveWidget(p, 'Arabidopsis_heatmap.html')

Candida_albicans_results <- read.csv("~/2023_self_regulatory_motifs/Candida_albicans/Candida_albicans_results.csv")

df <- Candida_albicans_results %>%
  group_by(Domain, Distance.Domain.Motif) %>%
  summarize(ocurrences = n())

#df$text <- df$Domain
## Add a variable for the display text with domain name, number of ocurrences and position

write.csv(df, "~/2023_self_regulatory_motifs/Candida_albicans/Candida_albicans_ocurrences.csv", row.names=FALSE)

p <- ggplot(df, aes(x=Distance.Domain.Motif, y=Domain, fill=ocurrences, text=ocurrences)) +
  geom_tile() +
  scale_fill_gradient(low = 'pink', high = 'darkred') +
  theme(legend.position = 'None',
        panel.background = element_rect(fill='white', color = 'black'))
p
ggplotly(p, tooltip='text')
# Check this:
#htmlwidgets::saveWidget(p, 'Arabidopsis_heatmap.html')
Danio_rerio_results <- read.csv("~/2023_self_regulatory_motifs/Danio_rerio/Danio_rerio_results.csv")
df <- Danio_rerio_results %>%
  group_by(Domain, Distance.Domain.Motif) %>%
  summarize(ocurrences = n())

#df$text <- df$Domain
## Add a variable for the display text with domain name, number of ocurrences and position

write.csv(df, "~/2023_self_regulatory_motifs/Danio_rerio/Danio_rerio_ocurrences.csv", row.names=FALSE)

p <- ggplot(df, aes(x=Distance.Domain.Motif, y=Domain, fill=ocurrences, text=ocurrences)) +
  geom_tile() +
  scale_fill_gradient(low = 'pink', high = 'darkred') +
  theme(legend.position = 'None',
        panel.background = element_rect(fill='white', color = 'black'))
p
ggplotly(p, tooltip='text')
# Check this:
#htmlwidgets::saveWidget(p, 'Arabidopsis_heatmap.html')
Dictyostelium_discoideum_results <- read.csv("~/2023_self_regulatory_motifs/Dictyostelium_discoideum/Dictyostelium_discoideum_results.csv")
df <- Dictyostelium_discoideum_results %>%
  group_by(Domain, Distance.Domain.Motif) %>%
  summarize(ocurrences = n())

#df$text <- df$Domain
## Add a variable for the display text with domain name, number of ocurrences and position

write.csv(df, "~/2023_self_regulatory_motifs/Dictyostelium_discoideum/Dictyostelium_discoideum_ocurrences.csv", row.names=FALSE)

p <- ggplot(df, aes(x=Distance.Domain.Motif, y=Domain, fill=ocurrences, text=ocurrences)) +
  geom_tile() +
  scale_fill_gradient(low = 'pink', high = 'darkred') +
  theme(legend.position = 'None',
        panel.background = element_rect(fill='white', color = 'black'))
p
ggplotly(p, tooltip='text')
# Check this:
#htmlwidgets::saveWidget(p, 'Arabidopsis_heatmap.html')
Drosophila_melanogaster_results <- read.csv("~/2023_self_regulatory_motifs/Drosophila_melanogaster/Drosophila_melanogaster_results.csv")
df <- Drosophila_melanogaster_results %>%
  group_by(Domain, Distance.Domain.Motif) %>%
  summarize(ocurrences = n())

#df$text <- df$Domain
## Add a variable for the display text with domain name, number of ocurrences and position

write.csv(df, "~/2023_self_regulatory_motifs/Drosophila_melanogaster/Drosophila_melanogaster_ocurrences.csv", row.names=FALSE)

p <- ggplot(df, aes(x=Distance.Domain.Motif, y=Domain, fill=ocurrences, text=ocurrences)) +
  geom_tile() +
  scale_fill_gradient(low = 'pink', high = 'darkred') +
  theme(legend.position = 'None',
        panel.background = element_rect(fill='white', color = 'black'))
p
ggplotly(p, tooltip='text')
# Check this:
#htmlwidgets::saveWidget(p, 'Arabidopsis_heatmap.html')
Escherichia_coli_results <- read.csv("~/2023_self_regulatory_motifs/Escherichia_coli/Escherichia_coli_results.csv")
df <- Escherichia_coli_results %>%
  group_by(Domain, Distance.Domain.Motif) %>%
  summarize(ocurrences = n())

#df$text <- df$Domain
## Add a variable for the display text with domain name, number of ocurrences and position

write.csv(df, "~/2023_self_regulatory_motifs/Escherichia_coli/Escherichia_coli_ocurrences.csv", row.names=FALSE)

p <- ggplot(df, aes(x=Distance.Domain.Motif, y=Domain, fill=ocurrences, text=ocurrences)) +
  geom_tile() +
  scale_fill_gradient(low = 'pink', high = 'darkred') +
  theme(legend.position = 'None',
        panel.background = element_rect(fill='white', color = 'black'))
p
ggplotly(p, tooltip='text')
# Check this:
#htmlwidgets::saveWidget(p, 'Arabidopsis_heatmap.html')
Glycine_max_results <- read.csv("~/2023_self_regulatory_motifs/Glycine_max/Glycine_max_results.csv")
df <- Glycine_max_results %>%
  group_by(Domain, Distance.Domain.Motif) %>%
  summarize(ocurrences = n())

#df$text <- df$Domain
## Add a variable for the display text with domain name, number of ocurrences and position

write.csv(df, "~/2023_self_regulatory_motifs/Glycine_max/Glycine_max_ocurrences.csv", row.names=FALSE)

p <- ggplot(df, aes(x=Distance.Domain.Motif, y=Domain, fill=ocurrences, text=ocurrences)) +
  geom_tile() +
  scale_fill_gradient(low = 'pink', high = 'darkred') +
  theme(legend.position = 'None',
        panel.background = element_rect(fill='white', color = 'black'))
p
ggplotly(p, tooltip='text')
# Check this:
#htmlwidgets::saveWidget(p, 'Arabidopsis_heatmap.html')
Homo_sapiens_results <- read.csv("~/2023_self_regulatory_motifs/Homo_sapiens/Homo_sapiens_results.csv")
df <- Homo_sapiens_results %>%
  group_by(Domain, Distance.Domain.Motif) %>%
  summarize(ocurrences = n())

#df$text <- df$Domain
## Add a variable for the display text with domain name, number of ocurrences and position

write.csv(df, "~/2023_self_regulatory_motifs/Homo_sapiens/Homo_sapiens_ocurrences.csv", row.names=FALSE)

p <- ggplot(df, aes(x=Distance.Domain.Motif, y=Domain, fill=ocurrences, text=ocurrences)) +
  geom_tile() +
  scale_fill_gradient(low = 'pink', high = 'darkred') +
  theme(legend.position = 'None',
        panel.background = element_rect(fill='white', color = 'black'))
p
ggplotly(p, tooltip='text')
# Check this:
#htmlwidgets::saveWidget(p, 'Arabidopsis_heatmap.html')
Methanocaldococcus_jannaschii_results <- read.csv("~/2023_self_regulatory_motifs/Methanocaldococcus_jannaschii/Methanocaldococcus_jannaschii_results.csv")
df <- Methanocaldococcus_jannaschii_results %>%
  group_by(Domain, Distance.Domain.Motif) %>%
  summarize(ocurrences = n())

#df$text <- df$Domain
## Add a variable for the display text with domain name, number of ocurrences and position

write.csv(df, "~/2023_self_regulatory_motifs/Methanocaldococcus_jannaschii/Methanocaldococcus_jannaschii_ocurrences.csv", row.names=FALSE)

p <- ggplot(df, aes(x=Distance.Domain.Motif, y=Domain, fill=ocurrences, text=ocurrences)) +
  geom_tile() +
  scale_fill_gradient(low = 'pink', high = 'darkred') +
  theme(legend.position = 'None',
        panel.background = element_rect(fill='white', color = 'black'))
p
ggplotly(p, tooltip='text')
# Check this:
#htmlwidgets::saveWidget(p, 'Arabidopsis_heatmap.html')
Mus_musculus_results <- read.csv("~/2023_self_regulatory_motifs/Mus_musculus/Mus_musculus_results.csv")
df <- Mus_musculus_results %>%
  group_by(Domain, Distance.Domain.Motif) %>%
  summarize(ocurrences = n())

#df$text <- df$Domain
## Add a variable for the display text with domain name, number of ocurrences and position

write.csv(df, "~/2023_self_regulatory_motifs/Mus_musculus/Mus_musculus_ocurrences.csv", row.names=FALSE)

p <- ggplot(df, aes(x=Distance.Domain.Motif, y=Domain, fill=ocurrences, text=ocurrences)) +
  geom_tile() +
  scale_fill_gradient(low = 'pink', high = 'darkred') +
  theme(legend.position = 'None',
        panel.background = element_rect(fill='white', color = 'black'))
p
ggplotly(p, tooltip='text')
# Check this:
#htmlwidgets::saveWidget(p, 'Arabidopsis_heatmap.html')
Oryza_sativa_results <- read.csv("~/2023_self_regulatory_motifs/Oryza_sativa/Oryza_sativa_results.csv")
df <- Oryza_sativa_results %>%
  group_by(Domain, Distance.Domain.Motif) %>%
  summarize(ocurrences = n())

#df$text <- df$Domain
## Add a variable for the display text with domain name, number of ocurrences and position

write.csv(df, "~/2023_self_regulatory_motifs/Oryza_sativa/Oryza_sativa_ocurrences.csv", row.names=FALSE)

p <- ggplot(df, aes(x=Distance.Domain.Motif, y=Domain, fill=ocurrences, text=ocurrences)) +
  geom_tile() +
  scale_fill_gradient(low = 'pink', high = 'darkred') +
  theme(legend.position = 'None',
        panel.background = element_rect(fill='white', color = 'black'))
p
ggplotly(p, tooltip='text')
# Check this:
#htmlwidgets::saveWidget(p, 'Arabidopsis_heatmap.html')
Rattus_norvegicus_results <- read.csv("~/2023_self_regulatory_motifs/Rattus_norvegicus/Rattus_norvegicus_results.csv")
df <- Rattus_norvegicus_results %>%
  group_by(Domain, Distance.Domain.Motif) %>%
  summarize(ocurrences = n())

#df$text <- df$Domain
## Add a variable for the display text with domain name, number of ocurrences and position

write.csv(df, "~/2023_self_regulatory_motifs/Rattus_norvegicus/Rattus_norvegicus_ocurrences.csv", row.names=FALSE)

p <- ggplot(df, aes(x=Distance.Domain.Motif, y=Domain, fill=ocurrences, text=ocurrences)) +
  geom_tile() +
  scale_fill_gradient(low = 'pink', high = 'darkred') +
  theme(legend.position = 'None',
        panel.background = element_rect(fill='white', color = 'black'))
p
ggplotly(p, tooltip='text')
# Check this:
#htmlwidgets::saveWidget(p, 'Arabidopsis_heatmap.html')
Saccharomyces_cerevisiae_results <- read.csv("~/2023_self_regulatory_motifs/Saccharomyces_cerevisiae/Saccharomyces_cerevisiae_results.csv")
df <- Saccharomyces_cerevisiae_results %>%
  group_by(Domain, Distance.Domain.Motif) %>%
  summarize(ocurrences = n())

#df$text <- df$Domain
## Add a variable for the display text with domain name, number of ocurrences and position

write.csv(df, "~/2023_self_regulatory_motifs/Saccharomyces_cerevisiae/Saccharomyces_cerevisiae_ocurrences.csv", row.names=FALSE)

p <- ggplot(df, aes(x=Distance.Domain.Motif, y=Domain, fill=ocurrences, text=ocurrences)) +
  geom_tile() +
  scale_fill_gradient(low = 'pink', high = 'darkred') +
  theme(legend.position = 'None',
        panel.background = element_rect(fill='white', color = 'black'))
p
ggplotly(p, tooltip='text')
# Check this:
#htmlwidgets::saveWidget(p, 'Arabidopsis_heatmap.html')
Schizosaccharomyces_pombe_results <- read.csv("~/2023_self_regulatory_motifs/Schizosaccharomyces_pombe/Schizosaccharomyces_pombe_results.csv")
df <- Schizosaccharomyces_pombe_results %>%
  group_by(Domain, Distance.Domain.Motif) %>%
  summarize(ocurrences = n())

#df$text <- df$Domain
## Add a variable for the display text with domain name, number of ocurrences and position

write.csv(df, "~/2023_self_regulatory_motifs/Schizosaccharomyces_pombe/Schizosaccharomyces_pombe_ocurrences.csv", row.names=FALSE)

p <- ggplot(df, aes(x=Distance.Domain.Motif, y=Domain, fill=ocurrences, text=ocurrences)) +
  geom_tile() +
  scale_fill_gradient(low = 'pink', high = 'darkred') +
  theme(legend.position = 'None',
        panel.background = element_rect(fill='white', color = 'black'))
p
ggplotly(p, tooltip='text')
# Check this:
#htmlwidgets::saveWidget(p, 'Arabidopsis_heatmap.html')
Zea_mays_results <- read.csv("~/2023_self_regulatory_motifs/Zea_mays/Zea_mays_results.csv")
df <- Zea_mays_results %>%
  group_by(Domain, Distance.Domain.Motif) %>%
  summarize(ocurrences = n())

#df$text <- df$Domain
## Add a variable for the display text with domain name, number of ocurrences and position

write.csv(df, "~/2023_self_regulatory_motifs/Zea_mays/Zea_mays_ocurrences.csv", row.names=FALSE)

p <- ggplot(df, aes(x=Distance.Domain.Motif, y=Domain, fill=ocurrences, text=ocurrences)) +
  geom_tile() +
  scale_fill_gradient(low = 'pink', high = 'darkred') +
  theme(legend.position = 'None',
        panel.background = element_rect(fill='white', color = 'black'))
p
ggplotly(p, tooltip='text')
# Check this:
#htmlwidgets::saveWidget(p, 'Arabidopsis_heatmap.html')