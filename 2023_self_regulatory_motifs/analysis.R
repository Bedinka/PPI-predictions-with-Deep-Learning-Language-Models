library(dplyr)
library(ggplot2)
library(plotly)
library(stringr)



for (i in SwissProt_results_1.6$Motif.seq) {
  if (nchar(i) > nchar(longest)){
    longest = i
    }
}
longest

SwissProt_results <- read.csv("~/2023_self_regulatory_motifs/SwissProt/SwissProt_results.csv")

SwissProt_results_1.5 <- read.csv("~/2023_self_regulatory_motifs/SwissProt/SwissProt_results_1.5.csv")

SwissProt_results_1.6 <- read.csv("~/2023_self_regulatory_motifs/SwissProt/SwissProt_results_1.7.csv")

autoinhibitory_sprot <- read.csv("~/2023_self_regulatory_motifs/autoinhibitory_sprot.csv")

for (r in 1:length(SwissProt_results_1.6$PDB)) {
  string = SwissProt_results_1.6$PDB[r]
  match <- str_extract(string, "(?<=AF-)[^-]+")
  if (match %in% autoinhibitory_sprot$Accession){
    print(match)
  }
}

length(setdiff(SwissProt_results_1.5$PDB, SwissProt_results_1.6$PDB))
length(setdiff(SwissProt_results_1.6$PDB, SwissProt_results_1.5$PDB))

df <- SwissProt_results %>%
  group_by(Domain, Distance.Domain.Motif) %>%
  summarize(Ocurrences = n())

#df <- df %>% mutate(text = paste('Peptide at', Distance.Domain.Motif, 'residues from', Domain, Ocurrences, 'times'))

df <- df[abs(df$Distance.Domain.Motif)>20,]

df2 <- SwissProt_results %>%
  group_by(Domain) %>%
  summarize(ocr = n())

df <- df[df$Domain %in% df2[df2$ocr > 2000, ]$Domain, ]

p <- ggplot(df, aes(x=Distance.Domain.Motif, y=Domain, fill=Ocurrences))+#, text=text)) +
  geom_tile() +
  geom_vline(aes(xintercept=0), color='black', linetype='longdash') +
  geom_vline(aes(xintercept=-20), color='black', linetype='longdash', alpha=0.4) +
  geom_vline(aes(xintercept=20), color='black', linetype='longdash', alpha=0.4) +
  scale_fill_gradient(low = 'pink', high = 'darkred') +
  theme(panel.background = element_rect(fill='white', color = 'black'),
        text = element_text(size = 18))
p


