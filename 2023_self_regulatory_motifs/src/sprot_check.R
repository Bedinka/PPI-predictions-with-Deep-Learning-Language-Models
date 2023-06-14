autoinhibitory_sprot <- read.csv("~/2023_self_regulatory_motifs/autoinhibitory_sprot.csv")
sprot_found <- read.csv("/media/storage1/2023_self_regulatory_motifs/sprot_found.csv", row.names=NULL)
sprot_not_found <- read.csv("/media/storage1/2023_self_regulatory_motifs/sprot_not_found.csv", row.names=NULL)


for (i in unique(sprot_found$Domain.Name)) {
  print(paste(sum(sprot_found$Domain.Name==i), i))
}
