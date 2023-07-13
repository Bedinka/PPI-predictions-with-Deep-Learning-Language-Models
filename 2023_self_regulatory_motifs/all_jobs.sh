#!/bin/bash

dirs=( '/home/roger/2023_self_regulatory_motifs/SwissProt' '/home/roger/2023_self_regulatory_motifs/Arabidopsis_thaliana/' '/home/roger/2023_self_regulatory_motifs/Caenorhabditis_elegans/' '/home/roger/2023_self_regulatory_motifs/Candida_albicans/' '/home/roger/2023_self_regulatory_motifs/Danio_rerio/' '/home/roger/2023_self_regulatory_motifs/Dictyostelium_discoideum/' '/home/roger/2023_self_regulatory_motifs/Drosophila_melanogaster/' '/home/roger/2023_self_regulatory_motifs/Escherichia_coli/' '/home/roger/2023_self_regulatory_motifs/Glycine_max/' '/home/roger/2023_self_regulatory_motifs/Homo_sapiens/' '/home/roger/2023_self_regulatory_motifs/Methanocaldococcus_jannaschii/' '/home/roger/2023_self_regulatory_motifs/Mus_musculus/' '/home/roger/2023_self_regulatory_motifs/Oryza_sativa/' '/home/roger/2023_self_regulatory_motifs/Rattus_norvegicus/' '/home/roger/2023_self_regulatory_motifs/Saccharomyces_cerevisiae/' '/home/roger/2023_self_regulatory_motifs/Schizosaccharomyces_pombe/' '/home/roger/2023_self_regulatory_motifs/Zea_mays/')


for dir in "${dirs[@]}"
do
    #cd $dir
    #pwd
    #cd ./compressed_D
    #find . -name '*cif.gz' -delete
    #cp ./compressed_D/* ./data/
    #cd data
    #gunzip -f *
    ./steps.sh $dir &
done
echo 'ALL FINISHED'
# cp /home/roger/2023_self_regulatory_motifs/SwissProt/SwissProt_results.csv /home/roger/2023_self_regulatory_motifs/SwissProt/SwissProt_results_1.8.csv
# cp /home/roger/2023_self_regulatory_motifs/SwissProt/SwissProt_results_1.8.csv ../ppi_web/SwissProt_results_1.8.csv
# ../ppi_web/mysite/manage.py import_csv ../SwissProt_results_1.8.csv ../uniprot_sprot.csv