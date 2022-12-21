module load R
module load python

pip install -r requirements.txt

halla -y species_halla.txt -x metab_halla.txt -o specmetab_bb3 -m spearman --fdr_alpha 0.05 

hallagram \
    -i specmetab_bb3/ \
    --cbar_label 'Pairwise Spearman' \
    --x_dataset_label 'Metabolites' \
    --y_dataset_label 'Species' \
    -n 20 \
    --block_border_width 2

hallagnostic -i specmetab_bb3/