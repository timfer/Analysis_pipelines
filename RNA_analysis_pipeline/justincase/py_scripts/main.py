#%%

import os
import platform

import pandas as pd

from functions.dir_utils import create_dir
from functions.preprocessing import annotate_doublets_enhanced
#%%
userID = "tferrari"
sampleIDs = ["NG037", "NG039", "NG042", "NG050", "TF049"]

if platform.system() == 'Windows':
    paths = {
        'iris_mPath': r'L:\projects\iris',
        'sig_mPath': r'L:\users\tferrari\Experiments\RNA_analysis\IRIS\CHUV_samples\signatures',
        'out_mPath': r'L:\users\tferrari\Experiments\RNA_analysis\IRIS\CHUV_samples\results\HD_controls\all_samples',
    }
else:
    paths = {
        'iris_mPath': '/home/tferrari/updepla/projects/iris',
        'sig_mPath': '/home/tferrari/updepla/users/tferrari/Experiments/RNA_analysis/IRIS/CHUV_samples/signatures',
        'out_mPath': '/home/tferrari/updepla/users/tferrari/Experiments/RNA_analysis/IRIS/CHUV_samples/results/HD_controls/all_samples',
    }
# Input pathways
paths['sig_inPath'] = os.path.join(paths['iris_mPath'], '1_scripts', 'iris_scRNAseq', '3_signature_genes')
paths['seq_inPath'] = os.path.join(paths['iris_mPath'], '4_sequencing_analysis')
paths['img_inPath'] = os.path.join(paths['iris_mPath'], '6_imaging_analysis')

# Output pathways
paths['fig_outPath'] = os.path.join(paths['out_mPath'], 'figures')
paths['obj_outPath'] = os.path.join(paths['out_mPath'], 'objects')
paths['table_outPath'] = os.path.join(paths['out_mPath'], 'tables')

l_outPaths = [paths['fig_outPath'], paths['obj_outPath'], paths['table_outPath']]
for outPath in l_outPaths:
    create_dir(outPath)

signature_PBMCs = pd.read_csv(os.path.join(paths['sig_mPath'], 'merged_indiv_sig.txt'),
                              sep='\t')
signature_CC_hum_mouse = pd.read_csv(os.path.join(paths['sig_mPath'], 'signature_cell_cycle_human_mouse.txt'),
                                      sep='\t')
#######
# Integration parameters to be modified
#######
species = 'human'  # 'human' or 'mouse'
integrate_over = 'expID' # 'expID' or 'CC'

# Threshold values
min_cells_percent = 0.005
min_gene_number = 500
mito_cutoff = 20
min_nUMI = 2500

#Set starting dimensions
dims_use = 50
#common variables to regress are: "n_UMI", "percent.mt", "Phase"
vars_to_regress =  []

set_seed = 300

#Set threshold values for marker DEGs
min_pct = 0.1
logfc_threshold = 0.2

#Set experiment IDs to be removed if necessary (only necessary if 
#want to maintain "uniqueness" of clusters for datasets with low relative cell counts)
cond2rmv = []

#Set whether to first split layers and then perform IntegrateLayers on the seurat.object
#Will split and integrate layers if "yes"
split_layers = "yes"
integrate_layers = "no"
# Integration method for integrate layers: CCAIntegration, RPCAIntegration, HarmonyIntegration
# FastMNNIntegration, scVIIntegration
#integration.method <- CCAIntegration

#Select the type of community clustering algorithm
cluster_algo = "Leiden"

#UMAP resolution
# resolution = 0.3
#Set UMAP reduction name
integration_name = "pca"
reduction_name = "umap"
#%%
#####
#Load Count matrix IRIS platform
#####
#Annotate doublets and store them in a list
l_image_doublet_ROI_annotation = annotate_doublets_enhanced(paths['img_inPath'], sampleIDs)

