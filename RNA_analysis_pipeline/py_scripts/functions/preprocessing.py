import os
import re

import pandas as pd
import numpy as np

import scanpy as sc
from anndata import AnnData

from pathlib import Path
from typing import List, Dict, Tuple, Optional
import warnings

def annotate_doublets_enhanced(
    PATH_input_IRIS_imaging: str, 
    sample_ids_IRIS: List[str],
    file_name: str = "complete_df.csv",
    verbose: bool = True
) -> Dict[str, pd.DataFrame]:
    """
    Enhanced version with better error handling and flexibility.
    
    Parameters:
    -----------
    PATH_input_IRIS_imaging : str
        Base path to IRIS imaging data directory
    sample_ids_IRIS : List[str]
        List of sample identifiers to process
    file_name : str, default "complete_df.csv"
        Name of the CSV file to read
    verbose : bool, default True
        Whether to print progress information
        
    Returns:
    --------
    Dict[str, pd.DataFrame]
        Dictionary with sample IDs as keys and corresponding annotation dataframes as values
    """
    
    l_image_doublet_ROI_annotation = {}
    failed_samples = []
    
    for sample_id in sample_ids_IRIS:
        # Use pathlib for more robust path handling
        file_path = Path(PATH_input_IRIS_imaging) / sample_id / "image_annotation" / "algo" / "cells_per_ROI" / file_name
        
        try:
            if file_path.exists():
                df = pd.read_csv(file_path)
                l_image_doublet_ROI_annotation[sample_id] = df
                
                if verbose:
                    print(f"✓ Loaded {sample_id}: {len(df)} rows, {len(df.columns)} columns")
            else:
                warnings.warn(f"File does not exist: {file_path}")
                failed_samples.append(sample_id)
                l_image_doublet_ROI_annotation[sample_id] = pd.DataFrame()
                
        except Exception as e:
            warnings.warn(f"Error loading {sample_id}: {str(e)}")
            failed_samples.append(sample_id)
            l_image_doublet_ROI_annotation[sample_id] = pd.DataFrame()
    
    if verbose and failed_samples:
        print(f"\nFailed to load {len(failed_samples)} samples: {failed_samples}")
    
    return l_image_doublet_ROI_annotation


