import os
os.environ['R_HOME'] = r'C:\PROGRA~1\R\R-4.5.0'  # Short path for "Program Files"
os.environ['PATH'] += r';C:\PROGRA~1\R\R-4.5.0\bin\x64'
import sys
sys.path.append(r'C:\PROGRA~1\R\R-4.5.0\bin\x64') # Ensure Python sees R


def create_dir(dirpath):
    '''
    Input:
    - dirpath: Directory path as a string.

    Action:
    - Checks if the directory exists at the specified path. If it exists, it attempts to create a new directory
      with a suffix `_r{x}`, incrementing `x` until a unique directory name is found.
    - If the directory does not exist, it attempts to create the directory.

    Output:
    - No return. Prints a message indicating whether the directory was created, already exists, or if an error occurred during creation.
    '''
    if not os.path.exists(dirpath):
        try:
            os.makedirs(dirpath)
            print(f"Directory '{dirpath}' created successfully.")
        except OSError as e:
            print(f"Error creating directory '{dirpath}': {e}")
    else:
        print(f"Directory '{dirpath}' already exists.")

        
def list_valid_directories(parent_path, pattern):
    '''
    Input:
    - parent_path: The path to the parent directory as a string.
    - pattern: A regular expression pattern to match directory names.

    Action:
    - Lists all subdirectories within the specified parent directory that match the given regular expression pattern.

    Output:
    - Returns a list of directory names that match the pattern.
    '''
    for item in os.listdir(parent_path):
        if os.path.isdir(os.path.join(parent_path, item)) and re.search(pattern, item):
            return item
        
def convert_rdata_to_h5ad(sample_path: str, species: str):
    """
    Convert RData objects to AnnData format - RUN THIS FIRST FOR EACH SAMPLE
    """
    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri
    from rpy2.robjects.conversion import localconverter
    
    # Construct RData path
    rdata_path = os.path.join(
        sample_path,
        "pipeline_output",f"{species}",
        "analysis", "objects","objects_plotting.RData"
    )
    
    if not os.path.exists(rdata_path):
        raise FileNotFoundError(f"RData file not found at {rdata_path}")
    
    # Load R environment
    ro.r(f'load("{rdata_path}")')
    
    with localconverter(ro.default_converter + pandas2ri.converter):
        pivot_table = ro.r('pivot_table')
        umi_matrix = ro.r('UMI_all_cells_SYMB')
    
    # Create AnnData object
    adata = AnnData(
        X=umi_matrix.T.astype(np.float32),  # Genes x cells -> cells x genes
        obs=pivot_table.set_index('Full_Cell_ID'),
        var=pd.DataFrame(index=umi_matrix.index.astype(str))
    )
    
    # Save as H5AD
    output_path = os.path.join(
        sample_path,
        f"pipeline_output/{species}/analysis/objects/adata.h5ad"
    )
    adata.write_h5ad(output_path)
    print(f"Converted {rdata_path} to {output_path}")
