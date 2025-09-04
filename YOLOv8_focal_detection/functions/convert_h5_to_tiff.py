import os
import numpy as np
import h5py
import tifffile
import re

from utils.dir_paths import paths
from functions.dir_utils import create_dir


def convert_hdf5_to_tiff():
    '''
    Input:
    - Assumes predefined global variables `dataPath` and `tiffPath`, where `dataPath` is the directory containing HDF5 files 
      and `tiffPath` is the target directory for the generated TIFF files.

    Action:
    - Calls `create_dir` to ensure the TIFF directory exists.
    - Iterates through directories in `dataPath`. Each directory name is processed to extract meaningful parts like plate codes and cell codes, 
      using regular expressions to filter out specific patterns (e.g., 'P' followed by numbers or 'CC' followed by numbers).
    - For each cell number (indicating a separate HDF5 file within a directory), it constructs a path for the corresponding TIFF files, 
      ensuring the directory structure mirrors that of the input but replaces HDF5 files with TIFFs.
    - Opens each HDF5 file, reads specific datasets corresponding to different imaging channels (e.g., bright field, 405 nm, 488 nm, 638 nm fluorescence),
      and saves each as a multi-page TIFF file, with one page per plane of the dataset. The datasets are transposed to have channels as the first dimension.
    - The function iterates through every cell (HDF5 file) within each directory and processes them accordingly, printing out the progress.

    Output:
    - No return value. This function creates TIFF files in a structured directory hierarchy mirroring that of the input HDF5 files, with console 
      output indicating progress and details of the conversion process.
    '''
    pattern = re.compile(r'((?<![A-Za-z])P\d+|CC\d+)') #Keep all elements where a 'P' is followed by a number, same for the 'CC' elements, remove all elements where a 'P' follows another letter (removes JP...)

    create_dir(paths['tiffPath'])

    i = 0
    for dir in os.listdir(paths['dataPath']):
        dir_names = (dir.replace('\\', '/').split('/')[-1]).split('_') #All directory name elements are separated by an underscore. Here we made it OS independent
        # Filtering elements that match the pattern
        
        filtered_elements = [element for element in dir_names if pattern.search(element)]
      
        if not filtered_elements:
            continue
        platecode, ccode = filtered_elements[:2]

        folder = os.path.join(paths['dataPath'], dir)
        print(folder)
        
        for cl_num in os.listdir(folder): #Cell number per cell code
            fpath = os.path.join(folder, cl_num)
            cl_num = cl_num.split('.hdf5')[0]

            print(f'pcode: {platecode}, ccode: {ccode}, cl_num: {cl_num}')

            cl_path = os.path.join(paths['tiffPath'], f'{platecode}_{ccode}', cl_num)
            os.makedirs(cl_path, exist_ok=True)
            
            with h5py.File(fpath, 'r') as file:
                channels = {
                    'bf': 'images/bf',
                    '405': 'images/405',
                    '488': 'images/488',
                    '638': 'images/638',
                }
                for channel_name, group_path in channels.items():
                    if group_path in file:
                        group = file[group_path]
                        planes = sorted(list(group.keys()), key=lambda x: int(x))
                        dset = [np.array(group[name]) for name in planes]
                        dset = np.stack(dset, axis=2).transpose((2, 0, 1))  # Transpose to (channels, height, width)
                        tifffile.imwrite(os.path.join(cl_path, f'{channel_name}.tif'), dset)
            i += 1