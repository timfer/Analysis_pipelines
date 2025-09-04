#%%
import time

from utils.dir_paths import paths
from utils.variables import intro_data

from functions.dir_utils import create_dir, remove_directory
from functions.convert_h5_to_tiff import convert_hdf5_to_tiff
from functions.convert_images_to_single_stack_tif import convert_images_to_single_stack_tif
from functions.yolo_detect_utils import run_predictions_and_save_results
from functions.centroid_analysis import create_centroid_dictionary
from functions.create_cell_dataframes import create_cell_dataframe

def main():
    '''Convert hdf5 output to TIFF stacks'''
    convert_hdf5_to_tiff()
    print(f".hdf5 files correctly converted to .tiff stacks for each imaged channel \n")

    '''Delete the HDF5_directory'''
    #remove_directory(paths['dataPath'])

    '''Create a directory to save results in (modify in config)'''
    create_dir(paths['resultPath'])
    print(f"Result directory created: {paths['resultPath']} \n")

    '''Convert files to '.tif' single stacks if necessary (YOLO can only take 2D images)'''
    convert_images_to_single_stack_tif(paths['tiffPath'], paths['yoloPath'])
    print(f".tiff stacks converted to single .tiff images \n")

    '''Run predictions on single '.tif' images'''
    run_predictions_and_save_results(paths['yoloPath'], paths['resultPath'])
    print(f"Yolov8-based prediction of focal cells complete \n")

    '''Create a centroid dictionary with all relevant information per cell'''
    create_centroid_dictionary()

    '''Create and save dataframes for relevant conditions: complete, single cell, no cell and multi cell dataframes'''
    create_cell_dataframe()
    print(f"Dataframes containing all relevant cell information have been created")
    print(f"Path to results: {paths['resultPath']} \n")

    '''Create summary .pdf with experimental information and image examples of images with no cells and multiple cells'''
    from functions.create_pdf import pdf_per_df
    pdf_per_df(intro_data)
    print(f"Summary .pdf file created in {paths['resultPath']}")

if __name__ == "__main__":
    start_time = time.time()  # Capture start time
    main()
    end_time = time.time()  # Capture end time
    total_time = end_time - start_time  # Calculate total duration
    print(f"Total execution time: {total_time:.2f} seconds")

# %%
