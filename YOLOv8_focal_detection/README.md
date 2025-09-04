# YOLOv8 cell detection pipeline

###
# Pre-requirement absolutes
###
#1: Imaging outputs HAVE to be of '.hdf5' format

#2: Cell '.hdf5' files must ONLY include the well NUMBER

#3: Directory names for .hdf5 files must include the information about the plate-code and cell-code following the P{x}_CC{x} format where {x} can be any numerical integer combination
 - Note1: The plate and cell codes MUST be separated by underscores 
 - Note2: Initials cannot be CC
 - Note3: Directory names can contain any other information but none of this will be included in the final dataframe other than the experimental ID, plate code, cell code and well code

###
# Pre-required inputs
###
#1: In the '/utils' folder open 'dir_path.py':
  - Modify the 'masterPath' to the master experimental directory (line 10 if WINDOWS or line 15 if Unix)
    - Can be absolute path to dir or relative if in related directories with no backtracking
  - Can change 'dataPath' for '.hdf5' master dir (currently '<experiment_dir>/HDF5_output')
  - Can change 'resultPath' for finalized dataframes (currently '<experiment_dir>/yolo_prediction_results')

###
# Optional inputs
###
#1: In the '/utils' folder open 'variables.py', you can modify the:
- 'eps': Modify depending on drift of cells (PBMCs ~32, NIH ~42, CTCs ~42-49 (depends on error of multi-cell counting)
- 'dropped_planes': Include planes in the 'ch{x}' with {x} being any numerical integer in range (0-indexed planes). Usually only for extreme planes
   - 'eps' & 'dropped_planes': are pixels (can check with FiJi straight line length)
   - IMPORTANT: If you want to fine-tune the drift and dropped_planes you only need to rerun 'create_cell_dataframe()' in 'main.py' as this is the step that determines centroid uniqueness
- 'intro_data': Modify for the printout on the first page of the '.pdf' outputs

###
# Notes
###
- Running convert_hdf5_to_tiff will include all channels imaged (bf, 405, 488, 561, 638) or any combination therein and outputs per cell will be named as such with '.tif'
- For the next step cells must be indexed as follows: /P{x}_CC{y} where {x} and {y} can be any numerical integer

#Run yolo preprocessing
- Since yolo can only run predictions on 2D images, .tiff stacks must be converted to single tiff images per imaging event
- These are saved to a folder called 'yolo_images' (if not modified in 'dir_path.py')

#Run yolo predict
- The predictions are saved to a subdirectory 'yolo_prediction_results'
- The prediction result directories are named and indexed the same as the directories containing the single tiff images
- The files are named 'bf_ch{x}_result.tif', where {x} is the numerical integer corresponding to the imaging plane
- These .tif images have drawn bounding boxes on them for inspection
- The directories per cell also contain pickle files which contain all relevant information about the cell (centroids and bounding box dimensions)

#Run create_centroid_dictionary() and create_cell_dataframe()
- Although these 2 functions are not compiled together they are part of a similar post-processing step
- Since a cell can be detected in multiple planes, one has to determine whether the centroids detected in different planes belong to the same cell or not
- The first creates a dictionary where centroids are analyzed and determined to be unique or not
	- Accounting for drift on our platform (IRIS), one has to take cell size and drift into account
	- Modify the 'eps' variable in the 'variables.py'. PBMCs: eps<35, Fucci: eps>=42
- Next, once the centroids are 'cleaned up', multiple dataframes are created containing all relevant information per imaging event
	- 'complete_df.csv' (all imaged cells), 'single_cell_df.csv' (only images determined to have a single cell), 'no_cell_df.csv' (all imaging events where no cell was detected), 'multi_cell_df.csv' (all imaging events where multiple cells are detected)
- All results are saved within the resultPath directory set in 'dir_path.py'
- All dataframes have the same columns: ['Full_Cell_ID', 'exp_ID', 'plate_code', 'well_code', 'cell_code', 'Cluster_names', 'Condition', 'nCount_RNA', 'nFeature_RNA', 'centroids', 'channels', 'cell_count'].
	- Here df['centroids]' contains lists of tuples and df['channels'] lists of strings
	- Here, the ['Cluster_names', 'Condition', 'nCount_RNA', 'nFeature_RNA'] will contain NaNs as these are updated after sequencing is available
