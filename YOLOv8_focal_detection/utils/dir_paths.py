import os, platform
'''
Of note: The directory paths following the "general experimental directory (masterPath)" 
overlap across different processes
'''

'''Modify masterPath to lead to experimental directory'''
if platform.system() == 'Windows':
    print(f"The prediction is run through Windows system")
    paths = {
        'masterPath': r'Windows_path_to_experimental_directory' #Absolute path to the experimental master directory
    }
else:
    print(f"The prediction is run through Unix system")
    paths = {
        'masterPath': 'Unix_path_to_experimental_directory' # Absolute path to the experimental master directory for Unix systems
    }
'''Modify according to result directory position. As set is a subdirectory of masterdir'''
paths['dataPath'] = os.path.join(paths['masterPath'], 'HDF5_output') #Path to HDF5 files
paths['resultPath'] = os.path.join(paths['masterPath'], 'yolo_prediction_results') #subdirectory of masterPath in which the prediction results will be saved

######################################################################
#---------------------------------------------------------------------
#####################################################################

#Unecessary to modify: Directories included or created automatically and are intermediary folders
paths['tiffPath'] = os.path.join(paths['masterPath'], 'TIFF_output') #Path to tiff stacks
paths['yoloPath'] = os.path.join(paths['masterPath'], 'yolo_images') #Path to single tif images

'''DO NOT TOUCH!!'''
paths['modelPath'] = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'models', 'runs', 'detect', 'train', 'weights') #Path to the best pretrained model (directory containing .pt file)