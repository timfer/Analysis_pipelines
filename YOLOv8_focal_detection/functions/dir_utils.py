import os
import re
import pickle
import shutil

#Directory names should include the {plate_code}_{cell_code} in the following format P{x}_CC{x} where {x} is any numerical integer
pattern = r'P\d+_CC\d+'

def list_valid_directories(parent_path):
    '''
    Input:
    - parent_path: The path to the parent directory as a string.
    - pattern: A regular expression pattern to match directory names.

    Action:
    - Lists all subdirectories within the specified parent directory that match the given regular expression pattern.

    Output:
    - Returns a list of directory names that match the pattern.
    '''
    valid_dirs = []
    for item in os.listdir(parent_path):
        if os.path.isdir(os.path.join(parent_path, item)) and re.search(pattern, item):
            valid_dirs.append(item)
    return valid_dirs

def list_tif_files(directory_path):
    '''
    Input:
    - directory_path: The path to the directory as a string.

    Action:
    - Lists all files in the specified directory that have a '.tif' extension.

    Output:
    - Returns a list of full paths to the files with a '.tif' extension.
    '''
    return [os.path.join(directory_path, f) for f in os.listdir(directory_path) if f.endswith('.tif')]

def load_pkl_file(pkl_path):
    '''
    Input:
    - pkl_path: The path to the .pkl (pickle) file as a string.

    Action:
    - Loads and returns the contents of the specified .pkl file.

    Output:
    - Returns the loaded object from the .pkl file.
    '''
    with open(pkl_path, 'rb') as pkl_file: #Load the corresponding .pkl file
        file = pickle.load(pkl_file)
        return file
    
def create_dir(dirpath):
    '''
    Input:
    - dirpath: Directory path as a string.

    Action:
    - Checks if the directory exists at the specified path. If it does not exist, it attempts to create the directory.

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

def remove_directory(dirpath):
    '''
    Input:
    - dirpath: A string specifying the path to the directory that needs to be removed.

    Action:
    - Attempts to remove the specified directory and all its contents using the `shutil.rmtree()` method.
    - Handles the process within a try-except block to manage exceptions, such as the directory not existing or being in use.

    Output:
    - No return value. Prints a message indicating whether the directory was successfully deleted or if an error occurred during the deletion process.
    '''
    # Remove the directory and all its contents
    try:
        shutil.rmtree(dirpath)
        dir_name = dirpath.replace('\\', '/').split('/')[-1]
        print(f"Directory {dir_name} and all its components has been successfully deleted")
    except OSError as error: print(f"Error: {error.strerror}")