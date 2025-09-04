import os
import tifffile

from functions.dir_utils import create_dir

def convert_images_to_single_stack_tif(imgPath, yoloPath):
    '''
    Input:
    - imgPath: The path to the directory containing folders for different sets of images. Each folder represents a different plate or experimental condition,
               and contains subfolders for individual cells. Each cell subfolder contains TIFF files for different imaging channels.
    - yoloPath: The path to the directory where the output single-image TIFF files will be stored. The directory structure will be mirrored from `imgPath`, 
                but each TIFF stack will be split into individual images.

    Action:
    - Iterates through each folder in `imgPath`, representing different plates or conditions.
    - For each cell in a folder, it finds the bright field (bf) TIFF file. This file is assumed to be a stack of images, one for each imaging channel or time point.
    - Creates a target directory within `yoloPath` that mirrors the structure in `imgPath` but is specifically for the output single-image files.
    - Reads the bright field TIFF stack and splits it into individual images. Each image is saved as a separate TIFF file within the appropriate target directory.
    - The naming convention for output files includes the channel index to differentiate between the images from the same stack.

    Output:
    - No return value. For each bright field image stack in the input directories, this function creates a series of single-image TIFF files in the output directory.
      This conversion is necessary for applications like YOLO (You Only Look Once), which require single images for object detection tasks.
    '''
    for folder in os.listdir(imgPath):
        target_path = os.path.join(imgPath, folder)
        for cell in os.listdir(target_path):
            target_dir = os.path.join(target_path, cell)
            bf_img = [file for file in os.listdir(target_dir) if file in ('bf.tif', 'bf.tiff')][0]
            bf_img_path = os.path.join(target_dir, bf_img)
            bf_target_path = os.path.join(yoloPath, folder, f'cell_{cell}') #YOLO can only run predictions on single images so we save the .tif hyperstack as single .tif images
            create_dir(bf_target_path)

            full_img = tifffile.imread(bf_img_path)
            
            for channel in range(full_img.shape[0]):
                img_name = f'bf_ch{channel}.tif'
                img = full_img[channel, :, :]
                tifffile.imwrite(os.path.join(bf_target_path, img_name), img)