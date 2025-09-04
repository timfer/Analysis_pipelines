import os
import cv2
from PIL import Image
import numpy as np
import pickle
from sklearn.cluster import DBSCAN

from utils.dir_paths import paths
from utils.variables import eps
from functions.dir_utils import list_valid_directories, list_tif_files, load_pkl_file

#Gaussian blur variables
kernel_s = (5,5)
kermel_m = (25,25)
kernel_l = (49,49)
sigma_blur = 20.0

def cluster_centroids(temp_centroids):
    '''
    Input:
    - temp_centroids: A dictionary with keys as channel names and values as another dictionary of detected objects' centroids and bounding box sizes.

    Action:
    - This function starts by aggregating all centroid coordinates, bounding box dimensions, and channel names from `temp_centroids` into a list.
    - It converts this list into a structured numpy array, facilitating easier manipulation and analysis of spatial data.
    - Utilizes DBSCAN clustering to group centroids based on their spatial proximity, aiming to identify distinct objects captured across multiple channels or frames.
    - The epsilon parameter of DBSCAN determines the maximum distance between two samples for one to be considered as in the neighborhood of the other, crucial for accurate clustering based on the specific spatial scale of the images.
    - After clustering, a new structured array is created to include additional fields for each centroid: the Gaussian blur difference (initialized with infinity) and the cluster label assigned by DBSCAN.
    - Cluster labels are then assigned back to the centroids, effectively grouping them into identifiable clusters corresponding to individual objects or focal points across the channels.

    Output:
    - temp_cent_arr: A structured numpy array containing the original centroid data along with the DBSCAN cluster labels and placeholder values for Gaussian blur differences.
    - unique_labels: An array of unique cluster labels (excluding any noise points identified by DBSCAN with a label of -1), representing the identified objects or focal points.
    '''
    #Now that we've collected the centroids for each channel for this given image, we will now
    #perform a DBSCAN with epsilon = 30-35 to regroup all related centroids into a single point
    temp_c_list = []
    for ch, data in temp_centroids.items():
        for clust_num, cent in data.items():
            temp_c_list.append((cent['centroid'], cent['bbox_wh'], ch))

    # Convert temp_c_list to a structured numpy array with named fields
    dtype = [('centroid', float, (2,)), ('bbox_wh', float, (2,)), ('channel', 'U10')]
    temp_cent_arr = np.array(temp_c_list, dtype=dtype)

    #Perform DBSCAN grouping of centroids to assign unique labels to unique centroids, minimal pixel distance between centroids (drift on IRIS)
    #can be modulated (take cell size into account) by changing the value of the 'eps' variable in the config.py file
    cent_for_clustering = np.stack(temp_cent_arr['centroid'])
    dbscan = DBSCAN(eps=eps, min_samples=1)
    clusters = dbscan.fit(cent_for_clustering)
    cluster_labels = clusters.labels_ #cluster_labels is an array of same length as temp_cent_arr
    #print(f'# of clusters: {len(np.unique(cluster_labels))}, labels: {cluster_labels}')

    # Define new dtype to include 'cluster_label' and the 'gauss_dif'
    new_dtype = [('centroid', '<f8', (2,)), ('bbox_wh', float, (2,)), ('channel', '<U10'), ('gauss_diff', '<f8'), ('cluster_label', '<i4')]
    
    # Create a new array with the new dtype and the same shape as temp_cent_arr
    new_array = np.empty(temp_cent_arr.shape, dtype=new_dtype)

    # Copy data from temp_cent_arr to the new array
    for field in temp_cent_arr.dtype.names:
        new_array[field] = temp_cent_arr[field]

    # Add the cluster_labels to the 'cluster_label' field of the new array
    new_array['cluster_label'] = cluster_labels
    new_array['gauss_diff'] = float('inf')
    temp_cent_arr = new_array
    
    unique_labels = np.unique(temp_cent_arr['cluster_label'][temp_cent_arr['cluster_label'] != -1]) #-1 corresponds to the background label

    return temp_cent_arr, unique_labels

def calculate_gaussian_difference(img):
    '''
    Input:
    - img: The image array on which Gaussian blur difference calculation is performed.

    Action:
    - Applies Gaussian blurring with two different kernel sizes to the image, then calculates the absolute difference between the two resulting images.
    - Computes the mean of the absolute difference to quantify the blur difference.

    Output:
    - gauss_diff: The calculated Gaussian blur difference.
    '''
    #Calculate the gaussian difference for each channel for optimal focal cell detection
    gauss_l = cv2.GaussianBlur(img, kernel_l, sigma_blur)
    gauss_s = cv2.GaussianBlur(img, kernel_s, sigma_blur)
    subtracted = np.absolute(gauss_s.astype('float32')-gauss_l.astype('float32'))
    gauss_diff = cv2.mean(subtracted)[0]

    return gauss_diff


def assign_focal_cell(cell_centroids, focal_channel):
    '''
    Input:
    - cell_centroids: A structured numpy array of cell centroids that includes Gaussian blur differences.
    - focal_channel: The channel determined to have the focal cell based on minimal Gaussian blur difference.

    Action:
    - Filters centroids by the focal channel.
    - If multiple centroids are present, selects the one with the lowest Gaussian blur difference as the focal cell.

    Output:
    - focal_cell: The centroid data for the identified focal cell.
    '''
    focal_cell =  cell_centroids[cell_centroids['channel'] == focal_channel]

    if len(focal_cell) > 1: #This means that there are ovelapping bounding boxes that have been attributed as 
        #Here we use the same gaussian difference techique as if the bounding box is centered on the cell than the difference will be minimized
        focal_cell = focal_cell[focal_cell == focal_cell[np.isfinite(focal_cell['gauss_diff'])][np.argmin(focal_cell[np.isfinite(focal_cell['gauss_diff'])]['gauss_diff'])]]

    return focal_cell


def crop_image_around_centroid(image_array, centroid, crop_shape):
    '''
    Crops an image array around a specified centroid within the given dimensions.

    This function calculates a cropping box centered around a given centroid (x, y) coordinate
    within the image and crops the image to the specified dimensions. If the calculated cropping
    box extends beyond the boundaries of the image, adjustments are made to fit within the image's
    dimensions.

    Inputs
    - image_array: np.array
        The image array to be cropped, with dimensions (channels, height, width).
    - centroid: tuple of int
        A tuple (x, y) representing the centroid around which the image will be cropped. 
        Coordinates are in the format (column_index, row_index) corresponding to (width, height).
    - crop_shape: tuple of int
        The dimensions (crop_height, crop_width) to which the image will be cropped.

    Returns
    - cropped_image: np.array
        The cropped image array with dimensions (channels, crop_height, crop_width). The cropping
        is performed such that the centroid is as centered as possible within the new dimensions, 
        given the constraints of the original image's size.

    Notes
    - If the centroid is near the edge of the image and the requested crop_shape would extend beyond
      the image's boundaries, the function automatically adjusts the cropping box to fit within the
      image. This may result in the centroid not being perfectly centered in the cropped image.
    '''
    height, width = image_array.shape
    crop_width, crop_height = crop_shape
    crop_width, crop_height = int(crop_width), int(crop_height)
    
    # Calculate the top-left corner of the crop box
    x, y = centroid
    top_left_x = max(0, int(x - crop_width // 2))
    top_left_y = max(0, int(y - crop_height // 2))
    
    # Adjust the crop box if it exceeds the image dimensions
    if top_left_x + crop_width > width:
        top_left_x = width - crop_width
    if top_left_y + crop_height > height:
        top_left_y = height - crop_height

    # Ensure adjustments are integers (in case of any floating point operations)
    top_left_x = int(top_left_x)
    top_left_y = int(top_left_y)
    
    # Crop the image
    cropped_image = image_array[top_left_y:top_left_y+crop_height, top_left_x:top_left_x+crop_width]
    
    return cropped_image

def assign_bbox_to_dict(centroid, temp_centroids):
    '''
    Input:
    - centroid: A data structure (likely from a YOLO model prediction) that contains detection information for a single image. It is expected to have a `.boxes` attribute with bounding box data in an XYWH format, where X and Y denote the centroid coordinates and W and H represent the width and height of the bounding box respectively.
    - temp_centroids: A dictionary that is intended to store information about each detected object within a particular image or channel. The keys are numeric identifiers for each detected object, and the values are dictionaries containing the 'centroid' coordinates and 'bbox_wh' (bounding box width and height) for each detection.

    Action:
    - The function first checks if there are any detections in the `centroid` input by examining `centroid[0].boxes`. If detections are present, it iterates over each detection.
    - For each detection, it extracts the centroid coordinates (X, Y) and the dimensions of the bounding box (width, height) from the `.boxes` attribute. These values are converted to floating-point numbers and then stored in the `temp_centroids` dictionary under a key corresponding to the detection number. The stored value is a dictionary containing the 'centroid' (a tuple of X and Y coordinates) and 'bbox_wh' (a tuple of width and height).
    - This process effectively translates the raw bounding box data from the detection results into a structured format that is easier to manipulate and analyze in subsequent steps.

    Output:
    - temp_centroids: The updated dictionary containing detailed information about each detection's centroid and bounding box dimensions. This output allows for further processing, such as spatial clustering of centroids or analysis of detected objects' dimensions and positions within the image.
    '''
    if centroid[0].boxes:
        for cent_num in range(int(centroid[0].boxes.shape[0])):
            cent = (float(centroid[0].boxes.xywh[cent_num][0].cpu()), float(centroid[0].boxes.xywh[cent_num][1].cpu()))
            temp_centroids[cent_num] = {'centroid': cent,
                                        'bbox_wh': tuple((float(centroid[0].boxes.xywh[cent_num][2].cpu()),
                                                          float(centroid[0].boxes.xywh[cent_num][3].cpu())))}
            
    return temp_centroids

def create_centroid_dictionary():
    '''
    Input:
    - resultPath: The directory path where YOLO detection results are stored, including .pkl files for detections and .tif files for images.
    - yoloPath: The directory path containing images processed by YOLO, used for further image analysis to identify focal cells.

    Action:
    - Traverses the directory structure within `resultPath` to gather detection results for each plate and cell combination, identified by directory names and file names respectively.
    - For each cell, loads the detection results (.pkl files) to reconstruct the centroids and associated data, then aggregates these across all relevant channels or imaging planes.
    - Calls `cluster_centroids` to group these centroids into clusters representing individual objects across channels or frames, utilizing spatial clustering via DBSCAN.
    - Iterates over these clusters to identify the focal cell or object within each, based on Gaussian blur differences calculated across the clustered centroids.
    - Specifically, for each cluster, the function extracts image regions corresponding to each centroid, applies Gaussian blurring to these regions, and calculates the blur difference. This step aims to identify the most in-focus object, presuming a minimized blur difference signifies focal clarity.
    - Organizes and compiles this focal cell data, along with the total count of identified objects, into a nested dictionary structure keyed by plate and cell codes, facilitating easy access to the focal cell data for each analyzed unit.

    Output:
    - centroid_dico: A comprehensive dictionary containing detailed centroid data for focal cells across all processed plates and cells, along with counts of identified objects or focal points.
    '''
    centroid_dico = {} #Track centroids
    for plate_cc in list_valid_directories(paths['resultPath']):
        #print(f'{plate_cc}')
        cc_path = os.path.join(paths['resultPath'], plate_cc)
        
        centroid_dico[plate_cc] = {}

        for cell in os.listdir(cc_path):
            #print(f"{cell}")
            centroid_dico[plate_cc][cell] = {}

            target_dir = os.path.join(cc_path, cell)

            temp_centroids = {}
            centroid_dico[plate_cc][cell] = {'data':{}, 'cell_count': None} #Track centroids

            tif_files = list_tif_files(target_dir)
            if tif_files:
                for tif_path in tif_files:
                    ch_name = tif_path.replace(('\\'), ('/')).split('/')[-1] .split('_')[1] #Make ch_name determination system independent
                    temp_centroids[ch_name] = {}

                    if tif_path.endswith('.tif'): #Simple check for actual tif files in list but not really necessary here
                        # Construct the corresponding .pkl filename
                        pkl_path = tif_path.replace('.tif', '.pkl')
                        img_centroid = load_pkl_file(pkl_path)

                        temp_centroids[ch_name] = assign_bbox_to_dict(img_centroid, temp_centroids[ch_name])
            else:
                print(f'Image {plate_cc}, cell: {cell} has no focal cells')

            if len(temp_centroids) > 0:
                #Now that we've collected the centroids for each channel for this given image, we will now
                #perform a DBSCAN with epsilon = 30-35 to regroup all related centroids into a single point
                temp_cent_arr, unique_labels = cluster_centroids(temp_centroids)

                for lbl in unique_labels:
                    focal_channel = None #Initiate focal_channel value
                    cell_centroids = temp_cent_arr[temp_cent_arr['cluster_label'] == lbl] #Select all centroids associated to a specific labels

                    if len(cell_centroids) > 0: #In this case we select the middle plane
                        for idx, cent in enumerate(cell_centroids):
                            #Extract the cropped cell image and perform gaussian blurring and calculate the difference to determine the focal plane between related centroids
                            img_path = os.path.join(paths['yoloPath'], plate_cc, cell, f'bf_{cent["channel"]}.tif')
                            img = Image.open(img_path) # Load the image with Pillow
                            img = np.array(img) #Convert to nd.array

                            cropped = crop_image_around_centroid(img, tuple(cent['centroid']), tuple(cent['bbox_wh']))

                            gauss_dif = calculate_gaussian_difference(cropped)
                            #print(f"gauss: {gauss_dif}")

                            cell_centroids['gauss_diff'][idx] = gauss_dif
                        #print(cell_centroids)
                        #Now that the gaussian difference has been calculated for each cell, then we select the lowest value and that is the focal channel for this specific cell
                        focal_channel = cell_centroids[np.isfinite(cell_centroids['gauss_diff'])]['channel'][np.argmin(cell_centroids[np.isfinite(cell_centroids['gauss_diff'])]['gauss_diff'])]

                    focal_cell = assign_focal_cell(cell_centroids, focal_channel)
                    centroid_dico[plate_cc][cell]['data'][lbl] = focal_cell

                centroid_dico[plate_cc][cell]['cell_count'] = len(unique_labels)
                print(f'Image {plate_cc}, cell: {cell} has {len(unique_labels)} focal cells')
    
    with open(os.path.join(paths['resultPath'], 'centroid_dico.pkl'), 'wb') as pkl_file:
        pickle.dump(centroid_dico, pkl_file)

    return centroid_dico