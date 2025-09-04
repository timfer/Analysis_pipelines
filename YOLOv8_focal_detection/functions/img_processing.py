################################
# Function for image cropping
################################
def crop_image_around_centroid(image_array, centroid, crop_shape):
    """
    Crops an image array around a specified centroid within the given dimensions.

    This function calculates a cropping box centered around a given centroid (x, y) coordinate
    within the image and crops the image to the specified dimensions. If the calculated cropping
    box extends beyond the boundaries of the image, adjustments are made to fit within the image's
    dimensions.

    Inputs
    ----------
    - image_array: np.array
        The image array to be cropped, with dimensions (channels, height, width).
    - centroid: tuple of int
        A tuple (x, y) representing the centroid around which the image will be cropped. 
        Coordinates are in the format (column_index, row_index) corresponding to (width, height).
    - crop_shape: tuple of int
        The dimensions (crop_height, crop_width) to which the image will be cropped.

    Returns
    ----------
    - cropped_image: np.array
        The cropped image array with dimensions (channels, crop_height, crop_width). The cropping
        is performed such that the centroid is as centered as possible within the new dimensions, 
        given the constraints of the original image's size.

    Notes
    ----------
    - If the centroid is near the edge of the image and the requested crop_shape would extend beyond
      the image's boundaries, the function automatically adjusts the cropping box to fit within the
      image. This may result in the centroid not being perfectly centered in the cropped image.
    """
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