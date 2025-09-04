import os
import pandas as pd
import numpy as np
import random
import ast  # For safely evaluating strings as Python expressions
from reportlab.lib.pagesizes import letter
from reportlab.pdfgen import canvas
import tempfile
import cv2

from utils.dir_paths import paths
from functions.convert_str_to_arrays import string_to_numpy_array, string_to_tuples_array

# Create and save the new dataframes based on 'cell_count'
df = pd.read_csv(os.path.join(paths['resultPath'], 'complete_df.csv'))
no_cell_df = pd.read_csv(os.path.join(paths['resultPath'], 'no_cell_df.csv'))
multi_cell_df = pd.read_csv(os.path.join(paths['resultPath'], 'multi_cell_df.csv'))

def safe_eval_centroid_str(centroid_str):

    try:
        # Use ast.literal_eval to safely evaluate the string back into Python objects
        tuples_list = ast.literal_eval(centroid_str)
        # Convert the inner representation into NumPy arrays
        centroids = [np.array(t[0]) for t in tuples_list]
        return centroids
    except (ValueError, SyntaxError, TypeError) as e:
        print(f"Error parsing centroids: {e}")
        return []
    
def create_pdf(dataframe, filename, imagePath, pdfPath, is_no_cell_df, intro_page_data=None):
    
    new_path = os.path.join(pdfPath, filename)
    if not os.path.exists(pdfPath):
        os.makedirs(pdfPath)
    c = canvas.Canvas(new_path, pagesize=letter)

    width, height = letter  # These are the dimensions of the page
    startY = 750  # Starting Y position for text on each page
    decrement = 180  # Height decrement for each cell image and text section
    imageHeight = 100  # Height for the images
    imageWidth = 100  # Width for the images, adjust as needed

    # Introduction Page
    if intro_page_data:
        c.drawString(72, height - 72, f"Total Cell Count: {intro_page_data['total_cell_count']}")
        c.drawString(72, height - 144, f"No Cells: {intro_page_data['no_cells_percentage']}%")
        c.drawString(72, height - 216, f"Single Cells: {intro_page_data['single_cells_percentage']}%")
        c.drawString(72, height - 288, f"Multi Cells: {intro_page_data['multi_cells_percentage']}%")
        c.drawString(72, height - 360, f"Date of Experiment: {intro_page_data['date_of_experiment']}")
        c.drawString(72, height - 432, f"Experimental ID: {intro_page_data['experimental_id']}")
        c.drawString(72, height - 504, f"Cell Type: {intro_page_data['cell_type']}")
        c.drawString(72, height - 576, f"Species: {intro_page_data['species']}")
        c.drawString(72, height - 648, f"Tissue: {intro_page_data['tissue']}")
        c.drawString(72, height - 720, f"Model Version: {intro_page_data['model_version']}")
        c.showPage()

    item_count = 0
    for index, row in dataframe.iterrows():
        # Calculate the current Y position for text and images
        currentY = startY - (item_count % 4) * decrement
        text = f"Full_Cell_ID: {row['Full_Cell_ID']}, Cell Count: {row['cell_count']}"
        c.drawString(72, currentY, text)
        
        # Construct image path
        img_name = f"{row['plate_code']}_{row['cell_code']}"
        cell_name = f"cell_{str(int(row['well_code'].split('WC')[-1]) - 1)}"

        if is_no_cell_df: channel_name = 'bf_ch9'
        else:
            channels = string_to_numpy_array(row['channels'])
            selected_channel_index = random.randint(0, len(channels) - 1)
            selected_channel = channels[selected_channel_index]
            channel_name = f"bf_ch{str(int(selected_channel))}"

            if 'centroids' in row and row['centroids']:
                centroids = string_to_tuples_array(row['centroids'])
                selected_centroid = centroids[selected_channel_index]
                selected_centroid = np.array(selected_centroid)
                #print(f"Name: {img_name}, cent: {selected_centroid}")
            else:
                selected_centroid = None

        image_path = os.path.join(imagePath, img_name, cell_name, f"{channel_name}.tif")
         
        if os.path.exists(image_path):
            img = cv2.imread(image_path)
        
            if not is_no_cell_df and selected_centroid.size > 0:
                for cent in centroids:
                    # Draw a rectangle around the selected centroid
                    x, y = int(cent[0]), int(cent[1])
                    cv2.rectangle(img, (x-64, y-64), (x+64, y+64), (0, 0, 255), 3)

            # Save the modified image to a temporary file
            with tempfile.NamedTemporaryFile(delete=False, suffix='.png') as tmpfile:
                temp_image_path = tmpfile.name
                cv2.imwrite(temp_image_path, img)

            # Draw the image on the PDF
            c.drawImage(temp_image_path, 72, currentY - imageHeight, width=400, height=imageHeight, preserveAspectRatio=True)

            # Remove the temporary file
            os.remove(temp_image_path)
        
        item_count += 1
        if item_count % 4 == 0:
            c.showPage()

    c.save()

def pdf_per_df(intro_data):

    intro_page_data = {
    "total_cell_count": len(df),
    "no_cells_percentage": (len(no_cell_df) / len(df)) * 100,
    "single_cells_percentage": (len(df[df['cell_count'] == 1])/len(df))*100,
    "multi_cells_percentage": (len(multi_cell_df)/len(df))*100,
    "date_of_experiment": intro_data['date_of_exp'],
    "experimental_id": df['exp_ID'].iloc[0],  # Assuming all rows have the same experimental ID,
    "cell_type": intro_data['cell_type'],
    "species": intro_data['species'],
    "tissue": intro_data['tissue'],
    "model_version": intro_data['model_version']
    }
    
    create_pdf(no_cell_df, 'no_cell_df.pdf', paths['yoloPath'], paths['resultPath'], is_no_cell_df=True, intro_page_data=intro_page_data)
    create_pdf(multi_cell_df, 'multi_cell_df.pdf', paths['yoloPath'], paths['resultPath'], is_no_cell_df=False, intro_page_data=intro_page_data)