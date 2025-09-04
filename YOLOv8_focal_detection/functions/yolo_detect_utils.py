# yolo_detect.py
from ultralytics import YOLO
import os
from PIL import Image
import numpy as np
import cv2
import pickle
import re

from utils.dir_paths import paths
from utils.variables import classes, conf, imgsz
from functions.dir_utils import pattern, create_dir

#Import best pretrained model
model = YOLO(os.path.join(paths['modelPath'], '20240331_yolov8m_best.pt')) #Load pre-trained model
#In this case (03.04.2024), the yolov8m model has been fine-tuned on a combination of images of different celly types and cell lines

classes = []
conf=0.8
imgsz=(512,2048)

'''The 2 functions below come from: 
https://medium.com/@kleve.2406/how-to-use-yolov8-and-yolo-nas-for-object-detection-8c5893939480
'''

def predict(model, img, classes=classes, imgsz=imgsz, conf=conf):
    '''
    Input:
    - model: The loaded YOLO model.
    - img: The image to predict on, which can be a path to an image file or a numpy array.
    - classes: Optional. Specific classes to detect.
    - imgsz: The image size to which the input images should be resized before prediction.
    - conf: The confidence threshold for detections.

    Action:
    - Runs the model's prediction method on the provided image, optionally filtering by classes, resizing images to `imgsz`, and applying a confidence threshold.
    
    Output:
    - results: The prediction results, including detected objects, their classes, confidence scores, and bounding boxes.
    '''
    if classes:
        results = model.predict(img, imgsz=imgsz, classes=classes, conf=conf)
    else:
        results = model.predict(img, imgsz=imgsz, conf=conf)
    return results



def predict_detect(model, img, classes=classes, conf=conf):
    '''
    Input:
    - model: The loaded YOLO model.
    - img: The image to run detection on, which can be a path to an image file or a numpy array.
    - classes: Specific classes to detect.
    - conf: The confidence threshold for detections.

    Action:
    - First, calls `predict` to run the model's prediction on the provided image.
    - Loads the image if a path is provided, or uses it directly if it's a numpy array.
    - Draws bounding boxes and class labels on the detected objects in the image.

    Output:
    - image: The input image annotated with bounding boxes and class labels for detected objects.
    - results: The prediction results, including detected objects, their classes, confidence scores, and bounding boxes.
    '''
    results = predict(model, img, classes=classes, imgsz=imgsz, conf=conf)

    if isinstance(img, str):  # If img is a path, load the image
        image = Image.open(img)
        image = np.array(image)
    elif isinstance(img, np.ndarray):  # If img is an array, use it directly
        image = img
    else:
        raise ValueError("img must be a file path or a numpy array")

    for result in results:
        for box in result.boxes:
            cv2.rectangle(image, (int(box.xyxy[0][0]), int(box.xyxy[0][1])),
                          (int(box.xyxy[0][2]), int(box.xyxy[0][3])), (255, 0, 0), 2)
            cv2.putText(image, f"{result.names[int(box.cls[0])]}",
                        (int(box.xyxy[0][0]), int(box.xyxy[0][1]) - 10),
                        cv2.FONT_HERSHEY_PLAIN, 1, (255, 0, 0), 1)
    return image, results



def run_predictions_and_save_results(yoloPath, resultPath):
    '''
    Input:
    - yoloPath: The directory path where the input images for YOLO predictions are located.
    - resultPath: The directory path where the prediction results (images with detections and pickle files) will be saved.

    Action:
    - Iterates through each image in `yoloPath`, runs object detection using `predict_detect`, and checks if any detections were made.
    - If detections are present, saves the annotated image and a pickle file containing the detection results.
    - Creates necessary directories for saving the results, maintaining the original directory structure.

    Output:
    - No return value. This function saves the annotated images and pickle files with detection results in `resultPath`.
    '''
    for folder in os.listdir(yoloPath):
        if not re.match(pattern, folder): #This enforces that the identified files have to have the plate_code and cell_code pattern in their name
            continue

        target_path = os.path.join(yoloPath, folder)
        for cell in os.listdir(target_path):
            target_dir = os.path.join(target_path, cell)
            result_dir = os.path.join(resultPath, folder, cell)
            
            create_dir(result_dir)
            
            for img in os.listdir(target_dir):
                name = img.split('.tif')[0]
                
                bf_img_path = os.path.join(target_dir, img)

                # Run detection and prediction
                image, results = predict_detect(model, bf_img_path, classes=classes, conf=conf)

                # Check if there are any detections before saving
                if results[0].boxes.data.shape[0] > 0:
                    cv2.imwrite(os.path.join(result_dir, f'{name}_result.tif'), image)
                    with open(os.path.join(result_dir, f'{name}_result.pkl'), 'wb') as file:
                        pickle.dump(results, file)