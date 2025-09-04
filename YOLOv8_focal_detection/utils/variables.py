'''Modify according to cell type and drift'''
#DBSCAN parameters
eps = 35 #Minimal pixel distance by which the centroids of cells must be separated (minimal drift): PBMC ~32-35, Fucci ~45-50, CTCs: TBD

'''
Modify according to channels to be dropped if frequent trash in extreme channels, 
or any other channel you wish to be dropped in post-processing (WARNING: Will not be
included in final dataframe)
'''
#Establish channels to be dropped for refined focal cell detection (variable depending on experiment)
dropped_planes = [] #['ch0', 'ch1', 'ch2', 'ch15', 'ch16', 'ch17', 'ch18', 'ch19'] #Usually only relevant for extreme planes

'''     
Modiy according to desired printout in final PDF, does not affect stats
'''
#Modify intro data to be in accord with experimental information
intro_data = {
    'date_of_exp': 'YYYY-MM-DD', # Add date of experiment or RNA-seq submission
    'cell_type': 'Cell Type',
    'species': 'Species of origin',
    'tissue': 'Tissue of origin',
    'model_version': 'V1' #DO NOT MODIFY: Relates to version of YOLOv8 fine-tuned model
}

######################################################################
#---------------------------------------------------------------------
######################################################################

#Unecessary to modify for IRIS datasets: Set YOLO prediction variables
imgsz=(512,2048) #(heightxwidth); IRIS image dimensions: (512, 2048). Can be changed to any value but must be real image dimensions
classes = [] #Will take the classes from the model, currently the model is only set up to detect the focal_cell class
conf=0.8 #Increasing confidence level will increase specificity but decrease sensitivity
