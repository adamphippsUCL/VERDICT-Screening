# Python script to convert .mat image volume files to mha so hey can be viewed in ITK SNAP

# Import relevant libraries
import SimpleITK as sitk
from scipy.io import loadmat
import numpy as np
import sys

# Read output folder
OutputFolder = str(open('output_folder.txt', 'r').read())

# Define Patient ID, model type, and volume name
PatNum = 'INN_221'
ModelNumbers = [999]#list(range(0,15))

VolumeName = 'fIC'

for ModelNum in ModelNumbers:
    
    try:
        # Load .mat file
        volume = loadmat(f'{OutputFolder}/VERDICT outputs/{PatNum}/Model {ModelNum}/{VolumeName}.mat')[VolumeName]
        # remove infinities
        volume[volume == np.inf] = 0
        # Remove nan
        volume[np.isnan(volume)] = 0
        
        # Change image orientation
        volume = np.moveaxis( volume , -1, 0)  

        # Save as mha file
        sitk.WriteImage( sitk.GetImageFromArray(volume), f'{OutputFolder}/VERDICT outputs/{PatNum}/Model {ModelNum}/{VolumeName}.mha' )
        
    except: None




