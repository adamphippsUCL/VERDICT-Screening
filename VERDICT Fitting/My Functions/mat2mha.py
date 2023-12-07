# Python script to convert .mat image volume files to mha so hey can be viewed in ITK SNAP

# Import relevant libraries
import SimpleITK as sitk
from scipy.io import loadmat
import numpy as np
import sys
import os
import glob

# Read output folder
OutputFolder = str(open('output_folder.txt', 'r').read())

# Find list of patient numbers
PatNums = [os.path.basename(path) for path in 
           glob.glob(r"C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\VERDICT Screening\Outputs\fIC ROIs\*")]
# PatNums = ['BAR_003']
print(PatNums)

ModelNames = ['Original VERDICT', 'RDIex903000']#, 'noVASCex903000']#list(range(0,15))

VolumeName = 'fIC'

normalise = True

for PatNum in PatNums:

    for ModelName in ModelNames:
        
        try:
            # Load .mat file
            volume = loadmat(f'{OutputFolder}/VERDICT outputs/{PatNum}/{ModelName}/{VolumeName}.mat')[VolumeName]
            # remove infinities
            volume[volume == np.inf] = 0
            # Remove nan
            volume[np.isnan(volume)] = 0
            
            # Change image orientation
            volume = np.moveaxis( volume , -1, 0)  
            
            if normalise:
                volume[volume>1] = 1

            # Save as mha file
            sitk.WriteImage( sitk.GetImageFromArray(volume), f'{OutputFolder}/VERDICT outputs/{PatNum}/{ModelName}/{VolumeName}.mha' )
            
            del volume
            
        except: None




