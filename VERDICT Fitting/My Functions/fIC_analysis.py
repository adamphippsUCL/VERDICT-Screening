# Python script to contain functions to analyse fIC values with lesion ROI

# Import relevnat libraries
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.path as mpltPath
import sys
import glob
import os
from scipy.io import loadmat, savemat
import pandas as pd
import pickle
from sklearn.metrics import roc_curve, roc_auc_score
import SimpleITK as sitk
from skimage.metrics import mean_squared_error
import scipy
import skimage

# Import DICOM from imgtools
sys.path.insert(0, r"C:\Users\adam\OneDrive - University College London\UCL PhD\Image-Processing")
import DICOM # type: ignore


# Function for saving ROI masks
def saveROImask(
                PatNum, 
                ROIName, 
                INNOVATE_path = r"D:\UCL PhD Imaging Data\INNOVATE VERDICT", 
                INNOVATE_ROIs_path = r"C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\INNOVATE\INNOVATE ROIs NT",
                output_path = r"C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\ISMRM Submission\Outputs\ROIs"):
    
    '''
    Function to save mask of ROI in RTStruct file.
    
    ROIs drawn on b3000 VERDICT scan as this is used as the registration target
    
    '''
    
    try:
        # Define path to image DICOMs
        Img_DICOM_path = rf"{INNOVATE_path}\{PatNum}\scans" 

        # Find filename for b3000 image
        b3000_DICOM_fnames = glob.glob(f'{Img_DICOM_path}/*b3000_80/DICOM/*.dcm')
        
        # Test to throw error if empty list
        test_valid = b3000_DICOM_fnames[0]
 
    except:
        
        # Define path to image DICOMs
        Img_DICOM_path = rf"{INNOVATE_path}\{PatNum}" 

        # Find filename for b3000 image
        b3000_DICOM_fnames = glob.glob(f'{Img_DICOM_path}\*b3000_80\DICOM\*.dcm')


    # Test if DICOM MF or SF
    if len(b3000_DICOM_fnames) == 1:
        # MF
        multiframe = True
        b3000_DICOM_fname = b3000_DICOM_fnames[0]
        b3000dcm = DICOM.MRIDICOM(b3000_DICOM_fname)
        
    elif len(b3000_DICOM_fnames) > 1:
        # SF
        multiframe = False
        b3000dcm = DICOM.MRIDICOM(DICOM_fnames = b3000_DICOM_fnames, multiframe = multiframe)
        
    else:
        print(f'No b3000 file for patient {PatNum}')
        sys.exit()
        
        
    # Create dcm object 
    b3000_ImageArray = b3000dcm.constructImageArray()
    b3000_bVals = b3000dcm.DICOM_df['Diffusion B Value'].values
    b0 = b3000_ImageArray[b3000_bVals == 0]
    
    # Define path to ROI DICOM
    RTStruct_path = f'{INNOVATE_ROIs_path}/{PatNum}'

    # Find RTstruct filename
    RTStruct_fname = glob.glob(f'{RTStruct_path}/*')[0]


    # === Create ROI mask

    # Instantiate contours object
    contours = DICOM.contours(RTStruct_fname)

    # Define lesion structure number (hardcoded here but should be found automatically in future)
    LesionStructNum = contours.Struct_Name_Num_dict[ROIName]
 

    LesionMask = contours.create_mask(Struct_Num = LesionStructNum, DICOM_dcm = b3000dcm)


    # Remove duplicate spatial slices
    LesionMask = LesionMask[b3000_bVals == 0]
    
    
    
    '''New code: accounting for z flips of fIC map'''
    
    # == First, load in b0 data from Matlab output
    
    # # Load .mat file
    # matb0 = loadmat(f'VERDICT outputs/{PatNum}/Model 1/b0from3000.mat')['b0from3000']
    
    # # Permute matlab b0 volume
    # matb0 = np.moveaxis( matb0 , -1, 0)  
    

    # # == Calculate MSE for different relative orientations
    # MSE0 = mean_squared_error(matb0, b0)
    # MSE1 = mean_squared_error(matb0, b0[::-1, : , :])

    # # If MSE0 > MSE1, python b3000 image and mask need to be flipped to match fIC
    # if MSE0 > MSE1:
    #     print('oof')
    #     LesionMask = LesionMask[::-1,:,:]
    #     b0 = b0[::-1,:,:]
    # else:
    #     None
        
    
    # == Save lesion mask
    try:
        os.makedirs(f'{output_path}/{PatNum}')
    except:
        None
      
    # Save as npy  
    np.save(f'{output_path}/{PatNum}/{ROIName}.npy', LesionMask)
    
    # Save as mha
    sitk.WriteImage( sitk.GetImageFromArray(LesionMask), f'{output_path}/{PatNum}/{ROIName}.mha' )
    
    # == Save b=0 from b3000 image
    
    # Save as npy  
    np.save(f'{output_path}/{PatNum}/b0.npy', b0)
    
    # Save as mha
    sitk.WriteImage( sitk.GetImageFromArray(b0), f'{output_path}/{PatNum}/b0.mha' )     
 
# Function for saving Nifti masks as npy arrays
def saveNiftimask(
    PatNum,
    Nifti_INNOVATE_ROIs_path = r"D:\UCL PhD Imaging Data\Nifti INNOVATE ROIs",
    output_path = r"C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\PhD_Year_1\ISMRM Submission\ROIs"):
    
    # Read Nifti image
    mask = sitk.GetArrayFromImage(
        sitk.ReadImage(f'{Nifti_INNOVATE_ROIs_path}/{PatNum}.nii.gz')
    )
    
    print(mask.shape)
    plt.figure()
    plt.imshow(mask[4])
    plt.show()
         
# Function for extracting fIC values from ROI
def extractROIfICs(
    PatNum,
    ROIName,
    ModelName,
    parameter = 'fIC',
    MaskType = 'numpy',
    VERDICT_output_path = r"C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\ISMRM Submission\Outputs\VERDICT outputs",
    ROI_path = r"C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\ISMRM Submission\Outputs\ROIs",
    output_path = r"C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\ISMRM Submission\Outputs\fIC results"
):
    
    '''
    Function to extract fIC values from within an ROI
    
    '''                
    
    
    # Load fIC array
    fIC = loadmat(f'{VERDICT_output_path}/{PatNum}/{ModelName}/{parameter}.mat')[parameter]
    
    # remove infinities
    fIC[fIC == np.inf] = 0
    
    # Remove nan
    fIC[np.isnan(fIC)] = 0
    
    # Permute array axes (account for MATLAB Python differences)
    fIC = np.moveaxis( fIC , -1, 0)   
    

    if MaskType == 'numpy':
        print('numpy')
        # Load ROI mask
        ROIMask = np.load(f'{ROI_path}/{PatNum}/{ROIName}.npy')
        
    elif MaskType == 'Analyze':
        
        ROIMask = sitk.GetArrayFromImage(sitk.ReadImage(f'{ROI_path}/{PatNum}/{ROIName}.img'))
        
    else:
        print('Incorrect mask type')
    
    
    # Make Bool
    ROIMask = (ROIMask != 0)
    
    # Extract fIC from ROI
    ROI_fIC = fIC[ROIMask]
    
    # Save as numpy array
    try:
        os.makedirs(f'{output_path}/{PatNum}/{ROIName}')
    except:
        None
        
    # Save as npy
    np.save(f'{output_path}/{PatNum}/{ROIName}/{ModelName}.npy', ROI_fIC)
    
    # Save fIC image as mha
    sitk.WriteImage(sitk.GetImageFromArray(fIC), f'{VERDICT_output_path}/{PatNum}/{ModelName}/{parameter}.mha')
    
# Function to read biopsy results     
def readBiopsyResults(
    biopsy_data_xlsx = r'C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\INNOVATE Data\INNOVATE patient groups 2.0.xlsx'
):
    
    '''
    Function to read biopsy data excel sheet and save results as binary dataframe
    
    '''
    
    BiopsyDF = pd.read_excel(biopsy_data_xlsx)
    
    # Clinically significant patients
    csPats = list( (BiopsyDF['Clinically significant (Gleason >=7)'].values)[1:] )
    
    
    # Ones (binary 1 for cs)
    ones = list(np.ones(len(csPats)))
    
    # Non-clincially significant patients 
    ncsPats = list( (BiopsyDF['Clinically insignificant (Gleason <7)'].values)[1:] )
    
    # Zeros (binary 0 for ncs)
    zeros = list(np.zeros(len(ncsPats)))
    
    
    # Make dataframe
    Patients = csPats + ncsPats
    Results = ones + zeros
    
    BiopsyResultsDF = pd.DataFrame({'Patient_ID': Patients, 'Biopsy_Result': Results})
    
    return BiopsyResultsDF
             
# Function to calculate avg fIC over ROIs
def avgROIfICs(
    PatNums, 
    ROIName,
    ModelName,
    avg_type = 'median',
    results_path = r"C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\ISMRM Submission\Outputs\fIC results"
    
):
    
    # Extract list of fIC results filenames
    fIC_fnames = glob.glob(f'{results_path}/*/{ROIName}/{ModelName}.npy')

    # For each file, extract patient number and calculate average fIC
    avgfICs = []
    thesePatNums = []
    
    for fname in fIC_fnames:
        
        # Find patient number
        PatNum = os.path.split( 
                               os.path.split(
                                   os.path.split(fname)[0]
                               )[0]
                               )[1]

        # Check if patient number in PatNums
        if PatNum not in PatNums:
            continue
        else:
            thesePatNums.append(PatNum)
            
        # Load fICs
        fICs = np.load(fname)
        
        # Calculate average
        if avg_type == 'mean':
            avg_fIC = np.mean(fICs)
        elif avg_type == 'median':
            avg_fIC = np.median(fICs)
        elif avg_type == 'UQ':
            avg_fIC = np.percentile(fICs, 75)
        elif avg_type == 'metric':
            avg_fIC = np.mean(fICs)/np.std(fICs)
        elif avg_type == 'max':
            avg_fIC = np.max(fICs)
        elif avg_type == 'contours':
            avg_fIC = contour_areas(PatNum, ModelName, ROIName)
            
        else:
            print('Incorrect input, default to mean')
            avg_fIC = np.mean(fICs)
            
        # Append to list
        avgfICs.append(avg_fIC)
        
     
       
    # Create dataframe
    fIC_DF = pd.DataFrame({'Patient_ID': thesePatNums, 'Avg_fIC': avgfICs})

    
    # Create directory
    try:
        os.makedirs(f'{results_path}/{avg_type} fIC Dataframes/{ModelName}')
    except:
        None
        
    # Save dataframe as pickle
    with open(f'{results_path}/{avg_type} fIC Dataframes/{ModelName}/{avg_type}_fIC_df.pickle', 'wb') as handle:
        pickle.dump(fIC_DF, handle, protocol=pickle.HIGHEST_PROTOCOL)
        
        
    # Save as matlab structure
    savemat(f'{results_path}/{avg_type} fIC Dataframes/{ModelName}/{avg_type}_fIC_df.mat', {'fICs': fIC_DF.to_dict(orient = 'list')})
    # # Save dataframe as txt
    # with open(f'{results_path}/Average fIC Dataframes/Model {ModelNum}/average_fIC_df.txt', 'w') as f:
    #     f.write(str(fIC_DF))
        
# Function to compute variance in fIC across ROI
def ROIvariance(
    PatNums,
    ROIName,
    ModelName,
    results_path = r"C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\ISMRM Submission\Outputs\fIC results"
):
    
    variances = []
    thesePatNums = []
    
    # Iterate over patient numbers
    for PatNum in PatNums:
        
        try:    
            # Load ROI fICs
            fICs = np.load(f'fIC results/{PatNum}/{ROIName}/{ModelName}.npy') 

            # Variance
            var = np.var(fICs)
            variances.append(var)
            
            # If success, append PatNum
            thesePatNums.append(PatNum)
            
        except:
            continue

    # Make array
    variances = np.asarray(variances)
    thesePatNums = np.asarray(thesePatNums)
    
    # Make dataframe
    fICvarDF = pd.DataFrame({'Patient_ID': thesePatNums, 'var_fIC': variances})
    
    # Create directory
    try:
        os.makedirs(f'{results_path}/variance fIC Dataframes/{ModelName}')
    except:
        None
        
    # Save dataframe as pickle
    with open(f'{results_path}/variance fIC Dataframes/{ModelName}/variance_fIC_df.pickle', 'wb') as handle:
        pickle.dump(fICvarDF, handle, protocol=pickle.HIGHEST_PROTOCOL)  
    
       
def fIC_ROC(
    ModelName,
    avg_type,
    results_path = r"C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\ISMRM Submission\Outputs\fIC results"
    ):
    
    '''
    Python function to generate ROC curve for lesion classification from 
    a specified model type
    
    General methods:
    
    1. Read in average fIC and biopsy results dataframes
    2. Create corresponding arrays of average fIC and biopsy outcomes (significant or insignificant)
    3. Use sklearn to generate ROC curve
    
    '''
    
    # Read in average fIC dataframe
    with open(f'{results_path}/{avg_type} fIC Dataframes/{ModelName}/{avg_type}_fIC_df.pickle', 'rb') as handle:
        fIC_DF = pickle.load(handle)
        
    # Read in biopsy results dataframe
    BiopsyResults_DF = readBiopsyResults()

    # Construct list of common patients (in Biopsy DF and fIC DF)
    fIC_PatList = fIC_DF['Patient_ID'].values
    Biopsy_PatList = BiopsyResults_DF['Patient_ID'].values

    PatList = fIC_PatList[ np.array([Pat in Biopsy_PatList for Pat in fIC_PatList]) ]
    
    # Construct arrays for biopsy results and fIC
    BiopsyResults = []
    fICs = []
    
    for PatNum in PatList:
        # Extract fIC
        fIC_Bools = (fIC_DF['Patient_ID'].values == PatNum)
        fICs.append(fIC_DF['Avg_fIC'].values[fIC_Bools][0])
        # Extract biopsy result
        Biopsy_Bools = (BiopsyResults_DF['Patient_ID'].values == PatNum)
        BiopsyResults.append(BiopsyResults_DF['Biopsy_Result'].values[Biopsy_Bools][0])
        
    # Make arrays
    fICs = np.asarray(fICs)    
    BiopsyResults = np.asarray(BiopsyResults)
    
    print(np.sum(BiopsyResults))
    
    # Make dataframe
    ResultsDF = pd.DataFrame({'Patient_ID': PatList, 'Biopsy_Result': BiopsyResults, 'Average_lesion_fIC': fICs })
        
    
    # Create ROC curve
    fpr, tpr, thresholds = roc_curve(y_true = BiopsyResults, y_score = fICs)
    
    # Calculate roc_auc_score
    roc_auc = roc_auc_score(y_true = BiopsyResults, y_score = fICs)
    
    return ResultsDF, fpr, tpr, thresholds, roc_auc


# Function to construct confusion table
def confusionTable(ROIdrawer, ModelNum, avg_type):
    
    # Load results dataframe
    with open(rf'ROC results\{ROIdrawer}\Model {ModelNum}\{avg_type}\Results_DF.pickle', 'rb') as handle:
        ResultDF= pickle.load(handle)

    # Load tpr
    tpr = np.load(rf'ROC results\{ROIdrawer}\Model {ModelNum}\{avg_type}\tpr.npy')
    fpr = np.load(rf'ROC results\{ROIdrawer}\Model {ModelNum}\{avg_type}\fpr.npy')

    # Load thresholds
    thresholds = np.load(rf'ROC results\{ROIdrawer}\Model {ModelNum}\{avg_type}\thresholds.npy')
    
    # Find Youden indx
    Yindx = np.where( (tpr-fpr)== np.max(tpr-fpr))[0][0]
    threshold = (thresholds[Yindx])
    sensitivity = (tpr[Yindx])
    specificity = 1 - (fpr[Yindx])

    # # Find threshold for 90% sensitivity
    # indx = np.sum(tpr<=0.9)
    # threshold = 0.39#0.5*(thresholds[indx-1]+thresholds[indx])
    # specificity = 1 - 0.5*(fpr[indx-1]+fpr[indx])
    
    print(threshold)
    print(sensitivity)
    print(specificity)
    
    # == At this threshold, find...
    
    # ... Number of true positives
    NTP = np.sum( (ResultDF['Average_lesion_fIC'].values >= threshold )*(ResultDF['Biopsy_Result'].values ) )
    
    # ... Number of false positives
    NFP = np.sum( (ResultDF['Average_lesion_fIC'].values >= threshold )*(1 - ResultDF['Biopsy_Result'].values ) )

    # ... Number of true negatives
    NTN = np.sum( (ResultDF['Average_lesion_fIC'].values <= threshold )*(1 - ResultDF['Biopsy_Result'].values ) )
    
    # ... Number of false negatives
    NFN = np.sum( (ResultDF['Average_lesion_fIC'].values <= threshold )*(ResultDF['Biopsy_Result'].values ) )

    print(NTP, NFP, NTN, NFN)

# confusionTable('NT', 2, 'median')


def contour_areas(
    PatNum, 
    ModelNum,
    ROIName = 'L1_b3000_NT',
    parameter = 'fIC',
    VERDICT_output_path = r"C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\PhD_Year_1\ISMRM Submission\VERDICT outputs",
    ROI_path = r"C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\PhD_Year_1\ISMRM Submission\ROIs"):
    
    # Function to experiment with contour drawing idea
    
    # First, read in fIC map and mask
   
    fIC = loadmat(f'{VERDICT_output_path}/{PatNum}/Model {ModelNum}/{parameter}.mat')[parameter]
    
    # remove infinities
    fIC[fIC == np.inf] = 0
    
    # Remove nan
    fIC[np.isnan(fIC)] = 0
    
    # Permute array axes (account for MATLAB Python differences)
    fIC = np.moveaxis( fIC , -1, 0)   
    
    ROIMask = np.load(f'{ROI_path}/{PatNum}/{ROIName}.npy')
    
    # Mulitply mask and fIC
    ROI = fIC*ROIMask
    
    sliceIndx = np.where(np.sum(ROI, axis = (1,2)) != 0)[0][0]
    
    ROIslice = ROI[sliceIndx]
    
    # Smooth
    ROIslice = scipy.ndimage.gaussian_filter(ROIslice, 1)
    
    # plt.figure
    # plt.imshow(ROI[sliceIndx])
    # plt.show()
    
    
    # Levels
    levels = np.arange(0,1, 0.02)
    Areas = []
    for level in levels:
        
        contours =  skimage.measure.find_contours(ROIslice, level)
        
        
        # Choose largest contour
        try:
            lens = np.array([len(cont) for cont in contours])
            maxindx = np.where(lens == np.max(lens))[0][0]
            contour = contours[maxindx]
        
        except:
            Area = 0
            Areas.append(Area)
            continue
        
        # Make path
        cont_path = mpltPath.Path(contour)
        
        
        # For efficiency, only query slice coordinates near contour
        cont_xmin = np.min( contour[...,0] ) - 1
        cont_xmax = np.max( contour[...,0] ) + 1
        cont_ymin = np.min( contour[...,1] ) - 1
        cont_ymax = np.max( contour[...,1] ) + 1
        
        # Construct coordinate array
        dx = 0.1
        xs = np.arange(cont_xmin, cont_xmax, dx)
        ys = np.arange(cont_ymin, cont_ymax, dx)
   
        Ys, Xs = np.meshgrid(ys, xs, indexing = 'ij')
        Coords = np.stack((Xs, Ys), axis = -1)

        Area = 0
        
        # Query each point
        for xindx in range(len(xs)):
            for yindx in range(len(ys)):
                
                point = Coords[yindx, xindx,:]
                
                if cont_path.contains_point(point):
                    Area += dx**2
             
        # print(Area)   
        Areas.append(Area)     

        # fig, ax = plt.subplots()
        # ax.imshow(ROIslice, cmap=plt.cm.gray)

        # for contour in contours:
        #     ax.plot(contour[:, 1], contour[:, 0], linewidth=2)


        # ax.axis('image')
        # ax.set_xticks([])
        # ax.set_yticks([])
        # plt.show(
    # Integrate area under graph
    if Areas[0]!=0:
        Areas = np.array(Areas)/Areas[0]   
    else:
        Areas = np.array(Areas)
        
    Score = scipy.integrate.trapezoid(Areas*levels, levels)
    
    return Score
    
    # plt.figure()
    # plt.plot(levels, Areas)
    # plt.show()
    
    

