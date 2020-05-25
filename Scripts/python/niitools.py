"""
Personalized functions to speed up fMRI processing.
"""

# %% Functions


def makeNiiCortexNetwork(basefilename, network, savename):
    """
    Change the ROI's IDvalue to the network's list values and save a .nii file
    with the new information.

    basefilename = the location of a nifti file, with ROI'S enumerated from
                    1 to n.
    network = a list | len(network)=n
    savename = a string with the save name of the new nifti file
    """

    import nibabel as nib

    img = nib.load(basefilename)

    new_data = changeROIsValues(img.get_fdata(), network)
    
    img = nib.Nifti1Image(new_data, img.affine, img.header)
    nib.save(img, savename)
    
    return img


def changeROIsValues(data, network):
    for i in range(1, int(data.max())+1):
        data[data == i] = network[i-1]
    return data

#def gifBrainNetwork(
# %% Tests

basefilename = 'D:\\GD_UNICAMP\\IC_NeuroFisica\\Projetos\\fMRI_TLE_2020\\Variables_and_Data_Info\\fMRI_Processing_Var\\Shen2013.nii'
network=[1 for i in range(1,238)]
makeBrainNetwork(basefilename,network,'test') 
