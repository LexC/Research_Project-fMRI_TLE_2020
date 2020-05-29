# %% Libraries ----------------------
import gl
import os

# %% General----------------------

gl.backcolor(255, 255, 255)

# %% background image ----------------------
gl.loadimage('mni152')

# %% overlay ----------------------
FILESPATH='D:\\GD_UNICAMP\\IC_NeuroFisica\\Projetos\\fMRI_TLE_2020\\Seed_based_analysis'
os.chdir(FILESPATH)

#networks_names=['DMN','SMOTOR','VISUAL','AUDIT','DAN']
networks_names=['DAN_R','DAN_L']
netfile_sufx= ['_Group1X','_Group2L','_Group3R','_Group4N']

for net in networks_names:
    SEEDFILENAME= net+'_seed'
    for sufx in netfile_sufx:

        FILENAME=net+sufx
            
        gl.overlaycloseall()
        gl.overlayload(FILENAME)
        gl.overlayload(SEEDFILENAME)
        gl.colorname (2,"3blue")
        
        gl.azimuthelevation(0,90)
        
        gl.savebmp(FILESPATH+'\\'+FILENAME)