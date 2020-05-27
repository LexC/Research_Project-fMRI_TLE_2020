# %% Libraries ----------------------
import gl
import os

# %% General----------------------
#gl.resetdefaults()
gl.backcolor(255, 255, 255)

# %% background image ----------------------
gl.loadimage('mni152')

# %% overlay ----------------------
FILESPATH='D:\\GD_UNICAMP\\IC_NeuroFisica\\Projetos\\fMRI_TLE_2020\\Seed_based_analysis'
os.chdir(FILESPATH)

NETFILE_NAME = 'DMN'
netfile_sufx= ['_Group1X','_Group2L','_Group3R','_Group4N']

for i in netfile_sufx:

    FILENAME=NETFILE_NAME+i
    SEEDFILENAME=NETFILE_NAME+"_seed"

    gl.overlaycloseall()
    gl.overlayload(FILENAME)
    gl.overlayload(SEEDFILENAME)
    gl.colorname (2,"3blue")
    
    #gl.opacity(1,50)
    
    # changing view ----------------------
    gl.azimuthelevation(140,20)
    ktime= 1
    ksteps= 72
    for x in range(0, ksteps):
        gl.azimuthelevation(140+(x*5),20)
        gl.wait(ktime)

        if x < 10:
            gl.savebmp(FILESPATH + '\\' + FILENAME + '_0' + str(x)+'.png')
        else:
            gl.savebmp(FILESPATH + '\\' + FILENAME + '_' + str(x)+'.png')

