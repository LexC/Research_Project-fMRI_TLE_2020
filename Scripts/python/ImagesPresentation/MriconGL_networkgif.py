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

networks_names=['DMN','SMOTOR','VISUAL','AUDIT','ECN_L','ECN_R']
netfile_sufx= ['_Group1X','_Group2L','_Group3R','_Group4N']

for net in networks_names:
    savefolder=FILESPATH + '\\' + net + '\\'
    for i in netfile_sufx:
    
        FILENAME=net+i
        SEEDFILENAME=net+"_seed"
    
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
                gl.savebmp(savefolder + FILENAME + '_0' + str(x)+'.png')
            else:
                gl.savebmp(savefolder + FILENAME + '_' + str(x)+'.png')

