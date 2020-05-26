# %% Libraries ----------------------
import gl
import os

# %% General----------------------
gl.resetdefaults()
gl.backcolor(255, 255, 255)
gl.colorbarvisible(0)

# %% background image ----------------------
gl.loadimage('mni152')

# %% overlay ----------------------
filepath='D:\\GD_UNICAMP\\IC_NeuroFisica\\Projetos\\fMRI_TLE_2020\\Seed_based_analysis'
os.chdir(filepath)
#gl.overlayload('Shen2013')
gl.overlayload('spmMotor')
#gl.opacity(1,50)

# %% #changing view ----------------------
gl.azimuthelevation(140,20)
ktime= 1
ksteps= 72
for x in range(0, ksteps):
  gl.azimuthelevation(140+(x*5),20)
  gl.wait(ktime)
  #gl.savebmp(str(x)+'.png')
