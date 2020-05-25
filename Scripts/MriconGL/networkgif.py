import gl
gl.resetdefaults()
gl.backcolor(255, 255, 255)

#open background image
gl.loadimage('mni152')

#open overlay
#filepath='D:\\GD_UNICAMP\\IC_NeuroFisica\\Projetos\\fMRI_TLE_2020\\Variables_and_Data_Info\\fMRI_Processing_Var\\Shen2013base.nii'
#gl.overlayload(filepah)
gl.overlayload('spmMotor')
#gl.opacity(1,50)


#changing view
gl.azimuthelevation(140,20)
ktime= 1
ksteps= 72
for x in range(0, ksteps):
  gl.azimuthelevation(140+(x*5),20)
  gl.wait(ktime)
  #gl.savebmp(str(x)+'.png')
