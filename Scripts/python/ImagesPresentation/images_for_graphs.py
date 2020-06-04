# -*- coding: utf-8 -*-
"""

to create a fade animation: https://note.nkmk.me/en/python-pillow-composite/

"""
# %% Libraries

import os
from PIL import Image


# %% Functions


def get_images(IMAGESDIR):
    """Return a list with the files location of .png files"""
    filesdir = []
    for root, subfolder, files in os.walk(IMAGESDIR):
        for filename in files:
            if filename[-4:] == '.png':
                filesdir += [os.path.join(root, filename)]
    return filesdir


def open_imgs_ingroups(imagesfiles,n):

    u=0
    imgs=[[]]

    for i, imgdir in enumerate(imagesfiles,start=1):
        imgs[u] += [Image.open(imgdir)]
        if i%(n)==0:
            u+=1
            imgs.append([])

    return imgs[:-1]


def make_gif(imgs_set,SAVE_NAME,time_per_frame):

    im=imgs_set[0]

    imgs_subset = [imgs_set[i] for i in range(1,len(imgs_set))]

    im.save(SAVE_NAME, save_all=True, append_images=imgs_subset,loop=0,
            duration=time_per_frame*len(imgs_subset)+1)


def make_all_gifs(imgs,imagesfiles,time_per_frame):

    for i,imgs_set in enumerate(imgs):

        SAVE_NAME = str()
        imgdir_aux = imagesfiles[i*len(imgs_set)].split('\\')

        for j in range(0,len(imgdir_aux)-1):
            if j == len(imgdir_aux)-2:
                SAVE_NAME += imgdir_aux[j] + '.gif'
            else:
                SAVE_NAME += imgdir_aux[j] + '\\'

        make_gif(imgs_set,SAVE_NAME,time_per_frame)




# %% Variables and Inputs

imgfolder = 'D:\\GD_UNICAMP\\IC_NeuroFisica\\Projetos\\fMRI_TLE_2020\\Graphs\\Step1'

time_per_frame = 100

# %% Main

imagesfiles = get_images(imgfolder)
imgs = open_imgs_ingroups(imagesfiles, 13)


make_all_gifs(imgs,imagesfiles,time_per_frame)

#%% testes

# for i,x in enumerate(imgs[20]):
#     x.show()