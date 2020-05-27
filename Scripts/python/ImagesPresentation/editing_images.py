# -*- coding: utf-8 -*-
"""
Created on Wed May 27 18:47:42 2020

@author: alexc
"""

# -*- coding: utf-8 -*-
"""

"""
# %% Libraries

import os
# import numpy as np
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

    GROUPNUMB=int(len(imagesfiles)/n)
    imgs=[[]]
    for i in range(GROUPNUMB):
        for j in range(n):
            imgs[i] += [Image.open(imagesfiles[i + j*GROUPNUMB])]
        imgs.append([])

    return imgs[:-1]


def comb_images_h(imgs):

    widths, heights = zip(*(i.size for i in imgs))
    
    total_width = sum(widths)
    max_height = max(heights)
    
    new_im = Image.new('RGB', (total_width, max_height))
    
    x_offset = 0
    for im in imgs:
        new_im.paste(im, (x_offset,0))
        x_offset += im.size[0]

    return new_im


def new72imgs(imgs):

    new_img=[]
    for x in imgs:
        new_img += [comb_images_h(x)]

    return new_img

def saving_imgs(new_imgs,imagesfiles):

    SAVENAME=imagesfiles[0][:-9]+'_'

    for i,x in enumerate(new_imgs):
        if i < 10:
            x.save(SAVENAME+'0'+str(i)+'.png')
        else:
            x.save(SAVENAME+str(i)+'.png')


# %% Variables and 

IMAGESDIR = "D:\\GD_UNICAMP\\IC_NeuroFisica\\Projetos\\fMRI_TLE_2020\\Seed_based_analysis"


# %% Main

imagesfiles = get_images(IMAGESDIR)
imgs = open_imgs_ingroups(imagesfiles,4)
new_imgs = new72imgs(imgs)
saving_imgs(new_imgs,imagesfiles)
