# -*- coding: utf-8 -*-
"""
Created on Mon May  2 11:17:29 2022

@author: chiec
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Apr  6 14:14:55 2021

@author: chiec
"""
import numpy as np
import cv2 as cv
#import matplotlib.pyplot as plt
#import time



runfile('E:/Chieco/Python/SessionClass.py', wdir='E:/Chieco/Python')

general_folder='E:/Chieco/BellybuttonBase_WetFoam/'
folder_for_images=general_folder+'train_images/'
folder_for_masks=general_folder+'masks/'

im_filename='DSC_0001_train'


Session_name=folder_for_images+im_filename+'.jpg'
make_mask = Session(Session_name)



mask_filename=im_filename
mask_image_filename=folder_for_masks+mask_filename+'.jpg'

mask_in_read = cv.imread(mask_image_filename)  
bw_img = cv.cvtColor(mask_in_read, cv.COLOR_BGR2GRAY)
mask_in=np.array(bw_img,dtype=float)


make_mask.mask = mask_in

make_mask.StartDraw()


#this is done periodically so we don't lose progress for a big file

mask_textfile=mask_filename+'.txt'
make_mask.SaveMask(mask_textfile)
mask_out = make_mask.mask

mask_imfile=mask_image_filename
cv.imwrite(mask_imfile, mask_out * 255)











