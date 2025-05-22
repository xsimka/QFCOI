close all
clear all

ref_image = imread('ref_Telescope.bmp');
test_image = imread('Telescope46AVIF.bmp');


result = QFCOI(ref_image,test_image)
