close all
clear all

ref_image = imread('ref_Flowers.bmp');
test_image = imread('Flowers51AVIF.bmp');


result = QFCOI(ref_image,test_image)
