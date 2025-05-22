close all
clear all

% Images are drawn from the OMNIQAD database available here: {https://zenodo.org/records/7607071}
ref_image = imread('ref_Telescope.bmp');
test_image = imread('Telescope46AVIF.bmp');


result = QFCOI(ref_image,test_image)
