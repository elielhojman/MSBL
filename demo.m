% USE EXAMPLE:
% Recovery of an object through the MSBL algorithm using 
%   noise photoacoustics measurements aquired over random 
%   speckle illumination. 
%   

L = 1000;  % Number of frames
SNR = 100; % signal-to-noise ratio
jump = 4;  % Ratio of image resize, used to speed the processing time

% Load the test image
load('leaf_image');
leaf = double(l2);

% Generate H matrix, using an acoustic PSF with full width at 70% of max of
% 11 pixels. 
H = createDicMat(116, 'psfSz', 11, 'type', 'Acoustic');

% Generate measurments
Y = generate_measurements(leaf, H, L, 1, SNR, jump);
% Generate dictionary to recover measurments
H4 = createDicMat(116, 'psfSz', 11, 'type', 'Acoustic', 'jump', jump);

[Xrec,g] = MSBL(H4, Y, 15, 200);
