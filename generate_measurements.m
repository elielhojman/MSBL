function [ Y, U, Yclean ] = generate_measurements(obj, H, frames, speckleSz, SNR, jump, jumpType )
%
% ===== INPUTS ======
%   obj        : The original object
%   H          : The Dictionary matrix
%   frames     : The number of frames to generate
%   speckleSz  : The speckle size (in pixels)
%   SNR        : SNR of the measurements
%   jump       : Reduce the images
%   jumpType   : Method to reduce images, 'bilinear', 'fft'
%
% ====== OUTPUT ======
%   Y          : Matrix N x L, L is number of frames. N Number of pixels in
%                image
%   U          : Matrix N x L, speckle illumination patterns
%   Yclean     : Y vector without the added noise
%

if nargin < 6
    jump = 1;    
end
if nargin < 7
    jumpType = 'bilinear';
end


Y = zeros(size(H,1), frames, 'single');
U = zeros(size(H,1), frames, 'single');
for frame = 1:frames
    I = genSpeckles(size(obj), speckleSz);
    U(:,frame) = I(:);
    x = obj.*I;    
    Y(:,frame) = H*x(:);
end

if jump > 1
    n = sqrt(size(Y,1));
    Y1 = [];
    for i = 1:frames;        
        if strcmp(jumpType, 'bilinear')
            imS = imresize(reshape(Y(:,i),n,n),1/jump,'bilinear');
        elseif strcmp(jumpType, 'skip')
            imS = interpft(reshape(Y(:,i),n,n), round(n/jump), round(n/jump));            
        else
            error('jumpType not recognized');
        end
        %imS = downsample(downsample(reshape(Y(:,i),n,n),jump)',jump)';
        Y1(:,i) = imS(:);        
    end
    Y = Y1;
end

S = max(mean(Y,2));
stdN = S / SNR;
Yclean = Y;
Y = Y + stdN * randn(size(Y));
Y = Y/mean(Y(:));
% Y = Y * size(Y,2) / mean(Y(:).^2);
% Yclean = Yclean * size(Yclean,2) / mean(Yclean(:).^2);

