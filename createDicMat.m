function [ H, PSF, Hclean ] = createDicMat( n, varargin )
%CREATEDICMAT Summary of this function goes here
%   n - is the size of the object we would like to recover
%   psfSz - the width of the PSF in pixels. Usually the number of pixels
%            which account for the standard deviation.
%   type - Acoustic or Gaussian
%   SNR - The signal to noise ratio that will be in the PSF
%   matrix. We have more columns than rows in H.

type = 'Gaussian';
SNR = Inf;
name = 'H';
oneD = 0;
jump = 1;
jumpType = 'bilinear';

for i=1:length(varargin)
    if mod(i,2) == 0; continue; end;
    switch varargin{i}
        case 'PSF'
            PSF = double(varargin{i+1});            
        case 'type'
            type = varargin{i+1};
        case 'psfSz'
            psfSz = varargin{i+1};            
        case 'SNR'
            SNR = varargin{i+1};
        case 'jump'
            jump = varargin{i+1};
        case 'jumpType'
            jumpType = varargin{i+1};
        case 'name'
            name = ['Sz' num2str(n(1)) '_' varargin{i}];
        case '1D'
            oneD = varargin{i+1};
        otherwise
            error('Parameter %s unknown', varargin{i});            
    end
end          

if oneD == 0 
    if isscalar(n) 
        n = [n n];
    end

    if ~(strcmp(type, 'Gaussian') || strcmp(type, 'Acoustic'))
        error('Invalid Dictionary type. Type must be "Gaussian" or "Acoustic"');
    end

    if length(n) > 2
        error('n should be a scalar or a 2D dimension size');
    end

    N = prod(n); 
else
    N = n;
end

H = zeros(N,'single');
%% PSF 
if ~exist('PSF', 'var')           
    [X, Y] = meshgrid(-(n(1)-1):n(1), -(n(2)-1):n(2));
    R = X.^2 + Y.^2;
    PSF = exp(-R/(psfSz.^2));    
    
    if strcmp('Acoustic', type) 
        PSF=PSF-circshift(PSF,[1 0]); % 1st axial derivative for PSF
        PSF=PSF-circshift(PSF,[1 0]); % 2nd axial derivative for PSF
        PSF=PSF/max(PSF(:));    
        PSF = -PSF;
    end
    % figure('Position',[100, 100, 200, 200]); imagesc(PSF(25:74,25:74)); axis image; title('PSF'); shg;
    %% Add noise
    stdN = max(PSF(:)) / SNR;
    PSF = PSF + stdN * randn(size(PSF)); 
    
    %% Check if the dictionary already exists
    name = ['Sz' num2str(n(1)) '_' num2str(n(2)) 'PSF' num2str(psfSz) type '_SNR' num2str(SNR)];
end
PSF = PSF ./ max(PSF(:));
if exist('dicMats.mat','file')
    load dicMats
    if isfield(dicMats,name)
        H = dicMats.(name);                
        return
    end
end

if jump > 0   
    % jump = ceil(psfSz/3);
    tmp = zeros(n);
    tmp = tmp(1:jump:end,1:jump:end);
    H = zeros(numel(tmp),N);
end

%H = single(H);
parfor i = 1:N;    
    if oneD == 0
        point = zeros(n,'single');
    else
        point = zeros(n,1,'single');
    end
    point(i) = 1;
    Htmp = imfilter(point, PSF,'same','conv');
    if jump > 0    
        if strcmp(jumpType, 'bilinear');            
            Htmp = imresize(Htmp,1/jump,'bilinear');
        elseif strcmp(jumpType, 'skip');            
            tmp = resample(double(Htmp),1,jump);
            Htmp = (resample(tmp',1,jump))';
        else
            error('jumpType not identified');
        end        
    % Unify close pixels according to the PSF. 
    % If the PSF has a width of 18 pixels, divide the 18x18 area into 9 regions.     
    end
    H(:,i) = single(Htmp(:));
end



if ~exist('dicMats','var')
    dicMats = {};
end

% H = H.*threshH;
dicMats.(name) = H;
% save('dicMats','dicMats');
H = single(H);
display('Dictionary created');
end

