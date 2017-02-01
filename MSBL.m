function [ Xrec, gammas, sigmaSqs, converged ] = MSBL_huji( H, Y, iters, SNR, varargin  )
%MYMSBL 
% ===== INPUTS =====:
%   H           : Dictionary matrix, size M x N
%   Y           : Measurements vector, matrix of size N x L, M = Total pixels in image, L = number of
%                 measurements
%   iters       : Number of iterations to perform
%   SNR         : Estimation of the signal-to-noise ratio
%   centerY     : boolean variable indicates if the measured data should be
%                 temporally averaged (default: True)
%   updateSigma : boolean variable indicates if the signal-to-noise ratio
%                 should be re-estimated after each iteration (default: False). 
%   showImage   : boolean variable, indicates if the results after each
%                 iteration should be displayed (default: True)
%   updateRule  : which update rule to use "fast" for equation (19) from
%                 Wipf's "An Empirical Bayesian strategy for solving the Simultaneous Sparse Approximation problem"
%                 or "slow" for equation (18)
%   gpu         : boolean variable, use GPU for calculations.
%   gamma       : Initial value for gamma
%   
% ===== OUTPUTS =====
%   Xrec        : The recovered object. Size N x L
%   gammas      : Value of gamma in each iteration. Size N x <iters> 
%   sigmaSqs    : Noise approximation values for each iteration
%   converged   : Values of gamma which converged in the last iteration. 
%

centerY = 1;
updateSigma = 0;
showImage = 1;
dataType = 'single';
updateRule = 'fast';
tol = 1e-14;
gpu = 0;


[M, N] = size(H);
gamma = ones(N,1,dataType);

for i=1:length(varargin)
    if mod(i,2) == 0; continue; end;
    switch varargin{i}
        case 'centerY'
            centerY = varargin{i+1};
        case 'updateSigma'
            updateSigma = varargin{i+1};
        case 'showImage'
            showImage =  varargin{i+1};
        case 'dataType'
            dataType = varargin{i+1};
        case 'updateRule'
            updateRule = varargin{i+1};
        case 'gpu'
            gpu = varargin{i+1};
        case 'gamma'
            gamma = varargin{i+1};
        otherwise
            error('Parameter %s unknown', varargin{i});
    end    
end

% Cast data type
eval(['H=',dataType,'(H);']);
eval(['Y=',dataType,'(Y);']);

S = max(mean(abs(Y),2));
sigmaSq = S/SNR * 1.1;

if centerY
    Y = Y - repmat(mean(Y,2),[1 size(Y,2)]);
end

tic;
if nargin <= 2
    iters = 15;
    sigmaSq = 0;
end
 
L = size(Y, 2);
converged = zeros(N, 1, dataType);
mu = zeros(N, L, dataType);
keep_list = [1:N]';
gammas = zeros(N, iters, dataType);
sigmaSqs = zeros(1, iters+1, dataType);
sigmaSqs(1) = sigmaSq;

if gpu
    H = gpuArray(H);
    gamma = gpuArray(gamma);
end

if showImage; figure('Position', [300, 500, 800, 600]); end;
 
for count = 1:iters
    try
        if showImage                        
            imagesc(reshape(gamma, sqrt(N), sqrt(N))); colorbar;axis image;
            title(['gammas iter:', num2str(count)]);            
            drawnow();            
        end
    catch
        display('Could not print image');
    end
    % ****** estimate the solution matrix *****
    Gamma = diag(gamma);
    Sigma_t = H*Gamma*H' + sigmaSqs(count)*eye(M);  
    Xi = Gamma*H'*pinv(Sigma_t, tol);
    
    mu = Xi * Y;
    Sigma_diag = gamma - (sum(Xi'.*(H*Gamma)))'; % sum(A.*B) = A*B        
    mu2_bar = sum(abs(mu).^2,2)/L;
    
    gamma_old = gamma; 
    if strcmp(updateRule, 'slow')
        gamma = mu2_bar + Sigma_diag;                  
    else
        gamma = mu2_bar ./(1 - Sigma_diag./gamma_old);                  
    end

    % ****** estimate the noise parameter *****
    if sigmaSqs(count) ~= 0
        if updateSigma
            sigmaSqs(count+1) = (norm(Y - H * mu,'fro')^2/L)/(M-N + sum(Sigma_diag./(gamma_old))) ;
            if isnan(sigmaSqs(count+1))
                sigmaSqs(count+1) = sigmaSqs(count);
            end
        else 
            sigmaSqs(count+1) = sigmaSqs(count);
        end
    end    
    
    if (sum(isinf(gamma(:))+isnan(gamma(:))) ...
            + sum(isinf(mu(:))+isnan(mu(:))))
        converged(isinf(gamma(:))) = 1;
        converged(isnan(gamma(:))) = 1;
         gamma(isinf(gamma(:))) = gamma_old(isinf(gamma(:)));
         gamma(isnan(gamma(:))) = gamma_old(isnan(gamma(:)));         
    end    
    
    gammas(:,count) = gather(gamma);   
end;
 
mu = mu .* (mu > 0);    
Xrec = zeros(N,L,'single');
Xrec(keep_list,:) = gather(mu);
Xrec = Xrec(1:N,:);
toc;
end
 

