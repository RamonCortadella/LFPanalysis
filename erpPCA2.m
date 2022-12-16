% erpPCA - Unrestricted, unstandardized covariance-based PCA with Varimax rota-
%          tion (cf. Kayser J, Tenke CE, Clin Neurophysiol, 2003, 114:2307-25)
%
% Usage: [LU, LR, FSr, VT] = erpPCA( X, nFac, maxIt, Tol, IfVerbose )
%
% Generic PCA-Varimax implementation emulating the PCA agorithms used by 
% BMDP-4M (Dixon, 1992) and SPSS 10.0 FACTOR (http://www.spss.com/tech/
% stat/Algorithms/11.5/factor.pdf). It expects a data matrix (X) of ERP 
% waveforms, with ERPs (cases) as rows and sample points (variables) as 
% columns. The routine returns the unrotated (LU) and Varimax-rotated (LR) 
% factor loadings as a variables-by-factors matrix, the rotated factor 
% scores (FSr) as a cases-by-factors matrix, and Eigenvalues and explained 
% variance as a variables-by-variance matrix (VT), with four columns 
% consisting of Eigenvalues and percentage of explained variance before 
% and after rotation.
%
% erpPCA employs Varimax4M (max. 100 iterations, 0.0001 convergence criterion,
% Kaiser's normalization; MatLab code by $jk available on request), which 
% emulates algorithms described by Harman (1967, pp. 304-308) as implemented 
% in BMDP-4M (Dixon, 1992, pp. 602-603).
%
% Copyright (C) 2003 by J�rgen Kayser (Email: kayserj@pi.cpmc.columbia.edu)
% GNU General Public License (http://www.gnu.org/licenses/gpl.txt)
% Updated: $Date: 2003/07/08 14:00:00 $ $Author: jk $
% Modfied and optimized by Anton. June 2007
function [LU,LR,FSr,VT] = erpPCA(X, varargin)

[nFac, maxIt, Tol, IfVerbose] = DefaultArgs(varargin,{[], 100,1e-3,1});

% get dimensions of input data matrix     
[cases, vars] = size(X);             
% compute covariance matrix               
D = cov(X);                          
% determine Eigenvectors and Eigenvalues  
[EM, EV] = eig(D);                   

% determine unrotated factor loadings     
UL = EM * sqrt(EV);                  

% sort initial Eigenvalues, keep indices  
[u, ux] = sort(diag(EV)','descend');        

% sort unrotated factor loadings          
LU = UL(:,ux);                       

% estimate the number of singular values  
try 
    rk = rank(corrcoef(X),1e-4);       % why not on D? would save on corrcoef(X) computation
catch
    rk = rank(D,1e-4);
end
% remove all linearly dependent components and their indices
LU = LU(:,1:rk);                     
u = u(1:rk); ux = ux(1:rk);          

% current sign of loading vectors        
s = ones(1,rk);                      

% determine direction of loading vectors 
s(abs(max(LU)) < abs(min(LU))) = -1; 

% redirect loading vectors if necessary  
LU = LU .* repmat(s,size(LU,1),1);   

% Varimax-rotate factor loadings         
RL = Varimax4M(LU,maxIt,Tol,1);       

% compute rotated Eigenvalues            
EVr = sum(RL .* RL);                 

% sort rotated Eigenvalues, keep indices 
[r, rx] = sort(EVr,'descend');                 

% sort rotated factor loadings           
LR = RL(:,rx);                       

% current sign of loading vectors        
s = ones(1,size(LR,2));              

% determine direction of loading vectors 
s(abs(max(LR)) < abs(min(LR))) = -1; 

% redirect loading vectors if necessary  
LR = LR .* repmat(s,vars,1);        

% compute total variance                 
tv = trace(EV);                      

% table explained variance for unrotated and Varimax-rotated components
VT = [u' 100*u'/tv r' 100*r'/tv ];              

% compute rotated FS coefficients        
FSCFr = LR * inv(LR' * LR);          % this is pseudo-inverse of LR

% rescale rotated FS coefficients by the corresponding SDs 
FSCFr = FSCFr .* repmat(sqrt(diag(D)),1,rk);       

if isempty(nFac) nFac = size(FSCFr,2); end
    
% compute rotated factor scores from the normalized raw data and  the
% corresponding rescaled factor score coefficients
FSr = zscore(X) * FSCFr(:,1:nFac);
LR = LR(:,1:nFac);
LU = LU(:,1:nFac);

if nargout==1
    out.LU=LU;
    out.LR=LR;
    out.FSr = FSr;
    out.VT= VT;
    LU=out;
end
   
