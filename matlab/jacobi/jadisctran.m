% The function s=jadisctran(n,alp,bet,x,w,p,iflag) performs the discrete Jacobi transforms 
% between the physical space (i.e., physical values) and frequency space
% (Jacobi expansion coefficients) at the Jacobi-Gauss-Lobatto points 
% Input:
%  n,alp,bet,x,w--- number of JGL points in x, where (x,w) can be computed by
%          [x,w]=jagslb(n,alp,bet). Note: x,w are column vectors  
%  iflag==0--- forward transform  
%    p--- (input) physical values at collocation points
%    s--- (output) expansion coefficients 
%  iflag not= 0--- backward transform  
%    p--- (input) expansion coefficients 
%    s--- (output) physical values at collocation points 
%
%  See Page 87 of the book: J. Shen, T. Tang and L. Wang, Spectral Methods:
%  Algorithms, Analysis and Applications, Springer Series in Compuational
%  Mathematics, 41, Springer, 2011. 
%  Use the function: lepolym()   
%  Last modified on August 31, 2011


function s=jadisctran(n,alp,bet,x,w,p,iflag)
 T=japolym(n-1,alp,bet,x); % compute the Legndre polynomials up to order n-1. Note: T(i,j)=J_{i-1}(x_j) 
 if iflag==0, 
     nv=[0:n-1]'; fgmn=2^(-(alp+bet+1))*exp(-gammaln(nv+alp+1)-gammaln(nv+bet+1)... 
         +gammaln(nv+1)+gammaln(nv+alp+bet+1)).*(2*nv+alp+bet+1);
     s=(T*(p.*w)).*[fgmn(1:end-1);fgmn(end)/((2+(alp+bet+1)/(n-1)))]; % see (3.156)
     return;
 end
 
 s=T'*p;  % see (3.155)
 return
 
 

