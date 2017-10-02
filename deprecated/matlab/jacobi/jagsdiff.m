% D=jagsdiff(n,x,alp,bet) computes the first-order differentiation matrix of size
% n by n, associated with the Jacobi-Gauss points x, which can be computed by 
% x=jags(n,alp,bet) 
% See Page 91 of the book: J. Shen, T. Tang and L. Wang, Spectral Methods:
%  Algorithms, Analysis and Applications, Springer Series in Compuational
%  Mathematics, 41, Springer, 2011. 
%  Use the function: japoly() 
% Last modified on September 4, 2011

function  D=jagsdiff(n,x,alp,bet)

if n==0, D=[]; return; end;
xx=x;[dy,y]=japoly(n,alp,bet,xx); nx=size(x); 
if nx(2)>nx(1), dy=dy'; y=y'; xx=x'; end;  %% dy is a column vector 
  D=(xx./dy)*dy'-(1./dy)*(xx.*dy)';  %% compute J_{n}^{alp,bet}'(x_j) (x_k-x_j)/J_{n}^{alp,bet}'(x_k);     
                                 % 1/d_{kj} for k not= j (see (3.164)) 
  D=D+eye(n);                    % add the identity matrix so that 1./D can be operated                                     
  D=1./D; 
  D=D-eye(n); D
  D=D+diag(0.5*((alp-bet)+(alp+bet+2)*xx)./(1-xx.^2));  % update the diagonal entries  
  return; 


