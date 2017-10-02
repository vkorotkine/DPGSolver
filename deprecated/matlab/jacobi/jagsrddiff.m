% D=jagsrddiff(n,x,alp,bet) computes the first-order differentiation matrix of size
% n by n, associated with the Jacobi-Gauss-Radau points x (with x(1)=-1), which can be computed by 
% x=jagsrd(n,alp,bet) 
% See Pages 89-90 of the book: J. Shen, T. Tang and L. Wang, Spectral Methods:
%  Algorithms, Analysis and Applications, Springer Series in Compuational
%  Mathematics, 41, Springer, 2011. 
%  Use the function: japoly() 
% Last modified on September 4, 2011

function  D=jagsrddiff(n,x,alp,bet)

if n==0, D=[]; return; end;
xx=x; [dy,y]=japoly(n-1,alp,bet+1,xx);  nx=size(x); 
if nx(2)>nx(1), dy=dy'; y=y'; xx=x'; end;  %% dy is a column vector

D=zeros(n);
%% Column 1 
D(:,1)=[-0.5*(n-1)*(n+alp+bet+1)/(bet+2); ...  
    (-1)^(n-1)*exp(gammaln(n)+gammaln(bet+2)-gammaln(n+bet+1)).*dy(2:end)];  
%% Column 2:n
 % row 1
 D(1,2:end)=(-1)^n*exp(gammaln(n+bet+1)-gammaln(n)-gammaln(bet+2))./...
    ((1+xx(2:end)).^2.*dy(2:end));
 % row 2:n-1
 Da=(xx(2:end)./((1+xx(2:end)).*dy(2:end)))*((1+xx(2:end)).*dy(2:end))'...
     -(1./((1+xx(2:end)).*dy(2:end)))*(xx(2:end).*(1+xx(2:end)).*dy(2:end))'; 
 Da=Da+eye(n-1);Da=1./Da; 
 Da=Da-diag(diag(Da));
 D(2:end,2:end)=Da+diag(0.5*((alp-bet+1)+(alp+bet+1)*xx(2:end))./(1-xx(2:end).^2));
    
return;