% D=jagslbdiff(n,x,alp,bet) computes the first-order differentiation matrix of size
% n by n, associated with the Jacobi-Gauss-Lobatto points x, which can be computed by 
% x=jagslb(n,alp,bet) 
% See Pages 89-90 of the book: J. Shen, T. Tang and L. Wang, Spectral Methods:
%  Algorithms, Analysis and Applications, Springer Series in Compuational
%  Mathematics, 41, Springer, 2011. 
%  Use the function: japoly() 
% Last modified on September 4, 2011

function  D=jagslbdiff(n,x,alp,bet)

if n==0, D=[]; return; end;
xx=x; [dy,y]=japoly(n-2,alp+1,bet+1,xx);  nx=size(x); 
if nx(2)>nx(1), dy=dy'; y=y'; xx=x'; end;  %% dy is a column vector

D=zeros(n);
%% Column 1 
D(:,1)=[0.5*(alp-(n-1)*(n+alp+bet))/(bet+2) ; ...  
    (-1)^n*0.5*exp(gammaln(n-1)+gammaln(bet+2)-gammaln(n+bet))*(1-xx(2:n-1)).*dy(2:n-1) ; ...
    (-1)^(n-1)*0.5*exp(gammaln(n+alp)+gammaln(bet+2)-gammaln(alp+2)-gammaln(n+bet))];  
%% Column 2:n-1
 % row 1
 D(1,2:n-1)=2*(-1)^(n-1)*exp(gammaln(n+bet)-gammaln(n-1)-gammaln(bet+2))./...
    ((1-xx(2:n-1)).*(1+xx(2:n-1)).^2.*dy(2:n-1));
 % row 2:n-1
 Da=(xx(2:n-1)./((1-xx(2:n-1).^2).*dy(2:n-1)))*((1-xx(2:n-1).^2).*dy(2:n-1))'...
     -(1./((1-xx(2:n-1).^2).*dy(2:n-1)))*(xx(2:n-1).*(1-xx(2:n-1).^2).*dy(2:n-1))'; 
 Da=Da+eye(n-2);Da=1./Da; 
 Da=Da-diag(diag(Da));
 D(2:n-1,2:n-1)=Da+diag(0.5*((alp-bet)+(alp+bet)*xx(2:n-1))./(1-xx(2:n-1).^2));
 % row n
 D(n,2:n-1)=-2*(exp(gammaln(n+alp)-gammaln(n-1)-gammaln(alp+2)))./...
     ((1+xx(2:n-1)).*(1-xx(2:n-1)).^2.*dy(2:n-1));
%% Column n 
D(:,n)=[0.5*(-1)^n*exp(gammaln(alp+2)+gammaln(n+bet)-gammaln(bet+2)-gammaln(n+alp));...
    0.5*exp(gammaln(n-1)+gammaln(alp+2)-gammaln(n+alp))*(1+xx(2:n-1)).*dy(2:n-1);
    0.5*((n-1)*(n+alp+bet)-bet)/(alp+2)];
   
return;