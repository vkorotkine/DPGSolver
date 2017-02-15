%y=gebauerpoly(n,lam,x) computes the Gegenbauer polynomial of degree n with parameter lam>-1/2 at a
%   vector-valued x
% [dy,y]=gebauerpoly(n, lam, x) also returns the values of the 1st-order 
%   derivative stored in dy
% See Page 81 of the book: G. Szego, Orthogonal Polynomials, Volumn 23,
% AMS, 1975.
%  Last modified on September 4, 2011    


function [varargout]=gebauerpoly(n,lam,x)
     
if nargout==1,
     if n==0, varargout{1}=ones(size(x));  return; end;
     if n==1, varargout{1}=2*lam*x; return; end;

     polylst=ones(size(x));	
     poly=2*lam*x;
     
   for k=2:n,
	  polyn=(2*(k+lam-1)*x.*poly-(k+2*lam-2)*polylst)/k;  
      polylst=poly; poly=polyn;	
   end;
      varargout{1}=polyn; return;
end;

   
if n==0, varargout{2}=ones(size(x)); vararout{1}=zeros(size(x)); return; end;
if n==1, varargout{2}=2*lam*x; varargout{1}=2*lam; return; end;

 polylst=ones(size(x)); pderlst=zeros(size(x)); poly=2*lam*x; pder=2*lam;
   
   for k=2:n,
	  polyn=(2*(k+lam-1)*x.*poly-(k+2*lam-2)*polylst)/k;
      pdern=(2*(k+lam-1)*x.*pder-(k+2*lam-2)*pderlst+2*(k+lam-1)*poly)/k;
      polylst=poly; poly=polyn;
	  pderlst=pder; pder=pdern;
   end;
      varargout{2}=polyn;
      varargout{1}=pdern;
   return;


    

  
     
	
