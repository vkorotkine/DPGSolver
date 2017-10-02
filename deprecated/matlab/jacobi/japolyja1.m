%y=japolyja1(n,alp,bet,x) computes the normalized Jacobi polynomial of degree n with parameter (alp,bet) at a
%  vector-valued x. 
%[dy,y]=japolyja1(n,alp, bet, x) also returns the values of the 1st-order 
%  derivative stored in dy
%  Note: here, the normalization is 
%   if alp>=bet, japolyja1(n,alp,bet,x)=japoly(n,alp,bet,x)/japoly(n,alp,bet,1); 
%   if alp<beta, japolyja1(n,alp,bet,x)=japoly(n,alp,bet,x)/japoly(n,alp,bet,-1);
% With this normalization, we have japolyja1(n,alp,bet,x)<=1, if alp,bet>=-1/2.
% See Page 78 of the book: J. Shen, T. Tang and L. Wang, Spectral Methods:
%  Algorithms, Analysis and Applications, Springer Series in Compuational
%  Mathematics, 41, Springer, 2011. 
% 
% Last modified on September 3, 2011    


function [varargout]=japolyja1(n,alp,bet,x)
  apb=alp+bet;     
  
% Case 1: alp>=bet
if alp>=bet, 
 if nargout==1,
     if n==0, varargout{1}=ones(size(x));  return; end;
     if n==1, varargout{1}=0.5*(alp-bet+(apb+2)*x)/(alp+1); return; end;

     polylst=ones(size(x));	
     poly=0.5*(alp-bet+(apb+2)*x)/(alp+1);
     
  for k=2:n,
	 a1=2.0*k*(k+apb)*(2.0*k+apb-2.);
	 a2=(2.0*k+apb-1.)*(alp^2-bet^2);
	 b3=(2.0*k+apb-2.); a3=b3*(b3+1.)*(b3+2.);
	 a4=2.0*(k+alp-1.)*(k+bet-1.)*(2.0*k+apb);
	 polyn=((a2+a3*x).*poly*(k/(k+alp))-a4*polylst*((k-1)/(k+alp-1))*(k/(k+alp)))/a1;  % See (3.110) and (3.111) 
     polylst=poly; poly=polyn;	
  end;
    varargout{1}=polyn; return;
 end;

   
 if n==0, varargout{2}=ones(size(x)); varargout{1}=zeros(size(x)); return; end;
 if n==1, varargout{2}=0.5*(alp-bet+(apb+2)*x)/(alp+1); varargout{1}=0.5*(apb+2)*ones(size(x))/(alp+1); return; end;

  polylst=ones(size(x));                  pderlst=zeros(size(x));
  poly=0.5*(alp-bet+(apb+2.)*x)/(alp+1);  pder=0.5*(apb+2.)*ones(size(x))/(alp+1);
   
  for k=2:n,
	 a1=2.0*k*(k+apb)*(2.0*k+apb-2.);
	 a2=(2.0*k+apb-1.)*(alp^2-bet^2);
	 b3=(2.0*k+apb-2.); a3=b3*(b3+1.)*(b3+2.);
	 a4=2.0*(k+alp-1.)*(k+bet-1.)*(2.0*k+apb);
	 polyn=((a2+a3*x).*poly*(k/(k+alp))-a4*polylst*((k-1)/(k+alp-1))*(k/(k+alp)))/a1;  % See (3.110) and (3.111)
	 pdern=((a2+a3*x).*pder*(k/(k+alp))-a4*pderlst*((k-1)/(k+alp-1))*(k/(k+alp))+a3*poly*(k/(k+alp)))/a1;	  	  
	 polylst=poly; poly=polyn;
	 pderlst=pder; pder=pdern;
  end;
      varargout{2}=polyn;
      varargout{1}=pdern;
  return;
end

% Case 2: alp<bet normalized by dividing J_n^{alp,bet}(-1)
if alp<bet, 
 if nargout==1,
     if n==0, varargout{1}=ones(size(x));  return; end;
     if n==1, varargout{1}=-0.5*(alp-bet+(apb+2)*x)/(bet+1); return; end;

     polylst=ones(size(x));	
     poly=-0.5*(alp-bet+(apb+2)*x)/(bet+1);
     
  for k=2:n,
	 a1=2.0*k*(k+apb)*(2.0*k+apb-2.);
	 a2=(2.0*k+apb-1.)*(alp^2-bet^2);
	 b3=(2.0*k+apb-2.); a3=b3*(b3+1.)*(b3+2.);
	 a4=2.0*(k+alp-1.)*(k+bet-1.)*(2.0*k+apb);
	 polyn=((a2+a3*x).*poly*(-k/(k+bet))-a4*polylst*((k-1)/(k+bet-1))*(k/(k+bet)))/a1;  % See (3.110) and (3.111) 
     polylst=poly; poly=polyn;	
  end;
    varargout{1}=polyn; return;
 end;

   
 if n==0, varargout{2}=ones(size(x)); varargout{1}=zeros(size(x)); return; end;
 if n==1, varargout{2}=-0.5*(alp-bet+(apb+2)*x)/(bet+1); varargout{1}=-0.5*(apb+2)*ones(size(x))/(bet+1); return; end;

  polylst=ones(size(x));                  pderlst=zeros(size(x));
  poly=-0.5*(alp-bet+(apb+2.)*x)/(bet+1);  pder=-0.5*(apb+2.)*ones(size(x))/(bet+1);
   
  for k=2:n,
	 a1=2.0*k*(k+apb)*(2.0*k+apb-2.);
	 a2=(2.0*k+apb-1.)*(alp^2-bet^2);
	 b3=(2.0*k+apb-2.); a3=b3*(b3+1.)*(b3+2.);
	 a4=2.0*(k+alp-1.)*(k+bet-1.)*(2.0*k+apb);
	 polyn=((a2+a3*x).*poly*(-k/(k+bet))-a4*polylst*((k-1)/(k+bet-1))*(k/(k+bet)))/a1;  % See (3.110) and (3.111)
	 pdern=((a2+a3*x).*pder*(-k/(k+bet))-a4*pderlst*((k-1)/(k+bet-1))*(k/(k+bet))+a3*poly*(-k/(k+bet)))/a1;	  	  
	 polylst=poly; poly=polyn;
	 pderlst=pder; pder=pdern;
  end;
     varargout{2}=polyn;
     varargout{1}=pdern;
  return;
end


    

  
     
	
