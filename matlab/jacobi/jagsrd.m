% x=jagsrd(n,alp,bet) computes n nodes of the Jacobi-Gauss-Radau quadrature with parameter (alp,bet)
% by using the eigen-method 
% [x,w]= jagsrd(n,alp,bet) also returns the weights (stored in w)
% Use the function japoly() and jags()
% Last modified on September 4, 2011


function [varargout]=jagsrd(n,alp,bet)

if n<1, disp('Input n >=1'); varargout{1}='Wrong input'; return; end;
if n==1, 
 varargout{1}=-1; 
 varargout{2}=exp((alp+bet+1)*log(2)+gammaln(bet+1)+ gammaln(alp+1)-gammaln(alp+bet+2)); 
 return;
end

x=jags(n-1,alp,bet+1);
varargout{1}=[-1;x]; 

if nargout==1, return; end;

  gn=(alp+bet+2)*log(2)+gammaln(n+alp)+gammaln(n+bet+1)-gammaln(n)-gammaln(n+alp+bet+1);
  gn=exp(gn);                           % Constant in the weight expression
  [dy,y]=japoly(n-1,alp,bet+1,x);        % Compute derivative of Jacobi polynomial of degree n-2 
                                        % at nodes 
  w=gn./((1+x).*(1-x.^2).*dy.^2);       % Compute the interior  weights 
  
  w0=(alp+bet+1)*log(2)+2*gammaln(bet+1)+gammaln(n)+gammaln(n+alp)...
    -gammaln(n+bet+1)-gammaln(n+alp+bet+1);
  w0=(bet+1)*exp(w0);                            % Weight corresponding to x=-1; 
  
  varargout{2} =[w0;w];
 return;

