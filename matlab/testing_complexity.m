%% TP

% Try with modal basis later using symmetry and anti-symmetry of the basis
% functions.
clc;

P = 7;
Pt = 6;
d = 1;

[rst,~,N] = CubatureTensorProd(GLOBAL,P,d,'GLL');

reorder = zeros(N,1);
reorder(1:2:N) = 1:ceil(N/2);
reorder(2:2:N) = N:-1:ceil(N/2)+1;

I       = eye(N);
Chi_rst = BasisTensorProd(P,rst);

T = I/Chi_rst; %nodal
% T = I; %modal

[rst,~,Nt] = CubatureTensorProd(GLOBAL,Pt,d,'GL');
rst

reorder_t = zeros(Nt,1);
reorder_t(1:2:Nt) = 1:ceil(Nt/2);
reorder_t(2:2:Nt) = Nt:-1:ceil(Nt/2)+1;

rst = rst(reorder_t,:)

% Chi_test = BasisTensorProd(P,rst)*T;
Chi_test = GradBasisTensorProd(P,rst)*T;

Chi_test = Chi_test(:,reorder);

% Fourier matrix
n = 2;
FN = zeros(n,n);
for i = 1:n
for j = 1:n
    FN(i,j) = exp(2*pi*1i*(i-1)*(j-1)/n);
end
end
if (n == 2); FN = real(FN); end
FN = 1/sqrt(n)*FN;

r = 1;
theta = [pi 0];

xi_r = ones(1,Nt);
for i = 0:floor(Nt/2)-1
    xi_r(i*2+(1:2)) = cos(theta);
end

Dr = diag(xi_r)*Chi_test

F = zeros(N,N);
for i = 0:floor(N/2)-1
    F(i*2+(1:2),i*2+(1:2)) = FN;
end
if (mod(N,2) == 1); F(N,N) = 1; end

Ft = zeros(Nt,Nt);
for i = 0:floor(Nt/2)-1
    Ft(i*2+(1:2),i*2+(1:2)) = FN;
end
if (mod(Nt,2) == 1); Ft(Nt,Nt) = 1; end

Dtmp = Ft*Dr*F';
Dtmp(abs(Dtmp) < 1e-14) = 0;
Drtilde = Dtmp

% break;

%% TRI
clc;

P = 3;
d = 2;

[rst,~,Nvn]      = CubatureTRI(GLOBAL,P,'alpha-opt');
if     (P == 2); rst = rst([1 3 6 2 5 4],:);
elseif (P == 3); rst = rst([1 4 10 2 7 8 3 9 5 6],:);
end
rst

I       = eye(Nvn);
Chi_rst = BasisTRI(P,rst);

T = I/Chi_rst;

GradChi_rst  = GradBasisTRI(P,rst);
for dim = 1:d
    GradChi_rst(:,:,dim) = GradChi_rst(:,:,dim)*T;
end

Dxi  = GradChi_rst(:,:,1);
Deta = GradChi_rst(:,:,2);

% Fourier matrix
N = 3;
FN = zeros(N,N);
for i = 1:N
for j = 1:N
    FN(i,j) = exp(2*pi*1i*(i-1)*(j-1)/N);
end
end
FN = 1/sqrt(N)*FN;

r = 1;
theta = [0 1/3*2*pi 2/3*2*pi];

if (P == 2)
    xi_r = [cos(theta) cos(theta)];
    xi_theta = -r*[sin(theta) sin(theta)];
    eta_r = [sin(theta) sin(theta)];
    eta_theta = r*[cos(theta) cos(theta)];
elseif (P == 3)
    xi_r = [cos(theta) cos(theta) cos(theta) 1];
    xi_theta = -r*[sin(theta) sin(theta) sin(theta) 0];
    eta_r = [sin(theta) sin(theta) sin(theta) 0];
    eta_theta = r*[cos(theta) cos(theta) cos(theta) 1];
end

Dr     = diag(xi_r)*Dxi     + diag(eta_r)*Deta;
Dr(abs(Dr) < 1e-14) = 0;
Dtheta = diag(xi_theta)*Dxi + diag(eta_theta)*Deta;
Dtheta(abs(Dtheta) < 1e-14) = 0
% Can see that it is block-circulant

if     (P == 2); F  = [FN FN*0; FN*0 FN];
elseif (P == 3); 
    F  = [FN zeros(3,7); ...
          zeros(3,3) FN zeros(3,4); ...
          zeros(3,6) FN zeros(3,1); ...
          zeros(1,9) 1];
end

Dtmp = F*Dr*F';
Dreal = real(Dtmp); Dreal(abs(Dreal) < 1e-14) = 0;
Dimag = imag(Dtmp); Dimag(abs(Dimag) < 1e-14) = 0;
Dreal;
Dimag;
Drtilde = Dreal+1i*Dimag;

Dtmp = F*Dtheta*F';
Dreal = real(Dtmp); Dreal(abs(Dreal) < 1e-14) = 0;
Dimag = imag(Dtmp); Dimag(abs(Dimag) < 1e-14) = 0;
Dreal
Dimag
Dthetatilde = Dreal+1i*Dimag;

real(F'*Dreal*F)


