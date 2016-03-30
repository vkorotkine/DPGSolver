%% Notes

% Note: Since there will be a FT followed by an inverse FT, it would be
%       beneficial to combine the two 1/sqrt(n) terms (to avoid the sqrt),
%       as in standard fft algorithms. Further, as the the 2x2 FN matrix is
%       filled with ones multiplied by a scaling constant, all of the
%       multiplication floating point operations can be eliminated in the
%       forward and backwards Fourier transforms (ensuring that the correct
%       addition/subtraction is used instead). Obviously, no complex
%       algebra is required in this case either, but this is no longer true
%       for n > 2.

%% TP

% Try with modal basis later using symmetry and anti-symmetry of the basis
% functions.
clc;

P = 7;
Pt = 7;
d = 1;

[rst,~,N] = CubatureTensorProd(GLOBAL,P,d,'GLL');

reorder = zeros(N,1);
reorder(1:2:N) = 1:ceil(N/2);
reorder(2:2:N) = N:-1:ceil(N/2)+1;
rst = rst(reorder,:)

I       = eye(N);
Chi_rst = BasisTensorProd(P,rst);

T = I/Chi_rst; %nodal
% T = I; %modal

[rst,~,Nt] = CubatureTensorProd(GLOBAL,Pt,d,'GL');

reorder_t = zeros(Nt,1);
reorder_t(1:2:Nt) = 1:ceil(Nt/2);
reorder_t(2:2:Nt) = Nt:-1:ceil(Nt/2)+1;

rst = rst(reorder_t,:)

Chi_test = BasisTensorProd(P,rst)*T;
% Chi_test = GradBasisTensorProd(P,rst)*T;

% Fourier matrix (Note the - sign)
n = 2;
FN = dftmtx(n);
% FN = zeros(n,n);
% for i = 1:n
% for j = 1:n
%     FN(i,j) = exp(-2*pi*1i*(i-1)*(j-1)/n);
% end
% end
% if (n == 2); FN = real(FN); end
% FN = 1/sqrt(n)*FN; Fscale = 1;
Fscale = 1/n;

r = 1;
theta = [pi 0];

xi_r = ones(1,Nt);
for i = 0:floor(Nt/2)-1
    xi_r(i*2+(1:2)) = cos(theta);
end

% Dr = diag(xi_r)*Chi_test %xi_r = +-1 cancels with r_x +-1
Dr = Chi_test

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

Dtmp = Fscale*Ft*Dr*F';
Dtmp(abs(Dtmp) < 1e-14) = 0;
Drtilde = Dtmp

% Application testing
Vec = rand([N 1]);

D_standard = Dr*Vec;
D_new1     = Fscale*F'*Drtilde*F*Vec;

[D_standard D_new1]
norm(D_standard-D_new1)

% break;

%% TRI
clc;

P = 2;
d = 2;

[rst,~,Nvn]      = CubatureTRI(GLOBAL,P,'alpha-opt');
if     (P == 2); rst = rst([2 5 4 1 3 6],:);
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
FN = dftmtx(N);
% FN = zeros(N,N);
% for i = 1:N
% for j = 1:N
%     FN(i,j) = exp(-2*pi*1i*(i-1)*(j-1)/N);
% end
% end
% FN = 1/sqrt(N)*FN;
Fscale = 1/N;

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
Dtheta = diag(xi_theta)*Dxi + diag(eta_theta)*Deta
% Dtheta(abs(Dtheta) < 1e-14) = 0
% Can see that it is block-circulant

if     (P == 2); F  = [FN FN*0; FN*0 FN];
elseif (P == 3); 
    F  = [FN zeros(3,7); ...
          zeros(3,3) FN zeros(3,4); ...
          zeros(3,6) FN zeros(3,1); ...
          zeros(1,9) 1];
end

Dtmp = Fscale*F*Dr*F';
Dreal = real(Dtmp); Dreal(abs(Dreal) < 1e-14) = 0;
Dimag = imag(Dtmp); Dimag(abs(Dimag) < 1e-14) = 0;
Dreal;
Dimag;
Drtilde = Dreal+1i*Dimag;

Dtmp = Fscale*F*Dtheta*F';
Dreal = real(Dtmp); Dreal(abs(Dreal) < 1e-14) = 0;
Dimag = imag(Dtmp); Dimag(abs(Dimag) < 1e-14) = 0;
Dreal;
Dimag;
Dthetatilde = Dreal+1i*Dimag

%% TRI 2

[rst,~,Nvn] = CubatureTRI(GLOBAL,P,'alpha-opt');
if     (P == 2); rst = rst([2 1 5 3 4 6],:);
elseif (P == 3); rst = rst([6 5 2 1 3 7 4 9 8 10],:);
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
FN = dftmtx(N);
if (P == 2)
    Fscale = ones(1,6)*1/N;
elseif (P == 3)
    Fscale = [1 ones(1,9)*1/N];
end
    
r = 1;
theta = [0 1/3*2*pi 2/3*2*pi];

if (P == 2)
    xi_r = [reshape(ones(2,1)*cos(theta),1,6)];
    xi_theta = -r*[reshape(ones(2,1)*sin(theta),1,6)];
    eta_r = [reshape(ones(2,1)*sin(theta),1,6)];
    eta_theta = r*[reshape(ones(2,1)*cos(theta),1,6)];
    
    reorder = [1 3 5 2 4 6];
elseif (P == 3)
    xi_r = [1 reshape(ones(3,1)*cos(theta),1,9)];
    xi_theta = -r*[0 reshape(ones(3,1)*sin(theta),1,9)];
    eta_r = [0 reshape(ones(3,1)*sin(theta),1,9)];
    eta_theta = r*[1 reshape(ones(3,1)*cos(theta),1,9)];
    
    reorder = [1 2 5 8 3 6 9 4 7 10];
end

Dr     = diag(xi_r)*Dxi     + diag(eta_r)*Deta;
Dtheta = diag(xi_theta)*Dxi + diag(eta_theta)*Deta

perm = I(reorder,:);

if     (P == 2); F  = [FN FN*0; FN*0 FN];
elseif (P == 3); 
    F  = [1 zeros(1,9); ...
          zeros(3,1) FN zeros(3,6); ...
          zeros(3,4) FN zeros(3,3); ...
          zeros(3,7) FN];
end

F = perm'*F*perm;

Dtmp = diag(Fscale)*F*Dr*F';
Dreal = real(Dtmp); Dreal(abs(Dreal) < 1e-14) = 0;
Dimag = imag(Dtmp); Dimag(abs(Dimag) < 1e-14) = 0;
Dreal;
Dimag;
Drtilde = Dreal+1i*Dimag;

Dtmp = diag(Fscale)*F*Dtheta*F';
Dreal = real(Dtmp); Dreal(abs(Dreal) < 1e-14) = 0;
Dimag = imag(Dtmp); Dimag(abs(Dimag) < 1e-14) = 0;
Dreal
Dimag
Dthetatilde = Dreal+1i*Dimag;

NCols = 1;
RHS_In = rand([Nvn NCols]);

O_std = Dtheta*RHS_In;
O_new2_slow = diag(Fscale)*F'*Dthetatilde*F*RHS_In;

if (P == 2)
    FT_In = [sum(RHS_In([1 3 5],:)); ...
             sum(RHS_In([2 4 6],:)); ...
             RHS_In(1,:)-1/2*sum(RHS_In([3 5],:)); ...
             RHS_In(2,:)-1/2*sum(RHS_In([4 6],:)); ...
             sqrt(3)/2*(RHS_In(3,:)-RHS_In(5,:)); ...
             sqrt(3)/2*(RHS_In(4,:)-RHS_In(6,:))];
         
    C1 = Dthetatilde(1:2,1:2);
    C2 = real(Dthetatilde(3:4,3:4));
    C3 = imag(Dthetatilde(5:6,5:6));
    
    FT_Out = [C1*FT_In(1:2,:); ...
              C2*FT_In(3:4,:)-C3*FT_In(5:6,:); ...
              C3*FT_In(3:4,:)+C2*FT_In(5:6,:)];
          
%     Dthetatilde*F*RHS_In
    O_new2 = diag(Fscale)*...
        [FT_Out([1 2],:)+2*FT_Out([3 4],:); ...
         FT_Out([1 2],:)-FT_Out([3 4],:)+sqrt(3)*FT_Out([5 6],:); ...
         FT_Out([1 2],:)-FT_Out([3 4],:)-sqrt(3)*FT_Out([5 6],:)]
           
    diag(Fscale)*F'*Dthetatilde*F*RHS_In
elseif (P == 3)
    FT_In = [RHS_In(1,:); ...
             sum(RHS_In([2 5 8] ,:)); ...
             sum(RHS_In([3 6 9] ,:)); ...
             sum(RHS_In([4 7 10],:)); ...
             RHS_In(2,:)-1/2*sum(RHS_In([5 8],:)); ...
             RHS_In(3,:)-1/2*sum(RHS_In([6 9],:)); ...
             RHS_In(4,:)-1/2*sum(RHS_In([7 10],:)); ...
             sqrt(3)/2*(RHS_In(5,:)-RHS_In(8,:)); ...
             sqrt(3)/2*(RHS_In(6,:)-RHS_In(9,:)); ...
             sqrt(3)/2*(RHS_In(7,:)-RHS_In(10,:))];
end

Dthetatilde

norm([O_std-O_new2_slow],'inf')

break;

%% Further testing
clc

for Nn = 1:30
% Nn = 12;
tmp = fft(rand([1 Nn]));

tmpR = real(tmp); I = find(tmpR == 0); tmpR(I) = [];
tmpI = imag(tmp); I = find(tmpI == 0); tmpI(I) = [];
tmpU = unique(abs([tmpR tmpI]));

size(tmpU)
end

clc;

n = 3;
FN = dftmtx(n);

Circ = circulant(rand([1 n]),1)

Circ_FT = FN*Circ*FN'
Circ_FT(abs(Circ_FT) < 1e-13) = 0

clc;

syms d1 d2 d3 d11 d21 d31;
j = 1i;

d11 = d1;
d21 = d2-d3*j;
d31 = d2+d3*j;

tmp = simplify(FN*[d11 d21 d31]') %should be fully real
conj(tmp)












