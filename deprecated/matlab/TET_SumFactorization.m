clc;
format short

%% 2D
n = 3;
% NOTE: SPARSITY AND VALUES OF FN AND FN*C*FN^T ARE MODIFIED BASED ON THE 
%       POSITION OF THE VERTICES!!!
% In 2D, FN2*C*FN2' is always the same for the cases in the loop below, but
% the sparsity of FN2 varies with varying operation counts in some cases.
% Thus, if the equilateral triangle vertices are used to define FN, then
% the choice of their position is no longer free. However, in the 2d case,
% FN can always be determined using theta (as below) and need not rely on
% the vertex positions. Note that this would be inconsistent with using
% vertex nodes for the specification of FN in 3D TETs.
tmp = rand([1 n]);
denom = 8;
for num = 0:denom

theta = 2*pi*(0:n-1)/n+2*pi/n*(num/denom);

f = @(theta) 1/n*[exp(1i*0*theta) ...
                  exp(1i*1*theta) ...
                  exp(1i*2*theta)];
       
% FN = dftmtx(n)
FN =  eye(n)/f(theta');

f2 = @(theta) 1/n*[0*theta+1 ...
                   sqrt(2)*cos(theta) ...
                   sqrt(2)*sin(theta)];

FN2 = eye(n)/f2(theta');


C = circulant(tmp,1);

% 3 unique values
num
% FN*C*FN'
FN2
% FN2*C*FN2'
% norm(FN*FN'-FN2*FN2')
end
% break;

%% 3D
% https://en.wikipedia.org/wiki/Spherical_harmonics
% https://en.wikipedia.org/wiki/Table_of_spherical_harmonics#Real_spherical_harmonics
clc;

n = 4;
% NOTE: SPARSITY AND VALUES OF FN AND FN*C*FN^T ARE MODIFIED BASED ON THE 
%       POSITION OF THE VERTICES!!!
% XYZ = [-1 0 -1/sqrt(2); 1 0 -1/sqrt(2); 0 -1 1/sqrt(2); 0 1 1/sqrt(2)];
XYZ = sqrt(2)*[-1 -1 -1; 1 1 -1; -1 1 1; 1 -1 1];
% XYZ = 1/2*[[-1 0 -1/sqrt(2); 1 0 -1/sqrt(2); 0 -1 1/sqrt(2); 0 1 1/sqrt(2)] + ...
%            sqrt(2)*[-1 -1 -1; 1 1 -1; -1 1 1; 1 -1 1]];
% XYZ = XYZ([1 2 3 4],:);
XYZ = XYZ([1 2 4 3],:);
% XYZ = XYZ([1 3 2 4],:);
% XYZ = XYZ([1 3 4 2],:);
% XYZ = XYZ([1 4 2 3],:);
% XYZ = XYZ([1 4 3 2],:);
% XYZ = XYZ([2 1 3 4],:);
% XYZ = XYZ([2 1 4 3],:);
% XYZ = XYZ([2 3 1 4],:);
% XYZ = XYZ([2 3 4 1],:);
% XYZ = XYZ([2 4 1 3],:);
% XYZ = XYZ([2 4 3 1],:);

alpha = pi*1/6;
beta = pi*1/6;
gamma = pi*1/6;
Rotation1 = [1 0 0; 0 cos(alpha) -sin(alpha); 0 sin(alpha) cos(alpha)];
Rotation2 = [cos(beta) 0 sin(beta); 0 1 0; -sin(beta) 0 cos(beta)];
Rotation3 = [cos(gamma) sin(gamma) 0; -sin(gamma) cos(gamma) 0; 0 0 1];
XYZ
% XYZ = XYZ*Rotation1*Rotation2*Rotation3

% break;

               
   
x = XYZ(:,1);
y = XYZ(:,2);
z = XYZ(:,3);
   
r = sqrt(sum(XYZ.^2,2));
theta = atan2(y,x);
phi = acos(z./r);

f = @(theta,phi) 1/2*sqrt(pi)*[0*theta + 1/2*sqrt(1/pi) ...
                  1/2*sqrt(3/(2*pi))*exp(-1i*theta).*sin(phi) ...
                  1/2*sqrt(3/pi)*cos(phi) ...
                  -1/2*sqrt(3/(2*pi))*exp(1i*theta).*sin(phi)];
              
FN =  eye(n)/f(theta,phi);

f2 = @(x,y,z,r) 1/2*sqrt(pi)*[0*theta + 1/2*sqrt(1/pi) ...
                 sqrt(3/(4*pi))*y./r ...
                 sqrt(3/(4*pi))*z./r ...
                 sqrt(3/(4*pi))*x./r];
             
FN2 =  eye(n)/f2(x,y,z,r);

C = circulant(rand([1 n]),1);

% 5 unique values
% 4 unique values using second XYZ above.
FN*C*FN'
FN2*C*FN2'

norm(FN*FN'-FN2*FN2')

FN2
round(FN*FN')

% break;

%% 3D testing
clc;
P = 2;

[xir_vGs,~,NvnGs] = CubatureTET(GLOBAL,1,'alpha-opt');
[rst,~,Nvn]       = CubatureTET(GLOBAL,P,'alpha-opt');

% Barycentric coordinates
% L_vGc = [rst ones(Nvn,1)]/[xir_vGs ones(d+1,1)]; %Wikipedia
% Note that the P1 basis is [1 x y z] =>
L_vGc = BasisTET(1,rst)/BasisTET(1,xir_vGs);

% Convert to vertices not collocated with the poles of the sphere (i.e. phi
% ~= 0 or pi)
xir_vGs = XYZ;
rst = L_vGc*XYZ;

if (P~= 2)
    Error only working for P2 atm
end

n = 4;
x = rst([1 3 6 10],1);
y = rst([1 3 6 10],2);
z = rst([1 3 6 10],3);

xyz = [x y z];

r = sqrt(sum(xyz.^2,2));
theta = atan2(y,x);
phi = acos(z./r);

f = @(theta,phi) ...
    1/2*sqrt(pi)*[0*theta + 1/2*sqrt(1/pi) ...
                  1/2*sqrt(3/(2*pi))*exp(-1i*theta).*sin(phi) ...
                  1/2*sqrt(3/pi)*cos(phi) ...
                  -1/2*sqrt(3/(2*pi))*exp(1i*theta).*sin(phi)];
              
FN =  eye(n)/f(theta,phi);
% F43 = eye(n)/f(theta,phi);

f2 = @(x,y,z,r) ...
    1/2*sqrt(pi)*[0*theta + 1/2*sqrt(1/pi) ...
                  sqrt(3/(4*pi))*y./r ...
                  sqrt(3/(4*pi))*z./r ...
                  sqrt(3/(4*pi))*x./r];
             
F43 =  eye(n)/f2(x,y,z,r);

C = circulant(rand([1 n]),1);

FN*C*FN'
F43*C*F43'

norm(FN*FN'-F43*F43')

F43
round(FN*FN')
% break;


%% n = 6 (4-rotational, 2-reflectional)

n = 6; n_r = 4; n_ref = 2;
x = rst([2 4 5 7 8 9],1);
y = rst([2 4 5 7 8 9],2);
z = rst([2 4 5 7 8 9],3);

xyz = [x y z];

Ind_r = [3 5 4 2];
x_r = xyz(Ind_r,1);
y_r = xyz(Ind_r,2);
xyz_r = [x_r y_r];

theta_r = atan2(y_r,x_r);

% FN = dftmtx(n)
f42 = @(theta) 1/n_r*[0*theta+1 ...
                     sqrt(2)*cos(theta) ...
                     sqrt(2)*sin(theta) ...
                     sqrt(1)*cos(2*theta)];

F42 = eye(n_r)/f42(theta_r)
% F42 = dftmtx(n_r);
F22 = dftmtx(n_ref);


C = circulant(rand([1 4]),1);

% 4 unique values
F42*C*F42'
% FN2*FN2'


clc

rst = rst([1 3 6 10 4 5 7 8 2 9],:);
L_vGc = BasisTET(1,rst)/BasisTET(1,xir_vGs)
x = rst(:,1);
y = rst(:,2);
z = rst(:,3);
xyz = [x y z];

r = sqrt(sum(xyz.^2,2));
theta = atan2(y,x);
phi = acos(z./r);

I       = eye(Nvn);
Chi_rst = BasisTET(P,rst);

T = I/Chi_rst;

GradChi_rst  = GradBasisTET(P,rst);
for dim = 1:d
    GradChi_rst(:,:,dim) = GradChi_rst(:,:,dim)*T;
end

Dxi  = GradChi_rst(:,:,1);
Deta = GradChi_rst(:,:,2);
Dzeta = GradChi_rst(:,:,3);

xi_r = cos(theta).*sin(phi);
xi_t = -r.*sin(theta).*sin(phi);
xi_p = r.*cos(theta).*cos(phi);
eta_r = sin(theta).*sin(phi);
eta_t = r.*cos(theta).*sin(phi);
eta_p = r.*sin(theta).*cos(phi);
zeta_r = cos(phi);
zeta_t = 0*r;
zeta_p = -r.*sin(phi);

Dr = diag(xi_r)*Dxi + diag(eta_r)*Deta + diag(zeta_r)*Dzeta;
Dt = diag(xi_t)*Dxi + diag(eta_t)*Deta + diag(zeta_t)*Dzeta;
Dp = diag(xi_p)*Dxi + diag(eta_p)*Deta + diag(zeta_p)*Dzeta;

OP = Dp;
% C = OP(5:8,5:8)
C = OP(1:4,1:4)
F42*C*F42'

F = [F43 zeros(4,6); zeros(4,4) F42 zeros(4,2); zeros(2,8) F22];
% F = [F42 zeros(4,6); zeros(4,4) F42 zeros(4,2); zeros(2,8) F22];

OP(abs(OP) < 1e-14) = 0; OP
tmp = F*OP*F'; tmp(abs(tmp) < 1e-14) = 0; tmp


clc

OP = [circulant(rand([1 4]),1) circulant(rand([1 4]),1) ...
      [circulant(rand([1 2]),1); circulant(rand([1 2]),1)]; ...
      circulant(rand([1 4]),1) circulant(rand([1 4]),1) ...
      [circulant(rand([1 2]),1); circulant(rand([1 2]),1)]; ...
      [circulant(rand([1 2]),1) circulant(rand([1 2]),1)] ...
      [circulant(rand([1 2]),1) circulant(rand([1 2]),1)] ...
      circulant(rand([1 2]),1)]
  
F = [F42 zeros(4,6); zeros(4,4) F42 zeros(4,2); zeros(2,8) F22];

tmp = F*OP*F'; tmp(abs(tmp) < 1e-14) = 0; tmp

range = 1:8;
tmp = F(range,range)*OP(range,range)*F(range,range)'; tmp(abs(tmp) < 1e-14) = 0; tmp
    


%% n = 6 (Failure: The 6-symmetry is not spherical (think of the octohedron
%%                 embedded in the sphere)

% 
% clc;
% 
% n = 6;
% x = rst([2 4 5 7 8 9],1);
% y = rst([2 4 5 7 8 9],2);
% z = rst([2 4 5 7 8 9],3);
% 
% xyz = [x y z];
% 
% r = sqrt(sum(xyz.^2,2));
% theta = atan2(y,x);
% phi = acos(z./r);
% 
% f = @(theta,phi) ...
%     [0*theta + 1/2*sqrt(1/pi) ...
%                   1/4*sqrt(15/(2*pi))*exp(-1i*2*theta).*sin(phi).^2 ...
%                   1/2*sqrt(15/(2*pi))*exp(-1i*theta).*sin(phi).*cos(phi) ...
%                   1/4*sqrt(5/pi)*(3*cos(phi).^2-1) ...
%                   -1/2*sqrt(15/(2*pi))*exp(1i*theta).*sin(phi).*cos(phi) ...
%                   1/4*sqrt(15/(2*pi))*exp(1i*2*theta).*sin(phi).^2];
%               
% FN =  eye(n)/f(theta,phi);
% 
% % f2 = @(x,y,z,r) ...
% %     [0*theta + 1/2*sqrt(1/pi) ...
% %                   1/2*sqrt(15/pi)*x.*y./r.^2 ...
% %                   1/2*sqrt(15/pi)*y.*z./r.^2 ...
% %                   1/4*sqrt(5/pi)*(-x.^2-y.^2+2*z.^2)./r.^2 ...
% %                   1/2*sqrt(15/pi)*z.*x./r.^2 ...
% %                   1/4*sqrt(15/pi)*(x.^2-y.^2)./r.^2];
%    
% % Working: Y[0,0], Y[1,-1], Y[1,0], Y[1,1], Y[2,0], Y[2,2]
% % Y[2,-1]: (Scale)*sqrt(1/72)*y.*z./r.^2
% % Y[2,1]: (Scale)*sqrt(1/72)*z.*x./r.^2
% % Y[2,2]: sqrt(3/72)*(x.^2-y.^2)./r.^2
% % Y[3,0]: (Scale)*z.*(2*z.^2-3*(x.^2-y.^2))./r.^3
% % Y[4,0]: (Scale)*(35*z.^4-30*(z.*r).^2+3*r.^4)./r.^4
% 
% % Y[0,0], Y[1,-1], Y[1,0], Y[1,1], Y[2,-1], Y[2,1]
% f2 = @(x,y,z,r) [1/6 + 0*r ...
%                  sqrt(1/12)*y./r ...
%                  sqrt(1/12)*z./r ...
%                  sqrt(1/12)*x./r ...
%                  sqrt(1/72)*(-(x.^2+y.^2)+2*z.^2)./r.^2 ...
%                  sqrt(3/72)*(x.^2-y.^2)./r.^2];
%              
% FN2 =  eye(n)/f2(x,y,z,r);
% 
% C = circulant(rand([1 n]),1);
% 
% % FN*C*FN'
% FN2*C*FN2'
% 
% % norm(FN*FN'-FN2*FN2')
% 
% FN2
% % format long
% FN2*FN2'
% 
% % xyz
% % f2(x,y,z,r)
% % rank(f2(x,y,z,r))

