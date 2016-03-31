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
ToBeDeleted = 0;

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

P = 3;
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
Dt = diag(xi_theta)*Dxi + diag(eta_theta)*Deta
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

Dtmp = Fscale*F*Dt*F';
Dreal = real(Dtmp); Dreal(abs(Dreal) < 1e-14) = 0;
Dimag = imag(Dtmp); Dimag(abs(Dimag) < 1e-14) = 0;
Dreal;
Dimag;
Dthetatilde = Dreal+1i*Dimag

%% TRI 2
clc
% For the next level of implementation, show that this works for PIn
% different from POut (such as for interpolation to a higher order) if both
% node sets are symmetric, similar to what was done for TP above.
P = 4;

[xir_vGs,~,NvnGs]   = CubatureTRI(GLOBAL,1,'alpha-opt');
[rst,~,Nvn] = CubatureTRI(GLOBAL,P,'alpha-opt');

% Barycentric coordinates
% L_vGc = [rst ones(Nvn,1)]/[xir_vGs ones(d+1,1)]; %Wikipedia
% Note that the P1 basis is [1 x y] =>
L_vGc = BasisTRI(1,rst)/BasisTRI(1,xir_vGs);

if (P == 1)
    Symms = 3;
    Symms_count = 1;
    rst = rst([1 2 3],:);
    reorder = [1 2 3];
elseif     (P == 2)
    Symms = [3 3];
    Symms_count = 2;
%     rst = rst([2 1 5 3 4 6],:);
%     reorder = [1 3 5 2 4 6];
    rst = rst([4 1 2 3 5 6],:);
    reorder = [1 3 5 2 4 6];
elseif (P == 3)
    Symms = [1 3 3 3];
    Symms_count = [1 3];
    rst = rst([6 5 2 1 3 7 4 9 8 10],:);
    reorder = [1 2 5 8 3 6 9 4 7 10];
elseif (P == 4)
    % All 3-symmetries
    Symms = [3 3 3 3 3];
    Symms_count = 5;
    rst = rst([7 10 6 2 1 8 3 4 9 5 11 12 14 13 15],:);
    reorder = [1 6 11 2 7 12 3 8 13 4 9 14 5 10 15];
elseif (P == 5)
    % All 3-symmetries
    Symms = [3 3 3 3 3 3 3];
    Symms_count = 7;
    rst = rst([13 8 12 3 7 2 1 9 10 4 15 5 11 6 14 17 18 16 20 19 21],:);
    reorder = [[1 1+7 1+2*7]+0 [1 1+7 1+2*7]+1 [1 1+7 1+2*7]+2 ...
               [1 1+7 1+2*7]+3 [1 1+7 1+2*7]+4 [1 1+7 1+2*7]+5 ...
               [1 1+7 1+2*7]+6];
elseif (P == 6)
    % First 6-symmetry (based on table 3 in Hesthaven(1998)), but the nodes
    % are not periodic in 60 degree intervals, slightly off... What does a
    % 6 symmetry provide with regards to the fft, if anything?
    Symms = [1 3 3 3 3 3 3 3 6];
    Symms_count = [1 7 1];
    rst = rst([16 9 19 14 3 8 2 1 15 10 12 4 5 18 6 13 7 11 17 ...
               24 22 25 23 27 26 28 21 20],:);
%     rst = rst([16 19 4 22 14 5 25 8 6 27 1 7 28 2 13 26 ...
%                   3 18 23 9 12 24 15 10 11 17 21 20],:);
    reorder = [1 [1 1+7 1+2*7]+1 [1 1+7 1+2*7]+2 [1 1+7 1+2*7]+3 ...
                 [1 1+7 1+2*7]+4 [1 1+7 1+2*7]+5 [1 1+7 1+2*7]+6 ...
                 [1 1+7 1+2*7]+7 [23:28]];
end
rst

clf; hold on;
x = rst(:,1); y = rst(:,2);
scatter(x,y);
a = [1:Nvn]'; b = num2str(a); c = cellstr(b);
dx = 0.025; dy = 0.025;
text(x+dx,y+dy,c);

[rst_P2,~,~] = CubatureTRI(GLOBAL,2,'alpha-opt');
scatter(rst_P2(:,1),rst_P2(:,2),'rs');
plot([0 rst_P2(2,1)],[0 rst_P2(2,2)],'r');
plot([0 rst_P2(4,1)],[0 rst_P2(4,2)],'r');
plot([0 rst_P2(5,1)],[0 rst_P2(5,2)],'r');
% break;

I       = eye(Nvn);
Chi_rst = BasisTRI(P,rst);

T = I/Chi_rst;

GradChi_rst  = GradBasisTRI(P,rst);
for dim = 1:d
    GradChi_rst(:,:,dim) = GradChi_rst(:,:,dim)*T;
end

Dxi  = GradChi_rst(:,:,1);
Deta = GradChi_rst(:,:,2);

% Fourier matrices
FN = zeros(max(Symms),max(Symms),max(Symms));
for N = unique(Symms)
    FN(1:N,1:N,N) = dftmtx(N);
end

Fscale = zeros(1,Nvn);
IndS = 0;
for i = Symms
    for j = IndS+(1:i)
        Fscale(j) = 1/i;
    end
    IndS = IndS+i;
end

if (ToBeDeleted)
% r = 1;
% 
% theta = zeros(1,max(Symms),max(Symms));
% for N = unique(Symms)
%     theta(1,1:N,N) = 0:2*pi*1/N:2*pi*(1-1/N);
%     if (N == 6)
%         theta(1,1:N,N) = pi/12+theta(1,1:N,N);
%     end
% end
% 
% xi_r = zeros(1,Nvn);
% xi_t = zeros(1,Nvn);
% eta_r = zeros(1,Nvn);
% eta_t = zeros(1,Nvn);
% 
% IndS = 1; jStart = 0;
% for i = unique(Symms)
%     Indj = 1;
%     for j = jStart + (1:i*Symms_count(IndS))
%         if (j > jStart+1 && mod(j-(jStart+1),Symms_count(IndS)) == 0)
%             Indj = Indj + 1;
%         end
%         xi_r(1,j)  = cos(theta(1,Indj,i));
%         xi_t(1,j)  = -r*sin(theta(1,Indj,i));
%         eta_r(1,j) = sin(theta(1,Indj,i));
%         eta_t(1,j) = r*cos(theta(1,Indj,i));
%     end
%     IndS = IndS+1;
%     jStart = jStart+j;
% end
end


r = sqrt(sum(rst.^2,2));
theta = atan2(rst(:,2),rst(:,1));

xi_r = cos(theta);
xi_t = -r.*sin(theta);
eta_r = sin(theta);
eta_t = r.*cos(theta);

Dr = diag(xi_r)*Dxi + diag(eta_r)*Deta;
Dt = diag(xi_t)*Dxi + diag(eta_t)*Deta
% Dt = BasisTRI(3,xir_vGs)*T

F = zeros(Nvn,Nvn);
IndS = 0;
Indi = 1;
for i = unique(Symms)
    range_l = 1:i;
    for j = 1:Symms_count(Indi)
        range_g = IndS+range_l;
        F(range_g,range_g) = FN(range_l,range_l,i);

        IndS = IndS+i;
    end
    Indi = Indi+1;
end
perm = I(reorder,:);

F = perm'*F*perm;

Dtmp = F*Dr*diag(Fscale)*F';
Dreal = real(Dtmp); Dreal(abs(Dreal) < 1e-14) = 0;
Dimag = imag(Dtmp); Dimag(abs(Dimag) < 1e-14) = 0;
Dreal;
Dimag;
Drtilde = Dreal+1i*Dimag;

Dtmp = F*Dt*diag(Fscale)*F';
Dreal = real(Dtmp); Dreal(abs(Dreal) < 1e-14) = 0;
Dimag = imag(Dtmp); Dimag(abs(Dimag) < 1e-14) = 0;
Dreal
Dimag
Dthetatilde = Dreal+1i*Dimag;

NCols = 1;
RHS_In = rand([Nvn NCols]);
OPtilde = Dthetatilde;
% OPtilde = Drtilde;

O_std = Dt*RHS_In;
O_new2_slow = diag(Fscale)*F'*OPtilde*F*RHS_In;

% Smart DFT
FT_In = zeros(Nvn,NCols);
Inds = 1;
Indg = 0;
for s = unique(Symms)
    sep = Symms_count(Inds);

    % may want to use the transpose of IndRHS as in 'Multiplication...'
    % below (to reduce memory stride)
    IndRHS = zeros(sep,s);
    for i = 0:sep-1
    for j = 0:s-1
        IndRHS(i+1,j+1) = i+1+j*sep;
    end
    end
    IndRHS = IndRHS+Indg;
        
    % Use a switch statement for this.
    if (s == 1) % 0 symmetry
        % 0
        ilocal = 0;
        for i = 0:sep-1
            FT_In(IndRHS(i+1,ilocal+1),:) = RHS_In(IndRHS(i+1,1),:);
        end
    elseif (s == 3) % 3 symmetry
        % 0
        ilocal = 0;
        for i = 0:sep-1
            FT_In(IndRHS(i+1,ilocal+1),:) = ...
                sum(RHS_In(IndRHS(i+1,:),:));
        end
        
        % 1
        ilocal = 1;
        for i = 0:sep-1
            FT_In(IndRHS(i+1,ilocal+1),:) = ...
                RHS_In(IndRHS(i+1,1),:) - ...
                1/2*sum(RHS_In(IndRHS(i+1,2:3),:));
        end
        
        % 2
        ilocal = 2;
        for i = 0:sep-1
            FT_In(IndRHS(i+1,ilocal+1),:) = sqrt(3)/2*(...
                RHS_In(IndRHS(i+1,2),:) - ...
                RHS_In(IndRHS(i+1,3),:));
        end
    elseif (s == 6) % 6 symmetry
        implement this.
    end
    
    Indg = Indg+s;
    Inds = Inds+1;
end

% [FT_In F*RHS_In]


% Multiplication in Fourier space
% Note: FT_Out can use same memory space as FT_In if Operator matrices are
%       square.
FT_Out = zeros(Nvn,NCols);
Inds = 1;
Indg = 0;
for s = unique(Symms)
    sep = Symms_count(Inds);
        
    % Transpose of IndRHS above (reduced memory stride)
    IndRHS = zeros(s,sep);
    for j = 0:sep-1
    for i = 0:s-1
        IndRHS(i+1,j+1) = j+1+i*sep;
    end
    end
    IndRHS = IndRHS+Indg;
       
    % switch statement
    if (s == 1)
        % Complicated by the fact that the operator is not diagonalized. It
        % was suggested in Fladrich(2008) that adding terms to make this
        % into a 3-symmetry ended up being faster (coming closer to peak
        % flops). Not sure if this is valid at low order though (P <= 4) =>
        % implement both starting without transformation to 3-symmetry.
        
        % Do s = 1 last as it affects previously computed components
    elseif (s == 3)
        OP = zeros(sep,sep,s);
        
        OP(:,:,1) = OPtilde(IndRHS(1,:),IndRHS(1,:));
        OP(:,:,2) = real(OPtilde(IndRHS(2,:),IndRHS(2,:)));
        OP(:,:,3) = imag(OPtilde(IndRHS(3,:),IndRHS(3,:)));
        
        % 0
        ilocal = 0;
        FT_Out(IndRHS(ilocal+1,:),:) = ...
            OP(:,:,1)*FT_In(IndRHS(1,:),:);
        
        % 1
        ilocal = 1;
        FT_Out(IndRHS(ilocal+1,:),:) = ...
            OP(:,:,2)*FT_In(IndRHS(2,:),:) - ...
            OP(:,:,3)*FT_In(IndRHS(3,:),:);
        
        % 2
        ilocal = 2;
        FT_Out(IndRHS(ilocal+1,:),:) = ...
            OP(:,:,3)*FT_In(IndRHS(2,:),:) + ...
            OP(:,:,2)*FT_In(IndRHS(3,:),:);
        
    elseif (s == 6)
        implement this.
        
    end

    Indg = Indg+s;
    Inds = Inds+1;
end

if (find(unique(Symms) == 1))
    Inds = 1;
    Indg = 0;
    
    % Compute FT_Out(1,:)
    for s = unique(Symms)
        if (s == 1)
            FT_Out(1,:) = OPtilde(1,1)*FT_In(1,:);
        elseif (s == 3)
            sep = Symms_count(Inds);
        
            IndRHS = zeros(s,sep);
            for j = 0:sep-1
            for i = 0:s-1
                IndRHS(i+1,j+1) = j+1+i*sep;
            end
            end
            IndRHS = IndRHS+Indg;
    
            % Need to convert to complex as the entries in the transformed
            % operator do not fit in the alternate format.
            
            % 0
            ilocal = 0;
            FT_InComplex = FT_In(IndRHS(1,:),:);
            
            FT_Out(1,:) = FT_Out(1,:) + ...
                    OPtilde(1,IndRHS(ilocal+1,:))*FT_InComplex;
                
            % 1
            ilocal = 1;
            FT_InComplex = FT_In(IndRHS(2,:),:)-1i*FT_In(IndRHS(3,:),:);
            
            FT_Out(1,:) = FT_Out(1,:) + ...
                    OPtilde(1,IndRHS(ilocal+1,:))*FT_InComplex;
                
            % 2
            ilocal = 2;
            FT_InComplex = FT_In(IndRHS(2,:),:)+1i*FT_In(IndRHS(3,:),:);
            
            FT_Out(1,:) = FT_Out(1,:) + ...
                    OPtilde(1,IndRHS(ilocal+1,:))*FT_InComplex;
        elseif (s == 6)
            implement this
        end
        
        Indg = Indg+s;
        Inds = Inds+1;
    end
    
    % Update FT_Out(2:end,:)
    % Investigate whether first column elements are always real.
    % For now I am assuming that they are real as this is what I have
    % observed so far.
    if (norm(imag(OPtilde(:,1)),'inf') ~= 0)
        Error need to change routine to account for complex.
    end
    
    % Outer product
    % FT_In(1,:) must be real
    % OPtilde(IndRHS(ilocal+1,1)) always real?
    Add_real = OPtilde(2:end,1)*FT_In(1,:);
    
    Inds = 1;
    Indg = 0;
    for s = unique(Symms)
        if (s == 1)
            % Do nothing (already added above).
        elseif (s == 3)
            sep = Symms_count(Inds);
            
            IndRHS = zeros(s,sep);
            for j = 0:sep-1
            for i = 0:s-1
                IndRHS(i+1,j+1) = j+1+i*sep;
            end
            end
            IndRHS = IndRHS+Indg;
            
            %Convert result to alternate format (assuming only real)
            FT_Out(IndRHS(1,:),:) = FT_Out(IndRHS(1,:),:) + ...
                Add_real(IndRHS(1,:)-1,:);
            
            FT_Out(IndRHS(2,:),:) = FT_Out(IndRHS(2,:),:) + ...
                Add_real(IndRHS(2,:)-1,:);
            
        elseif (s == 6)
            implement this
            
        end
        Indg = Indg+s;
        Inds = Inds+1;
    end
end

% [FT_Out OPtilde*F*RHS_In]

% Smart Inverse DFT
O_new2 = zeros(Nvn,NCols);
Inds = 1;
Indg = 0;
for s = unique(Symms)
    sep = Symms_count(Inds);

    % may want to use the transpose of IndRHS as in 'Multiplication...'
    % below (to reduce memory stride)
    IndRHS = zeros(sep,s);
    for i = 0:sep-1
    for j = 0:s-1
        IndRHS(i+1,j+1) = i+1+j*sep;
    end
    end
    IndRHS = IndRHS+Indg;
        
    % Use a switch statement for this.
    if (s == 1) % 0 symmetry
        % 0
        O_new2(1,:) = FT_Out(1,:);
    elseif (s == 3) % 3 symmetry
        
        % 0
        ilocal = 0;
        for i = 0:sep-1
            O_new2(IndRHS(i+1,ilocal+1),:) = ...
                FT_Out(IndRHS(i+1,1),:)+2*FT_Out(IndRHS(i+1,2),:);
        end
        
        % 1
        ilocal = 1;
        for i = 0:sep-1
            O_new2(IndRHS(i+1,ilocal+1),:) = ...
                FT_Out(IndRHS(i+1,1),:) - ...
                FT_Out(IndRHS(i+1,2),:) + ...
                sqrt(3)*FT_Out(IndRHS(i+1,3),:);
        end
        
        % 2
        % Note: repeated calculations here
        ilocal = 2;
        for i = 0:sep-1
            O_new2(IndRHS(i+1,ilocal+1),:) = ...
                FT_Out(IndRHS(i+1,1),:) - ...
                FT_Out(IndRHS(i+1,2),:) - ...
                sqrt(3)*FT_Out(IndRHS(i+1,3),:);
        end
    elseif (s == 6) % 6 symmetry
        implement this.
    end
    
    Indg = Indg+s;
    Inds = Inds+1;
end

O_new2 = diag(Fscale)*O_new2;


% [O_new2 diag(Fscale)*F'*Dthetatilde*F*RHS_In]

% Dthetatilde

norm(O_std-O_new2_slow,'inf')
norm(O_std-O_new2,'inf')

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












