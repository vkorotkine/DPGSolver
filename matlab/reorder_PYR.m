clc;
clf; hold on;
format long

% ToBeDeleted: Convert this to python!
% Note: Only need basisTET of order 1 for barycentric coordinates
plot_on = 0;

folder_name = 'Testing/Reordering/pyr/';
folder_nameU = [folder_name 'updated/'];
EType = 'PYR';
ETypeLower = 'pyr';
d = 3;
Nc = 5;
Nsymms = 3;
symms = [1 4 8];

NodeType = 'GLL';
% NodeType = 'GL';
% NodeType = 'WV';

if     (~isempty(strfind(NodeType,'GLL'))); PMax = 6;
elseif (~isempty(strfind(NodeType,'GL')));  PMax = 6;
elseif (~isempty(strfind(NodeType,'WV')));  PMax = 10;
end

rst_cEq = [          [-1  1  1 -1 0]' ...
                     [-1 -1  1  1 0]' ...
           sqrt(2)/5*[-1 -1 -1 -1 4]'];
rst_c   = [[-1  1  1 -1 0]' ...
           [-1 -1  1  1 0]' ...
           [-1 -1 -1 -1 1]'];

wPresent = 1;
wScale = 1/sqrt(2);

for P = 1:PMax
% for P = 8

% Read rst

if (~isempty(strfind(NodeType,'GL')))
    wPresent = 0;
    Nn = round(1/6*(P+1)*(P+2)*(2*P+3));
    
    if (~isempty(strfind(NodeType,'GLL')))
        fName = ['gauss-legendre-lobatto-n' num2str(Nn) '-su.txt'];
    elseif (~isempty(strfind(NodeType,'GL')))
        fName = ['gauss-legendre-n' num2str(Nn) '-su.txt'];
    end
    
    fID = fopen([folder_name fName],'r');

    formatSpec = '%f %f %f';
    sizeA = [3 Inf];
    A = fscanf(fID,formatSpec,sizeA);
    A = A';

    fclose(fID);

    rst = A(:,1:3);
elseif (~isempty(strfind(NodeType,'WV')))
    NnVec = [1 5 6 10 15 24 31 47 62 83];
    Nn = NnVec(P);
    
    fName = ['witherden-vincent-n' ...
             num2str(Nn) '-d' num2str(P) '-sp.txt'];

    fID = fopen([folder_name fName],'r');

    formatSpec = '%f %f %f %f';
    sizeA = [4 Inf];
    A = fscanf(fID,formatSpec,sizeA);
    A = A';

    fclose(fID);

    rst = A(:,1:3);
    w   = wScale*A(:,4);
end

Lv = [ones(Nn,1) rst prod(rst,2)]/[ones(Nc,1) rst_c prod(rst_c,2)];
% tmp_rst   = BasisTensorProd(1,rst);   tmp_rst   = tmp_rst(:,[1 2 3 5 8]);
% tmp_rst_c = BasisTensorProd(1,rst_c); tmp_rst_c = tmp_rst_c(:,[1 2 3 5 8]);
% Lv = tmp_rst/tmp_rst_c;
if (cond([ones(Nc,1) rst_c prod(rst_c,2)]) > 1e2)
    disp('High condition number in computation of BCoords');
    % Note: The condition number using basis_TP is actually higher than
    %       that of the non-normalized basis used above. THINK?
    cond([ones(Nc,1) rst_c sqrt(2)*prod(rst_c,2)])
    
    tmp_rst_c = BasisTensorProd(1,rst_c);
    tmp_rst_c = tmp_rst_c(:,[1 2 3 5 8]);
    cond(tmp_rst_c)
    break;
end

% Sort Lv such that barycentric coordinates are in descending order
for i = Nc:-1:1
    [~,I] = sortrows(1-Lv,i);
    Lv = Lv(I,:);
    if (wPresent)
        w = w(I);
    end
end

rst = Lv*rst_cEq;

if (plot_on)
    subplot(2,1,1); hold on; grid on;
    daspect([1 1 1]);
    
    x = rst_cEq(:,1); y = rst_cEq(:,2); z = rst_cEq(:,3);
    scatter3(x,y,z,'filled','b');
    
    x = rst(:,1); y = rst(:,2); z = rst(:,3);
    scatter3(x,y,z,'filled','g');
    view(73,22)
    
    a = (1:Nn)'; b = num2str(a); c = cellstr(b);
    dx = 0.025; dy = 0.025; dz = 0.025;
    text(x+dx,y+dy,z+dz,c);
end

% Find non-redundant barycentric coordinates
Lv_tmp = zeros(Nn,Nc);
w_tmp  = zeros(Nn,1);

symms_count = zeros(1,Nsymms);
symms_ind = zeros(1,Nn);

NsymmsP = 0;
for i = 1:Nn
    [Lv_unique,count_unique,place_unique] = find_symmetry3(Lv(i,:),Nc);
    Nuniq = size(count_unique,1);
    
    duplicate = 0;
    if (Nuniq == 1)
        % No chance of being repeated.
        NsymmsP = NsymmsP+1;
        Lv_tmp(NsymmsP,:) = Lv(i,:);
        
        symm_type = 1;
    else
        % Find if the Lv is a duplicate
        Lv_i_sorted = sort(Lv(i,1:min(Nc,4)));
        
        for j = 1:NsymmsP
            Lv_j_sorted = sort(Lv_tmp(j,1:min(Nc,4)));
            
            if (norm(Lv_i_sorted-Lv_j_sorted,'inf') < 1e1*eps)
                duplicate = 1;
                break;
            end
        end

        if (~duplicate)
            NsymmsP = NsymmsP+1;
            if (Nuniq == 2)
                if (max(count_unique) == 3); symm_type = 2;
                else                         symm_type = 2;
                end
            elseif (Nuniq == 3)
                symm_type = 2;
                tmp = [place_unique place_unique(1)];
                for j = 1:4
                    if (tmp(j) == tmp(j+1))
                        symm_type = 3;
                        break;
                    end
                end
            elseif (Nuniq == 4)
                symm_type = 3;
            end
            Lv_tmp(NsymmsP,:) = Lv(i,:);
        end
    end
    
    if (~duplicate)
        if (wPresent)
            w_tmp(NsymmsP,1) = w(i);
        end
        
        symms_count(symm_type) = symms_count(symm_type) + 1;
        symms_ind(NsymmsP) = symm_type;
    end
end

Lv_tmp(NsymmsP+1:Nn,:) = [];
symms_ind(NsymmsP+1:Nn) = [];
if (wPresent)
    w_tmp(NsymmsP+1:Nn,:) = [];
end

% symms_count
% 
% Lv_tmp
% symms_ind

% Sort Lv_tmp such that highest symmetries come first
[~,I] = sort(Nc-symms_ind);

Lv_reduced = Lv_tmp(I,:);
symms_ind_reduced = symms_ind(1,I);
if (wPresent)
    w_reduced = w_tmp(I,1);
end

% Check to make sure that the right number of nodes were found
if (sum(symms.*symms_count) ~= Nn)
    disp('Problem');
    break;
end

rst = BCoords_To_rst(Lv_reduced,symms_ind_reduced,symms_count,NsymmsP,Nc);

if (plot_on)
    Lv_reduced
    if (wPresent)
        w_reduced
    end

    rst
end

if (plot_on)
    subplot(2,1,2); hold on; grid on;
    daspect([1 1 1])
    x = rst(:,1); y = rst(:,2); z = rst(:,3);
    scatter3(x,y,z,'filled','g');

    r = sqrt(x.^2+y.^2);
    theta = 0:2*pi/100:2*pi;

    for i = 1:Nn
        xp = r(i)*cos(theta);
        yp = r(i)*sin(theta);
        zp = z(i)+theta*0;
        plot3(xp,yp,zp,'-r');
    end

    a = (1:Nn)'; b = num2str(a); c = cellstr(b);
    dx = 0.025; dy = 0.025; dz = 0.025;
    text(x+dx,y+dy,z+dz,c);
    
    x = rst_cEq(:,1); y = rst_cEq(:,2); z = rst_cEq(:,3);
    scatter3(x,y,z,'filled','bs');
end

% Output to file
fName = [NodeType num2str(P) '.txt'];

fID = fopen([folder_nameU fName],'w');
fprintf(fID,'%d %d\n',NsymmsP,Nsymms);
for i = Nsymms:-1:1
    fprintf(fID,'%d %d\n',symms_count(i),symms(i));
end
fprintf(fID,'\n');

fprintf(fID,['Barycentric Coordinates (Regular ' EType ')\n']);
for i = 1:NsymmsP
    fprintf(fID,' %18.15f',Lv_reduced(i,:));
    if (i ~= NsymmsP || wPresent)
        fprintf(fID,'\n');
    end
end

if (wPresent)
    fprintf(fID,'\n');
    fprintf(fID,'Weights\n');
    for i = 1:NsymmsP
        fprintf(fID,' %18.15f',w_reduced(i));
        if (i ~= NsymmsP)
            fprintf(fID,'\n');
        end
    end
end

fclose(fID);

end