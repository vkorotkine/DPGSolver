clc;
clf; hold on;
format long

% ToBeDeleted: Convert this to python!
% Note: Only need basisTET of order 1 for barycentric coordinates
plot_on = 1;

GLOBAL.TET.xir_vGs = [          [0 -1 0 1]' ...
                      1/sqrt(3)*[0 -1 2 -1]' ...
                      1/sqrt(6)*[3 -1 -1 -1]'];

folder_name = 'Testing/Reordering/tet/updated/';
EType = 'TET';
d = 3;
Nc = 4;
Nsymms = 5;

NodeType = 'AO';
% NodeType = 'SH';
% NodeType = 'WV';

if     (~isempty(strfind(NodeType,'AO'))); PMax = 15;
elseif (~isempty(strfind(NodeType,'SH'))); PMax = 6;
elseif (~isempty(strfind(NodeType,'WV'))); PMax = 10;
end

% for P = 1:PMax
for P = 2


% Read/Obtain rst
wPresent = 1;
[rst_cEq,~,Nn_c] = CubatureTET(GLOBAL,1,'alpha-opt');

if (~isempty(strfind(NodeType,'AO')))
    wPresent = 0;
    
    [rst_c,~,Nn_c] = CubatureTET(GLOBAL,1,'alpha-opt');
    [rst,w,Nn]     = CubatureTET(GLOBAL,P,'alpha-opt');
else
    wScale = 1/2*sqrt(2);
    rst_c = [[-1 -1 1 -1]' [1 -1 -1 -1]' [-1 -1 -1 1]'];
    if (~isempty(strfind(NodeType,'SH')))
        Nn = 1/6*(P+1)*(P+2)*(P+3);
        dP = [2 3 5 6 8 9];
    
        fName = ['shunn-ham-n' ...
                 num2str(Nn) '-d' num2str(dP(P)) '-spu.txt'];
    elseif (~isempty(strfind(NodeType,'WV')))
        NnVec = [1 4 8 14 14 24 35 46 59 81];
        Nn = NnVec(P);
        
        fName = ['witherden-vincent-n' ...
                 num2str(Nn) '-d' num2str(P) '-sp.txt'];
    end
	fID = fopen(['Testing/Reordering/tet/' fName],'r');
    
    formatSpec = '%f %f %f %f';
    sizeA = [4 Inf];
    A = fscanf(fID,formatSpec,sizeA);
    A = A';
    
    fclose(fID);
    
    rst = A(:,1:3);
    w   = wScale*A(:,4);
end

% break;

% barycentric coordinates
Lv = basis_TET_P1(rst)/basis_TET_P1(rst_c);
% [rst ones(Nn,1)]*inv([rst_c ones(Nc,1)])
% [rst ones(Nn,1)]/[rst_c ones(Nc,1)]

break;

% convert to regular tetrahedron
if (~isempty(strfind(NodeType,'WS')) || ...
    ~isempty(strfind(NodeType,'WV')))
        rst = Lv*rst_cEq;
end

if (plot_on)
    x = rst(:,1); y = rst(:,2); z = rst(:,3);
    plot3(x,y,z,'o');
    view(32,34)
end

% Find non-redundant barycentric coordinates
[~,I] = sortrows(1-Lv,1);
Lv = Lv(I,:);
if (wPresent)
    w = w(I);
end

[~,I] = sort(1-Lv,2);
for i = 1:Nn
    Lv(i,:) = Lv(i,I(i,:));
end

[~,I] = sortrows(1-Lv,1);
Lv = Lv(I,:);
if (wPresent)
    w = w(I);
end



Lv_tmp = zeros(Nn,Nc);
w_tmp  = zeros(Nn,1);
symms = [1 4 6 12 24]; %TET
symms_count = zeros(1,Nsymms);
symms_ind = zeros(1,Nn);

NsymmsP = 0;
for i = 1:Nn
    if (i == 1 || norm(Lv(i,1)-Lv(i-1,1),'inf') > 1e2*eps)
        NsymmsP = NsymmsP+1;
        
        if (wPresent)
            w_tmp(NsymmsP,1) = w(i);
        end
        
        % Find type of symmetry
        [Lv_unique,count_unique] = find_symmetry2(Lv(i,:),Nc);
        
        [~,I] = sort(count_unique);

        % Sort such that Lv is arranged with unique entries at the start of
        % the array
        IndLv = 0;
        for j = 1:length(count_unique)
            for k = 1:count_unique(I(j))
                IndLv = IndLv + 1;
                Lv_tmp(NsymmsP,IndLv) = Lv_unique(I(j));
            end
        end

        [~,I] = sort(Nc-count_unique);
        count_unique = count_unique(I);
        
        [ind] = find_symmetry_ind(count_unique,Nc);
        
        symms_count(ind) = symms_count(ind) + 1;
        
        symms_ind(NsymmsP) = symms(ind);
%         break;
    end
end

% symms
% symms_count
symms_ind(NsymmsP+1:Nn) = [];
Lv_tmp(NsymmsP+1:Nn,:) = [];
if (wPresent)
    w_tmp(NsymmsP+1:Nn,:) = [];
end

% Sort Lv_tmp such that highest symmetries come first
Lv_reduced = zeros(NsymmsP,Nc);
w_reduced = zeros(NsymmsP,1);
symms_ind_reduced = zeros(1,NsymmsP);

k = 0;
for j = Nsymms:-1:1
    symm = symms(j);
    for i = 1:NsymmsP
        if (symms_ind(i) == symm)
            k = k+1;
            Lv_reduced(k,:) = Lv_tmp(i,:);
            if (wPresent)
                w_reduced(k) = w_tmp(i);
            end
            symms_ind_reduced(k) = symms_ind(i);
        end
    end
end

if (k ~= NsymmsP)
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
end


% Output to file
fName = [NodeType num2str(P) '.txt'];

fID = fopen([folder_name fName],'w');
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