clc;
clf; hold on;
format short

% ToBeDeleted: Convert this to python!
% Note: Only need basisTRI of order 1 for barycentric coordinates
d = 2;
Nc = 3;
Nsymms = 3;

% NodeType = 'AO';
% NodeType = 'WS';
NodeType = 'WV';

if     (~isempty(strfind(NodeType,'AO'))) PMax = 16;
elseif (~isempty(strfind(NodeType,'WS'))) PMax = 8;
elseif (~isempty(strfind(NodeType,'WV'))) PMax = 20;
end

for P = 1:PMax
% for P = 6


% Read/Obtain rst
wPresent = 1;
[rst_cEq,~,Nn_c] = CubatureTRI(GLOBAL,1,'alpha-opt');
if (~isempty(strfind(NodeType,'AO')))
    wPresent = 0;
    
    [rst_c,~,Nn_c] = CubatureTRI(GLOBAL,1,'alpha-opt');
    [rst,w,Nn]     = CubatureTRI(GLOBAL,P,'alpha-opt');
else
    wScale = sqrt(3)/2;
    rst_c = [[-1 -1 1]' [1 -1 -1]'];
    if (~isempty(strfind(NodeType,'WS')))
        Nn = 1/2*(P+1)*(P+2);
        dP = [2 4 5 7 8 10 12 14];
    
        fName = ['williams-shunn-n' ...
                 num2str(Nn) '-d' num2str(dP(P)) '-spu.txt'];
    elseif (~isempty(strfind(NodeType,'WV')))
        NnVec = [1 3 6 6 7 12 15 16 19 25 28 33 37 42 49 55 60 67 73 79];
        Nn = NnVec(P);
        
        fName = ['witherden-vincent-n' ...
                 num2str(Nn) '-d' num2str(P) '-sp.txt'];
    end
	fID = fopen(['Testing/Reordering/tri/' fName],'r');
    
    formatSpec = '%f %f %f';
    sizeA = [3 Inf];
    A = fscanf(fID,formatSpec,sizeA);
    A = A';
    
    fclose(fID);
    
    rst = A(:,1:2);
    w   = wScale*A(:,3);
end

% break;

% barycentric coordinates
Lv = basis_TRI_P1(rst)/basis_TRI_P1(rst_c);

% convert to equilateral triangle
if (~isempty(strfind(NodeType,'WS')) || ...
    ~isempty(strfind(NodeType,'WV')))
        rst = Lv*rst_cEq;
end

hold on;
x = rst(:,1); y = rst(:,2);
scatter(x,y);


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
symms = [1 3 6]; %TRI
symms_count = zeros(1,Nc);
symms_ind = zeros(1,Nn);

NsymmsP = 0;
for i = 1:Nn
    if (i == 1 || norm(Lv(i,1)-Lv(i-1,1),'inf') > 1e2*eps)
        NsymmsP = NsymmsP+1;
        
        if (wPresent)
            w_tmp(NsymmsP,1) = w(i);
        end
        
        % Find type of symmetry
        [Lv_unique,count_unique,Nuniq] = find_symmetry(Lv(i,:),Nc);

%         Lv_tmp(NsymmsP,1:Nuniq) = Lv_unique;
        Lv_tmp(NsymmsP,:) = Lv(i,:);

        [~,I] = sort(Nc-count_unique);
        count_unique = count_unique(I);
        
        [ind] = find_symmetry_ind(count_unique,Nc);
        
        symms_count(ind) = symms_count(ind) + 1;
        
        symms_ind(NsymmsP) = symms(ind);
%         break;
    end
end

symms;
symms_count;
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
for j = Nc:-1:1
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

Lv_reduced
if (wPresent)
    w_reduced
end

rst = BCoords_To_rst(Lv_reduced,symms_ind_reduced,symms_count,NsymmsP,Nc,d);

daspect([1 1 0.0001])
x = rst(:,1); y = rst(:,2);
scatter(x,y,'filled');
a = (1:Nn)'; b = num2str(a); c = cellstr(b);
% c = cellstr([strcat(b,',') num2str(symms_ind_reduced')]);
dx = 0.025; dy = 0.025;
text(x+dx,y+dy,c);



% Output to file
fName = [NodeType num2str(P) '.txt'];

fID = fopen(['Testing/Reordering/tri/updated/' fName],'w');
fprintf(fID,'%d %d\n',NsymmsP,Nsymms);
for i = Nsymms:-1:1
    fprintf(fID,'%d %d \n',symms_count(i),symms(i));
end
fprintf(fID,'\n');

fprintf(fID,'Barycentric Coordinates (Regular TRI)\n');
for i = 1:NsymmsP
    fprintf(fID,'%18.15f ',Lv_reduced(i,:));
    fprintf(fID,'\n');
end

if (wPresent)
    fprintf(fID,'\n');
    fprintf(fID,'Weights\n');
    for i = 1:NsymmsP
        fprintf(fID,'%18.15f \n',w_reduced(i));
    end
end

fclose(fID);

end