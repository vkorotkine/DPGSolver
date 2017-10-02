clc;
clf; hold on;
format short

% ToBeDeleted: Convert this to python!
% Note: Only need basisTRI of order 1 for barycentric coordinates
d = 2;

% NodeType = 'AO';
% NodeType = 'WS';
NodeType = 'WV';

if     (~isempty(strfind(NodeType,'AO'))) PMax = 16;
elseif (~isempty(strfind(NodeType,'WS'))) PMax = 8;
elseif (~isempty(strfind(NodeType,'WV'))) PMax = 20;
end

for P = 1:PMax
% for P = 5:5


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

% barycentric coordinates
Lv = basis_TRI_P1(rst)/basis_TRI_P1(rst_c);

% convert to equilateral triangle
if (~isempty(strfind(NodeType,'WS')) || ...
    ~isempty(strfind(NodeType,'WV')))
        rst = Lv*rst_cEq;
end

subplot(2,1,1); hold on;
x = rst(:,1); y = rst(:,2);
scatter(x,y);

% sort barycentric coordinates by first column, reorder rst
[~,I] = sortrows(1-Lv,1);
rst = rst(I,:);
if (wPresent)
    w = w(I,:);
end

% find angles of nodes
t = mod(atan2(rst(:,2),rst(:,1))*360/(2*pi)+360,360);

% keep only those falling between [150, 270[ degrees
IndOver = zeros(1,Nn);
n1 = 0;
for i = 1:Nn
    if (t(i) >= (150-1e4*eps) && t(i) < (270-1e4*eps) || ...
        norm(rst(i,:),'inf') < 1e2*eps)
        n1 = n1 + 1;
        IndOver(n1) = i;
    end
end

rst_reduced = zeros(n1,d);
w_reduced   = zeros(n1,d);
for i = 1:n1
    rst_reduced(i,:) = rst(IndOver(i),:);
    if (wPresent)
        w_reduced(i,:) = w(IndOver(i),:);
    end
end

% find Lv, and t again for potentially reordered reduced basis
Lv_reduced = basis_TRI_P1(rst_reduced)/basis_TRI_P1(rst_cEq);
t = mod(atan2(rst_reduced(:,2),rst_reduced(:,1))*360/(2*pi)+360,360);

% sort by groups of equal Lv1 based on t
IndFinal = 1:n1;
for i = 1:n1
    count = 0;
    for j = i+1:n1
        if (norm(Lv_reduced(j,1)-Lv_reduced(i,1),'inf') < 100*eps)
            count = count + 1;
        else
            break; % break out of j loop
        end
    end
    
    if (count > 0)
        [~,I] = sort(t(i:i+count),1);
        IndFinal(i:i+count) = IndFinal(i-1+I);
    end
end

% Compute final Lv (To be stored)
rst_reduced = rst_reduced(IndFinal,:);
if (wPresent)
    w_reduced = w_reduced(IndFinal);
end
Lv_reduced = basis_TRI_P1(rst_reduced)/basis_TRI_P1(rst_cEq);

% Compute symmetries
% If the row in L_v does not have equal entries, this represents a 3
% symmetry, otherwise a 1 symmetry. Same for TETs except only looking at
% first 3 entries of the row; for PYRs first 4 entries
symms_count = zeros(1,2);
symms = [3 1];
Ncompare = 3;
for i = 1:n1
    symm1 = 0;
    norm_sum = 0;
    for j = 2:Ncompare
        norm_sum = norm_sum + norm(Lv_reduced(i,1)-Lv_reduced(i,j),'inf');
    end
    if (norm_sum > 100*eps)
        symms_count(1) = symms_count(1) + 1;
    else
        symms_count(2) = symms_count(2) + 1;
    end
end

disp(num2str([n1 2]));
disp(num2str([symms_count(1) symms(1)]));
disp(num2str([symms_count(2) symms(2)]));
disp(num2str(Lv_reduced,'%20.15f'));


subplot(2,1,1)
x = rst_reduced(:,1); y = rst_reduced(:,2);
scatter(x,y,'filled');
a = (1:n1)'; b = num2str(a); c = cellstr(b);
dx = 0.025; dy = 0.025;
text(x+dx,y+dy,c);

% Recover the nodes using permutations of the barycentric coordinates
Nn = sum(symms.*symms_count);

Lv = zeros(Nn,3);
w  = zeros(Nn,1);

IndLv = 0;
iStart = 0;
for i = 1:2
    if (i > 1)
        iStart = iStart + symms_count(i-1);
    end
    for j = 1:symms_count(i)
        IndLv = IndLv + 1;
        Lv(IndLv,:) = Lv_reduced(iStart+j,:);
        if (wPresent)
            w(IndLv) = w_reduced(iStart+j);
        end
        for k = 2:symms(i)
            IndLv = IndLv + 1;

            for l = 1:2
                Lv(IndLv,l+1) = Lv(IndLv-1,l);
            end
            Lv(IndLv,1) = Lv(IndLv-1,3);
            
            if (wPresent)
                w(IndLv) = w(IndLv-1);
            end
        end
    end
end

rst = Lv*rst_cEq;

subplot(2,1,2)
x = rst(:,1); y = rst(:,2);
scatter(x,y,'filled');
a = (1:Nn)'; b = num2str(a); c = cellstr(b);
dx = 0.025; dy = 0.025;
text(x+dx,y+dy,c);

disp(' ');
disp(' ');
disp(num2str(rst,'%20.15f'));
if (wPresent)
    disp(' ');
    disp(num2str(w,'%20.15f'));
end

% Output to file
fName = [NodeType num2str(P) '.txt'];

fID = fopen(['Testing/Reordering/tri/updated/' fName],'w');
fprintf(fID,'%d %d\n',n1,2);
for i = 1:2
    fprintf(fID,'%d %d \n',symms_count(i),symms(i));
end
fprintf(fID,'\n');

fprintf(fID,'Barycentric Coordinates (Regular TRI)\n');
for i = 1:n1
    fprintf(fID,'%18.15f ',Lv_reduced(i,:));
    fprintf(fID,'\n');
end

if (wPresent)
    fprintf(fID,'\n');
    fprintf(fID,'Weights\n');
    for i = 1:n1
        fprintf(fID,'%18.15f \n',w_reduced(i));
    end
end

fclose(fID);

end