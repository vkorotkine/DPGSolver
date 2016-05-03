clc;
clf; hold on;
format long

% Convert to python.

folder_name = 'Testing/Reordering/pyr/';
d = 3;

% NodeType = 'GLL';
% NodeType = 'GL';
NodeType = 'WVHToP';

if     (~isempty(strfind(NodeType,'WVHToP'))); PMax = 11;
elseif (~isempty(strfind(NodeType,'GLL')));    PMax = 10;
elseif (~isempty(strfind(NodeType,'GL')));     PMax = 10;
end

for P = 1:PMax
% for P = 4

if (~isempty(strfind(NodeType,'GL')))
    [abc,wTP,Nn]    = CubatureTensorProd(GLOBAL,P,d,NodeType);
else
    NnVec = [1 6 6 14 14 34 34 58 58 90 90];
    Nn = NnVec(P);
    
    fName = ['witherden-vincent-n' ...
             num2str(Nn) '-d' num2str(P) '-sp.txt'];
    fID = fopen([folder_name 'hex/' fName],'r');

    formatSpec = '%f %f %f %f';
    sizeA = [4 Inf];
    A = fscanf(fID,formatSpec,sizeA);
    A = A';

    fclose(fID);

    abc = A(:,1:3);
    wTP = A(:,4);
end

rst = [0.5*(1-abc(:,3)).*abc(:,1) ...
       0.5*(1-abc(:,3)).*abc(:,2) ...
       sqrt(2.0)/2.0*(0.6+abc(:,3))];
   
if (P > 1 || P > 0 && ~isempty(strfind(NodeType,'GL')))
    w = wTP*2.0^(-2.5).*(1.0-abc(:,3)).^2;
else
    w = wTP*1/(3*sqrt(2)).*(1.0-abc(:,3)).^2;
end
sum(w)

% x = rst(:,1); y = rst(:,2); z = rst(:,3);
% scatter3(x,y,z,'filled','b');

% Output to file
if (~isempty(strfind(NodeType,'GL')))
    fName = [NodeType 'W' num2str(P) '.txt'];
    fID = fopen([folder_name fName],'w');
else
    fName = [NodeType num2str(P) '.txt'];
    fID = fopen([folder_name 'HEX_To_PYR/' fName],'w');
end

for i = 1:Nn
    fprintf(fID,' %18.15f',[rst(i,:) w(i)]);
    if (i ~= Nn)
        fprintf(fID,'\n');
    end
end

fclose(fID);

end