clc;
clf; hold on;
format long

% Convert to python.

folder_name = 'Testing/Reordering/pyr/';
d = 3;

% NodeType = 'WVHToP';
% NodeType = 'GLL';
NodeType = 'GL';
NodeType = 'GJ';

if     (~isempty(strfind(NodeType,'WVHToP'))); PMax = 11;
elseif (~isempty(strfind(NodeType,'GLL')));    PMax = 10;
elseif (~isempty(strfind(NodeType,'GL')));     PMax = 10;
elseif (~isempty(strfind(NodeType,'GJ')));     PMax = 10;
end

% for P = 1:PMax
for P = 1

if (~isempty(strfind(NodeType,'GL')))
    [abc,wTP,Nn]    = CubatureTensorProd(GLOBAL,P,d,NodeType);
elseif (~isempty(strfind(NodeType,'WVHToP')))
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
elseif (~isempty(strfind(NodeType,'GJ')))
    [abc2D,wTP2D,Nn2D]    = CubatureTensorProd(GLOBAL,P,d-1,'GL');
    [abc1D,wTP1D]=jags(P+1,2.0,0.0);
    
    wTP1D = wTP1D/sum(wTP1D)*2;
    
    Nn1D = Nn2D^(1/(d-1));
    abc = zeros(Nn1D^d,d);
    wTP = zeros(Nn1D^d,1);
    
    for k = 1:Nn1D
    for j = 1:Nn1D
    for i = 1:Nn1D
        abc((k-1)*Nn1D^2+(j-1)*Nn1D+(i),:) = [abc2D((j-1)*Nn1D+(i),:) abc1D(k,:)];
        wTP((k-1)*Nn1D^2+(j-1)*Nn1D+(i),:) = wTP2D((j-1)*Nn1D+(i),:)*wTP1D(k,:);
    end
    end
    end
    
%     abc
%     wTP = wTP*3/4
end

sum(wTP)

rst = [0.5*(1-abc(:,3)).*abc(:,1) ...
       0.5*(1-abc(:,3)).*abc(:,2) ...
       sqrt(2.0)/2.0*(0.6+abc(:,3))];
   
if (P > 1 || P > 0 && ~isempty(strfind(NodeType,'GL')))
    w = wTP*2.0^(-2.5).*(1.0-abc(:,3)).^2;
elseif (~isempty(strfind(NodeType,'WVHToP')))
    w = wTP*1/(3*sqrt(2)).*(1.0-abc(:,3)).^2;
elseif (P > 0 && ~isempty(strfind(NodeType,'GJ')))
     w = wTP*2.0^(-2.5).*(1.0-abc(:,3)).^2;
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