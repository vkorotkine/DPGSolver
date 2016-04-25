function [rst] = BCoords_To_rst(Lv,symms_ind,symms_count,NLv,Nc,d)
% Lv = Lv_reduced;
% symms_ind = symms_ind_reduced;
% symms_count;
% NLv = NsymmsP;
% Nc;
% d;
    
	if (Nc == 3) %TRI
        rst_c = [ -1 -1/sqrt(3); ...
                   1 -1/sqrt(3); ...
                   0  2/sqrt(3)];
        symms = [1 3 6];
        Nn = sum(symms.*symms_count);
        
        LvC = zeros(Nn,Nc);
        
        % static array
        perms = [ 1 2 3; 3 1 2; 2 3 1; ...
                  1 3 2; 2 1 3; 3 2 1];
        
        IndLvC = 0;
        for i = 1:NLv
            Lv_tmp = Lv(i,:);
            Nperm = symms_ind(i);
            
            LvC(IndLvC+(1:Nperm),:) = Lv_tmp(perms(1:Nperm,:));
            IndLvC = IndLvC + Nperm;
        end
        rst = LvC*rst_c;
	end



return