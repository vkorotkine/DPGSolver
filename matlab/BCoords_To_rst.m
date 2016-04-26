% function [rst] = BCoords_To_rst(Lv,symms_ind,symms_count,NLv,Nc)
Lv = Lv_reduced;
symms_ind = symms_ind_reduced;
symms_count;
NLv = NsymmsP;
Nc;
    
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
    elseif (Nc == 4) %TET
        rst_c = [ -1 -1/sqrt(3) -1/sqrt(6); ...
                   1 -1/sqrt(3) -1/sqrt(6); ...
                   0  2/sqrt(3) -1/sqrt(6); ...
                   0  0          3/sqrt(6)];
        symms = [1 4 6 12 24];
        Nn = sum(symms.*symms_count);
        
        % count number of 1 symmetries
        N1 = sum(symms_count(1:2));
        
        LvC = zeros(Nn,Nc);
        
        %static array
        perms_TRI = [ 1 2 3; 3 1 2; 2 3 1; ...
                      1 3 2; 2 1 3; 3 2 1];
              
        IndLvC = 0;
        Ind1 = 0;
        for i = 1:NLv
            Lv_tmp = Lv(i,:);
            
            if (symms_ind(i) == 1)
                Ind1 = Ind1 + 1;
                % 1 symmetry
                LvC(Nn-N1+Ind1,:) = Lv_tmp([1 2 3 4]);
            elseif (symms_ind(i) == 4)
                Ind1 = Ind1 + 1;
                % 1 symmetry
                LvC(Nn-N1+Ind1,:) = Lv_tmp([2 3 4 1]); % or [2 2 2 1]
                
                % 3 symmetry
                Lv_tmpTRI = Lv_tmp;
                Nperm = 3;
                LvC(IndLvC+(1:Nperm),:) = ...
                    Lv_tmpTRI([perms_TRI(1:Nperm,:) ones(Nperm,1)*4]);
                
                IndLvC = IndLvC + Nperm;
            elseif (symms_ind(i) == 6)
                % 3 symmetries (2x)
                Lv_tmpTRI = Lv_tmp;
                Nperm = 3;
                LvC(IndLvC+(1:Nperm),:) = ...
                    Lv_tmpTRI([perms_TRI(1:Nperm,:) ones(Nperm,1)*4]);
                
                IndLvC = IndLvC + Nperm;
                
                Lv_tmpTRI = Lv_tmp([1 3 4 2]);
                Nperm = 3;
                LvC(IndLvC+(1:Nperm),:) = ...
                    Lv_tmpTRI([perms_TRI(1:Nperm,:) ones(Nperm,1)*4]);
                
                IndLvC = IndLvC + Nperm;
                
            end
            
        end
                  
    end
    
    LvC



return
