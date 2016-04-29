function [rst] = BCoords_To_rst(Lv,symms_ind,symms_count,NLv,Nc)
% Lv = Lv_reduced;
% symms_ind = symms_ind_reduced;
% symms_count;
% NLv = NsymmsP;
% Nc;

% static array
perms_TRI = [ 1 2 3; 3 1 2; 2 3 1; 1 3 2; 2 1 3; 3 2 1];
if (Nc == 3) %TRI
    rst_c = [ -1 -1/sqrt(3); ...
               1 -1/sqrt(3); ...
               0  2/sqrt(3)];
    symms = [1 3 6];
    Nn = sum(symms.*symms_count);

    LvC = zeros(Nn,Nc);

    IndLvC = 0;
    for i = 1:NLv
        Lv_tmp = Lv(i,:);
        Nperm = symms_ind(i);

        LvC(IndLvC+(1:Nperm),:) = Lv_tmp(perms_TRI(1:Nperm,:));
        IndLvC = IndLvC + Nperm;
    end
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
    perms_TET = [ 1 2 3 4; 4 1 2 3; 3 4 1 2; 2 3 4 1];

    IndLvC = 0;
    Ind1 = 0;
    for i = 1:NLv
        Lv_tmp = Lv(i,:);

        if (symms_ind(i) == 1)
            NTRIsymms = [1 0 0];
            TETperms = 1;
        elseif (symms_ind(i) == 4)
            NTRIsymms = [1 1 0];
            TETperms = [4 1];
        elseif (symms_ind(i) == 6)
            NTRIsymms = [0 2 0];
            TETperms = [1 4];
        elseif (symms_ind(i) == 12)
            NTRIsymms = [0 2 1];
            TETperms = [4 3 1];
        elseif (symms_ind(i) == 24)
            NTRIsymms = [0 0 4];
            TETperms = [1 2 3 4];
        end

        Indperm = 0;

        % Can only ever have one 1-symmetry in a 4-tuple of BCoords
        if (NTRIsymms(1) == 1)
            Indperm = Indperm + 1;

            Ind1 = Ind1 + 1;
            LvC(Nn-N1+Ind1,:) = Lv_tmp(perms_TET(TETperms(Indperm),:));
        end

        for TRIsymm = 2:3
            for j = 1:NTRIsymms(TRIsymm)
                Indperm = Indperm + 1;

                Nperm = sum(1:TRIsymm);
                Lv_tmpTRI = Lv_tmp(perms_TET(TETperms(Indperm),:));

                LvC(IndLvC+(1:Nperm),:) = ...
                    Lv_tmpTRI([perms_TRI(1:Nperm,:) ones(Nperm,1)*4]);

                IndLvC = IndLvC + Nperm;
            end
        end
    end
elseif (Nc == 5) %PYR
    perms_QUAD = [ 1 2 3 4; 4 1 2 3; 3 4 1 2; 2 3 4 1];
    rst_c = [          [-1  1  1 -1 0]' ...
                       [-1 -1  1  1 0]' ...
             sqrt(2)/5*[-1 -1 -1 -1 4]'];
    symms = [1 4 8];
    Nn = sum(symms.*symms_count);
    
    % count number of 1 symmetries
    N1 = symms_count(1);
    
    LvC = zeros(Nn,Nc);
    
    IndLvC = 0;
    Ind1 = 0;
    
    for i = 1:NLv
        Lv_tmp = Lv(i,:);
        
        if (symms_ind(i) == 1)
            Nperm = 1;
        else
            Nperm = 4;
            if (symms_ind(i) == 2); NQUADsymms = 1;
            else                    NQUADsymms = 2;
            end
        end
        
        if (Nperm == 1)
            Ind1 = Ind1+1;
            LvC(Nn-N1+Ind1,:) = Lv_tmp;
        else
            for j = 1:NQUADsymms
                if (j == 1); Lv_tmpQUAD = Lv_tmp;
                else         Lv_tmpQUAD = Lv_tmp([2 1 4 3 5]);
                end
                LvC(IndLvC+(1:Nperm),:) = ...
                        Lv_tmpQUAD([perms_QUAD ones(Nperm,1)*Nc]);
                IndLvC = IndLvC + Nperm;
            end
        end
    end
end

% LvC
rst = LvC*rst_c;



return
