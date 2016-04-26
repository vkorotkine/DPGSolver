function [ind] = find_symmetry_ind(count_unique,Nc)
    
	if (Nc == 3) %TRI
        if     (count_unique(1) == 3); ind = 1;
        elseif (count_unique(1) == 2); ind = 2;
        elseif (count_unique(1) == 1); ind = 3;
        end
    elseif (Nc == 4) %TET
        if     (count_unique(1) == 4); ind = 1;
        elseif (count_unique(1) == 3); ind = 2;
        elseif (count_unique(1) == 2)
            if     (count_unique(2) == 2); ind = 3;
            elseif (count_unique(2) == 1); ind = 4;
            end
        elseif (count_unique(1) == 1); ind = 5;
        end
	end



return