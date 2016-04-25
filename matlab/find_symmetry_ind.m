function [ind] = find_symmetry_ind(count_unique,Nc)
    
	if (Nc == 3) %TRI
        if     (count_unique(1) == 3); ind = 1;
        elseif (count_unique(1) == 2); ind = 2;
        elseif (count_unique(1) == 1); ind = 3;
        end
	end



return