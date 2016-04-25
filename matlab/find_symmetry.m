function [Lv_unique,count_unique,Nuniq] = find_symmetry(Lv,Nc)

Lv_unique = zeros(1,Nc);
count_unique = zeros(Nc,1);
Nuniq = 0;
if (Nc == 3) % TRI
    for i = 1:Nc
        duplicate = 0;
        for j = 1:i-1
            tmp = Lv_unique(j);
            
            if (norm(tmp-Lv(1,i),'inf') < 1e2*eps)
                duplicate = 1;
                count_unique(j) = count_unique(j) + 1;
                break;
            end
        end
        if (~duplicate)
            Nuniq = Nuniq + 1;
            Lv_unique(Nuniq) = Lv(i);
            count_unique(i) = count_unique(i) + 1;
        end
    end
end

Lv_unique(Nuniq+1:Nc) = [];
count_unique(Nuniq+1:Nc) = [];


return