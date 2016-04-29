function [Lv_unique,count_unique,place_unique] = find_symmetry3(Lv,Nc)
% i = 12; tmp2 = Lv; Lv = tmp2(i,:);

Lv_unique = zeros(1,Nc);
count_unique = zeros(Nc,1);
place_unique = zeros(1,min(Nc,4));
Nuniq = 0;

for i = 1:min(Nc,4)
    duplicate = 0;
    for j = 1:Nuniq
        tmp = Lv_unique(j);

        if (norm(tmp-Lv(1,i),'inf') < 1e2*eps)
            duplicate = 1;
            count_unique(j) = count_unique(j) + 1;
            place_unique(i) = place_unique(j);
            break;
        end
    end
    if (~duplicate)
        Nuniq = Nuniq + 1;
        Lv_unique(Nuniq) = Lv(i);
        count_unique(Nuniq) = count_unique(Nuniq) + 1;
        place_unique(i) = Nuniq;
    end
end

count_unique(Nuniq+1:Nc) = [];
Lv_unique(Nuniq+1:Nc) = [];

return