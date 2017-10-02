function [Chi_v] = basis_TRI_P1(rs)
% rs = rst;

Nb = 3;
Nn = size(rs,1);

a = zeros(Nn,1);
b = zeros(Nn,1);

con = zeros(Nn,Nb);

for node = 1:Nn
    r = rs(node,1);
    s = rs(node,2);
    
    if (norm(2-sqrt(3)*s) > 1e-12)
        a(node,1) = 3*r/(2-sqrt(3)*s);
    else
        a(node,1) = 0;
    end
    b(node,1) = 1/3*(2*sqrt(3)*s-1);
    
    BF = 1;
    for i = 0:1
    for j = 0:1-i
        con(node,BF) = (1-b(node,1))^i;
        
        BF = BF + 1;
    end
    end
end

Chi_v = (2/3^0.25)*[1/2+0*a sqrt(1/8)*(3*b+1) sqrt(3/8)*(1-b).*a];

return