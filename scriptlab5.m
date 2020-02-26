A =[0 , 0, 0; 1, 0 , -1; 0,1,2];
B = [0,-1;1,0;0,1];
C = [1,0,1];
lambdad = [-1, -1, -1];
n = length(A);
m = size(B,2);
calc = B;
R = B;
u = [1;-2]
for i = 1:n-1
    calc = A*calc;
    R = [R calc];
end
R
r = rank(R)
if (r == n)
    g = rand([m, 1]);
    F0 = rand([m,n]);
    A0 = A + B*F0
    b0 = B*g
    R0 = b0
    calc = b0;
    for i = 1:n-1
        calc = A0*calc;
        R0 = [R0 calc];
    end
    if (rank(R0) == n)
        R0_inv = inv(R0)
        qt = R0_inv(n,:)
        pol_carac = 1;
        for i = 1:length(lambdad)
            pol_carac = pol_carac * (A0+eye(n)*(-lambdad(i)))
        end
        ft = (-qt) * pol_carac
        F = F0 + g*ft
    end 
end
G = eye(m)




