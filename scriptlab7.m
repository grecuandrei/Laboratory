A = [0,0,0;1,0,-1;0,1,2]
B = [0,-1;1,0;0,1]
C = [1,0,1]
D = [0, 0]
lambdad = [-1,-1,-1]
lambdae = [-2,-2,-2]
n = length(A);
calc = B;
R = B;
for i = 1:n-1
    calc = A*calc;
    R = [R calc];
end
r = rank(R)
if (r == n)
    calc = C;
    Q = C;
    for i = 1:n-1
        calc = calc*A;
        Q = [Q ;calc];
    end
    r = rank(Q)
    if (r == n)
        g = [1;1];
        F0 = [0,0,0;0,0,0];
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
        Astar = A.'
        Bstar = C.'
        m = size(Bstar,2);
        Rstar = Bstar;
        calc = Bstar;
        for i = 1:n-1
            calc = Astar*calc;
            Rstar = [Rstar calc];
        end
        r = rank(Rstar)
        if (r == n)
            R0_inv = inv(Rstar)
            qt = R0_inv(n,:)
            pol_carac = 1;
            for i = 1:length(lambdae)
                pol_carac = pol_carac * (Astar+eye(n)*(-lambdae(i)))
            end
            Lt = (-qt) * pol_carac
            L = Lt.'
            J = A + L*C
            K = -L
            H = B
            M = eye(n)
            N = zeros(n,m)
            ne = n
            Ac = J + H*F*M
            Bc = K + H*F*N
            Fc = F*M
            Gc = F*N
            nc = n
        end
    end
end


