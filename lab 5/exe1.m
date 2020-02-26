A = [0, 0, 0; 1, 0, -1; 0, 1, 2];
B = [0 -1;1 0;0 1];
C = [1 0 1];

lambda = [-20 -20 -20];
R = ctrb(A, B);
[n, m] = size(R);
if(rank(R) == n)
    g = [1; 1];
    F0 = [0 0 0; 0 0 0];
    A0 = A + B*F0;
    b0 = B*g;
    R0 = ctrb(A0, b0);
    [j, p] = size(R0);
    if(rank(R0) == j)
        R0_inv = inv(R0);
        q_t = R0_inv(p, :);
        psi_d = (A0 - (-20)*eye(3))*(A0 - (-20)*eye(3))*(A0 - (-20)*eye(3));
        f_t = (-1)*q_t*psi_d;
    end
    F = F0 + g*f_t;
    G = eye(2);
end
