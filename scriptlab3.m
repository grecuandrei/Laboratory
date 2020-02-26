A = [-1, 0, 0; -1, -2 ,-101; 1, 1, 0];
B = [0,1;1,-1;0,0];
C = [1,0,101];
D = [0,0];
syms s;
syms h;
syms tau;
phi(h) = ilaplace(inv(s*eye(3)-A),h);
Ad=phi(h);
Ad1 = double(phi(0.01))
integr(h) = int(phi(tau)*B*1,tau,[0;h]);
Bd = integr(h);
Bd1 = double(integr(0.01))
Cd = C;
Dd = D;
[num1s, den1s] = ss2tf(A,B,C,D,1)
[num2s, den2s] = ss2tf(A,B,C,D,2)
H11s = tf(num1s, den1s)
H21s = tf(num2s, den2s)
H11z = c2d(H11s,0.01)
H21z = c2d(H21s,0.01)
[num1z, den1z] = tfdata(H11z,'v')
[num2z, den2z] = tfdata(H21z,'v')



