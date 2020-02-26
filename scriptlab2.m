A = [-4, -4, 0; 1, 0 ,0; 0, 1, 0];
B = [1;0;0];
C = [0,1,2];
D = 0;
x0 = [1;-1;1];
syms s;
syms t;
u = 1/s;
phi(t) = ilaplace(inv(s*eye(3)-A))
xl(t) = phi*x0
T(s) = C*inv(s*eye(3)-A)*B
yf(s) = u*T(s)
yf(t) = ilaplace(yf(s))
xf(t) = int(phi(t-tau)*B*1,tau,[0;t])
y(t) = C*(xl(t) + xf(t))