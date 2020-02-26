syms s;
T(s) = [1/((s+1)*(s+3)), (s-2)/((s+1)*(s+2))]
m = 2;
p = 1;
cmmmc = (s+1)*(s+2)*(s+3)
r = polynomialDegree(cmmmc)
coef = coeffs(cmmmc)
R = [s+2, s^2+s-6]
rest1 = coeffs(R(1))
rest2 = coeffs(R(2))
R0 = [rest1(1), rest2(1)]
R1 = [rest1(2), rest2(2)]
R2 = [0, rest2(3)]
z = zeros(m);
e = eye(m);
Ac = [z, e, z; z, z, e; -coef(1)*e, -coef(2)*e, -coef(3)*e]
Bc = [z;z;e]
Cc = [R0, R1, R2]
Ac = double(Ac);
Bc = double(Bc);
Cc = double(Cc);
zo = zeros(p);
eo = eye(p)
A0 = [0, 0, -coef(1)*eo; eo ,zo, -coef(2)*eo; zo, eo, -coef(3)*eo]
B0 = [R0; R1; R2]
C0 = [zo, zo, eo]
A0 = double(A0);
B0 = double(B0);
C0 = double(C0);

syms lambda;
pol_car = det(lambda*eye(3) - A0)
coefic = sym2poly(pol_car)
ok = 0;
for k = 1:size(coefic)
    if coefic(k) < 0
        disp('Nu se poate aplica hurwitz')
        ok = 1;
        break;
    end
end
if ok == 0
    [H, delta] = hurwitz(coefic);
end
ok = 0;
for k = 1:size(delta)
    if delta(k) < 0
        disp('Nu e intern stabil')
        ok = 1;
        break;
    end
end
if ok == 0
    disp('Sistemul e intern stabil')
end
syms a;
c = real(roots(sym2poly(cmmmc)))
ok = 0;
for k = 1:size(c)
    if c(k) > 0
        disp('Nu e extern stabil')
        ok = 1;
        break;
    end
end
if ok == 0
    disp('Sistemul e extern stabil')
end

