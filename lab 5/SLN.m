A = [0 0 0; 1 0 -1; 0 1 2]
[m n]=size(A)
B=[0 -1; 1 0; 0 1]
[q p]=size(B)
C=[1 0 1]
u=[1; -2]
lambda_d=[-1 -1 -1]
[e r]=size(lambda_d)
R=ctrb(A,B)
if (rank(R)<n)
    exit
end
F0=rand(p,n)
g0=rand(p,1)

A0=A+B*F0
b0=B*g0
[l z]=size(A0)
R0=ctrb(A0,b0)
while (rank(R0)~= z)
    F0=rand(p,n)
    g0=rand(p,1)
    A0=A+B*F0
    b0=B*g0
    [l z]=size(A0)
    R0=ctrb(A0,b0)
end
R0_invers=inv(R0)
[j h]= size(R0_invers)
q_t=R0_invers(j,1:h)
polinom=eye(l,z)
while (r>0)
    polinom=polinom*(A0-eye(l,z)*lambda_d(e,r))
    r=r-1;
end
  f_t=-q_t*polinom
   
  F=F0+g0*f_t
  v=[1; -2]
  G=eye(p)