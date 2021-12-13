function x = HHLinearLeastSquare(A,b)
a = A(:,1);
[rA,cA] = size(A);
if(a(1,1)>0)
    alpha = -norm(a,2);
else
    alpha = norm(a,2);
end
[row,column] = size(a);
I = eye(row);
e1 = I(:,1);
v = a-alpha*e1;
vt = transpose(v);
H1 = I-2*((v*vt)/(vt*v));
H1A = H1*A;
R = H1A;
c = H1*b;
A1 = H1A(2:end,2:end);
for i = 2:cA
    a = A1(:,1);
    if(a(1)>0)
        alpha = -norm(a,2);
    else
        alpha = norm(a,2);
    end
    [rowi,columni] = size(a);
    I = eye(rowi);
    e1 = I(:,1);
    v = a- alpha*e1;
    vt = transpose(v);
    Hkd = I-2*((v*vt)/(vt*v));
    HkdA = Hkd*A1;
    Ik = eye(i-1);
    [r,c1] = size(Hkd);
    zr = zeros(i-1,c1);
    zl = zeros(r,i-1);
    Hk = [Ik,zr;zl,Hkd];
    c = Hk*c;
    R = Hk*R;
    A1 = HkdA(2:end,2:end);
end
R1 = R(1:cA,:);
c1 = c(1:cA,:);
x = R1\c1;

