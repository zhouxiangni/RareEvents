function [h11, h12, h21, h22] = CalHessian (x, y)

x0 = [1 0 -0.5 -1];
y0 = [0 0.5 1.5 1];

a = [-1 -1 -6.5 0.7];
b = [0 0 11 0.6];
c = [-10 -10 -6.5 0.7];
A = [-200 -100 -170 15];

h11 = 0; h12 = 0; h21 = 0; h22 = 0;

for k = 1:4
  a1(k) = a(k)*(x-x0(k))*(x-x0(k));
  b1(k) = b(k)*(x-x0(k))*(y-y0(k));
  c1(k) = c(k)*(y-y0(k))*(y-y0(k));
  V(k) = A(k)*exp(a1(k) + b1(k) + c1(k));

  Var1 = 2*a(k)*(x-x0(k)) + b(k)*(y-y0(k));
  Var2 = 2*c(k)*(y-y0(k)) + b(k)*(x-x0(k));
  Var3 = V(k);

  h11 = h11 + 2*Var3*a1(k) + Var3*Var1*Var1;
  h12 = h12 +   Var3*b1(k) + Var3*Var1*Var2;
  h21 = h21 +   Var3*b1(k) + Var3*Var1*Var2;
  h22 = h22 + 2*Var3*c1(k) + Var3*Var2*Var2;
end
