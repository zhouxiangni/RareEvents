function output = EqnVar(t, input, dim, dynfun_handle,Jacobi_handle)

 


n = dim;

% Define the variational equation with original dynamical system 
% see P306 "Practical Numerical algorithms for Chaotic Systems" by T.S.
% Parker and L.O. Chua 

% x' = f(x,t)
% \Phi' = D_x f(x,t)\Phi

%initial condition : x(0)=x_0
%\Phi(0) = Identity matrix 

%n  is the dimension of x 
%input(1:n)  for x
%input(n+1:n*(n+1))  for \Phi 
%input(n+(i-1)*n+j)  for  \Phi(i,j), with i=1,,n; j= 1..n 

output = zeros(n*(n+1),1);

% original dynamics 
output(1:n) = dynfun_handle(1,input(1:n));



% dynamics for Phi(x) 
Phi   = zeros(n,n);
for i = 1: n
    for j = 1: n 
        Phi(i,j) = input(n+(i-1)*n+j);
    end
end
 
 

% call outside Jacobi function 
Jacob = Jacobi_handle(input(1:n));

% % right hand side 
Phi = Jacob*Phi;
for i = 1: n
    for j = 1: n 
        output(n+(i-1)*n+j) = Phi(i,j);
    end
end

%output=output;

