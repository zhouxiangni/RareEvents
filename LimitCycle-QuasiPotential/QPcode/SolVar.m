

function output = SolVar(input, T,dim,dynfun_handle,Jacobi_handle )

%input are row vector with size n 
%output are row vector with size n+n*n

% solve the variational queation from time [0,T] using ode45 subroutine 
% initial value of x is input : row vector, dim=n 

init = zeros(1,dim*(dim+1));

init(1:dim)= input(1:dim);
%initial vaule for variational is identity matrix 
for i = 1: dim ;
    init(dim+(i-1)*dim+i)= 1; 
end
 


options = odeset('RelTol', 1.e-14, 'AbsTol',1.e-22);
[~,y] = ode45(@(t,input) EqnVar(t, input,dim, dynfun_handle,Jacobi_handle), ...
    [0,T/2,T], init, options); % T/2 just technical for matlab 
output = y(3,:); %  row vector




