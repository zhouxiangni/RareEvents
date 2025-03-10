
function [X, T] = LimitCycleShooting (x0,T0,dynfun_handle,Jacobi_handle, RelTol)
% find the limit cycle by Newton's shooting method
% H(x,T):=\phi(x,T) -x 
% where \phi(x,T) is the solution of ODE \dot{x}=b(x) at time T with initial x
% with the constraint  b(x) is perpendicular to  the dispalcement dx

%\Phi(x,t) = derivative_x  \phi(x,t) satisyfing matrix-valued ODE  

%Practical Numerical algorithms for Chaotic Systems" by T.S.
% Parker and L.O. Chua 


% x0 is row-vector 
% T0 is scalar 

dim = length(x0);
sz = size(x0);
if sz(1) > 1
    x0=x0';
end 

%dim = 2; %dimension of the dynamical system

v = [x0, T0]; 
dv = v'; 
b = zeros(dim+1,1) ; 


err = 1; 

 
IteMatrix = zeros(dim+1,dim+1); %Newton-method interation matrix 

tmp = zeros(1,dim*(dim+1)); % temperal use for transformation between matrix and vector 
 
while (err > RelTol)
  % transform matrix into vector 
    
    % integrate until  time v(n+1), i.e. T
    
    tmp(1:dim+dim*dim) = SolVar(v(1:dim), v(dim+1),dim,dynfun_handle,Jacobi_handle); %row vector
    
    % construct right-hand-side for the Newton scheme
    b(1:dim) = v(1:dim)-tmp(1:dim);    b(dim+1) = 0 ; 
    
    % back to matrix form 
    for i = 1: dim;
      for j = 1: dim; 
          IteMatrix(i,j) = tmp(dim+(i-1)*dim+j); 
      end
    end   
    
    
    IteMatrix(1:dim,1:dim) = IteMatrix(1:dim,1:dim) -eye(dim,dim);
    IteMatrix(1:dim, dim+1) =  dynfun_handle(0,tmp(1:dim)') ; %column vector
   
    IteMatrix(dim+1, 1:dim) = (dynfun_handle(0,  v(1:dim)') )' ;  %row vector  
    IteMatrix(dim+1,dim+1) = 0 ; 
   
   
   
    dv = IteMatrix\b ;    %coln vector
    v = v + dv' ;
    
    T=v(dim+1)
    err = norm(dv,inf)/norm(v,inf)
% pause
end 

X=v(1:dim);
T=v(dim+1)


 
options = odeset('RelTol', 1.e-5, 'AbsTol',1.e-8);
[t,y] = ode45('dynfun', [0,T], X, options);
figure(2); hold off;
if (dim==2) 
 plot(y(:,1),y(:,2),'-o');grid on ; title('LC', 'FontSize',18)
end  
 if (dim>=3)
  plot3(y(:,1),y(:,2),y(:,3),'-o');grid on ;title('LC', 'FontSize',18)
 end 
 

figure(1); hold off;
title('LC', 'FontSize',18); xlabel('time'); ylabel('x, y')
plot(t,y); grid on;


end 

 

