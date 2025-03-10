function out = HamiltonODE(t,xp,dynfun_handle,Jacobi_handle,diffun_handle)

%H(x,p) = dot(b(x),p) + (p, a(x)p)./2;

%ODE23tb is recommended 

sxp=size(xp);
if (sxp(2) >1 ) 
    xp =xp';
end 

n =length(xp);  
x= xp(1:2:n-1); 
p = xp(2:2:n);
Hp = dynfun_handle(0,x)+ diffun_handle(x) *p;
Hx =  Jacobi_handle(x)' * p  ;
out=[Hp;-Hx];
out=reshape(out,[n/2 2])';
out=out(:);
% 
% 
% function out=DHamiltonODE(t,xp)
% n=length(xp);
% x=xp(1:2:end-1);
% p = xp(2:2:n);
% Hqq = cos(q);
% Hqp = Jacobi_handle(x);
% Hpp = diffun_handle(x);
% Hpq = Hqp;
% out=[ (Hqp)  (Hpp); - (Hqq) - (Hqp)];