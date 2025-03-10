function J = Jacobi_PRDE1D(t,G,T,tmesh,m,a)
%all are scalars
%input:  t, G   
%output: dG/dt
% parameters:  T -period 
%              tmesh:  tmesh(1)=0, tmesh(last) = T;
%              m and a :  vectors with same size of tmesh
t = mode(t,T);
m = interp1(tmesh,m,t); % Interpolate the data set (ft,f) at time t
a = interp1(tmesh,a,t); % Interpolate the data set (gt,g) at time t
J = -2.*m -2.*a.*G; % Evaluate ODE at time t
    