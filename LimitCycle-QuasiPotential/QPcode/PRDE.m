
function dGdt = PRDE(t,G,Period,tmesh,m,a)
% Model the  Riccati Differential Equation 
% dG/dt = - m' * G - G' * m - G' * a * G;
%  G' are applied to enforce the formal symmetry 
% This ODE model is used for Int_PRDE_Period ODE Solver 


%all are scalars
%input:  t, G   
%output: dG/dt
% parameters:  T -period 
%              tmesh:  tmesh(1)=0, tmesh(last) = T;
%              m and a :  vectors with same size of tmesh


t = mod(t,Period);

if (length(G) == 1)
    
    m = interp1(tmesh,m,t,'spline'); % Interpolate the data set (ft,f) at time t
    
    a = interp1(tmesh,a,t,'spline' ); % Interpolate the data set (gt,g) at time t
  
    
    dGdt = -2.*m.*G -a.*G.*G; % Evaluate ODE at time t
    
else
    sG  = size(G);
    if sG(2) > 1
        G=G';
    end
    % G: col vec with dim d*d;
    % m(it,:),a: matrix with size dG * dG
    d= sqrt(length(G));
    if (d*d ~= length(G))
        disp(' ERROR'); return
    end
    
    Gm = reshape(G,d,d);
    Mm = zeros(size(Gm));
    Am = eye(size(Gm));
    for i = 1: d
        for j = 1: d
            Mm(i,j) = interp1(tmesh,m(:,i,j),t,'spline'); % Interpolate the data set (ft,f) at time t
            Am(i,j)= interp1(tmesh,a(:,i,j),t,'spline' ); % Interpolate the data set (gt,g) at time t
        end 
    end
        dGdt = -Mm' * Gm - Gm' * Mm - Gm' * Am * Gm;
        dGdt = reshape(dGdt,[],1);
        
end