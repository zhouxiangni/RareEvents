
function dGdt = TransitJacobi(t,G,Period,tmesh,LC,Jacobi)
% Model the  Riccati Differential Equation 
% dG/dt = - J*G;


t = mod(t,Period);
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
X=zeros(d,1);
for i = 1: d
%      X(i)=interp1(tmesh,LC(:,i),t,'spline');
     X(i) = ppval(LC{i},t);
end 
dGdt =  Jacobi(X)*Gm;
dGdt = reshape(dGdt,[],1);

