function  LC_arc = MakeArcLengthPara(tmesh, LC, ns)
% 

sLC = size(LC);
n=sLC(1);
d=sLC(2);

if length(tmesh) ~= n
    error('length does not match in MakeArcLengthPara')
end 


LC_arc = zeros(ns,d);

arc = zeros(1,n); arc(1) = 0 ;  
for it = 2: n
    arc(it) = norm( LC(it,:)-LC(it-1,:));
end
arc = cumsum(arc);
arcLength = arc(n);
arc = arc./arcLength; 


% interpolate LC on strictly uniform mesh with new size 'ns'
for id = 1: d  % parallel in each dimension
    PP_spline = spline(arc,LC(:,id));
    LC_arc(1:ns,id) = ppval(PP_spline,linspace(0,1,ns));
end


