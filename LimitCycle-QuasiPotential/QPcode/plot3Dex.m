
arc = zeros(1,nt+1); arc(1) = 0 ;  
for it = 2: nt+1
    arc(it) = norm( LC(it,:)-LC(it-1,:));
end
arc = cumsum(arc);
arcLength = arc(nt+1);
arc = arc./arcLength; 


ns = 160;
LC_arc = zeros(ns+1,d);
% interpolate LC on strictly uniform mesh with new size 'ns'
for id = 1: d  % parallel in each dimensiona
    LC_arc(1:ns+1,id) = interp1(arc,LC(:,id),linspace(0,1,ns+1),'spline');
end

for i = 1: d 
    for j = 1: d
         E_arc(1:ns+1,i,j) = interp1(arc,E(:,i,j),linspace(0,1,ns+1),'spline');
    end 
end

for i = 1:2
        G_eig_arc(1:ns+1,i ) = interp1(arc,G_eig(:,i),linspace(0,1,ns+1),'spline');
for j = 1 :2 
        G_eigV_arc(1:ns+1,i,j) = interp1(arc,G_eigV(:,i,j),linspace(0,1,ns+1),'spline');
end 
end 


figure(5); hold off;
plot3(LC_arc(:,1), LC_arc(:,2),LC_arc(:,3),'k','LineWidth',3); hold on
title('The 3D example', 'FontSize',20);
xlabel('x_1');    ylabel('x_2'); zlabel('x_3');
 set(gca,'fontsize',20); box on ; grid on ;

ntheta=30;
theta=linspace(0,2*pi,ntheta+1);

tmpm = zeros(3,2);
vec = zeros(3,1);
tmp = zeros(2,1);
X = zeros(ns, ntheta);
Y = zeros(size(X));
Z = zeros(size(X));
delta = 0.00002;

uniform_arc=linspace(0,1,ns+1);
for i =1 :ns+1
for j = 1:ntheta+1
        tmpm(1:3,1:2) = E_arc(i,1:3,2:3);
      
        tmp(1:2,1) = cos(theta(j)) * sqrt(2*delta/G_eig_arc(i,1)).* G_eigV_arc(i,1:2,1) ...
            +        sin(theta(j)) * sqrt(2*delta/G_eig_arc(i,2)).* G_eigV_arc(i,1:2,2); 
        
%         tmp(1:2,1) = cos(theta(j)) * sqrt(2*delta/G_eig_arc(1)).* G_eigV_arc(i,1:2,1) ...
%             +        sin(theta(j)) * sqrt(2*delta/G_eig_arc(1)).* G_eigV_arc(i,1:2,2); 
%        tmp(1,1) = cos(theta(j)) * sqrt(2*delta);
%        tmp(2,1) = sin(theta(j)) * sqrt(2*delta); 
        
        vec(1:d,1) = LC_arc(i,1:d)' + tmpm* tmp;  
         
                
     X(i,j) = vec(1);
    Y(i,j) = vec(2);
    Z(i,j) = vec(3);
    C(i,j) = log(G_eig_arc(i,1) / G_eig_arc(i,2));
end 
end 

hold on;
plot3(1,1,1,'o','MarkerSize',15,'LineWidth',4);
plot3(0,0,0,'x','MarkerSize',20,'LineWidth',4);
 mesh(X,Y,Z,C); colormap winter;grid on ;
    
%% quiver3



[X, Y] =meshgrid(linspace(0.5,2.5,3),linspace(0.5,1.5,3));
zz = linspace(0,1.8,4);

Z=zeros(3,3);
U=zeros(size(Z));
V=zeros(size(U));
W=zeros(size(U));


case_id = 3;
for k = 1: 4
for i = 1: 3
for j = 1:3
        Z(i,j)=zz(k);
        tmp = dynfun(1,[X(i,j) ;Y(i,j);  Z(i,j)]);
         U(i,j)=tmp(1)./norm(tmp);
         V(i,j)=tmp(2)./norm(tmp);
         W(i,j)=tmp(3)./norm(tmp);
         
end
end
% quiver3(X,Y,Z,U,V,W,0.3,'k','MaxHeadSize',0.3,'LineWidth',2); hold on 
end 

tmp = linspace(0,5,200);
eps = 1.0;
[tmp,traj] = ode45('dynfun', tmp,   eps.*[ 1 2 1  ], options);
traj= MakeArcLengthPara(tmp, traj, 100);
plot3(traj(:,1), traj(:,2), traj(:,3), 'r','LineWidth',2); hold on
view([68 35])


