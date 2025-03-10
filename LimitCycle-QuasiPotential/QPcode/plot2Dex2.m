clear all 

global dim 


global case_id 


%dim of the dyn sys
dim = 2;
d=dim;


% compute which example
% refer to the functions 'dynfun' and 'Jacoi'

case_id = 22 ;

 
T =     5.796622979461920;
X = [-0.274756695273899   2.338574878703913 ] ;
 

dynfun_handle = @(t,x) dynfun(t,x);
Jacobi_handle = @(x) Jacobi(x);
diffun_handle = @(x) diffun(x);


[X, T] = LimitCycleShooting (X,T,dynfun_handle,Jacobi_handle, 1e-14)
sprintf('The Limit Cycle  is Found. Period = %1.10g',T) 

options = odeset('RelTol', 1.e-13, 'AbsTol',1.e-15);
nt = 200;
tmesh = linspace(0,T,nt+1);
[tmesh,LC1] = ode45('dynfun', tmesh, X, options);

T =       4.448614936286606;
X = [ 0.013962075678092  -2.998457193321529 ] ;


[X, T] = LimitCycleShooting (X,T,dynfun_handle,Jacobi_handle, 1e-14)
sprintf('The Limit Cycle  is Found. Period = %1.10g',T) 

tmesh = linspace(0,T,nt+1); 
[tmesh,LC2] = ode45('dynfun', tmesh, X, options);


ns = 40;
LC1= MakeArcLengthPara(tmesh, LC1, ns);
LC2= MakeArcLengthPara(tmesh, LC2, ns);



figure(5); hold off;
plot(LC1(:,1), LC1(:,2),'k','LineWidth',3); hold on
plot(LC2(:,1), LC2(:,2),'k','LineWidth',3); hold on
title('The example with two limit cycles', 'FontSize',20);
xlabel('x');    ylabel('y');
axis([-4 2 -4 5]);
set(gca,'fontsize',20); box on ;

xsaddle = [-0.9; 0.694166171569431];
plot(xsaddle(1), xsaddle(2),'x', 'MarkerSize',20,'LineWidth',6)
J = Jacobi_handle(xsaddle);

[ V, D]=eig(J); eps =0.0001;
tmesh = linspace(0,13,200);
[tmesh,traj] = ode45('dynfun', tmesh, xsaddle + eps.* V(:,1), options);
traj= MakeArcLengthPara(tmesh, traj, 100);
plot(traj(:,1), traj(:,2),'r','LineWidth',2); hold on
tmesh = linspace(0,13,200);
[tmesh,traj] = ode45('dynfun', tmesh, xsaddle - eps.* V(:,1), options);
traj= MakeArcLengthPara(tmesh, traj, 100);
plot(traj(:,1), traj(:,2),'r','LineWidth',2); hold on

case_id = 221 ; 
eps =0.00001;
tmesh = linspace(0,20,200);
[tmesh,traj] = ode45('dynfun', tmesh, xsaddle + eps.* V(:,2), options);
traj= MakeArcLengthPara(tmesh, traj, 100);
plot(traj(:,1), traj(:,2),'g','LineWidth',2); hold on
tmesh = linspace(0,20,200);
[tmesh,traj] = ode45('dynfun', tmesh, xsaddle - eps.* V(:,2), options);
traj= MakeArcLengthPara(tmesh, traj, 100);
plot(traj(:,1), traj(:,2),'g','LineWidth',2); hold on
case_id =22;


