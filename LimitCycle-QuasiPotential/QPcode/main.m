clear all 
global dim 


 global case_id 
 global diff_case_id


%dim of the dyn sys
dim = 5;
d=dim;


% compute which example
% refer to the functions 'dynfun' and 'Jacoi'

case_id = 22 ;

diff_case_id = 1;

if  dim >  2
    diff_case_id = 1; 
    case_id = dim;
end 

%% FIND LC
switch case_id
    case 2  %Van del Pol
        T =   6.663286859323122; 
        X = [2.006332678322606   0.100501421512104];
    case 22 % two 2D limit cycle;
        T =  5.796622979461945;
        X = [-0.274756695273898   2.338574878703912 ] ;
    case 3   % 3D case Lokta-Volterra
        T = 8.353353138196825; %1.49
        X = [2.708602981596172   0.230781048538252   0.968203002422098];
          T=         7.9782;%1.486
         X=[   2.6680    0.2628    0.9289];
         T=          6.7965;%1.48
         X=[   2.4380    0.4260    0.7475];
    case 4
        X=[  0.616379917269709   2.018202595822800   0.918988109831200   1.328148584084487];
        T= 4.665719993;
    case 50  % abanded do not use !
        T=  40.965603257866853;
        X= [0.277019130840132   0.325439424442765   0.168193371137287 ...
            0.391929893431921   0.189676992242229];
    case 5   % 5D example
        T = 8.116492229;
        X=[  -0.238681653608981   0.740973231274708  -0.595576139666053 ...
            0.296597987646723  -0.838180009885698 ];
        X=[   0.297649663880028   0.246207302670090  -0.642191694312627  ...
            0.61133696767963     -0.911663731284832];


        
end

 


dynfun_handle = @(t,x) dynfun(t,x);
Jacobi_handle = @(x) Jacobi(x);
diffun_handle = @(x) diffun(x);

 


    [X, T] = LimitCycleShooting (X,T,dynfun_handle,Jacobi_handle, 1e-12);
   sprintf('The Limit Cycle  is Found. Period = %1.10g',T) 
   X
  
 
options = odeset('RelTol', 1.e-13, 'AbsTol',1.e-15);
nt = 20000;

tmesh = linspace(0,T,nt+1);

[tmesh,LC] = ode45('dynfun', tmesh, X, options);


figure(1)
plot(tmesh,LC); grid on 


%%

% # of images in output 
if (dim <= 3) 
    ns = 20000;
else 
    ns = 2000; 
end 


if d >= 5
%     
%     tmesh(nt+1:2*nt+1) = tmesh(nt+1) + tmesh(1:nt+1);
%     LC(nt+2:2*nt+1,1:d) = LC(2:nt+1,1:d);


    [tmesh, LC, E, Edt, Einv]= FindE2(tmesh,LC,ns, ... 
                        dynfun_handle,Jacobi_handle,diffun_handle);
else
    [tmesh, LC, E, Edt, Einv]= FindE(tmesh,LC,ns);
end


  
[M,A,Omega] = construct_MA(tmesh,LC,E,Edt,Einv);
  
nt = ns;

 
 
disp(' The basis and coefficients matrices in PRDE are constructed .')
 

  
%%  
disp(' ... START TO SOLVE PRDE ... ')

     

%% find initial guess of interval of G0 at d = 2 
if (d == 2)
    M=reshape(M,[],1);
    Check_M =  dot(M(1:nt), diff(tmesh(1:nt+1)));
    if (Check_M >0 )
        sprintf('The stability of LC has problem: int M(s)ds= %0.5g <0', Check_M)
        return
    end
    
    
    
      
    

% %       %bisection search -- 
%     G0tmp=linspace(0.01,0.08,10); tmp=zeros(size(G0tmp));
%     for i = 1:10
%         G0 = G0tmp(i);  
%         options = odeset('RelTol', 1.e-3, 'AbsTol',1.e-5);
%         [t,y] = ode45(@(t,G) PRDE(t,G,T,tmesh,M,A), [0,T], G0, options);
%         tmp(i)=y(length(t))-G0;
%         i
%     end
%     figure(1); hold off; title('GT-G0 vs G0', 'FontSize', 22);
%     plot(G0tmp,tmp); grid on ; 
%     disp('Check the Figure 1 and select the region containg zero')
%     
    if (case_id==2)
        switch diff_case_id 
            case 1
                G0guess = [1.20 1.90]; 
                G0guess = 1.565;% Van del Pol + a==1
            case 2
                G0guess = [1 3]; 
                %G0guess=2.042
            case 3
                G0guess = [5 8]; 
                %G0guess=  7.;
        end 
    end
    if (case_id == 22)
        if  (diff_case_id  == 1)
          G0guess = 0.30133; 
         end 
    end 
    
    
    G0guess = mean(G0guess);
    options = odeset('RelTol', 1.e-5, 'AbsTol',1.e-7);
    rel_err = 1;
    figure(1); hold off; title('Trajectory of PRDE', 'FontSize', 20);
    for i = 1: 50
        if (rel_err < 1e-6) && (i > 2)
            break;
        end 
        if (rel_err> 100)
            break;
        end 
        [tmp,ytmp] = ode45(@(t,G) PRDE(t,G,T,tmesh,M,A), [0,T], G0guess, options);
        rel_err = abs(G0guess- ytmp(length(tmp))) / abs(G0guess);
        G0guess = ytmp(length(tmp));
        plot(tmp,ytmp); hold on
        i 
        rel_err
    end
    if ( i == 50)  
        disp('WARNING:  check the Figure 1 and See if PRDE converges to a stable periodic solution')
         pause;
    else 
         sprintf(' * The Guess of G0 for PRDE is found = %0.16g', G0guess)
    end 
    
    
     % solve the nonlinear equation of G0 satisfies
%      optfzero = optimset('TolX',1e-20);
%     G0 = fzero(@(G0) Int_PRDE1D_Period(G0,T, tmesh, M,A), G0guess,optfzero);
    % check options here 
    G0 = G0guess;
    sprintf(' *** The initial of G for PRDE is found = %0.16g', G0)
      
    options = odeset('RelTol', 1.e-6, 'AbsTol',1.e-12);
    [tmesh,G] = ode45(@(t,G) PRDE(t,G,T,tmesh,M,A), tmesh, G0, options);
%     abs(G(1)-G(nt+1))/abs(G(1))
    
     if (case_id == 22 )
         b_vec=zeros(nt+1,1);
        for it = 1: nt+1 
            b_vec(it) = norm( dynfun_handle(0,LC(it,1:d) ) ); 
        end 
         figure(3); subplot(2,1,1);
         plot(tmesh, LC, tmesh, b_vec,'LineWidth',3);
         xlabel('time', 'FontSize',24); ylabel('limit cycle','FontSize',24)
          set(gca,'fontsize',20); box on ;
          subplot(2,1,2);
         plot(tmesh, 1./G, '-','LineWidth',3);  
     xlabel('time', 'FontSize',24); ylabel('$G$','FontSize',24);

        set(gca,'fontsize',20); box on ;
     end 


    %% draw contour plot 
    LC_lb = zeros(nt+1,d);    LC_ub = zeros(nt+1,d);
    
    figure(5); hold off;
    plot(LC(:,1), LC(:,2),'k','LineWidth',6);
    
    
    switch case_id 
        case 2
            delta = 0.02; % the contoure of z:  0.5*G*z*z = delta 
        case 22 
            delta = 2e-5;
    end     
    % the contoure of z:  0.5*G*z*z = delta 
    for it = 1: nt+1 
        LC_lb(it,1:d) = LC(it,1:d) - E(it,1:d,2).* sqrt(2*delta./G(it)) ;
        LC_ub(it,1:d) = LC(it,1:d) + E(it,1:d,2).* sqrt(2*delta./G(it)) ;
    end 
    
    figure(5); hold on;
    plot(LC_lb(:,1), LC_lb(:,2),':b','LineWidth',3); hold on 
    plot(LC_ub(:,1), LC_ub(:,2),':b','LineWidth',3);
    
    switch case_id 
        case 2
            delta = 0.1; % the contoure of z:  0.5*G*z*z = delta 
        case 22 
            delta = 1e-4;
    end  
    for it = 1: nt+1 
        LC_lb(it,1:d) = LC(it,1:d) - E(it,1:d,2).* sqrt(2*delta./G(it)) ;
        LC_ub(it,1:d) = LC(it,1:d) + E(it,1:d,2).* sqrt(2*delta./G(it)) ;
    end
    %inner
    plot(LC_lb(:,1), LC_lb(:,2),':r','LineWidth',3); hold on
    %outside
    plot(LC_ub(:,1), LC_ub(:,2),':r','LineWidth',3);
    
    title('The Contour of Quasi-potential of Van Del Pol example', 'FontSize',20);
    xlabel('x');    ylabel('y');
    axis([-3 3 -4 4]);
    set(gca,'fontsize',20); box on ;
    
    shootindx =  [1:floor(nt/2.8):nt];
    delta_indx = 1e-5;
    XP0=zeros(2*d,1);
    
   if case_id == 2 
    for inout = [ -1]
    for kk =  1: length(shootindx) 
        it =          shootindx(kk)
        delta = delta_indx(1);
        z =   inout.*sqrt(2*delta./G(it));
        
        XP0(1:2:2*d-1) = LC(it,1:d) +  E(it,1:d,2).* z;
        XP0(2:2:2*d) = getMomentum(z,it,nt,d,tmesh,LC, Omega,E, Einv, G) ;
        
        options = odeset('Events',@Exitevents,'RelTol', 1.e-12, 'AbsTol',1.e-18);
        fhandle = @(t,xp) HamiltonODE(t,xp,dynfun_handle,Jacobi_handle,diffun_handle);
        [t, y] = ode45(fhandle, [0, 100],XP0,options);
        x= y(:,1:2:2*d-1); p = y(:,2:2:d*2); 
        t(end)
        Ha=zeros(length(t),1);
        for i=1:length(t);
            Ha(i) = dot(dynfun_handle(0,x(i,:)), p(i,:)) + ...
                       .5.*dot(p(i,:), diffun_handle(x(i,:))* (p(i,:)') ) ;
        end 
        figure(1); plot(t,Ha)
        if ( max(abs(Ha))> 1.e-1) 
            error('ERROR: too large H >1e-4');
        end 
         figure(2); hold on;  plot(x(:,1),x(:,2),'LineWidth',2);   
%         XP0(2:2:2*d) = XP0(2:2:2*d) ./ norm(XP0(2:2:2*d));
%         quiver(XP0(1),XP0(3),XP0(2),XP0(4), 0.3);
        for it = 1:length(t);
            Ha(it) = dot(p(it), diffun_handle(x(it))*p(it))./2 ;
        end     
        ntmp = max(20000, length(t));
        int_mesh=linspace(0,t(end),ntmp+1);
        Ha = interp1(t,Ha,int_mesh,'spline') .* (t(end)./ntmp);
        x = interp1(t,x,int_mesh,'spline');
        x = [LC(it,1) LC(it,2) ;x]; 
        Ha = [0  Ha]; Ha = cumsum(Ha);%action
        
%         nx = 50; ny= 100;
%         QPGridX=linspace(0,3,nx);  QPGridY=linspace(-3,3,ny);
%         nx = length(QPGridX);ny = length(QPGridY);
%         dx = QPGridX(2)-QPGridX(1);  dy = QPGridY(2)-QPGridY(1);
%         QPGrid=ones(nx,ny).*(900000);
%         for i = 1: nx
%         for j = 1: ny
%             wh =[ QPGridX(i)-dx/2  QPGridX(i)+dx/2 ...
%                    QPGridY(j)-dy/2  QPGridY(j)+dy/2   ];
%             INDX = ( x(:,1) > wh(1) ) & (x(:,1) < wh(2) ) & ...
%                    ( x(:,2) > wh(3) ) & (x(:,2) < wh(4) );
%              if  ~ isempty(Ha(INDX)) 
%                 QPGrid(i,j) = min( QPGrid(i,j), mean(Ha(INDX)));
%              end 
%         end 
%         end 
%         [GX, GY] = meshgrid(QPGridX,QPGridY);
%         QPGrid (  QPGrid  >= 100000) = NaN;
%         figure(200); surf(mesh(GX,GY,QPGrid));
    end 
    end 
       
   end 
    
    
   
    
   
    
    
    
elseif (d>2) 
    
  
        %% high dim cases
        
    % long time integral to find good initial guess
    
    if (case_id ==3 ) 
          G0guess = [  11.253284528769219  -0.939026434888983; ...
                       -0.939026434888976   0.078444754755502 ];
%           G0guess =[   ];
    end 
    
  
    
    if case_id ==5
        G0guess = [  1.271547034945875   0.716853530816152  -0.598924434501274   0.042011086691818; ...
   0.716853530816145   3.637677632992334  -3.395810451456156   0.190044588297110; ...
  -0.598924434501266  -3.395810451456181   3.864016936637070   0.400708857420006; ...
   0.042011086691814   0.190044588297112   0.400708857419989   2.297208983165200 ];
    end 
    
    % start with a large initial positive definite matrix 
    % in col vec form
     rel_err = 1;
     Gnew=zeros(d-1,d-1);
 
     
    figure(1); hold off; title('Trajectory of PRDE', 'FontSize', 20);
     M2 = M;  A2 = A;  FF = eye(d-1,d-1); FF(3,3)=-1; F(4,4)=-1; 
        for it = 1: nt+1
            M2(it,1:2,3:4) = -M(it,1:2,3:4);
            M2(it,3:4,1:2) = -M(it,3:4,1:2);
            A2(it,1:2,3:4) = -A(it,1:2,3:4);
            A2(it,3:4,1:2) = -A(it,3:4,1:2);
        end 
        
        
    anti_period = 0;
    if (anti_period == 0)
    for i = 1: 50
        if (rel_err < 3.e-4) &&  (i > 3)
            break;
        end 
        options = odeset('RelTol', max(1e-6, 1.e-2/(2.^i)), ...
                         'AbsTol', max(1e-12,1.e-4/(2.^i)) );
        [tmp,ytmp] = ode45(@(t,G) PRDE(t,G,T,tmesh,M,A),  [0:T/200:T], ...
                            reshape(G0guess,[],1), options);
        Gnew = reshape( ytmp(length(tmp),:), d-1,d-1);
        
         plot(tmp, ytmp(:,3),'x'); hold on

        [tmp,ytmp] = ode45(@(t,G) PRDE(t,G,T,tmesh,M2,A2), [0:T/200:T], ...
                            reshape(Gnew,[],1), options);
        Gnew = reshape( ytmp(length(tmp),:), d-1,d-1);
        rel_err = norm(G0guess- Gnew) / norm(G0guess);
        G0guess = Gnew;
        
        plot(tmp+T, ytmp(:,3),'o'); hold on
        i
        rel_err
        if (rel_err>10) 
%             break
        end 
    end
    end 
    
    if (anti_period == 1 ) 
       
        
        for i = 1: 1
        if (rel_err < 2.e-4) &&  (i > 3)
            break;
        end 
        options = odeset('RelTol', max(1e-6, 1.e-2/(2.^i)), ...
                         'AbsTol', max(1e-12,1.e-4/(2.^i)) );
        [tmp,ytmp] = ode45(@(t,G) PRDE(t,G,T,tmesh,M,A), [0:T/200:T], ...
                            reshape(G0guess,[],1), options);
        Gnew = reshape( ytmp(length(tmp),:), d-1,d-1);
        Gnew = F*Gnew*F;
        
        plot(tmp, ytmp(:,2:3)); hold on

        
        [tmp,ytmp] = ode45(@(t,G) PRDE(t,G,T,tmesh,M2,A2), [0:T/200:T], ...
                            reshape(Gnew,[],1), options);
        Gnew = reshape( ytmp(length(tmp),:), d-1,d-1);
        Gnew = F*Gnew*F;
        
        
        rel_err = norm(G0guess-Gnew ) / norm(G0guess);
        G0guess = Gnew;
        
        
        plot(T+tmp, -ytmp(:,3)); hold on
        
        i
        rel_err
        if (rel_err>10) 
%             break
        end 
    end
    end 
    
    if ( i == 50)  || rel_err > 10
        disp('WARNING:  check the Figure 1 and See if PRDE converges to a stable periodic solution')
         pause;
    else 
         sprintf(' * The Guess of G0 for PRDE is found ');
         disp(G0guess)
        % check symmegry 
        for j = 1:d-1
            for k = j+1: d-1
                if (j ~= k )
                    if abs(G0guess(j,k)-G0guess(k,j))/ abs(G0guess(j,k)) > 1e-6
                        disp('WARNING: the initial of G found is not symmetric !!' ); pause
                    end 
                end 
            end 
        end 
        if  min(eig(G0guess)) <= 0 
            disp('WARNING: the initial of G found is not positive definite !!' ); pause
        end 
    end 
      
 
        
    G_eig =zeros(nt+1,1);
    G_eigV=zeros(nt+1,d-1,d-1);
    G=zeros(nt+1,d-1,d-1);
    
    G0=reshape(G0guess,[],1);
    options = odeset('RelTol', 1.e-6, 'AbsTol',1.e-12);
    ytmp=[];
    [tmesh,ytmp] = ode45(@(t,G) PRDE(t,G,T,tmesh,M,A), tmesh, G0, options);
     ytmp(nt+1,:)=ytmp(1,:);
     
    for it =1 : nt 
        G(it,1:d-1,1:d-1) = reshape(ytmp(it,:),d-1,d-1);
        for j = 1:d-1
            for k = j+1: d-1
                if abs(G(it,j,k)-G(it,k,j))/ abs(G(it,j,k)) > 1e-6
                    disp('WARNING: the solutionf of PRDE G is not symmetric !!' ); pause
                end
                tmp = G(it,j,k)./2  + G(it,k,j)./2;
                 G(it,j,k) = tmp;  G(it,k,j)=tmp;          % strictly symmetry
            end
        end
        tmp=zeros(d-1,d-1); tmp(1:d-1,1:d-1)=G(it,1:d-1,1:d-1);
        [V, D] = eig(tmp); D=diag(D);
        [D, I]=sort(D,'descend');
        V(1:d-1,1:d-1) = V(1:d-1,I);
        
        G_eigV(it,1:d-1,1:d-1)=V;
        G_eig(it,1:d-1) =  D;
        if  min( G_eig(it,1:d-1)) <= 0 
            disp('WARNING: the initial of G found is not positive definite !!' ); pause
        end 
    end 
    
    
    
   
    figure(10);hold off;
%     subplot(2,1,1);
%     b_vec=zeros(nt+1,1);
%     for it = 1: nt+1
%         b_vec(it) = norm( dynfun_handle(0,LC(it,1:d) ) );
%     end
%     plot(tmesh, b_vec,'LineWidth',3);
%     xlabel('time', 'FontSize',24); ylabel('limit cycle','FontSize',24)
%     set(gca,'fontsize',20); box on ;
%     subplot(2,1,2);
%    
    hold off; plot(tmesh(1:ns),G_eig(1:ns,:),'.-','LineWidth',4); grid on
    xlabel('time', 'FontSize',24); ylabel('$G$','FontSize',24);
    set(gca,'fontsize',20); box on ;
    title('Eigenvalues of G', 'FontSize',20);
    
    % fix the sign of the eigenvectors for continuity in time
    tmp=zeros(d-1,1); tmp1=zeros(d-1,1);
    for it = 2:nt
        for ic =1:d-1
            tmp = G_eigV(it,1:d-1,ic);
            tmp1 = G_eigV(it-1,1:d-1,ic);
            if ( norm(tmp-tmp1)> norm (tmp+tmp1))
                G_eigV(it,1:d-1,ic)= - G_eigV(it,1:d-1,ic);
            end
        end
    end
    
    G(nt+1,:,:)=G(1,:,:); 
    G_eig(nt+1,:)= G_eig(1,:); 
    G_eigV(nt+1,:,:)=G_eigV(1,:,:);
    
    
    
    
figure(11)
 % uniform_arc(1:ns) matches E(1:ns)
 for ic = 1: d-1
     subplot(d-1,1,ic); 
    hold off;
     for id = 1:d-1
        plot(tmesh(1:nt),G_eigV(1:nt,id,ic),'o-','MarkerSize',6); hold on ;
     end  
     title (' eigenvectors of G','FontSize',18)
 end 
 
 
 
 
 
 end
    
    
    
        

