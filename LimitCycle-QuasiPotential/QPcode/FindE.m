function  [tmesh, LC, E, Edt, invE] = FindE(tmesh,LC,ns)

%% comment 
%tmesh (input):  col vector with size nt+1 : tmesh(1)=0, tmesh(nt+1)=Period
%      (output): col vector with size ns+1 : tmesh(1)=0, tmesh(nt+1)=Period
%                uniform grided 

%LC:  size nt+1 by d ;  each row is a point in R^d space ; match 'tmesh'
%       
% ns:  the number of discrete points on  output LC 

% output:   E, Edt,    size(1:ns+1,d,d),  with 'ns+1' points on LC
%           M and A,   size(1:ns+1,d-1,d-1). The first row/col deleted already 


%% old testing code; abanded
if 1 == 0
N=20;
T=pi*2; dt = T/(N); a= 0.0;
b = a+T;
x=a+[0:N-1].*dt;  x=x';
xm = (x(2:N)+x(1:N-1))./2;
f= exp(cos(x).*sin(2*x));
df = f.*( -sin(x).*sin(2.*x) + 2.*cos(x).*cos(2.*x));
% figure;
% plot(2:N/2, imag(FT(2:N/2)), 2:N/2, -imag(FT(N:-1:N+2-N/2)),'o');
 
[F, dfdx] = fourierderivative(f,a,b,N);%%%
dfdx2 =  fourierdiff(f,1,[a, b],N);
fdfdx =  ( f(2:N)-f(1:N-1) ) ./ ( x(2:N)-x(1:N-1));
% figure(1);plot(x,f-F);
figure(2);

%  plot(x(1:N-1), (abs(fdfdx-df(1:N-1))))
 plot(x, (abs(dfdx-df)),'s-',x, (abs(dfdx2-df)),'x-' );
 legend('my','out-fun')
end 


%% start the code here ....

global dim

d = dim; %dim of the dynamics

 

E=zeros(ns+1,d,d);
Edt=zeros(ns+1,d,d);
invE=zeros(ns+1,d,d);



%input  nt, and LC, T 
nt = length(tmesh)-1;
%LC = zeros(nt,d);
%definition:  LC(nt+1) := LC(1); 
Period = tmesh(nt+1);
LC(nt+1,:)=LC(1,:);  % make strict the period condition


% check the generacy of LC
% sprintf('SVD of the Limit Cycle is ')
% S=svd(LC)'




%% arclength 's' parametrized Limit Cycle

 LC_arc = zeros(ns+1,d);
arc = zeros(1,nt+1); arc(1) = 0 ;  
for it = 2: nt+1
    arc(it) = norm( LC(it,:)-LC(it-1,:));
end %dist p_nt+1=p_1 to p_n closed loop.  
arc = cumsum(arc);
arcLength = arc(nt+1);
 



%% interpolate LC on strictly uniform mesh with new size 'ns'
%augmented the periodic interval to have two periods
X = zeros(2*nt,1); 
tmesh_extended = extend_period(tmesh(1:nt),nt,Period); 
for id = 1: d  % parallel in each dimension
    X(1:2*nt) = extend_period(LC(1:nt,id), nt, 0);
    PP_spline = spline(tmesh_extended(1:2*nt),X(1:2*nt));
   
    LC(1:ns+1,id) = ppval(PP_spline,linspace(0,Period,ns+1));
    %spline derivative 
    pp2 = fnder(PP_spline);
    Edt(1:ns+1,id,1) = ppval(pp2,linspace(0,Period,ns+1));
end 

LC(ns+2:nt+1,:)=[];

tmesh  = linspace(0,Period,ns+1)'; 
nt = ns; %update LC, tmesh, nt.

%% 
for id = 1:d 
%debug
    E(1:nt,id,2) = diff(LC(1:nt+1,id))./ diff(tmesh(1:nt+1));
    [~,E(1:nt,id,1)] = fourierderivative(LC(1:nt,id),0,Period,nt) ;
end         
E(nt+1,:,:)=E(1,:,:);

b_vec=zeros(nt+1,d);
for it = 1: nt+1 
   b_vec(it,1:d) = dynfun(0,LC(it,1:d) ); 
end 

%% test tangent direction 
%debug
if  ( 1==2) 
id = 1; 
figure(1);
   plot( tmesh(1:nt+1), log(abs(E(1:nt+1,id,2)./b_vec(1:nt+1,id)-1)),...
         tmesh(1:nt+1), log(abs(E(1:nt+1,id,1)./b_vec(1:nt+1,id)-1)), '-or',...
         tmesh(1:nt+1), log(abs(Edt(1:nt+1,id,1)./b_vec(1:nt+1,id)-1)), '-xb');
   legend('finite-diff','fourier-diff','spline-diff')
   title('tangent direction','FontSize',22)
   pause
plot_Fourier(LC(1:nt,id)); title('Fourier mode of LC','FontSize',22);
 pause;  
plot_Fourier(E(1:nt,id,1)); title('Fourier mode of tangent','FontSize',22)
 pause; close(gcf)
end
%conclusion :  spline derivative is good.
    
%% Build The First Basis Vector 


tmpv = zeros(d,1);  tmpv1 = zeros(d,1);
norm_monitor = zeros(ns,1);

%for id=1:d
    %E(1:nt,id,1) = Edt(1:nt,id,1);%spline derivative  or Fourier are both OK
%end 
E(1:nt,1:d,1) = b_vec(1:nt,1:d);  %or simply use the b_vec

% first basis vector (normalizes)
for it = 1:nt
    tmpv(1:d) = E(it,1:d,1); E(it,1:d,1) = E(it,1:d,1) ./ norm(tmpv);
end 

E(nt+1,:,:)=E(1,:,:);

%% Build Other  Basis Vector 

if ( d ==2 ) 
    for it = 1:nt+1 %(x,y)-> (-y,x)
         E(it,1,2) = -E(it,2,1) ;
         E(it,2,2) = E(it,1,1) ;
    end 
elseif (d==3) %d ==3
    % second basis
    % use the trick as in 2D for the first two-componets
    ic =2;% second basis vector 
    for it = 1: ns+1
        E(it,2,ic) = -E(it,3,1); % project to (y,z) plane and then rotate 
        E(it,3,ic) = E(it,2,1) ;
        E(it,1,ic) = E(it,1,1); 
        tmpv  = E(it,1:d,1);
        E(it, 1:d, ic) = E(it, 1:d, ic) - tmpv .* ...
            dot(E(it, 1:d, ic),  tmpv ) / dot(tmpv,tmpv); %project out
        
        tmpv(1:d) = E(it, 1:d, ic);
        norm_monitor(it)= norm(tmpv); %normalizatioin
        E(it,1:d,ic) = E(it,1:d,ic) ./ norm( tmpv );
    end
    if  min(norm_monitor(1:ns)) < 0.001 
       sprintf('Warning in normalization of Second Basis Vector for dim=3: ic=%d,  norm=%0.5g',...
                    ic, min(norm_monitor(1:ns)))
          pause
    end 
    ic= 3;% the third basis  
    for it =1 : ns+1 
        tmpv(1:d) = E(it, 1:d, 1); tmpv1(1:d) = E(it, 1:d, 2);
        E(it,1:d,ic) =  cross(tmpv,tmpv1);
        tmpv(1:d) = E(it, 1:d, ic);
        norm_monitor(it)= norm(tmpv);
        E(it,1:d,ic) = E(it,1:d,ic) ./ norm( tmpv );
    end 
    %orthonormal basis
elseif (d == 4) 
    GS_d = 1;  
    
    rng(4889);
    tmesh_extended = extend_period(tmesh(1:nt),nt,Period);
    for ic = 2:d % basis vector RECURSIVELY
        for id = 1:d
            X(1:2*nt) = extend_period(E(1:nt,id,ic-1), nt, 0);
            PP_spline = spline(tmesh_extended(1:2*nt),X(1:2*nt));
            pp2 = fnder(PP_spline);
            E(1:nt+1,id,ic)= ppval(pp2,tmesh(1:nt+1));
            %spline derivative
                     [E(1:nt,id,ic),~]=fourierderivative(E(1:nt,id,ic),0,Period,50-ic*5);
            %low-pass filter
            
            %         [~,E(1:ns,id,ic)]=fourierderivative(E(1:ns,id,ic-1),0,arcLength,ns);
            
        end % derivative by using spline interpolation or Fourier
        
        if  (ic >=4 ) %replace by Gaussian Random Function % no need
            sprintf('Gaussian Random Function is used for ic= %d', ic)
            for id = 1:d
                RK = randn(nt,1)+ 1i.*randn(nt,1);
                E_mode= randi([2 8],1,1);
                RK(E_mode +1 :Nx-E_mode)=0;
                E(1:nt,id,ic) = real(fft(RK(1:nt))) ;
            end
        end
        
        
        
        %%% build next basis vector, starting from the derivative of previous
        %      % followed by the Gram-Schmildt procedure.
        for it = 1:nt
            for ir = 1: min(ic-1,GS_d) %project out the old vectors
                tmpv = E(it,1:d,ir);
                E(it, 1:d, ic) = E(it, 1:d, ic) - tmpv(1:d)  .* ...
                    dot(E(it, 1:d, ic),  tmpv) / dot(tmpv,tmpv);
            end
            
            %normalization
            tmpv(1:d) = E(it, 1:d, ic);
            norm_monitor(it)= norm(tmpv);
            E(it,1:d,ic) = E(it,1:d,ic) ./ norm( tmpv );
            
            %debug only
            %           if (ic == 2)
            %                E(it,1:d,ic) =  norm( tmpv );
            %           end
            %         if ic == 2
            %             tmpv = E(it,1:d,ic-1); tmpv1(1:d) = E(it, 1:d, ic);
            %             E(it,1:d,ic) = ( dot(tmpv1(1:d),  tmpv(1:d))./norm(tmpv(1:d))./norm(tmpv1(1:d)));
            %         end
            
            
        end
        
        
        if  min(norm_monitor(1:ns)) < 0.01
            sprintf('Warning in normalization: ic=%d,  norm=%0.5g',...
                ic, min(norm_monitor(1:ns)))
            pause
        end
        %
        E(nt+1,:,:)=E(1,:,:);
        
    end
    
    
 

end %end if d==?


%% inverse of E(t)

tmp = zeros(d,d);
for id = 1:d
    for it=1:nt+1
        tmp(1:d,1:d) = E(it,1:d,1:d);
        invE(it,:,:) = inv(tmp);
    end 
end 
invE(nt+1,:,:)=invE(1,:,:);
 

figure(21); grid on ;  
for ic = 1: d
    subplot(d,1,ic); hold off; 
    for id = 1: d 
        plot_Fourier(E(1:ns,id,ic));  hold on ;
        axis([0 3 -20 4]);   
        title(['Fourier mode of E' ,num2str(ic) ], 'FontSize',20);   
    end 
end 


figure(31)
for ic = 1: d
    subplot(d,1,ic); 
    hold off;
     for id = 1:d 
        plot(tmesh(1:nt+1),E(1:nt+1,id,ic),'.-'); hold on ; grid on ;
     end  
     title ([' E' ,num2str(ic)],'FontSize',18)
 end 
 
 
 
figure(32)
 % uniform_arc(1:ns) matches E(1:ns)
 for ic = 1: d
     subplot(d,1,ic); 
    hold off;
     for id = 1:d 
        plot(tmesh(1:nt+1),invE(1:nt+1,id,ic),'.-'); hold on ; grid on ;
     end  
     title ([' inverse of E' ,num2str(ic)],'FontSize',18)
 end 
 
%  diagonistic plot
figure(44); plot_SVD_check(E(1:ns,1:d,1:d)); 
title('SVD (log10) of E: check singularity   ','FontSize',24);
axis([0 ns -3 2]);   



  
 %% %time derivative of E(t)
E(ns+1,:,:)=E(1,:,:);
tmesh_extended = extend_period(tmesh(1:nt),nt,Period); 
for ic = 1:d
for ir = 1:d 
     X(1:2*ns) = extend_period(E(1:ns,ir,ic), ns, 0);
     PP_spline= spline(tmesh_extended(1:2*ns),X(1:2*ns));
     pp2 =  fnder(PP_spline);
     Edt(1:ns,ir,ic)= ppval(pp2,tmesh(1:ns));
end
end
Edt(1+ns,:,:)=Edt(1,:,:);


 figure(41)
 % uniform_arc(1:ns) matches E(1:ns)
 for ic = 1: d
      subplot(d,1,ic);
      plot_Fourier(Edt(1:ns,1,ic));
      title (['Fourier Mode of Edt'  num2str(ic)],'FontSize',18);
 end

 
  figure(51)
 for ic = 1:d 
      subplot(d,1,ic);
   hold off
      for id= 1:d  
          plot(tmesh(1:ns+1), Edt(1:ns+1,id,ic),'-'); hold on ;
      end 
      title ([' Edt'  num2str(ic)], 'FontSize',18)
 end 

 
 

%%%%%%%%% 
%  END of CONSTRUCTION OF BASIS 
%%%%%%%







%% functions


function  XX = extend_period(X,nt,T)
    % extend the function 'X' from one period T
    % to two periods T;
if (nt ~= length(X)) 
       error('erro in extend_period')
end
XX= zeros(2*nt,1);
nt1 = floor(nt/2); nt2 = nt - nt1;
XX(1:nt1) = X(nt-nt1+1:nt)-T;
XX(nt1+1:nt1+nt) = X(1:nt);
XX(nt1+nt+1:nt1+nt+nt2) = X(1:nt2)+T;
end 


function [Q,R] =  mgs(X)
    % Modified Gram-Schmidt.  [Q,R] = mgs(X);
    % G. W. Stewart, "Matrix Algorithms, Volume 1", SIAM, 1998.
    [n,p] = size(X);
    Q = zeros(n,p);
    R = zeros(p,p);
    for k = 1:p
        Q(:,k) = X(:,k);
        for i = 1:k-1
            R(i,k) = Q(:,i)'*Q(:,k);
            Q(:,k) = Q(:,k) - R(i,k)*Q(:,i);
        end
        R(k,k) = norm(Q(:,k))';
        Q(:,k) = Q(:,k)/R(k,k);
    end
end


    function [] = plot_Fourier(f)
        Fk= fft(f); Fk=Fk(1: floor(length(f)/2));
        plot(log10(1: floor(length(f)/2)), log10(abs(Fk)./length(f)),'o-');
        grid on ;
        ylabel('log10');
    end

    


function [F, dfdx] = fourierderivative(f,a,b,E_mode)
 %https://www.mathworks.com/matlabcentral/fileexchange/39700-fourier-derivative?s_tid=gn_loc_drop
% FOURIERDERIVATIVE Fourier derivative
%     dfdx = FOURIERDERIVATIVE(f,a,b) approximates the derivative a
%     discrete function f over the domain (a,b).  f, a vector, must be
%     uniformly sampled, periodic, and contain an even number of samples.
%     For best results, f should be periodic such that f(x + a) = f(x + b).
%     As an example,
%
%          x = linspace(0,pi);
%          f = exp(cos(x).*sin(2*x));
%          dfdx = fourierderivative(f,0,pi);
%
%     Results for nonperiodic f are dubious.
% Warning :: the derivative is not trustable near boundaries 'a' or 'b'.

sf =size(f);
Nx = max(sf);

% the original choice:
%   k= [0:Nx/2-1 0 -Nx/2+1:-1]; 
% is equivalent to
   k = [0:Nx/2-1  -Nx/2:-1];
%

%better one 
%  k=[0:Nx/2-2 0 Nx/2 -Nx/2+1:-1  ] ; 
%best choice with lowest erro 
% refer to https://www.mathworks.com/matlabcentral/fileexchange/52642-fourier-differentiation
 
k = 2*pi/(b-a).* k;
if(sf(1) > 1 );
    k = k';
end
Fk=fft(f);

E_mode = min(Nx,E_mode);
Fk(E_mode +1 :Nx-E_mode)=0;

dfdx = real( ifft(1i*k.*Fk) );
F = real(ifft(Fk));
% 
% 
% figure(2);
% plot(f); hold on
% plot(log(abs(fft(f))),'o');
% 
% %plot(abs(ifft(Fk)),'r') ;
% 
% %dfdx2 = real( ifft(1i*k.* Fk) );
% 
% plot(dfdx,'r') ;
% 
% hold off
% pause

end

    function  pp2 = fnder(PP_spline)
        [breaks,coefs,l,k,dpp] = unmkpp(PP_spline);
        pp2 = mkpp(breaks,repmat(k-1:-1:1,dpp*l,1).*coefs(:,1:k-1),dpp);
    end 
end 

