function  [tmesh, LC, E, Edt, invE] = FindE2(tmesh,LC,ns,... 
                            dynfun_handle,Jacobi_handle,diffun_handle)


%tmesh (input):  col vector with size nt+1 : tmesh(1)=0, tmesh(nt+1)=Period
%      (output): col vector with size ns+1 : tmesh(1)=0, tmesh(nt+1)=Period
%                uniform grided 

%LC:  size nt+1 by d ;  each row is a point in R^d space ; match 'tmesh'
%       
% ns:  the number of discrete points on LC after reparametrization 

% output:   E, Edt,    size(1:ns+1,d,d),  with 'ns+1' points on LC
%           M and A,   size(1:ns+1,d-1,d-1). The first row/col deleted already 


%% start the code here ....

global dim

d = dim; %dim of the dynamics

 

E=zeros(ns+1,d,d); 
invE= zeros(size(E));  
Edt= zeros(size(E));  



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

    LC(1:ns+1,id,1) = ppval(PP_spline,linspace(0,Period,ns+1));
    %spline derivative 
%     pp2 = fnder(PP_spline);
%     Edt(1:ns+1,id,1) = ppval(pp2,linspace(0,Period,ns+1));
end 


LC(ns+2:nt+1,:)=[];
tmesh  = linspace(0,Period,ns+1)'; 
nt = ns; %update LC, tmesh, nt.





%% 
% for id = 1:d 
% %debug
%     E(1:nt,id,2) = diff(LC(1:nt+1,id))./ diff(tmesh(1:nt+1));
%     [~,E(1:nt,id,1)] = fourierderivative(LC(1:nt,id),0,Period,nt) ;
% end         
% E(nt+1,:,:)=E(1,:,:);

b_vec=zeros(nt+1,d);
for it = 1: nt+1 
   b_vec(it,1:d) = dynfun_handle(0,LC(it,1:d) ); 
end 



%% Build The First Basis Vector 

tmpm = zeros(d,d);
tmpv = zeros(d,1);  tmpv1 = zeros(d,1);
norm_monitor = zeros(ns+1,d+3);

 

E(1:nt,1:d,1) = b_vec(1:nt,1:d);  %or simply use the b_vec

% first basis vector (normalizes)
for it = 1:nt
    tmpv(1:d) = E(it,1:d,1); E(it,1:d,1) = E(it,1:d,1) ./ norm(tmpv);
end 

E(nt+1,:,:)=E(1,:,:);


 
 
%% Build Other  Basis Vector 

    
    options = odeset('RelTol', 1.e-5, 'AbsTol',1.e-8);
    
    disp(' ... solving Monodromy Matrices (be patient)...');
    
    for id = 1: d  % parallel in each dimension
        LC_spline{id}=spline(tmesh,LC(:,id));
    end
    
    indx= [1:5:nt+1];
    Lyp_Exp = zeros(length(indx),d);

    for kk = 1: length(indx)
        it=indx(kk);
%       fhandle=@(t,phi) TransitJacobi(t,phi,Period,tmesh(1:nt+1),LC(1:nt+1,1:d),Jacobi_handle);
        fhandle=@(t,phi) TransitJacobi(t,phi,Period,tmesh,LC_spline,Jacobi_handle);

        [tmp,ytmp] = ode45(fhandle,[tmesh(it),Period+tmesh(it)], reshape(eye(d,d),[],1), options);
        tmpm =reshape(ytmp(length(tmp),:),d,d);
       
        [V,D,W] = eig(tmpm); D=diag(D);
        D= (real(D)); % it SHOULD be real, originally.        
        V=real(V); 
        W=real(W);
        [D, I]=sort(D,'descend');% must be the same for all `it'
      

        V(1:d,1:d) = V(1:d,I);
  
        W(1:d,1:d) = W(1:d,I);  

        E(it,1:d,2:d) =  W(1:d,2:d); 
        %left eigenvector , not orthogonal           
%          norm_monitor(it,1:d)= D(1:d);
%          tmpv(1:d) = E(it,1:d,1);
%          norm_monitor(it,6) = norm(tmpm*tmpv-tmpv)./norm(tmpv);
    end
%     
% K = 2;
% for kk = 1:K
%     tmesh(nt*kk+2:(kk+1)*nt+1) = Period*kk  + tmesh(2:nt+1);
%     LC(nt*kk+2:(kk+1)*nt+1,1:d) = LC(2:nt+1,1:d);
% end 
% for id = 1: d  % parallel in each dimension
%     LC_spline{id}=spline(tmesh,LC(:,id));
% end 
% Period=Period*K;
% indx= [1:25:nt+1];
%     for kk = 1: length(indx)
%         it=indx(kk);
%         fhandle=@(t,phi) TransitJacobi(t,phi,Period,tmesh,LC_spline,Jacobi_handle);
% 
%         [tmp,ytmp] = ode45(fhandle,[tmesh(it),Period+tmesh(it)], ...
%                                    reshape(eye(d,d),[],1), options);
%         tmpm =reshape(ytmp(length(tmp),:),d,d);
%        
%         [V,D,W] = eig(tmpm); D=diag(real(D)); % it SHOULD be real, originally.        
%         V=real(V); W=real(W);   [D, I]=sort(D,'descend');
%         V(1:d,1:d) = V(1:d,I);  W(1:d,1:d) = W(1:d,I);  
%         E(it,1:d,4) =  W(1:d,4);%  eigenvector        
%         E(it,1:d,5) =  W(1:d,5);%  eigenvector   
%            
%          Lyp_Exp(kk,1:d)=D;
% 
%     end
%     
%     figure(66); plot(Lyp_Exp,'o-')

   %  first basis vector (normalizes)
     E(1:nt+1,1:d,1) = b_vec(1:nt+1,1:d);  %or simply use the b_vec
     for it = 1:nt+1
         tmpv(1:d) = E(it,1:d,1); E(it,1:d,1) = E(it,1:d,1) ./ norm(tmpv);
     end 
   
   
    % fix the sign probblem
    tmpv = zeros(d,1);  tmpv1 = zeros(d,1);
    
    for ic = 2:d
        tmpv(1:d) = E(1,1:d,ic);
        tmpv1(1:d) = ones(size(tmpv));
        if ( norm(tmpv1-tmpv)> norm(tmpv1+tmpv)) %change sign 
            E(1,1:d,ic) = -E(1,1:d,ic);
        end
        for kk =2: length(indx)
            it =indx(kk);
            tmpv(1:d) = E(indx(kk),1:d,ic);
            tmpv1(1:d) = E(indx(kk-1),1:d,ic);
            if ( norm(tmpv1-tmpv)> norm(tmpv1+tmpv)) %change sign 
                E(it,1:d,ic) = -E(it,1:d,ic);
            end
            norm_monitor(indx(kk),d+2) = dot( tmpv,tmpv1);
            tmpv(1:d) = E(indx(kk),1:d,ic);
            tmpv1(1:d) = E(indx(kk-1),1:d,ic);
            norm_monitor(indx(kk),d+3) = dot( tmpv,tmpv1);
        end
    end
    
 %last vector 
%     for it = 1: nt +1
%         for ic = 1: 4
%             tmpm(ic,1:d) = E(it,1:d,ic);
%         end 
%         tmpm(5,1:d)=ones(1,d);
%         rhs = [0 ;-1 ;1 ; 0; 1];
%         tmpv = tmpm \ rhs;
%         E(it,1:d,5)= tmpv';
%     end 
        
 
 
 %% interpolte all tmesh points  
 
for ic = 1:d
for ir = 1:d 
      PP_spline= spline(tmesh(indx),E(indx,ir,ic));
      E(1:ns+1,ir,ic)= ppval(PP_spline,tmesh(1:ns+1));
      pp2 =  fnder(PP_spline);
      Edt(1:ns+1,ir,ic)= ppval(pp2,tmesh(1:ns+1));
end
end

 


% ALL time it are obtained

%% normalization 
for it = 1:nt +1 
for ic = 1:d
    tmpv = E(it,1:d,ic);
    E(it,1:d,ic) = E(it,1:d,ic) ./ norm(tmpv);
end 
end 

 


% check period 
for ic = 1:d 
    if ( norm(E(nt+1,:,ic)-E(1,:,ic)) > 1.e-1)
        disp(['NOT period for E' num2str(ic) ]) ; 
    end 
end 

if length(tmesh)>2*nt
tmesh(nt+2:end) = [];
LC(nt+2:end,:) = [];
Period = Period /K;
end 

 

%% calculate Edt
indx = 1:nt+1;
for ic = 1:d
for ir = 1:d 
      PP_spline= spline(tmesh(indx),E(indx,ir,ic));
      pp2 =  fnder(PP_spline);
      Edt(indx,ir,ic)= ppval(pp2,tmesh(indx));
end
end

%%  calcualte the inverse 
for it = 1:length(tmesh);
    tmpm(1:d,1:d) = E(it,1:d,1:d);
    invE(it,1:d,1:d)=inv(tmpm);
end

 


 
 

%%  plot result 
   indx = 1: nt+1;
   
   %  diagonistic plot
figure(44); plot_SVD_check(E(indx,1:d,1:d)); 
title('SVD (log10) of E: check singularity   ','FontSize',24);
% axis([0 ns -3 2]); 
 
 % uniform_arc(1:ns) matches E(1:ns)
 for ic = 1:d
       
figure(10+ic)
%      subplot(d,1,ic); 
    hold off;
     for id = 1:d 
        plot(tmesh(indx),E(indx,id,ic),'o-','MarkerSize',6); hold on ; grid on ;
     end  
     title ([' E' ,num2str(ic)],'FontSize',18);
     axis([0,Period, -1 1]);
 end 


 figure(21); grid on ;  
for ic = 1: d
    subplot(d,1,ic); hold off; 
    for id = 1: d 
        plot_Fourier(E(1:ns,id,ic));  hold on ;
        axis([0 3 -20 4]);   
        title(['Fourier mode of E' ,num2str(ic) ], 'FontSize',20);   
    end 
end 

 
 % uniform_arc(1:ns) matches E(1:ns)
 for ic = 1: d
     figure(30+ic)
%      subplot(d,1,ic); 
     hold off;
     for id = 1:d 
        plot(tmesh(indx),invE(indx,id,ic),'o-','MarkerSize',6); hold on ; grid on ;
     end  
     title ([' inverse of E' ,num2str(ic)],'FontSize',18)
 end 
 

 figure(41)
 % uniform_arc(1:ns) matches E(1:ns)
 for ic = 1: d
      subplot(d,1,ic);
      plot_Fourier(Edt(1:ns,1,ic));
      title (['Fourier Mode of Edt'  num2str(ic)],'FontSize',18);
 end

  
 for ic = 1:d 
%       subplot(d,1,ic);
figure(50+ic)
   hold off
      for id= 1:d  
          plot(tmesh(indx), Edt(indx,id,ic),'o-','MarkerSize',6); hold on ;
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

