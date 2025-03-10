
function [M,A,Omega] = construct_MA(tmesh,LC,E,Edt,invE)

global dim 
d = dim ;
ns = length(tmesh)-1;
M = zeros(ns+1,d,d); 
Omega=zeros(size(M)); 
A = zeros(size(M));



%% contruct of Omega
% input:: :C, E, Edt, Jacobi, difffun
for it = 1: length(tmesh)
    tmpE(1:d,1:d) = E (it,1:d,1:d); 
    tmpEdt(1:d,1:d) = Edt (it,1:d,1:d);
    
    if rcond(tmpE)<1e-5
        sprintf('E is bad-conditioned at it=%d', it)
        rcond(tmpE)
        pause
    end
    
    tmpEinv(1:d,1:d) = invE(it,1:d,1:d);
    
    tmp1 = tmpEinv * tmpEdt;
    Omega(it,1:d,1:d) =  tmp1(1:d,1:d);
    
     
    M(it,1:d,1:d) = tmpEinv * Jacobi(reshape(LC(it, 1:d), [],1)) * tmpE ;
    A(it,1:d,1:d) = tmpEinv * diffun(reshape(LC(it, 1:d), [],1)) * (tmpEinv)' ;
    
end

%check Omega 
if (d == 2) 
    if max(abs(Omega)) > 1e-8
        sprintf('d=2, the Omega is non zero: %0.g ', max(abs(Omega)))
    end 
end 
if (d==3)
    if (max(abs(Omega(:,1,1))) + max(abs(Omega(:,2,2))) > 1e-8)
                disp('d=2 or 3, diagonal Omega is non zero!!')
    end 
    if (max(abs(Omega(:,1,2)+Omega(:,2,1))) > 1e-10)
        disp('d=3, the Omega is not antisymmetric' )
    end 
end 

M = M - Omega;

 % check M and Omega for dim =2,3

%% check Jacobi's stability for d= 2
% if ( d == 2) 
%      figure(1); hold off; plot(tmesh,M(:,2,2),tmesh,Omega(:,2,2));
% end 


 

 
     indx = 1: ns+1;
 for ic = 1:d
     figure(60+ic);    hold off
      for id= 1:d
          plot(tmesh(indx), Omega(indx,id,ic),'o-','MarkerSize',6); hold on ;
      end 
      title ([' Omega ' num2str(ic)] , 'FontSize',18)
 end 
 
 for ic = 1:d
      figure(70+ic);    hold off
       for id= 1:d  
          plot(tmesh(indx), M(indx,id,ic),'o-','MarkerSize',6); hold on ;
      end 
      title ([' M ' num2str(ic)] , 'FontSize',18)
 end 
 
 
  for ic = 1:d 
      figure(90+ic);    hold off
      for id= 1:d  
          plot(tmesh(indx), A(indx,id,ic), 'o-','MarkerSize',6); hold on ;
      end 
      title ([' A' num2str(ic) ], 'FontSize',18)
 end 
 

 
 
M(:,:,1)=[];
M(:,1,:)=[];
A(:,:,1)=[];
A(:,1,:)=[];
Omega(:,:,1)=[];
Omega(:,1,:)=[];



  