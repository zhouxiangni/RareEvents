

 
function Jacob = Jacobi(x)

global case_id

 % x is column vector 

 

sx= size(x);
if sx(2) > 1
    x=x';
end
dim = length(x);

Jacob = zeros(dim,dim);
%i,j-entry is the d b_i / dx_j

 
if case_id == 2 
  
        %Van der Pol
        Jacob(1,1) = 0;
        Jacob(1,2) = 1;
        Jacob(2,1) = -2.*x(1).*x(2)-1 ;
        Jacob(2,2) = (1-x(1).^2);
        return;
end 

if case_id == 22 
        Jacob(1,1) = 1 - x(1)*x(1);
        Jacob(1,2) = 1 - x(2)*x(2)./3;
        Jacob(2,1) = 1;
        Jacob(2,2) = 0;
        return;
end
% n=3;
if case_id == 3     
        mu = 1.48;
        A=[2  5 0.5; 0.5 1 mu; 1 0.5 1];
        r = sum(A,2);
        tmp = A*x;
        for i = 1: 3
            for j = 1: 3
                Jacob(i,j) = - x(i) * A(i,j);
            end
            Jacob(i,i) =  Jacob(i,i) + ( r(i) - tmp(i));
        end
end 

if case_id == 4
       %PHYSICAL REVIEW E 79, 051131  2009 
       b=3.3; c=2; cp = 1;
       Jacob(1,1) =   -  ( 1+b) + 2*c*x(1)*x(2);
       Jacob(1,2) =   c*x(1)*x(1) ;
       Jacob(2,1) =  -2*c*x(1)*x(2);
       Jacob(2,2) =  -c*x(1)*x(1);
       Jacob(3,2) =  -x(3);
       Jacob(3,3) = -(1+x(2)) +2*cp*x(3)*x(4);
       Jacob(3,4) = cp*x(3)*x(3);
       Jacob(4,2) = x(3);
       Jacob(4,3) = x(2) -2*cp*x(3)*x(4);
       Jacob(4,4) = -cp*x(3)*x(3);
end 


%n=5;
if    case_id == 50
        a = 0.50 ;  b = 1.5 ;
        A = ones(dim,dim).*a;
        for i = 1:dim
            A(i,i) = 1;
        end
        A(dim,1)=b;
        for i = 1: dim-1
            A(i,i+1)=b;
        end
        tmp = A*x;
        for i = 1: dim
            for j = 1: dim
                Jacob(i,j) = - x(i) * A(i,j);
            end
            Jacob(i,i) =  Jacob(i,i) + ( 1 - tmp(i));
        end
        return
end 

if    case_id == 5 
        Jacob(1,1) = -1; 
        Jacob(1,2) = -2.*dh(x(2))+2.*dh(x(2))*h(x(5)) ...
                   +2.*dh(x(2)).*h(x(4)) - 2.*dh(x(2))*h(x(4))*h(x(5));
        Jacob(1,4) = 2.*h(x(2)).*dh(x(4))- 2.*h(x(2))*dh(x(4))*h(x(5)) ;
        Jacob(1,5) = -2.*dh(x(5))+2.*h(x(2))*dh(x(5)) -2.*h(x(2))*h(x(4))*dh(x(5));
        
        Jacob(2,1) =-2.*dh(x(1)); 
        Jacob(2,2) = -1;
        
        Jacob(3,2) = -2.*dh(x(2))+2.*dh(x(2))*h(x(5));
        Jacob(3,3) = -1;
        Jacob(3,5) = 2.*h(x(2))*dh(x(5));
        
        Jacob(4,3) = -2.*dh(x(3));
        Jacob(4,4) = -1;
        
        Jacob(5,2) = -2*dh(x(2))+2.*dh(x(2)).*h(x(4));
        Jacob(5,4) = -2*dh(x(4))+2.*h(x(2)).*dh(x(4));
        Jacob(5,5) = -1;
        return 
end 


    function s=h(x)
        eta = 2.0;
        s = 0.5.*(1+tanh(x.*eta));
    end

    function d= dh(x)
        eta = 2.0;
        d = 0.5*eta.*(sech(x.*eta).^2 );
    end

end 

