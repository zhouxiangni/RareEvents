
 
function result=dynfun(t,x)

global case_id
%input and output are both column vectors

if nargin == 1
    t=0;
end 

 
 sx= size(x);
 if sx(2) > 1
     x=x';
 end
 
dim = length(x);
result=zeros(dim,1);

        
   
        
 switch case_id 
    
    case 2  
%n=2; %Van der Pol
    result = [x(2) ;  (1-x(1).^2).* x(2) - x(1)] ;
    
    case 22 
    result = [x(1) - x(1)^3 ./3 + x(2) - x(2)^3./9 ; ...
              x(1)+0.9];
          
     case 221 
           result = [x(1) - x(1)^3 ./3 + x(2) - x(2)^3./9 ; ...
              x(1)+0.9];
          result = - result ;  % time inverse of 22 case
         
        

% 3D Lotka?Volterra 
%n=3;
    case 3 
     mu = 1.480;  A=[2 5 0.5; 0.5 1 mu; 1 0.5 1];
    r = sum(A,2);
    result = x.*(r-A*x);

 % refer to Limit Cycles in Competition Communities
 % by Michael E Gilpin, The American Naturalist. vol 109, No. 965, pp51-60,
 % 1975
 
     case 4
         %PHYSICAL REVIEW E 79, 051131  2009 
       b=3.3; c=2; cp = 1;
       result(1) = 1 - x(1)*( 1+b-c*x(1)*x(2));
       result(2) =     x(1)*(   b-c*x(1)*x(2));
       result(3) = 1 - x(3)*( 1+x(2) -cp*x(3)*x(4));
       result(4) =     x(3)*(   x(2) -cp *x(3)*x(4));
% n=5;
    case 50 
        a = 0.50 ;  b = 1.5 ; c = 0.02;
        A = ones(dim,dim).*a;
        for i = 1:dim
            A(i,i) = 1;
        end
        A(dim,1)=b;
        for i = 1: dim-1
            A(i,i+1)=b;
        end
        result = x.* ( ones(dim,1) - A*x)+c;
 
     case 5 
        result(1:dim) = 1-x(1:dim); 
        result(1) = result(1) -2.*h(x(5)) -2.* h(x(2)) + 2.* h(x(2)).*h(x(5)) ...
                    +2.* h(x(2)).*h(x(4)) - 2.* h(x(2)).*h(x(4)).*h(x(5));
        result(2:5) = result(2:5)-2.*h(x(1:4));
        result(3) = result(3) + 2.* h(x(2)).*h(x(5));
        result(5) = result(5) - 2.*h(x(2)) +  2.* h(x(2)).*h(x(4));
 end 
 
 

function s=h(x)
    eta = 2.0;
    s = 0.5.*(1+tanh(x.*eta));
end 

end 


 
 
 
