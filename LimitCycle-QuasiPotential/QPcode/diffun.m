 

 function dm = diffun(x)
 
 global diff_case_id

%input : x , column vector with size d
%output :  dm , matrix with size d*d

  

sx=size(x);
if (sx(1)==1)
    x=x';
end

d = length(x);
dm = eye(d,d);

switch diff_case_id
    case 1  % constant 1
        dm = eye(d,d);
        return ; 
    case 2 % for van del poll 
        if (x(1) < 0.0)
            dm = zeros(d,d);
        end
    case 3 % langevian
        dm(1,1)=0;
    case 4 
        if (x(2)<2.8)   
            dm = 50.*eye(d,d);
        end 
end 