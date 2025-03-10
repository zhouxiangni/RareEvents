  function p=getMomentum(z,itin,nt,d,tmesh,LC, Omega,E, Einv, G)
       
        sz = size(z);
        if sz(2) >1 
            z = z';
        end 
        if length(z) ~= d-1
            error('getMomentum: dim of z is wrong')
        end 
        
        Gz = G(itin)*z;
        p=zeros(d,1);  tmpm=zeros(d-1,d-1);
        p(1:d) = E(itin,1:d,1);
        p = p ./ norm( dynfun(0,LC(itin,1:d)));
        
        tmp = zeros(nt+1,1);
        for it = 1:nt+1
            tmp(it) = dot(z, G(it)*z);
        end 
        PP_spline = spline(tmesh(1:nt+1),tmp(1:nt+1));
        pp2 = fnder(PP_spline);
        tmp = ppval(pp2,tmesh(itin));
        tmpm(1:d-1,1:d-1) = Omega(itin,1:d-1,1:d-1);
        tmp = tmp/2 - dot(z, tmpm' *Gz);
        
        p = p.*tmp;

        for i = 1:d-1
            for j = 1:d
            p(j) = p(j) + Gz(i).* Einv(itin,i+1,j);
            end 
        end 
    end    
