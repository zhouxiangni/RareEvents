  function MD=Find_Monodromy(tmesh,J)
  
%     nt = length(tmesh)-1;
%     T =  tmesh(nt+1);
%     Jsize = size(J); 
%     d = Jsize(2);
%     
%     for it = 1: nt
%         tspan =[tmesh(it) tmesh(it)+T/2 tmesh(it) + T];
%         [~,MD(it,)]=ode45(@PODE, tspan,)
%     end 
%     
%     MD(nt+1,:,:) = MD(1,:,:);
  end  