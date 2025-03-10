function GT = Int_PRDE1D_Period(G0,T, tmesh, M,A)
    options = odeset('RelTol', 1.e-4, 'AbsTol',1.e-8);
    [t,y] = ode45(@(t,G) PRDE(t,G,T,tmesh,M,A), [0,T/2,T], G0, options);
    GT = y(length(t))-G0;
end 