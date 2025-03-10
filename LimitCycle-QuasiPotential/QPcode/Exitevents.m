 
function [value,isterminal,direction] = Exitevents(t,y)
% Locate the time when height passes through zero in a decreasing direction
% and stop integration.
value = max(abs(y))-3.0;     % detect height = 0
isterminal = 1;   % stop the integration
direction = 1;   % positive direction