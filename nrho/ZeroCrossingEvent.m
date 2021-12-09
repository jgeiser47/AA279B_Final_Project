function [value,isterminal,direction] = ZeroCrossingEvent(t,y)
value = y(2); % The value that we want to be zero
isterminal = 1;  % Don't halt integration 
direction = 0;   % The zero can be approached from either direction
end

