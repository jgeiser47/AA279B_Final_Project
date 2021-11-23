function plot_contour_blank(s)
% Purpose: Iterate over 2D array of departure dates and arrival dates to
%          create a porkchop plot. 

% Create contour plot
figure(); hold on; grid on; axis equal;
contour_helper(s);
title('C_3 vs Departure/Arrival Date');
hcb = colorbar; 
hcb.Title.String = 'C_3';
hcb.Title.FontWeight = 'bold';
end