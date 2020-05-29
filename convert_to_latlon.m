load hawaii_8_shift_val.mat
R = 6378206.4;    
% Convert meters to lat/lon
XNODES = ((XNODES + shift)./(R*cos(CPPLAT))+CPPLON)*180/pi;
YNODES = (180/pi).*YNODES./R;
