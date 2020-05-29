%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% write_input_parameter.m
% Creates parameter input file for the 2D DI lava flow model
% Written by: Colton J. Conroy
%                11.26.18
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Boundary condition

wall_bc_type = 1;        % 0 = no slip, 1 = no normal flow

% Dynamic parameters

rho  = 1350;        % density
hint = 10.0;        % steady thickness h 
n    = 1.750;       % power law exponent in nonlinear viscosity relation
A    = -4.55d0;     % for n = 1
B    = 5574.2;      % for n = 1
C    = 609.4;       % for n = 1
b0   = 2.73e4;      % experimental reference temperature (for n neq 1)
mu0  = 6.025e-7;    % viscosity @ infinite temperature   (for n neq 1)
tau_yield = 0;      % yield strength 

% Thermal parameters

T_int   = 1425;        % inlet temperature
T_wall  = 1283;        % wall temperature (NOT farfield boundary)
T_crust = 750;         % ground temperature
T_air   = 363;         % air temperature
eff     = 0.850;       % effusivity constant
kt      = 0.10;        % thermal conductivity of basal boundary

% Miscellaneous

vid = 0;               % output for video? (0 = no, 1 = yes)
vid_frame_rate = 10;   % video frame rate 
limiter = 0;           % slope limiting? (0 = no, 1 = yes) Currently not implemented.

% Write parameters to text file

write_input_file

