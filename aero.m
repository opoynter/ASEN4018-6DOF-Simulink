function [F_total, M_total] = aero(vel, omega, euler_angles, controls, ...
    A_lon_SI, A_lat_SI, vel0, F_aero_trim, T_body_trim, mass, g)
% Author: Darpan Singh - Fall 2026, Albatross Senior Design
% CALCULATES AERODYNAMIC FORCES AND MOMENTS in continuous loop

%% Calculate perts
dVel = vel - vel0;

dOmega = omega; % assumption that initial angular rates are 0

dControls = controls; % assumption that initial controls are 0

% Create lat and lon pert vectors

% du, dw, q, wdot, delElevator (0 for wdot is placeholder)
dX_lon = [dVel(1); dVel(3); dOmega(2); 0; dControls(1)];

% dv, p, r, delAileron, delRudder
dX_lat = [dVel(2); dOmega(1); dOmega(3); dControls(2); dControls(3)]; 

%% Calculate change in F and M

% dX, dZ, dM
dFM_lon = A_lon_SI' * dX_lon;

% dY, dL, dN
dFM_lat = A_lat_SI' * dX_lat;

%% Calculate wdot from Longintudinal FM

% gravitational force in Z direction
g_z = mass * g * cos(euler_angles(1)) * cos(euler_angles(2));

% other forces in Z direction
Z0 = F_aero_trim(3) + dFM_lon(2);

Z_wdot = A_lon_SI(4,2);

wdot = (Z0 + T_body_trim(3) + g_z + ...
    mass*(omega(2)*vel(1) - omega(1)*vel(2))) / (mass - Z_wdot);

%% Update dFM_lon

dFM_lon(4) = wdot;

dFM_lon = A_lon_SI' * dX_lon;

end