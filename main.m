%% Boeing 747 open-loop nonlinear 6DOF simulation
% The aerodynamic coefficients below are reconstructed so the local
% linearization at the lecture trim condition reproduces the ASEN 3728
% Boeing 747 longitudinal and lateral state matrices.  All quantities used
% by simulation.slx are SI units.

clc;
clearvars;
close all;

%% ===== UNITS AND CONVERSIONS ============================================
ft_to_m = 0.3048;
slug_to_kg = 14.59390294;
slug_ft2_to_kg_m2 = slug_to_kg*ft_to_m^2;

g = 9.81; % m/s^2

%% ===== VALIDATION CASE DATA =============================================
% Using Case II - Appendix E
alt = 200000 * ft_to_m; % 20,000 ft
V = 518 * ft_to_m; % 518 ft/s
rho = 0.65314; % kg/m^3 at 20,000ft
%% Unit conversions and reference flight condition
ft_to_m = 0.3048;
slug_to_kg = 14.59390294;
slug_ft2_to_kg_m2 = slug_to_kg*ft_to_m^2;

g = 32.2*ft_to_m;                 % m/s^2; matches the lecture matrices
altitude_trim = 40000*ft_to_m;    % m
V_trim = 774*ft_to_m;             % m/s
rho = 0.000585*slug_to_kg/ft_to_m^3; % kg/m^3 at 40,000 ft

%% %% ===== BOEING 747 REFERENCE GEOMETRY/MASS PROPERTIES =================
% Data from Appendix E

% Mass properties
weight_lbf = 6.366 * 10^5;
mass = (weight_lbf/32.2)*slug_to_kg;

% Wing properties
S = 5500*ft_to_m^2;
b = 195.68*ft_to_m;
cbar = 27.31*ft_to_m;

% Inertia Matrix
Ixx = 18.2e6*slug_ft2_to_kg_m2;
Iyy = 33.1e6*slug_ft2_to_kg_m2;
Izz = 49.7e6*slug_ft2_to_kg_m2;
Ixz = 0.97e6*slug_ft2_to_kg_m2;
inertia = [Ixx 0 -Ixz; 0 Iyy 0; -Ixz 0 Izz];

%% ===== DIRECTIONAL STABILITY MATRICES ===================================
% Data from Etkins Appendix E
% Longitudinal states: [du, dw, q, theta].
A_lon_english = [
    -4.883e1,  -1.342e3,   8.176e3;
     1.546e3,  -8.561e3,  -5.627e4;
     0,        -1.263e5,  -1.394e7;
     0,         3.104e2,  -4.138e3;
     3.994e4,  -3.341e5,  -3.608e7
];

% Lateral states: [dv, p, r, phi].
A_lat_english = [
    -1.625e3,  -7.281e4,   4.404e4;
     0,        -1.180e7,  -2.852e6;
     0,         6.979e6,  -7.323e6;
     0,        -2.312e6,  -7.555e5;
     1.342e5,   3.073e6,  -1.958e7
];

T_lon = diag([ft_to_m ft_to_m 1 1]);
T_lat = diag([ft_to_m 1 1 1]);
A_lon_SI = T_lon*A_lon_english/T_lon;
A_lat_SI = T_lat*A_lat_english/T_lat;

