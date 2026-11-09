
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
lbf_to_N = 4.448221615;

g = 9.81; % m/s^2

%% ===== VALIDATION CASE DATA =============================================
% Using Case II - Appendix E
alt = 20000 * ft_to_m; % 20,000 ft
V_trim = 518 * ft_to_m; % 518 ft/s
rho = 0.65314; % kg/m^3 at 20,000ft

% Case II trim orientation
xi_trim = deg2rad(-6.8);
alpha_trim = -xi_trim;
theta_trim = alpha_trim; % Level flight

%% %% ===== BOEING 747 REFERENCE GEOMETRY/MASS PROPERTIES =================
% Data from Appendix E

% Mass properties
weight_lbf = 6.366 * 10^5;
mass = (weight_lbf/32.2)*slug_to_kg;

% Wing properties
S = 5500*ft_to_m^2;
b = 195.68*ft_to_m;
cbar = 27.31*ft_to_m;
AR = b^2 / S;

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

% Conversion factors for input variables
lon_inputs = [ft_to_m ft_to_m 1 ft_to_m 1];
lat_inputs = [ft_to_m 1 1 1 1];

% Conversion factors for output forces and moments
outputs = [lbf_to_N lbf_to_N lbf_to_N*ft_to_m];

% Convert English derivatives to SI
A_lon_SI = diag(1./lon_inputs)*A_lon_english*diag(outputs);
A_lat_SI = diag(1./lat_inputs)*A_lat_english*diag(outputs);

%% ===== INITIAL CONDITIONS ===============================================

% setting initial state vector
pos0 = [0 0 -alt]; % Position vector (N E D), m
vel0 = [V_trim*cos(alpha_trim) 0 V_trim*sin(alpha_trim)]; % (u v w), m/s
euler0 = [0 theta_trim 0]; % Euler angles (phi theta psi), rad
omega0 = [0 0 0]; % [p q r], rad/s

%% ===== AERO COEFFICIENTS ================================================

% Find q and CL
qbar = 0.5 * rho * V_trim^2;
L = mass * g;

% Drag coefficients
CD_trim = 0.040; % appendix e (case 2)

% LIFT coefficients
CL_trim = L/(qbar*S) - CD_trim*tan(alpha_trim);

% -- FINDING CL_ALPHA --

% Lift stability derivative from Appendix E, Table E.3

% Convert force derivatives to SI units
lbf_to_N = 4.448221615;

% Extract dimensional force derivatives (SI)
Xu = A_lon_SI(1,1);
Xw = A_lon_SI(2,1);
Zu = A_lon_SI(1,2);
Zw = A_lon_SI(2,2);

% Force derivatives with respect to alpha at constant airspeed
X_alpha = -vel0(3)*Xu + vel0(1)*Xw;
Z_alpha = -vel0(3)*Zu + vel0(1)*Zw;

% Lift-curve slope (1/rad)
CL_alpha = (sin(alpha_trim)*X_alpha - ...
    cos(alpha_trim)*Z_alpha)/(qbar*S) - CD_trim;

% Final Lift Coefficient
alpha = atan2(vel0(3), vel0(1));
CL = CL_trim + CL_alpha*(alpha - alpha_trim);


%% ====== DIMENSIONAL DERIVATIVES =========================================

% initial conditions
d_lat = zeros(5,1);
d_lon = zeros(5,1);

L_trim = qbar * S * CL_trim;
D_trim = qbar * S * CD_trim;

F_aero_trim = [ -D_trim * cos(alpha_trim) + L_trim * sin(alpha_trim)
    0
    -D_trim * sin(alpha_trim) - L_trim * cos(alpha_trim)];

% Simplified thrust through CG, aligned with body X
T_trim = D_trim/cos(alpha_trim);
T_body_trim = [T_trim; 0; 0];

M_aero_trim = [0; 0; 0];
M_thrust_trim = [0; 0; 0];

%% ===== PROPULSION (linear throttle) =========================
% T_max: 4x JT9D-3 maximum thrust at Case II (20,000 ft, M ~0.50).
% Thrust data taken from NASA CR-114494 Vol. II p.12.2-4: Fn/delta 41,000 lbf at max EPR, M 0.5
% delta @ 20,000 ft = 0.4599  leading to 19,000 lbf per engine (rough estimate from graph reading)
n_eng     = 4;
T_max_eng = 19000*lbf_to_N;
T_max     = n_eng*T_max_eng; % N total around 338 kN

% Throttle setting that gives the trim thrust (using linear model)
throttle_trim = T_trim/T_max;        % ~0.50

%% ===== VERIFY TRIM EQUILIBRIUM ===========================================

% Gravity in body axes
F_gravity_trim = mass*g * [
    -sin(theta_trim);
     0;
     cos(theta_trim)
];

% Total forces and moments
F_total_trim = F_aero_trim + T_body_trim + F_gravity_trim;
M_total_trim = M_aero_trim + M_thrust_trim;

% Display results
disp('Total Forces at Trim (N):');
disp(F_total_trim);

disp('Total Moments at Trim (N*m):');
disp(M_total_trim);

%% ===== AERO FUNCTION VALIDATION =========================================

controls0 = zeros(3,1);

[F_test, M_test] = aero( ...
    vel0(:), omega0(:), euler0(:), controls0, ...
    A_lon_SI, A_lat_SI, vel0(:), ...
    F_aero_trim, T_body_trim, mass, g);

fprintf('Forces at trim (N):\n');
disp(F_test);

fprintf('Moments at trim (N*m):\n');
disp(M_test);

controls_test = [deg2rad(-1); 0; 0];

[F_test, M_test] = aero( ...
    vel0(:), omega0(:), euler0(:), controls_test, ...
    A_lon_SI, A_lat_SI, vel0(:), ...
    F_aero_trim, T_body_trim, mass, g);

disp('Forces with elevator perturbation:');
disp(F_test);

disp('Moments with elevator perturbation:');
disp(M_test);

