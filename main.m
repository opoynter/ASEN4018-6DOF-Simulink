%% Boeing 747 open-loop nonlinear 6DOF simulation
% The aerodynamic coefficients below are reconstructed so the local
% linearization at the lecture trim condition reproduces the ASEN 3728
% Boeing 747 longitudinal and lateral state matrices.  All quantities used
% by simulation.slx are SI units.

clc;
clearvars;
close all;

%% Unit conversions and reference flight condition
ft_to_m = 0.3048;
slug_to_kg = 14.59390294;
slug_ft2_to_kg_m2 = slug_to_kg*ft_to_m^2;

g = 32.2*ft_to_m;                 % m/s^2; matches the lecture matrices
altitude_trim = 40000*ft_to_m;    % m
V_trim = 774*ft_to_m;             % m/s
rho = 0.000585*slug_to_kg/ft_to_m^3; % kg/m^3 at 40,000 ft

%% Boeing 747 mass properties and reference geometry (Etkin example)
weight_lbf = 636636;
mass = (weight_lbf/32.2)*slug_to_kg;
S = 5500*ft_to_m^2;
b = 195.68*ft_to_m;
cbar = 27.31*ft_to_m;

Ixx = 18.2e6*slug_ft2_to_kg_m2;
Iyy = 33.1e6*slug_ft2_to_kg_m2;
Izz = 49.7e6*slug_ft2_to_kg_m2;
Ixz = 0.97e6*slug_ft2_to_kg_m2;
inertia = [Ixx 0 -Ixz; 0 Iyy 0; -Ixz 0 Izz];

%% Lecture linear models and SI similarity transformation
% Longitudinal states: [du, dw, q, theta].
A_lon_english = [-0.006868   0.01395     0       -32.2;
                 -0.09055   -0.3151    773.98     0;
                  0.0001187 -0.001026  -0.4285     0;
                  0           0          1          0];

% Lateral states: [dv, p, r, phi].
A_lat_english = [-0.0558      0        -774       32.2;
                 -0.003865   -0.4342     0.4136    0;
                  0.001086   -0.006112  -0.1458    0;
                  0           1          0          0];

T_lon = diag([ft_to_m ft_to_m 1 1]);
T_lat = diag([ft_to_m 1 1 1]);
A_lon_SI = T_lon*A_lon_english/T_lon;
A_lat_SI = T_lat*A_lat_english/T_lat;

%% Reconstruct a coefficient model at the trim point
qbar_trim = 0.5*rho*V_trim^2;
CL0 = mass*g/(qbar_trim*S);

% A representative induced-drag polar is retained; the additional speed
% and alpha derivatives make its local dimensional derivatives agree with
% the lecture A matrix.
CD0 = 0.020;
k_drag = 0.045;
CD_trim = CD0 + k_drag*CL0^2;
T_trim = qbar_trim*S*CD_trim;

CL_V = ((-mass*A_lon_SI(2,1)) - rho*V_trim*S*CL0) ...
       *V_trim/(qbar_trim*S);
CL_alpha = (-mass*V_trim*A_lon_SI(2,2) - T_trim)/(qbar_trim*S);
CL_q = 0;

CD_V = (-mass*A_lon_SI(1,1) - rho*V_trim*S*CD_trim) ...
       *V_trim/(qbar_trim*S) - 2*k_drag*CL0*CL_V;
CD_alpha = (mass*g - mass*V_trim*A_lon_SI(1,2))/(qbar_trim*S) ...
           - 2*k_drag*CL0*CL_alpha;

Cm0 = 0;
Cm_V = A_lon_SI(3,1)*Iyy*V_trim/(qbar_trim*S*cbar);
Cm_alpha = A_lon_SI(3,2)*Iyy*V_trim/(qbar_trim*S*cbar);
Cm_q = A_lon_SI(3,3)*Iyy*(2*V_trim)/(qbar_trim*S*cbar^2);

% Lateral force derivatives.
CY_beta = A_lat_SI(1,1)*mass*V_trim/(qbar_trim*S);
CY_p = A_lat_SI(1,2)*mass*(2*V_trim)/(qbar_trim*S*b);
CY_r = (A_lat_SI(1,3) + V_trim)*mass*(2*V_trim)/(qbar_trim*S*b);

% Account for Ixz when converting desired p-dot/r-dot derivatives into
% rolling- and yawing-moment coefficient derivatives.
I_pr = [Ixx -Ixz; -Ixz Izz];
LN_v = I_pr*A_lat_SI([2 3],1);
LN_p = I_pr*A_lat_SI([2 3],2);
LN_r = I_pr*A_lat_SI([2 3],3);

Cl_beta = LN_v(1)*V_trim/(qbar_trim*S*b);
Cn_beta = LN_v(2)*V_trim/(qbar_trim*S*b);
Cl_p = LN_p(1)*(2*V_trim)/(qbar_trim*S*b^2);
Cn_p = LN_p(2)*(2*V_trim)/(qbar_trim*S*b^2);
Cl_r = LN_r(1)*(2*V_trim)/(qbar_trim*S*b^2);
Cn_r = LN_r(2)*(2*V_trim)/(qbar_trim*S*b^2);

% Parameter vector consumed by the Aerodynamics MATLAB Function block.
aero_params = [rho; S; b; cbar; V_trim; ...
               CL0; CL_V; CL_alpha; CL_q; ...
               CD0; k_drag; CD_V; CD_alpha; ...
               Cm0; Cm_V; Cm_alpha; Cm_q; ...
               CY_beta; CY_p; CY_r; ...
               Cl_beta; Cl_p; Cl_r; ...
               Cn_beta; Cn_p; Cn_r];

gravity_earth = [0; 0; mass*g];
thrust_body = [T_trim; 0; 0];

%% Select an initial-condition experiment
% Options: "trim", "short_period", "phugoid", "roll", "dutch_roll",
% or "spiral".  Each input is a perturbation from trimmed controls and
% thrust; the absolute trim thrust remains applied by the model.
scenario = "short_period";

Vb_0 = [V_trim 0 0];
pqr_0 = [0 0 0];
eul_0 = [0 0 0];
pos_0 = [0 0 -altitude_trim];

switch scenario
    case "trim"
        % No state perturbation.
    case "short_period"
        alpha0 = deg2rad(1);
        Vb_0 = V_trim*[cos(alpha0) 0 sin(alpha0)];
    case "phugoid"
        Vb_0(1) = V_trim + 2;
    case "roll"
        pqr_0(1) = deg2rad(1);
    case "dutch_roll"
        beta0 = deg2rad(1);
        Vb_0 = V_trim*[cos(beta0) sin(beta0) 0];
    case "spiral"
        eul_0(1) = deg2rad(1);
    otherwise
        error("Unknown scenario: %s", scenario);
end

%% Simulation settings and execution
stop_time = 120;                   % s; long enough to see the phugoid
model_name = "simulation";
run_simulation = true;

if run_simulation
    sim_out = sim(model_name, "StopTime", string(stop_time));

    Vb_ts = sim_out.Vb_ts;
    eul_ts = sim_out.euler_ts;
    pqr_ts = sim_out.pqr_ts;
    pos_ts = sim_out.position_ts;
    air_ts = sim_out.airdata_ts;

    t = Vb_ts.Time;
    Vb = Vb_ts.Data;
    eul = eul_ts.Data;
    pqr = pqr_ts.Data;
    pos = pos_ts.Data;
    % MATLAB Function vector outputs are logged as [signal, 1, time].
    air = squeeze(air_ts.Data).';

    figure("Name", "Boeing 747 nonlinear 6DOF response");
    tiledlayout(3,2);
    nexttile; plot(t,Vb); grid on; ylabel("V_b (m/s)"); legend("u","v","w");
    nexttile; plot(t,rad2deg(eul)); grid on; ylabel("Euler angles (deg)"); legend("phi","theta","psi");
    nexttile; plot(t,rad2deg(pqr)); grid on; ylabel("Body rates (deg/s)"); legend("p","q","r");
    nexttile; plot(t,-pos(:,3)); grid on; ylabel("Altitude (m)");
    nexttile; plot(air_ts.Time,air(:,1)); grid on; ylabel("V_a (m/s)"); xlabel("Time (s)");
    nexttile; plot(air_ts.Time,rad2deg(air(:,2:3))); grid on; ylabel("Air angles (deg)"); xlabel("Time (s)"); legend("alpha","beta");
    sgtitle("Boeing 747: " + replace(scenario,"_"," "));

    % Compare the nonlinear perturbation states with the lecture model at
    % the same initial condition.  Agreement is expected only while the
    % perturbations remain small because the Simulink model is nonlinear.
    if any(scenario == ["trim" "short_period" "phugoid"])
        x0_ref = [Vb_0(1)-V_trim; Vb_0(3); pqr_0(2); eul_0(2)];
        A_reference = A_lon_SI;
        x_nonlinear = [Vb(:,1)-V_trim Vb(:,3) pqr(:,2) eul(:,2)];
        labels = ["Delta u (m/s)" "Delta w (m/s)" "q (rad/s)" "theta (rad)"];
        comparison_name = "Longitudinal";
    else
        x0_ref = [Vb_0(2); pqr_0(1); pqr_0(3); eul_0(1)];
        A_reference = A_lat_SI;
        x_nonlinear = [Vb(:,2) pqr(:,1) pqr(:,3) eul(:,1)];
        labels = ["Delta v (m/s)" "p (rad/s)" "r (rad/s)" "phi (rad)"];
        comparison_name = "Lateral-directional";
    end

    % The nonlinear solver uses a nonuniform time grid, so evaluate the
    % continuous-time linear solution directly at those sample times.
    x_ref = zeros(numel(t),4);
    for sample_index = 1:numel(t)
        x_ref(sample_index,:) = ...
            (expm(A_reference*(t(sample_index)-t(1)))*x0_ref).';
    end

    figure("Name", comparison_name + " validation");
    tiledlayout(2,2);
    for state_index = 1:4
        nexttile;
        plot(t,x_nonlinear(:,state_index),"LineWidth",1.2);
        hold on;
        plot(t,x_ref(:,state_index),"--","LineWidth",1.2);
        grid on;
        xlabel("Time (s)");
        ylabel(labels(state_index));
        legend("Nonlinear 6DOF","Lecture linear model");
    end
    sgtitle(comparison_name + " small-disturbance comparison");

    fprintf("\nReference longitudinal eigenvalues (1/s):\n");
    disp(eig(A_lon_SI));
    fprintf("Reference lateral eigenvalues (1/s):\n");
    disp(eig(A_lat_SI));
end
