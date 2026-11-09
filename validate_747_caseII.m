%% Validate Etkin/Reid Boeing 747-100, Appendix E Case II
% Place this file beside main.m and aero.m, then run:
%   validate_747_caseII
%
% PURPOSE
%  1) Check exact steady-flight trim.
%  2) Independently assemble 4x4 longitudinal and lateral linear equations
%     from the Appendix E.3 dimensional derivatives and rigid-body physics.
%  3) Finite-difference the ACTUAL aero.m implementation coupled with
%     rigid-body equations to check it produces the same local dynamics.
%  4) Report predicted natural modes (eigenvalues).
%
% IMPORTANT: This script does NOT execute or linearize simulation.slx.
% For full validation, compare these A/B matrices and modes with the
% linearization or time response of the Simulink model.

run('main.m');

%% Trim equilibrium
x0 = [vel0(:); omega0(:); euler0(:)]; % [u v w p q r phi theta psi]
ctl0 = zeros(3,1);                  % elevator / aileron / rudder PERTURBATIONS
[F_trim_check, M_trim_check] = aero(vel0, omega0, euler0, ctl0, ...
    A_lon_SI, A_lat_SI, vel0, F_aero_trim, T_body_trim, mass, g);
trimF = norm(F_trim_check);
trimM = norm(M_trim_check);
fprintf('\n=== CASE II TRIM ===\n');
fprintf('F_total norm: %.8g N\nM_total norm: %.8g N*m\n',trimF,trimM);
fprintf('u0 = %.6f m/s, w0 = %.6f m/s, theta0 = %.4f deg\n', ...
    vel0(1),vel0(3),rad2deg(euler0(2)));
assert(trimF < 1e-6 && trimM < 1e-6,'FAIL: force/moment trim');

%% Reference longitudinal A from dimensional Etkin derivatives + EOM
% State order [du dw q dtheta]
m=mass; u0=vel0(1); w0=vel0(3); theta0=euler0(2);
Iyy=inertia(2,2);
Xu=A_lon_SI(1,1); Xw=A_lon_SI(2,1); Xq=A_lon_SI(3,1);
Zu=A_lon_SI(1,2); Zw=A_lon_SI(2,2); Zq=A_lon_SI(3,2);
Zwdot=A_lon_SI(4,2);
Mu=A_lon_SI(1,3); Mw=A_lon_SI(2,3); Mq=A_lon_SI(3,3);
Mwdot=A_lon_SI(4,3);

A_lon_ref=zeros(4);
A_lon_ref(1,:)=[Xu/m, Xw/m, Xq/m-w0, -g*cos(theta0)];
A_lon_ref(2,:)=[Zu, Zw, Zq+m*u0, -m*g*sin(theta0)]/(m-Zwdot);
A_lon_ref(3,:)=[Mu,Mw,Mq,0]/Iyy + (Mwdot/Iyy)*A_lon_ref(2,:);
A_lon_ref(4,:)=[0 0 1 0];

%% Reference lateral A from E.3 + x-z coupled inertias
% State order [dv p r dphi]
Yv=A_lat_SI(1,1); Yp=A_lat_SI(2,1); Yr=A_lat_SI(3,1);
Lv=A_lat_SI(1,2); Lp=A_lat_SI(2,2); Lr=A_lat_SI(3,2);
Nv=A_lat_SI(1,3); Np=A_lat_SI(2,3); Nr=A_lat_SI(3,3);
Ipr=inertia([1 3],[1 3]); % [Ixx -Ixz; -Ixz Izz]
A_lat_ref=zeros(4);
A_lat_ref(1,:)=[Yv/m, Yp/m+w0, Yr/m-u0, g*cos(theta0)];
A_lat_ref(2:3,:)=Ipr\[Lv,Lp,Lr,0;Nv,Np,Nr,0];
A_lat_ref(4,:)=[0,1,tan(theta0),0];

%% Reference control-effect matrices, no throttle input in aero.m
% B_lon: 4x1 elevator; B_lat: 4x2 [aileron rudder]
Xe=A_lon_SI(5,1); Ze=A_lon_SI(5,2); Me=A_lon_SI(5,3);
B_lon_ref=[Xe/m; Ze/(m-Zwdot); ...
    (Me+Mwdot*Ze/(m-Zwdot))/Iyy; 0];
Yda=A_lat_SI(4,1); Ydr=A_lat_SI(5,1);
Lda=A_lat_SI(4,2); Ldr=A_lat_SI(5,2);
Nda=A_lat_SI(4,3); Ndr=A_lat_SI(5,3);
B_lat_ref=zeros(4,2);
B_lat_ref(1,:)=[Yda,Ydr]/m;
B_lat_ref(2:3,:)=Ipr\[Lda,Ldr;Nda,Ndr];

%% Independently finite-difference the actual aero.m implementation
% Body 6DOF EOM: Vdot=F/m-omega x V;
% I*omegadot = M-omega x (I*omega); Euler angle kinematics.
h=1e-5;
A_num=zeros(9,9);
for k=1:9
    dx=zeros(9,1);dx(k)=h;
    A_num(:,k)=(caseII_rhs(x0+dx,ctl0,A_lon_SI,A_lat_SI,vel0, ...
           F_aero_trim,T_body_trim,inertia,m,g) ...
                  -caseII_rhs(x0-dx,ctl0,A_lon_SI,A_lat_SI,vel0, ...
           F_aero_trim,T_body_trim,inertia,m,g))/(2*h);
end
B_num=zeros(9,3);
for k=1:3
    du=zeros(3,1);du(k)=h;
    B_num(:,k)=(caseII_rhs(x0,du,A_lon_SI,A_lat_SI,vel0, ...
           F_aero_trim,T_body_trim,inertia,m,g) ...
                  -caseII_rhs(x0,-du,A_lon_SI,A_lat_SI,vel0, ...
           F_aero_trim,T_body_trim,inertia,m,g))/(2*h);
end
lonIdx=[1 3 5 8]; latIdx=[2 4 6 7];
A_lon_num=A_num(lonIdx,lonIdx);
A_lat_num=A_num(latIdx,latIdx);
B_lon_num=B_num(lonIdx,1);
B_lat_num=B_num(latIdx,2:3);

fprintf('\n=== REFERENCE LONGITUDINAL A, [du dw q theta] ===\n');disp(A_lon_ref);
fprintf('Eigenvalues, 1/s:\n');disp(eig(A_lon_ref));
fprintf('\n=== REFERENCE LATERAL A, [dv p r phi] ===\n');disp(A_lat_ref);
fprintf('Eigenvalues, 1/s:\n');disp(eig(A_lat_ref));

errs=[max(abs(A_lon_num(:)-A_lon_ref(:))),...
      max(abs(A_lat_num(:)-A_lat_ref(:))),...
      max(abs(B_lon_num(:)-B_lon_ref(:))),...
      max(abs(B_lat_num(:)-B_lat_ref(:)))];
fprintf('\n=== MAX ABSOLUTE FINITE-DIFFERENCE ERRORS ===\n');
fprintf('A_lon: %.4g   A_lat: %.4g   B_lon: %.4g   B_lat: %.4g\n',errs);
if all(errs<1e-4)
    fprintf('PASS: aero.m + rigid-body equations reproduce Case II local dynamics.\n');
else
    warning('Local implementation differs from the reference EOM; inspect A/B matrices.');
end
fprintf('\nNEXT: Compare these matrices/modes with an actual SIMULINK linearization.\n');

%% Local function, intentionally independent of the Simulink 6DOF block.
function xdot=caseII_rhs(x,controls,A_lon_SI,A_lat_SI,vel0,...
                      F_aero_trim,T_body_trim,inertia,m,g)
    Vb=x(1:3);omega=x(4:6);eta=x(7:9);
    [F,M]=aero(Vb,omega,eta,controls,A_lon_SI,A_lat_SI,vel0,...
               F_aero_trim,T_body_trim,m,g);
    Vdot=F/m-cross(omega,Vb);
    omegadot=inertia\(M-cross(omega,inertia*omega));
    phi=eta(1);theta=eta(2);
    % Body angular rates [p q r] -> [phi_dot theta_dot psi_dot]
    eulerdot=[omega(1)+tan(theta)*(omega(2)*sin(phi)+omega(3)*cos(phi));
              omega(2)*cos(phi)-omega(3)*sin(phi);
              (omega(2)*sin(phi)+omega(3)*cos(phi))/cos(theta)];
    xdot=[Vdot;omegadot;eulerdot];
end
