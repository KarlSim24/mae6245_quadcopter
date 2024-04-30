syms x y z real
syms vx vy vz real
syms ax ay az real
syms phi theta psi real %euler angles
syms phidot thetadot psidot real
syms phiddot thetaddot psiddot real 
syms m g Ixx Iyy Izz real 
syms K l b Ax Ay Az real
syms omega1 omega2 omega3 omega4 real

i = [1 0 0]';
j = [0 1 0]';
k = [0 0 1]';

% rotation derivations
R_x = [1    0       0; ...
       0  cos(phi) -sin(phi); ...
       0  sin(phi) cos(phi)];
   
R_y = [cos(theta)  0   sin(theta); ...
       0           1         0; ...
      -sin(theta) 0   cos(theta)]; 
   
R_z = [cos(psi) -sin(psi)  0; ...
       sin(psi)  cos(psi)  0; ...
       0           0       1]; 
%%%% ZYX Euler angle representation
R = R_z*R_y*R_x;

% inertial, energy, and vels
om_b = phidot*i +  R_x'*(thetadot*j) + R_x'*R_y'*(psidot*k);
I = [Ixx 0 0; 0 Iyy 0; 0 0 Izz];%body frame inertia, because quadcopter is symmetrical
v = [vx; vy; vz];
T = 0.5*m*(v')*v + 0.5*om_b'*I*om_b; % kinetic energy
V = m*g*z; % potential energy
L = T-V; %lagrangian

% derive equations of motion using Euler-Lagrange method
q = [x y z  phi theta psi ];
qdot = [vx vy vz phidot thetadot psidot ];
qddot = [ax ay az phiddot thetaddot psiddot];


% external forces and torques
Thrust = [0; 0; K*(omega1^2+omega2^2+omega3^2+omega4^2)];
Drag = [Ax*vx; Ay*vy; Az*vz];
F_ext = R*Thrust-Drag;
tau_phi = K*l*(omega4^2 - omega2^2);
tau_theta = K*l*(omega3^2 - omega1^2);
tau_psi = b*(omega1^2-omega2^2+omega3^2-omega4^2);
tau_ext = [tau_phi; tau_theta; tau_psi];

T_ext = [F_ext; tau_ext];

% solving the euler-lagrange functions for each DOF
for ii=1:6
    dLdqdot(ii) = diff(L,qdot(ii));
    ddt_dLdqdot(ii) = diff(dLdqdot(ii),q(1))*qdot(1) + diff(dLdqdot(ii),qdot(1))*qddot(1)+...
                      diff(dLdqdot(ii),q(2))*qdot(2) + diff(dLdqdot(ii),qdot(2))*qddot(2)+...
                      diff(dLdqdot(ii),q(3))*qdot(3) + diff(dLdqdot(ii),qdot(3))*qddot(3)+...
                      diff(dLdqdot(ii),q(4))*qdot(4) + diff(dLdqdot(ii),qdot(4))*qddot(4)+...
                      diff(dLdqdot(ii),q(5))*qdot(5) + diff(dLdqdot(ii),qdot(5))*qddot(5)+...
                      diff(dLdqdot(ii),q(6))*qdot(6) + diff(dLdqdot(ii),qdot(6))*qddot(6);
    dLdq(ii) = diff(L,q(ii));

    EOM(ii) = ddt_dLdqdot(ii) - dLdq(ii) - T_ext(ii);
end


% post process equations (from the UIC code)
A = jacobian(EOM,qddot);
for ii=1:6
    B(ii,1) = -subs(EOM(ii),qddot,[0 0 0 0 0 0]);
end

for ii=1:6
    for jj=1:6
        string = [ 'A(',num2str(ii),',',num2str(jj),')='];
        disp([string, char(simplify(A(ii,jj))), ';']);
    end
end

disp(' ');
for ii=1:6
    string = [ 'B(',num2str(ii),',1)='];
    disp([string, char(simplify(B(ii,1))), ';']);
end

disp('X = A\B;');
disp(' ');

% world frame velocity (used for attitude control (which we didn't get working))
om = (psidot*k) + R_z*(thetadot*j)+ R_z*R_y*(phidot*i);
R_we = jacobian(om,[phidot,thetadot,psidot]);
R_we = jacobian(om,[phidot,thetadot,psidot]);
for ii=1:3
    for jj=1:3
        string = [ 'R_we(',num2str(ii),',',num2str(jj),')='];
        disp([string, char(simplify(R_we(ii,jj))), ';']);
    end
end
disp('omega = R_we*[phidot; thetadot; psidot];');
disp('omega_x(i) = omega(1); omega_y(i) = omega(2); omega_z(i) = omega(3);');
disp(' ');

% body frame velocity (used in attitude control)
%om_b = phidot*i +  R_x'*(thetadot*j) + R_x'*R_y'*(psidot*k);
R_be = jacobian(om_b,[phidot,thetadot,psidot]);
for ii=1:3
    for jj=1:3
        string = [ 'R_be(',num2str(ii),',',num2str(jj),')='];
        disp([string, char(simplify(R_be(ii,jj))), ';']);
    end
end
disp('omega_body = R_ube*[phidot; thetadot; psidot];');
disp('omega_body_x(i) = omega_body(1); omega_body_y(i) = omega_body(2); omega_body_z(i) = omega_body(3);');

 