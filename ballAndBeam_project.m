%% 
clear all
close all
%% Plant of the SS

% Define the physical values
g = 9.8;            %Gravity acceleration
m_B = 0.064;        %Ball mass
m_b = 0.65;         %Beam mass
R = 0.0254;         %Ball radius
L = 0.425;          %Beam Length
d = 0.12;           %Lever length
delta_2 = 0.2;      %Equilibrium Point of ball position
K_m = 0.00767;      %Back EMF constant
K_i = 0.00767;      %Torque Constant
K_g = 14;           %Gear ratio
R_m = 2.6;          %Motor resistance
n_motor = 0.69;     %Motor efficiency
n_gearbox = 0.85;   %Gearbox efficiency

n_total = n_motor + n_gearbox;  %Total efficiency
J_b = (1/2)*m_b*L^2;            %Beam moment of inertia 

% Define the A matrix for SS
A = [0 0 1 0 ;
    0 0 0 1 ;
    0 -(((m_B*J_b*g) + (m_B^2*g*delta_2^2))/(J_b+m_B*delta_2^2)^2) -((K_g^2*K_i*K_m*n_total)/(R_m*(J_b+m_B*delta_2^2))*(L^2/d^2)) 0;
    -(5*g/7) 0 0 0];

% Define the B matrix for SS
B = [0 ;
    0;
    ((K_g*K_i*n_total)/(R_m*(J_b+m_B*delta_2^2)))*(L/d);
    0];

% Define the C matrix for SS
C = [0 1 0 0];

% Define the D matrix for SS
D = 0;

%Output perturbance
OP = 0.1;
%Input perturbance
IP = 0.8;
%parameter uncertainty
S = 1.1;

% Define the system order
n = size(A,1);
%% Plant Without Control
PWC = ss(A,B,C,D);
figure(1), step(PWC,10);
%% Pole Placement Control
% Design control for OS=5% and Ts=1
OS=4;
Ts=0.7;
% Define the dominant poles
E = abs(log(OS/100))/sqrt(pi^2+(log(OS/100))^2);
Wn = 4/(E*Ts);

polescl1 = -E*Wn+Wn*sqrt(E^2-1);
polescl2 = -E*Wn-Wn*sqrt(E^2-1);
polescl3 = real(polescl1)*10;
polescl4 = polescl3+1;

Polescl = [polescl1 polescl2 polescl3 polescl4];
K1 = place(A,B,Polescl);

% Determine the precompensation gain
Kpre1 = inv(-(C-D*K1)*inv(A-B*K1)*B + D);
%% ITAE Control
% Design control for Wn = 1 rad/s
% s^4 + 2.1*Wn*s^3 + 3.4*Wn^2*s^2 + 2.7*Wn^3*s + Wn^4
% Define the dominant poles
Polescl = (roots([1 2.1*Wn 3.4*Wn^2 2.7*Wn^3 Wn^4]))';
K2 = place(A,B,Polescl);

% Determine the precompensation gain
Kpre2 = inv(-(C-D*K2)*inv(A-B*K2)*B + D);
%% LQR Control

Q =  [10   0     0     0;
      0     1     0    0;
      0     0     1    0;
      0     0     0    1.23];
R = 1.50;

[K3, S, e] = lqr(A,B,Q,R);

% Determine the precompensation gain
Kpre3 = inv(-(C-D*K3)*inv(A-B*K3)*B + D);

% PUC3
Acl3 = (A-B*K3);
PUC3 = ss(Acl3,B,C,D);
figure(3) 
[y1,t] = step(PUC3*Kpre3,10); grid on
stepResults = stepinfo(y1,t);
stepResults.SettlingTime
stepResults.Overshoot

%% State Observers
% define the observer poles 10 times the placed poles
PolesObs = Polescl(1:n)*10;
% find the L-values
L = place(A',C', PolesObs)';
% Define the Observer Space State
Aobs = A - L*C;
Bobs = [B-L*D L];
Cobs = eye(n);
Dobs = zeros(n,2);
% Observer + Controller Space State
Aoc = A-L*C-B*K1+L*D*K1;
Boc = [B-L*D L];
Coc = -K1;
Doc = [1 0];
%% Static Error
% Feedback state control + Static Error
Aext = [A zeros(n,1) ; -C 0];
Bext = [B ; -D];
Cext = [C 0];
Dext = 0;
% Poles placement
polescl5 = polescl4+1;
Polescl = [polescl1 polescl2 polescl3 polescl4 polescl5];
Kext = place(Aext,Bext,Polescl);
K4 = Kext(1:n);
Kpre4 = Kext(n+1);
% Observer + Controller Space State + Static Control
Aocc = [A-L*C-B*K4+L*D*K4 -B*Kpre4+L*D*Kpre4; zeros(1,n) 0];
Bocc = [zeros(n,1) L; 1 -1];
Cocc = [-K4 -Kpre4];
Docc = zeros(1,2);