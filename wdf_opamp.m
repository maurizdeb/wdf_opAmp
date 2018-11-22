clear all
close all
clc

%%

addpath class_definitions;

syms Ra Rb Rc Rd Re Rf;

Fs = 44100;
G = 10e-3;
N = Fs/10;
t=0:N-1;
t_label = (1:length(t))./Fs;

trigger = zeros(length(t),1);
trigger(1:round(Fs/1000)) = 1;
trigger(round(Fs/1000)+1:end) = 0;
y = G*trigger;

%R1 = R(500);
%R1.Conductance;
V1 = V(0,0);
V1.Conductance;
R2 = R(10e+6);
R2.Conductance;
RL = R(10e+3);
RL.Conductance;
C1 = C(1/(2*Fs*(1e-9)));
C1.Conductance;
C2 = C(1/(2*Fs*(1e-9)));
C2.Conductance;


output = zeros(length(y), 1);

X =[0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0;
    0, 1/Rb + 1/Rc, 0, 0, 0,0, -1/Rb, -1/Rc, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 1/Rd+(1/Re), 0,0, 0, 0, 0, -1/Rd, -1/Re, 0, 0, 0, -1, 0, 0, 0, 0;
    0, 0, 0, 1/Rf, 0, 0, 0, 0, 0, 0, -1/Rf, 0, -1, 0, -1, 0, 0, 1;
    0, 0, 0, 0, 1/Ra, -1/Ra, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, -1;
    0, 0, 0, 0,-1/Ra, 1/Ra, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0;
    0, -1/Rb, 0, 0, 0, 0, 1/Rb, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0;
    0, -1/Rc, 0, 0, 0, 0, 0, 1/Rc, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;
    0, 0, -1/Rd, 0, 0, 0, 0, 1/Rd, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0;
    0, 0, -1/Re, 0, 0, 0, 0, 0, 0, 1/Re, 0, 0, 0, 0, 0, 1, 0, 0;
    0, 0, 0, -1/Rf,0, 0, 0, 0, 0, 0, 1/Rf, 0, 0, 0, 0, 0, 1, 0;
    -1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, -1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, -1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, -1,0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, -1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0;
    1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];


X = X(2:end, 2:end);
X_inv = inv(X);



left_term = 2*[zeros(6, 10), diag([Ra, Rb, Rc, Rd, Re, Rf]), zeros(6,1)];
right_term = [zeros(6,10), eye(6), zeros(6,1)];

S = eye(6) + ((left_term)*X_inv)*(right_term.');

R_PortRes = solve(S(5,5) == 0, Re);


%%
double(subs(S, R_PortRes));

for i=1:N
    WU_B = WaveUp(R2);
    WU_C = WaveUp(C1);
    WU_D = WaveUp(C2);
    WU_E = WaveUp(R1);
    WU_F = WaveUp(RL);
    
    WU_R = S(1,:)*[0, WU_B, WU_C, WU_D, WU_E, WU_F]';
    
    
end

