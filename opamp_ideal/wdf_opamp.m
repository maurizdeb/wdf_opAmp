clear all
close all
clc

%%

addpath ../class_definitions;

syms Ra Rb Rc Rd Re Rf;

Fs = 44100;
G = 10e-3;
f0 = 100;
N = Fs/10;
t=0:N-1;
t_label = (1:length(t))./Fs;

trigger = zeros(length(t),1);
trigger(1:round(Fs/10000)) = 1;
trigger(round(Fs/10000):end) = 0;
y = G.*trigger;

R1 = R(500);
R1.Conductance;
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

%WERNER'S MATRIX
% X =[1/Rc, 0, 0, 0, 0, 0, 0, -1/Rc, 0, 0, 0, -1, 0, 0, 0, 0, -1, -1;
%     0, 1/Ra, 0, 0, 0, -1/Ra, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
%     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, -1, 0, 0, 0;
%     0, 0, 0, 1/Rb, 0, 0, -1/Rb, 0, 0, 0, 0, 0, 0, -1, 0, -1, 0, 0;
%     0, 0, 0, 0, 1/Rd + 1/Re + 1/Rf, 0, 0, 0, -1/Rd, -1/Re, -1/Rf, 0, 0, 0, 0, 0, 0, 1;
%     0, -1/Ra, 0, 0, 0, 1/Ra, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0;
%     0, 0, 0, -1/Rb, 0, 0, 1/Rb, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0;
%     -1/Rc, 0, 0, 0, 0, 0, 0, 1/Rc, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;
%     0, 0, 0, 0, -1/Rd, 0, 0, 0, 1/Rd, 0, 0, 0, 0, 0, 1, 0, 0, 0;
%     0, 0, 0, 0, -1/Re, 0, 0, 0, 0, 1/Re, 0, 0, 0, 0, 0, 1, 0, 0;
%     0, 0, 0, 0, -1/Rf, 0, 0, 0, 0, 0, 1/Rf, 0, 0, 0, 0, 0, 1, 0;
%     -1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
%     0, 0, -1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
%     0, 0, 0, -1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
%     0, 0, -1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0;
%     0, 0, 0, -1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0;
%     -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0;
%     0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];

%ALBERTINI'S MATRIX
X =[1/Ra + 1/Rf, 0, 0, 0, 0, -1/Ra, 0, 0, 0, 0, -1/Rf, 0, 0, 0, 0, -1, 0, -1;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, -1, 0, -1, 1;
    0, 0, 1/Rb + 1/Rc, 0, 0, 0, -1/Rb, -1/Rc, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 1/Rd + 1/Re, 0, 0, 0, -1/Rd, -1/Re, 0, 0, 0, -1, 0, 0, 0, 0;
    -1/Ra, 0, 0, 0, 0, 1/Ra, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0;
    0, 0, -1/Rb, 0, 0, 0, 1/Rb, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0;
    0, 0, -1/Rc, 0, 0, 0, 0, 1/Rc, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;
    0, 0, 0, 0, -1/Rd, 0, 0, 0, 1/Rd, 0, 0, 0, 0, 0, 1, 0, 0, 0;
    0, 0, 0, 0, -1/Re, 0, 0, 0, 0, 1/Re, 0, 0, 0, 0, 0, 1, 0, 0;
    -1/Rf, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1/Rf, 0, 0, 0, 0, 0, 1, 0;
    0, 0, 0, -1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, -1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, -1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, -1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    -1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0;
    0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0;
    0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];


X = X(2:end, 2:end);
X_inv = inv(X);

left_term = 2*[zeros(6, 10), diag([Ra, Rb, Rc, Rd, Re, Rf]), zeros(6,1)];
right_term = [zeros(6,10), eye(6), zeros(6,1)];

S = eye(6) + ((left_term)*X_inv)*(right_term.');

[R_PortRes, param, cond] = solve(S(5,5) == 0, Re, 'ReturnConditions', true);

assume(cond);
R_PortRes = double(subs(R_PortRes, [Ra, Rb, Rc, Rd, Rf],[V1.PortRes, R2.PortRes, C1.PortRes, C2.PortRes, RL.PortRes]));
S = double(subs(S, [Ra, Rb, Rc, Rd, Re, Rf], [V1.PortRes, R2.PortRes, C1.PortRes, C2.PortRes, R_PortRes, RL.PortRes]));


%%
r = (R1.PortRes - R_PortRes)/(R1.PortRes+R_PortRes);

for i=1:N
    
    V1.E = y(i);
    WU_A = WaveUp(V1);
    WU_B = WaveUp(R2);
    WU_C = WaveUp(C1);
    WU_D = WaveUp(C2);
    
    WU_F = WaveUp(RL);
    
    WU_R = S(5,:)*[WU_A, WU_B, WU_C, WU_D, 0, WU_F]';
    WD_R = r*WU_R;
    
    b = S*[WU_A, WU_B, WU_C, WU_D, WD_R, WU_F]';
    
    V1.WD = b(1);
    R2.WD = b(2);
    C1.WD = b(3);
    C2.WD = b(4);
    %R1.WD = b(3);
    RL.WD = b(6);
    
    output(i) = Voltage(RL);
end

NFFT = 2^nextpow2(N);
Y = abs(fft(output, NFFT));
Y = Y(1:NFFT/2+1);
%OUT(2:end-1) = OUT(2:end-1);

subplot(2,1,1);
plot(t_label, y, '--'); 
hold on;
plot(t_label, output);
hold off;
grid on;
xlabel('Time(s)');
ylabel('Voltage(V)');

subplot(2,1,2);
%plot(Fs/NFFT*((0:(NFFT/2))),pow2db(abs(OUT(1:NFFT/2+1))));
plot(Fs/NFFT*((0:(NFFT/2))), 20*log10(Y));
grid on;
xlabel('Frequency [Hz]');
ylabel('Magnitude [dB]');
