clear all;
close all;
clc;

%% 

addpath ../class_definitions

syms Ra;

Fs = 44100;
N = Fs/10;
t = 0:N-1;
f0 = 100;
input = 20.*sin((2*pi*f0/Fs).*t);
R1 = R(10);
R2 = R(10e+3);
R3 = R(10e+3);
R4 = R(4.75e+3);
C1 = C(1/(2*Fs*(1e-6)));
C2 = C(1/(2*Fs*(100e-12)));
L1 = L(2*Fs*(18e-3));

s2 = ser(C2, R2);

X = [L1.Conductance, 0, 0, 0, 0, 0, 0, -L1.Conductance, 0, 0, 0, 0, 0, 0, 1;
    0, R4.Conductance, 0, 0, 0, 0, 0, 0, -R4.Conductance, 0, 0, 0, 0, 1, 0;
    0, 0, R3.Conductance, 0, 0, 0, 0, 0, -R3.Conductance, 0, 0, 0, 1, 0, 0;
    0, 0, 0, C1.Conductance, 0, 0, 0, -C1.Conductance, 0, 0, 0, 1, 0, 0, 0;
    0, 0, 0, 0, s2.Conductance, 0, -s2.Conductance, 0, 0, 0, 1, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 1/Ra, -1/Ra, 0, 0, 1, 0, 0, 0, 0, 0;
    0, 0, 0, 0, -s2.Conductance, -1/Ra, 1/Ra + s2.Conductance, 0, 0, 0, 0, 0, -1, 0, 0;
    -L1.Conductance, 0, 0, -C1.Conductance, 0, 0, 0, L1.Conductance + C1.Conductance, 0, 0, -1, 0, 0, 0, 0;
    0, -R4.Conductance, -R3.Conductance, 0, 0, 0, 0, 0, R3.Conductance+R4.Conductance, 0, 0, -1, 0, 0, 0;
    0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 1, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0;
    0, 0, 1, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];

X_inv = inv(X);

R_diag=diag([Ra, s2.PortRes, C1.PortRes, R3.PortRes, R4.PortRes, L1.PortRes]);
first_term = 2*[zeros(6, 9), R_diag ];
second_term = [zeros(6,9), eye(6)].';
S = eye(6) + (first_term*X_inv)*second_term;

R_PortRes = double(solve(S(1,1) == 0));

S = double(subs(S, R_PortRes));
S(1,1) = 0;

output = zeros(length(input), 1);

for i=1:N
    %Up_s1 = -input(i) - WaveUp(R1);
    WU_B = WaveUp(s2);
    WU_C = WaveUp(C1);
    WU_D = WaveUp(R3);
    WU_E = WaveUp(R4);
    WU_F = WaveUp(L1);
    
    WU_R = S(1,:)*[0 WU_B WU_C WU_D WU_E WU_F]';
        
    WU_s1 = -WU_R - WaveUp(R1);
    WD_s1 = 2*input(i) - WU_s1;
    s1_PortRes = R_PortRes + R1.PortRes;
    s1_refCoeff = R1.PortRes/s1_PortRes;
    R1.WD = (R1.WU+WU_R+WD_s1)*(-s1_refCoeff)+R1.WU;
    WD_R = -R1.WU-WD_s1 + s1_refCoeff*(WD_s1 + R1.WU+WU_R);
    
    
    b = S*[WD_R WU_B WU_C WU_D WU_E WU_F]';
    
    %WU_R = b(1);
    
    s2.WD = b(2);
    C1.WD = b(3);
    R3.WD = b(4);
    R4.WD = b(5);
    L1.WD = b(6);
    
    output(i) = -Voltage(R4);

end

t = (1:length(input))./Fs;
hi = plot(t, input, '--'); hold on;
ho = plot(t, output); hold off;
grid on;
xlabel('Time (s)');
ylabel('Voltage (V)');
legend([hi, ho], 'Source Voltage E', 'Voltage Over R4');