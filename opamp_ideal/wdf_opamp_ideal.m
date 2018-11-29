clear all
close all
clc

%% modeling of MNA matrix

addpath ../class_definitions

syms Ra Rb Rc Rd Re Rf;

Xmat = MNAMatrix(11, 6, 1);
addResistorStamp(Xmat, Ra, 4, 8);
addResistorStamp(Xmat, Rb, 3, 9);
addResistorStamp(Xmat, Rc, 4, 10);
addResistorStamp(Xmat, Rd, 5, 11);
addResistorStamp(Xmat, Re, 2, 6);
addResistorStamp(Xmat, Rf, 4, 7);
addVoltageSourceStamp(Xmat, 8,1,1);
addVoltageSourceStamp(Xmat, 9,5,2);
addVoltageSourceStamp(Xmat, 10,5,3);
addVoltageSourceStamp(Xmat, 11,1,4);
addVoltageSourceStamp(Xmat, 6,1,5);
addVoltageSourceStamp(Xmat, 7,3,6);
addNullorStamp(Xmat, 5, 1, 2, 3, 7);

%% modeling of circuit

Fs = 44100;
R1 = R(500);
V1 = V(0,0);
R2 = R(10e+6);
RL = R(10e+3);
C1 = C(1/(2*Fs*(1e-9)));
C2 = C(1/(2*Fs*(1e-9)));

%% adaptation of R

R_vect = [Ra, Rb, Rc, Rd, Re, Rf];
R_value = [R2.PortRes, C2.PortRes, RL.PortRes, V1.PortRes, C1.PortRes];

Rjunc = RJunction(Xmat, R_vect, R_value, Ra);

S = Rjunc.S;
R_PortRes = Rjunc.PortRes;

%% WDF circuit simulation

G = 10e-3;
N = Fs/10;
t=0:N-1;

trigger = zeros(length(t),1);
trigger(1:round(Fs/10000)) = 1;
trigger(round(Fs/10000):end) = 0;
y = G.*trigger;

output = zeros(length(y), 1);

r = (R1.PortRes - R_PortRes)/(R1.PortRes+R_PortRes);

for i=1:N
    
    V1.E = y(i);
    WU_B = WaveUp(R2);
    WU_C = WaveUp(C2);
    WU_D = WaveUp(RL);
    WU_E = WaveUp(V1);
    WU_F = WaveUp(C1);
    
    WU_R = S(1,:)*[0, WU_B, WU_C, WU_D, WU_E, WU_F]';
    WD_R = r*WU_R;
    
    b = S*[WD_R, WU_B, WU_C, WU_D, WU_E, WU_F]';
    
    R2.WD = b(2);
    C2.WD = b(3);
    RL.WD = b(4);
    V1.WD = b(5);
    C1.WD = b(6);
    
    output(i) = Voltage(RL);
    
end

t_label =(1:length(t))./Fs;
NFFT = 2^nextpow2(N);
%OUT1 = fft(output, NFFT)/NFFT;
OUT1 = abs(fft(output, NFFT));
OUT = OUT1(1:NFFT/2+1);
OUT(2:end-1) = 2*OUT(2:end-1);
subplot(2,1,1);
plot(t_label, y, '--'); 
hold on;
plot(t_label, output);
hold off;
grid on;
xlabel Time(s);
ylabel Voltage(V);
subplot(2,1,2);
%plot(Fs/NFFT*((0:(NFFT/2))),pow2db(abs(OUT(1:NFFT/2+1))));
f = Fs/NFFT*((0:(NFFT/2)));
f1 = (NFFT/Fs)*1000;
f2 = (NFFT/Fs)*3000;
OUT_db = 20*log10(OUT);
plot(f(floor(f1):floor(f2)), OUT_db(floor(f1):floor(f2)));
grid on;
xlabel('Frequency (Hz)');
ylabel('magnitude (dB)');