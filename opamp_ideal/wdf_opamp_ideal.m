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
ConnectedPorts = [R2, C2, RL, V1, C1];

Rjunc = RJunction(Xmat, R_vect, ConnectedPorts, Ra);

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
    
    WU_R = WaveUp(Rjunc);
    WD_R = r*WU_R;
    Rjunc.WD = WD_R;
    
    output(i) = Voltage(RL);
    
end

%plot signal in time
t_label =(1:length(t))./Fs;
NFFT = 2^nextpow2(N);
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

%plot fft of signal
subplot(2,1,2);
f = Fs/NFFT*((0:(NFFT/2)));
f1 = (NFFT/Fs)*1000;
f2 = (NFFT/Fs)*3000;
OUT_db = 20*log10(OUT);
plot(f(floor(f1):floor(f2)), OUT_db(floor(f1):floor(f2)));
axis([1000 3000 0 max(OUT_db(floor(f1):floor(f2)))+10]);
grid on;
xlabel('Frequency (Hz)');
ylabel('magnitude (dB)');