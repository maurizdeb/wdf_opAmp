clear all
close all
clc

%% creation of MNA X matrix

addpath ../class_definitions;
syms Ra Rb Rc Rd Re Rf Rg Rh Ri Rj Rk;

Xmat = MNAMatrix(21, 11, 4);
%adding resistor stamps
addResistorStamp(Xmat, Ra, 2, 11);
addResistorStamp(Xmat, Rb, 3, 12);
addResistorStamp(Xmat, Rc, 4, 13);
addResistorStamp(Xmat, Rd, 4, 14);
addResistorStamp(Xmat, Re, 6, 15);
addResistorStamp(Xmat, Rf, 1, 16);
addResistorStamp(Xmat, Rg, 8, 17);
addResistorStamp(Xmat, Rh, 7, 18);
addResistorStamp(Xmat, Ri, 7, 19);
addResistorStamp(Xmat, Rj, 7, 20);
addResistorStamp(Xmat, Rk, 4, 21);
%adding voltage sources stamps
addVoltageSourceStamp(Xmat, 11, 1, 1);
addVoltageSourceStamp(Xmat, 12, 1, 2);
addVoltageSourceStamp(Xmat, 13, 1, 3);
addVoltageSourceStamp(Xmat, 14, 3, 4);
addVoltageSourceStamp(Xmat, 15, 5, 5);
addVoltageSourceStamp(Xmat, 16, 5, 6);
addVoltageSourceStamp(Xmat, 17, 7, 7);
addVoltageSourceStamp(Xmat, 18, 1, 8);
addVoltageSourceStamp(Xmat, 19, 2, 9);
addVoltageSourceStamp(Xmat, 20, 4, 10);
addVoltageSourceStamp(Xmat, 21, 2, 11);
%adding VCVS stamps
addVCVSStamp(Xmat, 3, 1, 6, 10, 12, 3.1625);
addVCVSStamp(Xmat, 3, 4, 10, 9, 13, (200e+3));
addVCVSStamp(Xmat, 4, 1, 9, 1, 14, 3.1625);
addVCVSStamp(Xmat, 5, 1, 8, 1, 15, 1);

%% modeling of junctions and ports of circuit
Fs = 44100;

R1 = R(500);
Ib1 = I(90e-9, 2*(5e+6));
Ccm1 = C(1/(4*Fs*(2e-12)));
Ccm2 = C(1/(4*Fs*(2e-12)));
Vin = V(0,0);
Voff = 1e-3;
p1 = par(Vin, par(Ib1, Ccm1));
Ib2 = I(70e-9, 2*(5e+6));
p2 = par(Ib2, Ccm2);
Rid = R(5e+6);
Cid = C(1/(2*Fs*(1.4e-12)));
p3 = par(Rid,Cid);
Rbw = R(100e+3);
Cbw = C(1/(2*Fs*(0.3183e-6)));
Rout = R(75);
RL = R(10e+3);
C2 = C(1/(2*Fs*(1e-9)));
R2 = R(10e+6);
C1 = C(1/(2*Fs*(1e-9)));

%% adaptation of R junction

R_vect = [Ra, Rb, Rc, Rd, Re, Rf, Rg, Rh, Ri, Rj, Rk];
ConnectedPorts = [p1, p2, p3, Rbw, Cbw, Rout, RL, C2, R2, C1];

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
    
    Vin.E = y(i)-Voff;
    
    WU_R = WaveUp(Rjunc);
    Rjunc.WD = r*WU_R;
    
    
    output(i) = Voltage(RL);
end

t_label =(1:length(t))./Fs;
NFFT = 2^nextpow2(N);
OUT1 = abs(fft(output, NFFT));
OUT = OUT1(1:NFFT/2+1);
OUT(2:end-1) = 2*OUT(2:end-1);

%plot time signal
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