clear all;
close all;
clc;

%% 

addpath ../class_definitions

syms Ra Rb Rc Rd Re Rf;

Fs = 44100;
N = Fs/10;
t = 0:N-1;
f0 = 100;
trigger = zeros(length(t),1);
trigger(1) = 1;
G = 20;
%input = 20.*sin((2*pi*f0/Fs).*t);
input = trigger;

%% creation of MNA X matrix

Xmat = MNAMatrix(10,6, 0);

addResistorStamp(Xmat, Ra, 7, 8);
addResistorStamp(Xmat, Rb, 6, 8);
addResistorStamp(Xmat, Rc, 5, 9);
addResistorStamp(Xmat, Rd, 4, 10);
addResistorStamp(Xmat, Re, 3, 10);
addResistorStamp(Xmat, Rf, 2, 9);
addVoltageSourceStamp(Xmat, 7,1,1);
addVoltageSourceStamp(Xmat, 6,9,2);
addVoltageSourceStamp(Xmat, 5,10,3);
addVoltageSourceStamp(Xmat, 4,8,4);
addVoltageSourceStamp(Xmat, 3,1,5);
addVoltageSourceStamp(Xmat, 2,1,6);

%% creation of R junction and circuit

R1 = R(10);
R2 = R(10e+3);
R3 = R(10e+3);
R4 = R(4.75e+3);
C1 = C(1/(2*Fs*(1e-6)));
C2 = C(1/(2*Fs*(100e-12)));
L1 = L(2*Fs*(18e-3));
s2 = ser(C2, R2);

R_vect = [Ra, Rb, Rc, Rd, Re, Rf];
ConnectedPorts = [s2, C1, R3, R4, L1];

Rjunc = RJunction(Xmat, R_vect, ConnectedPorts, Ra);

s1 = ser(R1, Rjunc);
output = zeros(length(input), 1);

%% evolution of circuit

for i=1:N
        
    WU_s1 = WaveUp(s1);
    
    s1.WD = 2*input(i) - WU_s1;
    
    output(i) = -Voltage(R4);

end

ltspice_vect = LTspice2Matlab('wdf_Rnode.raw', 5);
t = t./Fs;
NFFT = 2^nextpow2(N);
out_fft = fft(output, NFFT);
OUT1 = abs(out_fft);
OUT1_phase = (180/pi)*angle(out_fft);
OUT = OUT1(1:NFFT/2+1);
OUT1_phase = OUT1_phase(1:NFFT/2+1);
%hi = plot(t, input, '--'); 
fig = figure();
f = Fs/NFFT*((1:(NFFT/2+1)));
f1 = (NFFT/Fs)*100;
f2 = (NFFT/Fs)*7000;
OUT_db = 20*log10(OUT);
yyaxis left;
hold on;
semilogx(f(f1:f2), OUT_db(f1:f2), 'r','DisplayName','WDF Magnitude');
semilogx(ltspice_vect.freq_vect, 20*log10(abs(ltspice_vect.variable_mat)), 'b', 'DisplayName', 'LTSpice Magnitude');
ylabel('magnitude (dB)');
%axis([1000 3000 -6 68]);
yyaxis right;
semilogx(f(f1:f2), OUT1_phase(f1:f2), 'r','DisplayName', 'WDF Phase');
%ylim([-100, 100]);
ylabel('phase');
semilogx(ltspice_vect.freq_vect, (180/pi)*angle(ltspice_vect.variable_mat),'b', 'DisplayName', 'LTSpice Phase');
hold off
%axis([1000 3000 -100 100]);
grid on;
xlabel('Frequency (Hz)');
legend;
% hold on;
% plot(t, output, 'r','DisplayName', ' WDF Voltage Over R4'); 
% plot(ltspice_vect.time_vect, ltspice_vect.variable_mat, '--b', 'DisplayName', 'LTspice Voltage Over R4');
% hold off;
% grid on;
% xlabel('Time (s)');
% ylabel('Voltage (V)');
% legend;

hgexport(fig, 'Generic Rnode 44100Hz');