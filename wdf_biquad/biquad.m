clear all
close all
clc

%% modeling of MNA matrix (adaptation on Rb)

addpath ../class_definitions

syms Ra Rb Rc Rd Re Re Rf Rg Rh Ri;

Xmat = MNAMatrix(17, 9, 3);
addResistorStamp(Xmat, Ra, 9, 1);
addResistorStamp(Xmat, Rb, 10, 7);
addResistorStamp(Xmat, Rc, 11, 2);
addResistorStamp(Xmat, Rd, 12, 3);
addResistorStamp(Xmat, Re, 13, 4);
addResistorStamp(Xmat, Rf, 14, 5);
addResistorStamp(Xmat, Rg, 15, 7);
addResistorStamp(Xmat, Rh, 16, 4);
addResistorStamp(Xmat, Ri, 17, 1);
addVoltageSourceStamp(Xmat, 9, 6, 1);
addVoltageSourceStamp(Xmat, 10, 6, 2);
addVoltageSourceStamp(Xmat, 11, 6, 3);
addVoltageSourceStamp(Xmat, 12, 2, 4);
addVoltageSourceStamp(Xmat, 13, 3, 5);
addVoltageSourceStamp(Xmat, 14, 4, 6);
addVoltageSourceStamp(Xmat, 15, 5, 7);
addVoltageSourceStamp(Xmat, 16, 8, 8);
addVoltageSourceStamp(Xmat, 17, 8, 9);
addNullorStamp(Xmat, 2, 1, 8, 6, 10);
addNullorStamp(Xmat, 4, 1, 1, 3, 11);
addNullorStamp(Xmat, 7, 1, 1, 5, 12);

%% modeling of circuit

Fs = 48000;
V1 = V(0, 10000);
R7 = R(10000);
R4 = R(10000);
R5 = R(16000);
C1 = C(1/(2*Fs*(10e-9)));
R6 = R(16000);
C2 = C(1/(2*Fs*(10e-9)));
R3 = R(10000);
R2 = R(500);

%% adaptation of R

R_vect = [Ra, Rb, Rc, Rd, Re, Rf, Rg, Rh, Ri];
ConnectedPorts = [V1, R7, R5, C1, R6, C2, R3, R2]; %all connected ports without R7 that is the adaptation port

Rjunc = RJunction(Xmat, R_vect, ConnectedPorts, Rc);

S = Rjunc.S; %S scattering matrix
R_PortRes = Rjunc.PortRes; % adapted resistance
%% circuit simulation

N = Fs/10;
t = 0:N-1;
input = zeros(length(t), 1);
input(1) = 1;
r = (10000-R_PortRes)/(10000+R_PortRes);

output_LP = zeros(length(t),1);
output_BP = zeros(length(t),1);
output_HP = zeros(length(t),1);

for i=1:N
    
    V1.E = input(i);
    
    WU_R = WaveUp(Rjunc);
    WD_R = r*WU_R;
    Rjunc.WD = WD_R;
    
    output_BP(i) = Voltage(R2) - Voltage(R3);
    output_LP(i) = output_BP(i) - Voltage(R6) - Voltage(C2);
    output_HP(i) = output_BP(i) + Voltage(C1) + Voltage(R5);
    
end

% fft calculation
NFFT = 2^nextpow2(N);
OUT_BP = fft(output_BP, NFFT);
OUT_LP = fft(output_LP, NFFT);
OUT_HP = fft(output_HP, NFFT);
OUT_BP_phase = (180/pi)*angle(OUT_BP);
OUT_LP_phase = (180/pi)*angle(OUT_LP);
OUT_HP_phase = (180/pi)*angle(OUT_HP);
OUT_BP = abs(OUT_BP);
OUT_LP = abs(OUT_LP);
OUT_HP = abs(OUT_HP);

OUT_BP = OUT_BP(1:NFFT/2);
OUT_BP_phase = OUT_BP_phase(1:NFFT/2);
OUT_LP = OUT_LP(1:NFFT/2);
OUT_LP_phase = OUT_LP_phase(1:NFFT/2);
OUT_HP = OUT_HP(1:NFFT/2);
OUT_HP_phase = OUT_HP_phase(1:NFFT/2);

BP_LTspice = LTspice2Matlab('biquad.raw', 4);
LP_LTspice = LTspice2Matlab('biquad.raw', 6);
HP_LTspice = LTspice2Matlab('biquad.raw', 1);

%% plots

fig = figure();
subplot(3,1,1);
f = Fs/NFFT*((0:(NFFT/2-1)));
yyaxis left;
semilogx(f, 20*log10(OUT_BP), 'r','DisplayName','WDF Bandpass Output');
hold on;
semilogx(BP_LTspice.freq_vect, 20*log10(abs(BP_LTspice.variable_mat)), 'b', 'DisplayName', 'LTSpice Bandpass Output');
hold off;
ylabel('magnitude (dB)');
axis([10 10000 -40 20]);
yyaxis right;
semilogx(f, OUT_BP_phase, 'm','DisplayName', 'WDF Bandpass Phase');
%ylim([-100, 100]);
ylabel('phase');
hold on;
semilogx(BP_LTspice.freq_vect, (180/pi)*angle(BP_LTspice.variable_mat),'g', 'DisplayName', 'LTSpice Bandpass Phase');
hold off;
axis([10 10000 -200 200]);
grid on;
xlabel('Frequency (Hz)');
legend;

subplot(3,1,2);
yyaxis left;
semilogx(f, 20*log10(OUT_LP), 'r','DisplayName','WDF Lowpass Output');
hold on;
semilogx(LP_LTspice.freq_vect, 20*log10(abs(LP_LTspice.variable_mat)), 'b', 'DisplayName', 'LTSpice Lowpass Output');
hold off;
ylabel('magnitude (dB)');
axis([10 10000 -40 20]);
yyaxis right;
semilogx(f, OUT_LP_phase, 'm','DisplayName', 'WDF Lowpass Phase');
%ylim([-100, 100]);
ylabel('phase');
hold on;
semilogx(LP_LTspice.freq_vect, (180/pi)*angle(LP_LTspice.variable_mat),'g', 'DisplayName', 'LTSpice Lowpass Phase');
hold off;
axis([10 10000 -200 200]);
grid on;
xlabel('Frequency (Hz)');
legend;

subplot(3,1,3);
yyaxis left;
semilogx(f, 20*log10(OUT_HP), 'r','DisplayName','WDF Highpass Output');
hold on;
semilogx(HP_LTspice.freq_vect, 20*log10(abs(HP_LTspice.variable_mat)), 'b', 'DisplayName', 'LTSpice Lowpass Output');
hold off;
ylabel('magnitude (dB)');
axis([10 10000 -40 20]);
yyaxis right;
semilogx(f, OUT_HP_phase, 'm','DisplayName', 'WDF Highpass Phase');
%ylim([-100, 100]);
ylabel('phase');
hold on;
semilogx(HP_LTspice.freq_vect, (180/pi)*angle(HP_LTspice.variable_mat),'g', 'DisplayName', 'LTSpice Lowpass Phase');
hold off;
axis([10 10000 -200 200]);
grid on;
xlabel('Frequency (Hz)');
legend;

%hgexport(fig, 'Biquad State Variable');