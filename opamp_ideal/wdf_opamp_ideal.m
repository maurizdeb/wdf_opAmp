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

Fs = 48000;
R1 = R(500);
V1 = V(0,1);
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

G = 100e-3;
N = Fs/10;
t=0:N-1;

trigger = zeros(length(t),1);
trigger(1) = 1;
trigger(2:end) = 0;
y = G.*trigger;

output = zeros(length(y), 1);

r = (500-R_PortRes)/(500+R_PortRes);

for i=1:N
    
    V1.E = y(i);
    
    WU_R = WaveUp(Rjunc);
    WD_R = r*WU_R;
    Rjunc.WD = WD_R;
    
    output(i) = Voltage(RL);
    
end

trigger = zeros(length(t),1);
trigger(1:round(Fs/10000)) = 1;
trigger(round(Fs/10000):end) = 0;
y = G.*trigger;

for i=1:N
    
    V1.E = y(i);
    
    WU_R = WaveUp(Rjunc);
    WD_R = r*WU_R;
    Rjunc.WD = WD_R;
    
    output_time(i) = Voltage(RL);
    
end

x_time = LTspice2Matlab('opamp_ideal.raw', 2);
x = LTspice2Matlab('opamp_ideal_freq.raw', 2);
%plot signal in time
t_label =(1:length(t))./Fs;
NFFT = 2^nextpow2(N);
out_fft = fft(output, NFFT);
OUT1 = abs(out_fft);
OUT1_phase = (180/pi)*angle(out_fft);
OUT = OUT1(1:NFFT/2+1);
OUT1_phase = OUT1_phase(1:NFFT/2+1);
%OUT(2:end-1) = 2*OUT(2:end-1);

fig = figure();
subplot(2,1,1);
% plot(t_label, y, '--'); 
hold on;
plot(t_label, output_time, 'r', 'DisplayName', 'WDF');
plot(x_time.time_vect, x_time.variable_mat,'--b', 'DisplayName', 'LTspice');
hold off;
grid on;
xlabel Time(s);
ylabel Voltage(V);
legend;

%plot fft of signal
subplot(2,1,2);
f = Fs/NFFT*((0:(NFFT/2)));
f1 = (NFFT/Fs)*1000;
f2 = (NFFT/Fs)*3000;
OUT_db = 20*log10(OUT);
yyaxis left;
hold on;
semilogx(f(floor(f1):floor(f2)), OUT_db(floor(f1):floor(f2)), 'r','DisplayName','WDF Magnitude');
semilogx(x.freq_vect, 20*log10(abs(x.variable_mat)), 'b', 'DisplayName', 'LTSpice Magnitude');
ylabel('magnitude (dB)');
axis([1000 3000 -6 68]);
yyaxis right;
semilogx(f(floor(f1):floor(f2)), OUT1_phase(floor(f1):floor(f2)), 'r','DisplayName', 'WDF Phase');
%ylim([-100, 100]);
ylabel('phase');
semilogx(x.freq_vect, (180/pi)*angle(x.variable_mat),'b', 'DisplayName', 'LTSpice Phase');
hold off
axis([1000 3000 -100 100]);
grid on;
xlabel('Frequency (Hz)');
legend;

hgexport(fig, 'Bridged T result 100000Hz');