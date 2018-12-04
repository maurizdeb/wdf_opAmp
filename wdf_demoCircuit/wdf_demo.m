%%
clear all
close all
clc
%%
% WDF Linear Circuit Example

addpath ../class_definitions;
Fs = 48000;
N = Fs/10;
gain = 30;
f0 = 100;
t = 0:N-1;
input = gain.*sin(2*pi*f0/Fs.*t);
output = zeros(1, length(input));
V1 = V(0,1); % voltage source with 1 Ohm series resistance
R2 = R(100e+3); % 1000 kOhm resistor
R3 = R(80);
CapVal = 3.5e-5;
C1 = C(1/(2*Fs*CapVal));
IndVal = 18000e-6;
L1 = L(2*Fs*IndVal);

p1 = par(R2, L1);
s1 = ser(V1, p1);
s2 = ser(s1,R3);
%Vdiode = 0; initial value for the value of voltage over the diode

cStateIn = 0;
cStateOut = 0;
cIn = 0;
cOut = 0;
p = ((1/Fs)-(2*s2.PortRes*CapVal))/((1/Fs)+(2*s2.PortRes*CapVal));

for i=1:N
    V1.E = input(i);
    WaveUp(s2);
    cIn = s2.WU;
    cOut = p*(cIn-cStateOut)+cStateIn;
    s2.WD = cOut;
    cStateIn = cIn;
    cStateOut = cOut;
    output(i) = Voltage(s2);
end

% plots
t = (1:length(input))./Fs;
hi = plot(t, input, '--'); hold on;
ho = plot(t, output); hold off;
grid on;
xlabel('Time (s)');
ylabel('Voltage (V)');
legend([hi, ho], 'Source Voltage E', 'Voltage Over R3');