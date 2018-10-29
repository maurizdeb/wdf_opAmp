%%
clear all
close all
clc
%%
addpath class_definitions

Fs = 48000;
N = Fs/10;
gain = 30;
f0 = 100;
t = 0:N-1;
input = gain.*sin(2*pi*f0/Fs.*t);
output = zeros(1, length(input));
R1 = R(953);
CapVal = 2.2e-6;
C1 = C(1/(2*Fs*CapVal));
s1 = ser(R1,C1);
V1 = V(0,0);

s2 = ser(V1,C1);

for i=1:N
    
    V1.E = input(i);
    WaveUp(s2);
    r = (R1.PortRes-s2.PortRes)/(R1.PortRes+s2.PortRes);
    %s1.PortRes;
    
    s2.WD = r*s2.WU;
    
    
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