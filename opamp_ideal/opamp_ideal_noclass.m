clear all
close all
clc
%%

syms Rc;
%%
Fs = 48000;
V1_Res = 1;
C1 = 1/(2*Fs*(1e-9));
R1 = 500;
R2 = 10e+6;
C2 = 1/(2*Fs*(1e-9));
RL = 10e+3;

%%

X = [1/Rc, 0, 0, 0, 0, 0, 0, -1/Rc, 0, 0, 0, -1, 0, 0, 0, 0, -1, -1;
    0, (1/V1_Res), 0, 0, 0, -(1/V1_Res), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, -1, 0, 0, 0;
    0, 0, 0, (1/C1), 0, 0, -(1/C1), 0, 0, 0, 0, 0, 0, -1, 0, -1, 0, 0;
    0, 0, 0, 0, (1/R2)+(1/C2)+(1/RL), 0, 0, 0, -1/R2, -1/C2, -1/RL, 0, 0, 0, 0, 0, 0, 1;
    0, -(1/V1_Res), 0, 0, 0, (1/V1_Res), 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0;
    0, 0, 0, -(1/C1), 0, 0, (1/C1), 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0;
    -(1/Rc), 0, 0, 0, 0, 0, 0, (1/Rc), 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;
    0, 0, 0, 0, -(1/R2), 0, 0, 0, (1/R2), 0, 0, 0, 0, 0, 1, 0, 0, 0;
    0, 0, 0, 0, -(1/C2), 0, 0, 0, 0, (1/C2), 0, 0, 0, 0, 0, 1, 0, 0;
    0, 0, 0, 0, -(1/RL), 0, 0, 0, 0, 0, (1/RL), 0, 0, 0, 0, 0, 1, 0;
    -1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, -1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, -1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, -1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, -1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0;
    -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0;
    0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];

%%

X = X(2:end,2:end);
X_inv = inv(X);
left_term = 2*[zeros(6, 10), diag([V1_Res, C1, Rc, R2, C2, RL]), zeros(6,1)];
right_term = [zeros(6, 10), eye(6), zeros(6,1)];
S = eye(6) + left_term*(X_inv*(right_term'));
R_PortRes = double(solve(S(3,3)==0));
S = double(subs(S, Rc, R_PortRes));

%% Frequency analysis

G = 100e-3;
N = Fs/10;
t=0:N-1;

trigger = zeros(length(t),1);
trigger(1) = 1;
trigger(2:end) = 0;
y = G.*trigger;

output = zeros(length(y), 1);

r = (R1-R_PortRes)/(R1+R_PortRes);
C1_WD = 0;
C2_WD = 0;


for i=1:N
    
    Wu_V1 = y(i);
    Wu_C1 = C1_WD;
    Wu_R2 = 0;
    Wu_C2 = C2_WD;
    Wu_RL = 0;
    WU_R = S(3,:)*[Wu_V1, Wu_C1, 0, Wu_R2, Wu_C2, Wu_RL]';
    WD_R = r*WU_R;
    b = S*[Wu_V1, Wu_C1, WD_R, Wu_R2, Wu_C2, Wu_RL]';
    
    V1_WD = b(1);
    C1_WD = b(2);
    R2_WD = b(4);
    C2_WD = b(5);
    RL_WD = b(6);
    
    output(i) = (Wu_RL+RL_WD)/2;
    
end

%% Time analysis

trigger = zeros(length(t),1);
trigger(1:end) = 0;
trigger(1:floor(Fs/10000)+1) = 1;
y = G.*trigger;

C1_WD = 0;
C2_WD = 0;
output_time = zeros(length(y), 1);
for i=1:N
    
    Wu_V1 = y(i);
    Wu_C1 = C1_WD;
    Wu_R2 = 0;
    Wu_C2 = C2_WD;
    Wu_RL = 0;
    WU_R = S(3,:)*[Wu_V1, Wu_C1, 0, Wu_R2, Wu_C2, Wu_RL]';
    WD_R = r*WU_R;
    b = S*[Wu_V1, Wu_C1, WD_R, Wu_R2, Wu_C2, Wu_RL]';
   
    V1_WD = b(1);
    C1_WD = b(2);
    R2_WD = b(4);
    C2_WD = b(5);
    RL_WD = b(6);
    
    output_time(i) = (RL_WD+Wu_RL)/2;
end

%% Plots

x_time = LTspice2Matlab('opamp_ideal_time.raw', 2);
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
[maxOUT, max_i] = max(OUT_db);
max_i = (Fs/NFFT)*(max_i-1);
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

hgexport(fig, 'Bridged T result 48000Hz without class implem');