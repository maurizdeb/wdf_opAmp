clear all
close all
clc

%% creation of MNA X matrix

addpath ../class_definitions;
syms Ra Rb Rc Rd Re Rf Rg Rh Ri Rj Rk

XA = sym(zeros(21,21));
XB = sym(zeros(15,21));
XC = sym(zeros(21,15));
XD = sym(zeros(15,15));
XA([2,11],[2,11]) = XA([2,11],[2,11]) + [1/Ra, -1/Ra; -1/Ra, 1/Ra];
XA([3, 12],[3, 12]) = XA([3,12],[3,12]) + [1/Rb, -1/Rb; -1/Rb, 1/Rb];
XA([4, 13], [4, 13]) = XA([4, 13], [4, 13]) + [1/Rc, -1/Rc; -1/Rc, 1/Rc];
XA([4, 14], [4, 14]) = XA([4, 14], [4, 14]) + [1/Rd, -1/Rd; -1/Rd, 1/Rd];
XA([6, 15], [6, 15]) = XA([6, 15], [6, 15]) + [1/Re, -1/Re; -1/Re, 1/Re];
XA([1, 16], [1, 16]) = XA([1, 15], [1, 16]) + [1/Rf, -1/Rf; -1/Rf, 1/Rf];
XA([8, 17], [8, 17]) = XA([8, 17], [8, 17]) + [1/Rg, -1/Rg; -1/Rg, 1/Rg];
XA([7, 18], [7, 18]) = XA([7, 18], [7, 18]) + [1/Rh, -1/Rh; -1/Rh, 1/Rh];
XA([7, 19], [7, 19]) = XA([7, 19], [7, 19]) + [1/Ri, -1/Ri; -1/Ri, 1/Ri];
XA([7, 20], [7, 20]) = XA([7, 20], [7, 20]) + [1/Rj, -1/Rj; -1/Rj, 1/Rj];
XA([4, 21], [4, 21]) = XA([4, 21], [4, 21]) + [1/Rk, -1/Rk; -1/Rk, 1/Rk];
XB(1,11) = XB(1,11)+1;
XB(1,1) = XB(1,1) -1;
XB(2,12) = XB(2,12)+1;
XB(2,1) = XB(2,1)-1;
XB(3,13) = XB(3,13)+1;
XB(3,1) = XB(3,1)-1;
XB(4,14)=XB(4,14)+1;
XB(4,3) = XB(4,3)-1;
XB(5,15) = XB(5,15)+1;
XB(5,5) = XB(5,5)-1;
XB(6,16) = XB(6,16)+1;
XB(6,5) = XB(6,5)-1;
XB(7,17)=XB(7,17)+1;
XB(7,7)=XB(7,7)-1;
XB(8,18)=XB(8,18)+1;
XB(8,1)=XB(8,1)-1;
XB(9,19)=XB(9,19)+1;
XB(9,2)=XB(9,2)-1;
XB(10,20) = XB(10,20)+1;
XB(10,4)=XB(10,4)-1;
XB(11,21) = XB(11,21)+1;
XB(11,1)=XB(11,1)-1;
XB(12, 6) = XB(14, 6)+1;
XB(12, 10) = XB(12,10)-1;
XB(13,10)=XB(13,10)+1;
XB(13,9)=XB(13,9)-1;
XB(14,9)=XB(14,9)+1;
XB(14,1)=XB(14,1)-1;
XB(15,8)=XB(15,8)+1;
XB(15, 1)=XB(15,1)-1;
XC = XC + XB';
XB(12,3) = XB(12,3)-1;
XB(12,1)=XB(12,1)+1;
XB(13,3) = XB(13,3)-1;
XB(13,4)=XB(13,4)+1;
XB(14,4)=XB(14,4)-1;
XB(14,1)=XB(14,1)+1;
XB(15,5)=XB(15,5)-1;
XB(15,1)=XB(15,1)+1;

X=[XA, XC; XB, XD];

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

X = X(2:end, 2:end);
X_inv = inv(X);

R_vect = [Ra, Rb, Rc, Rd, Re, Rf, Rg, Rh, Ri, Rj, Rk];
R_value = [p1.PortRes, p2.PortRes, p3.PortRes, Rbw.PortRes, Cbw.PortRes, Rout.PortRes, RL.PortRes, C2.PortRes, R2.PortRes, C1.PortRes];
R_diag = diag(R_vect);
left_term = 2*[zeros(11,20), R_diag, zeros(11, 4)];
right_term = [zeros(11,20), eye(11), zeros(11,4)];
S = eye(11) + left_term*(X_inv*(right_term'));

[R_PortRes, param, cond] = solve(S(1,1)==0, Ra, 'ReturnConditions', true);
assume(cond);
R_PortRes = double(subs(R_PortRes, R_vect(2:end), R_value));
R_value = [R_PortRes, R_value];
S = double(subs(S, R_vect, R_value));

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
    WUPs = [WaveUp(p1), WaveUp(p2), WaveUp(p3), WaveUp(Rbw), WaveUp(Cbw), WaveUp(Rout), WaveUp(RL), WaveUp(C2), WaveUp(R2), WaveUp(C1)];
    
    WU_R = S(1,:)*([0, WUPs]');
    WD_R = r*WU_R;
    
    b = S*([WD_R, WUPs]');
    
    p1.WD = b(2);
    p2.WD=b(3);
    p3.WD=b(4);
    Rbw.WD=b(5);
    Cbw.WD=b(6);
    Rout.WD=b(7);
    RL.WD=b(8);
    C2.WD=b(9);
    R2.WD=b(10);
    C1.WD=b(11);
    
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
plot(Fs/NFFT*((0:(NFFT/2))), 20*log10(OUT));
grid on;
xlabel('Frequency (Hz)');
ylabel('magnitude (dB)');