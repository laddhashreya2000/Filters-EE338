f_samp = 330e3;
%Band Edge speifications
fp1 = 53.5;
fs1 = 49.5;
fs2 = 77.5;
fp2 = 73.5;

Wc1 = fp1*2*pi/f_samp;
Wc2  = fp2*2*pi/f_samp;

%Kaiser paramters
A = -20*log10(0.15);
if(A < 21)
    beta = 0;
elseif(A <51)
    beta = 0.5842*(A-21)^0.4 + 0.07886*(A-21);
else
    beta = 0.1102*(A-8.7);
end
N_min = ceil((A-8) / (2.285*0.02424*pi));     %empirical formula for N_min

%Window length for Kaiser Window
n=N_min+20;

%Ideal bandpass impulse response of length "n"
bp_ideal = ideal_lp(0.4577*pi,n) - ideal_lp(0.312*pi,n);

%Kaiser Window of length "n" with shape paramter beta calculated above
kaiser_win = (kaiser(n,beta))';

FIR_BandPass = bp_ideal .* kaiser_win;
fvtool(FIR_BandPass);         %frequency response
disp(bp_ideal);
figure;
plot(bp_ideal);
%magnitude response
[H,f] = freqz(FIR_BandPass,1,1024, f_samp);
figure;
plot(f,abs(H))
grid
line([0;18e4], [0.85;0.85], 'Color', 'red');
line([0;18e4], [1.15;1.15], 'Color', 'red');
line([0;18e4], [0.15;0.15], 'Color', 'red');