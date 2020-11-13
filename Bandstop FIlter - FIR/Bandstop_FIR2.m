f_samp = 260e3;
%Band Edge speifications
fs1 = 48.7e3;
fp1 = 44.7e3;
fp2 = 72.7e3;
fs2 = 68.7e3;

%Kaiser paramters
A = -20*log10(0.15);
if(A < 21)
    beta = 0;
elseif(A <51)
    beta = 0.5842*(A-21)^0.4 + 0.07886*(A-21);
else
    beta = 0.1102*(A-8.7);
end

Wn = [(fs1+fp1)/2 (fs2+fp2)/2]*2/f_samp;        %average value of the two paramters
N_min = ceil((A-7.95) / (2.285*0.02424*pi));       %empirical formula for N_min

%Window length for Kaiser Window
n=N_min + 13;

%Ideal bandstop impulse response of length "n"

bs_ideal =  ideal_lp(pi,n) -ideal_lp(0.5435*pi,n) + ideal_lp(0.3595*pi,n);

%Kaiser Window of length "n" with shape paramter beta calculated above
kaiser_win = (kaiser(n,beta))';

FIR_BandStop = bs_ideal .* kaiser_win;
fvtool(FIR_BandStop);         %frequency response
disp(bs_ideal);
figure;
plot(bs_ideal);
%magnitude response
[H,f] = freqz(FIR_BandStop,1,1024, f_samp);
figure;
plot(f,abs(H))
grid
line([0;18e4], [0.85;0.85], 'Color', 'red');
line([0;18e4], [1.15;1.15], 'Color', 'red');
line([0;18e4], [0.15;0.15], 'Color', 'red');