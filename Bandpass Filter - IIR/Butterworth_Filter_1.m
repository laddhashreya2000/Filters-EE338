%Butterworth Analog LPF parameters
Wc = 1.0788;              %cut-off frequency
N = 8;                  %order 

%poles of Butterworth polynomial of degree 8 in the open CLHP 
p1 = Wc*cos(pi/2 + pi/16) + Wc*sin(pi/2 + pi/16)*1i;
p2 = Wc*cos(pi/2 + pi/16) - i*Wc*sin(pi/2 + pi/16);
p3 = Wc*cos(pi/2 + pi/16+pi/8) + i*Wc*sin(pi/2 + pi/16+pi/8);
p4 = Wc*cos(pi/2 + pi/16+pi/8) - i*Wc*sin(pi/2 + pi/16+pi/8);
p5 = Wc*cos(pi/2 + pi/16+2*pi/8) + i*Wc*sin(pi/2 + pi/16+2*pi/8);
p6 = Wc*cos(pi/2 + pi/16+2*pi/8) - i*Wc*sin(pi/2 + pi/16+2*pi/8);
p7 = Wc*cos(pi/2 + pi/16+3*pi/8) + i*Wc*sin(pi/2 + pi/16+3*pi/8);
p8 = Wc*cos(pi/2 + pi/16+3*pi/8) - i*Wc*sin(pi/2 + pi/16+3*pi/8);

fprintf("p1 = %d", real(p1));
fprintf("+ i %d", imag(p1));
disp(" ");
fprintf("p2 = %d", real(p2));
fprintf("+ i %d", imag(p2));
disp(" ");
fprintf("p3 = %d", real(p3));
fprintf("+ i %d", imag(p3));
disp(" ");
fprintf("p4 = %d", real(p4));
fprintf("+ i %d", imag(p4));
disp(" ");
fprintf("p5 = %d", real(p5));
fprintf("+ i %d", imag(p5));
disp(" ");
fprintf("p6 = %d", real(p6));
fprintf("+ i %d", imag(p6));
disp(" ");
fprintf("p7 = %d", real(p7));
fprintf("+ i %d", imag(p7));
disp(" ");
fprintf("p8 = %d", real(p8));
fprintf("+ i %d", imag(p8));

%Band Edge speifications
fp1 = 53.5;
fs1 = 49.5;
fs2 = 77.5;
fp2 = 73.5;

%Transformed Band Edge specs using Bilinear Transformation
f_samp = 330;         
wp1 = tan(fp1/f_samp*pi);
ws1 = tan(fs1/f_samp*pi); 
ws2 = tan(fs2/f_samp*pi);
wp2 = tan(fp2/f_samp*pi);

%Parameters for Bandpass Transformation
W0 = sqrt(wp1*wp2);
B = wp2-wp1;

[num,den] = zp2tf([],[p1 p2 p3 p4 p5 p6 p7 p8],Wc^N);   %TF with poles p1-p8 and numerator Wc^N and no zeroes
                                                        %numerator chosen to make the DC Gain = 1

%Evaluating Frequency Response of Final Filter
syms s z;
analog_lpf(s) = poly2sym(num,s)/poly2sym(den,s);        %analog LPF Transfer Function
analog_bsf(s) = analog_lpf((s*s + W0*W0)/(B*s));        %bandpass transformation
discrete_bsf(z) = analog_bsf((z-1)/(z+1));              %bilinear transformation

%coeffs of analog bsf
[ns, ds] = numden(analog_bsf(s));                   %numerical simplification to collect coeffs
ns = sym2poly(expand(ns));                          
ds = sym2poly(expand(ds));                          %collect coeffs into matrix form
k = ds(1);    
ds = ds/k;
ns = ns/k;
disp("Coefficient of analog bpf");
disp(ns);
disp(ds);
%coeffs of discrete bsf
[nz, dz] = numden(discrete_bsf(z));                     %numerical simplification to collect coeffs                    
nz = sym2poly(expand(nz));
dz = sym2poly(expand(dz));                              %coeffs to matrix form
k = dz(1);                                              %normalisation factor
dz = dz/k;
nz = nz/k;
disp("Coefficient of discrete bpf");
disp(nz);
disp(dz);
fvtool(nz,dz)                                           %frequency response

%magnitude plot (not in log scale) 
[H,f] = freqz(nz,dz,1024*1024, 330e3);
figure;
plot(f,abs(H))
grid
line([0;18e4], [0.85;0.85], 'Color', 'red');
line([0;18e4], [0.15;0.15], 'Color', 'red');