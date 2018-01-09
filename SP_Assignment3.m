%% Q1
Fs = 5000; % [Hz]
T  = 2;    % [sec]
t  = linspace(-T/2, T/2, T*Fs);

x  = cos(120*pi*t).*sinc(50*t); 

% downsample by fds
fds = [250 100];
x1  = downsample(x, Fs./fds(1)); t1 = downsample(t, Fs./fds(1));
x2  = downsample(x, Fs./fds(2)); t2 = downsample(t, Fs./fds(2));

figure(); plot(t, x);
hold on
stem(t1, x1);
stem(t2, x2);
axis([-0.05 0.05 -1 1]);
xlabel('Time [sec]'); ylabel('x(t)');
hold off

% schematic plot of fourier transform
A  = 1/100;
f2 = [-85-fds(2) -85-fds(2) -35-fds(2) -35-fds(2) -85 -85 ...,
    -fds(2)+35 -fds(2)+35 -35 -35 -fds(2)+85 -fds(2)+85 fds(2)-85 fds(2)-85 ...
    35 35 fds(2)-35 fds(2)-35 85 85 fds(2)+35 fds(2)+35 fds(2)+85 fds(2)+85];
A2 = [0 A A 0 0 A A 2*A 2*A A A 0 0 A A 2*A 2*A A A 0 0 A A 0];
figure(); plot(f2, A2); hold on
plot([-fds(2) -fds(2)], [0 2*A],'--k', [fds(2) fds(2)], [0 2*A], '--k');
hold off
ax = gca;
ax.XTick = [-85-fds(2) -35-fds(2) -fds(2) -85 ...
    -fds(2)+35 -35 -fds(2)+85 fds(2)-85 ...
    35 fds(2)-35  85 fds(2) fds(2)+35 fds(2)+85];
ax.XTickLabel = {'185','-135','-f_{s}','-85','-65','-35','-15','15','35',...
    '65','85','f_{s}','135','185'};
xlabel('f [Hz]'); ylabel('|X(f)|');
f1 = [-85-fds(1) -85-fds(1) -35-fds(1) -35-fds(1) -fds(1)+35 -fds(1)+35 -fds(1)+85 -fds(1)+85 ...
    -85 -85 -35 -35 35 35 85 85 fds(1)-85 fds(1)-85 fds(1)-35 ...
    fds(1)-35 fds(1)+35 fds(1)+35 fds(1)+85 fds(1)+85];
A1 = [0 A A 0 0 A A 0 0 A A 0 0 A A 0 0 A A 0 0 A A 0];
figure(); plot(f1,A1);
hold on
plot([-fds(1) -fds(1)], [0 A],'--k', [fds(1) fds(1)], [0 A], '--k');
hold off
ax = gca;
ax.XTick = [-85-fds(1) -35-fds(1) -fds(1)  ...
    -fds(1)+35 -fds(1)+85 -85 -35   ...
    35 85 fds(1)-85 fds(1)-35 fds(1) fds(1)+35 fds(1)+85];
ax.XTickLabel = {'-335','-285','-f_{s}','-215','-165','-85','-35','35','85',...
    '165','215','f_{s}','285','335'};
xlabel('f [Hz]'); ylabel('|X(f)|');

% fft 
L = T*Fs;
X = fft(x)/L*T; % L = T/D
f = (-L/2:L/2-1)*Fs/L; 

Xshift = fftshift(abs(X));
figure(); plot(f, Xshift);
xlabel('f [Hz]'); ylabel('|X(f)|');

L1 = T*fds(1);
X1 = fft(x1)/L1*T;
f1 = (-L1/2:L1/2-1)*fds(1)/L1;

X1shift = fftshift(abs(X1));
figure(); plot(f1 , X1shift);
xlabel('f [Hz]'); ylabel('|X_{1}(f)|');

L2 = T*fds(2);
X2 = fft(x2)/L2*T;
f2 = (-L2/2:L2/2-1)*fds(2)/L2;


X2shift = fftshift(abs(X2));
figure(); plot(f2 , X2shift);
xlabel('f [Hz]'); ylabel('|X_{2}(f)|');

%% Q2
clear all %#ok
close all 

% create transient signal
Fs = 300; % [Hz]
T  = 2;  % [sec]
t  = linspace(0, T, Fs*T);
x  = exp(-5*t).*sin(2*pi*5*t);

figure(); plot(t, x); xlabel('Time [sec]'); ylabel('x(t)');

% fft calculation
L = T*Fs;
X = abs(fft(x))/L;
k = 0:L-1;
figure(); stem(k, X); xlabel('k'); ylabel('|X(f)|');
xlim([0 25]);

% copy x(t)
x_rep = repmat(x,1,2);
t_rep = linspace(0, 2*T, 2*Fs*T);

figure(); plot(t_rep, x_rep); xlabel('Time [sec]'); ylabel('x_{rep}(t)');

% fft of repeating x(t)
L_rep = T*Fs; 
X_rep = abs(fft(x_rep))/L_rep;
k_rep = 0:L_rep*2-1;
figure(); stem(k_rep,X_rep); xlabel('k'); ylabel('|X_{rep}(f)|');
xlim([0 50]);

% compare
f     = (-length(X)/2:length(X)/2-1)*Fs/length(X);
f_rep = (-length(X_rep)/2:length(X_rep)/2-1)*Fs/length(X_rep);

Xs     = fftshift(X);
Xs_rep = fftshift(X_rep);
figure(); stem(f, Xs); hold on; stem(f_rep, Xs_rep);
xlabel('f [Hz]'); ylabel('|X(f)|');
legend('Original','Repeated');
xlim([-25 25]);
hold off

%% Q3
clear all %#ok
close all
f  = [1536, 1578.667, 1621.333]; % [Hz]
fc = f(2);
fm = f(2) - f(1); % we could also use f(3) - f(2)

% plot
T  = 0.1;
Fs = 1e4;
t  = linspace(0, T, Fs*T);
x  = 0.5*(cos(2*pi*t*(fc + fm)) + cos(2*pi*t*(fc - fm))) + cos(2*pi*t*fc);
figure(); plot(t, x);
xlabel('Time [sec]'); ylabel('x(t)');

minT = 2/fm;

% plot fft with given requirements
Fs  = 5*2*(fc + fm);
N   = 1024;
T   = N/Fs;
t   = linspace(0, T, N);
Kam = 0.5;
x   = Kam*0.5*(cos(2*pi*t*(fc + fm)) + cos(2*pi*t*(fc - fm))) + cos(2*pi*t*fc);
X   = abs(fft(x))/N;
X   = X(1:N/2 + 1);

X(2:end-1) = 2*X(2:end-1);
f = Fs*(0:N/2)/N;
figure(); stem(f, X);
xlabel('f [Hz]'); ylabel('|X(f)|');
xlim([0, 1.5*(fm + fc)]);

% plot padded fft
Fs  = 5*2*(fc + fm);
N   = 1024;
T   = N/Fs;
t   = linspace(0, T, N);
Kam = 0.5;
x   = Kam*0.5*(cos(2*pi*t*(fc + fm)) + cos(2*pi*t*(fc - fm))) + cos(2*pi*t*fc);
x   = [x, zeros(1,1024)]; % pad with 1024 zeros
X   = abs(fft(x))/length(x);
X   = X(1:length(x)/2 + 1);

X(2:end-1) = 2*X(2:end-1);
f = Fs*(0:length(x)/2)/length(x);
figure(); stem(f, X);
xlabel('f [Hz]'); ylabel('|X(f)|');
xlim([0, 1.5*(fm + fc)]);

% downsample
Fs  = 5*2*(fc + fm);
N   = 1024;
T   = N/Fs;
t   = linspace(0, T, N);
Kam = 0.5;
x   = Kam*0.5*(cos(2*pi*t*(fc + fm)) + cos(2*pi*t*(fc - fm))) + cos(2*pi*t*fc);
x   = downsample(x, 8);
N   = N/8;
X   = abs(fft(x))/N;
X   = X(1:N/2 + 1);

X(2:end-1) = 2*X(2:end-1);
f = Fs/8*(0:length(x)/2)/length(x);
figure(); stem(f, X);
xlabel('f [Hz]'); ylabel('|X(f)|');
%% Q4
clear all %#ok
close all
M1 = 16; M2 = 32;
f1 = 20;       % [Hz]
f2 = 20*M1/M2; % [Hz]
fm = 20*M1;    % [Hz]
T  = 5;        % [sec]
k = 1:20;
Fs = 1e5;     % the minimum is k(end)*fm*2

% plot fourier series approximation 
a_n = 4*(1 - (-1).^k)./(pi^2*k.^2);
t = linspace(-T/2, T/2, T*Fs);
h_n = [];
for ii = 1:k(end)
    h_n = [h_n, cos(2*pi*ii*t.'*fm)]; %#ok
end
xApprox = h_n*a_n.';
figure(); plot(t, xApprox); xlabel('Time [sec]'); ylabel('x(t)');


x = xApprox + xApprox.*cos(2*pi*f1*t.') + xApprox.*cos(2*pi*f2*t.');
f = (-length(x)/2:length(x)/2 -1)*Fs/(length(x));
X = fft(x);
X = fftshift(abs(X)/length(x));

figure(); plot(f,X);
xlim([-5e3 5e3]);
xlabel('f [Hz]'); ylabel('|X(f)|');

%% Q5

Fs = 2^13; % [Hz]
f0 = 440;  % [Hz]
n = [N.E N.F N.G N.A N.G N.E N.F N.E N.D N.F N.E N.E N.F N.G N.A N.G N.E...
    N.F N.E N.D N.C N.E N.A N.B N.C5 N.C5 N.B N.B N.A N.G N.G N.A N.G N.E...
    N.F N.E N.D N.C];
T = [0.2 .2 .2 .2 .4 .4 .2 .2 .2 .2 .8 .2 .2 .2 .2 .4 .4 .4 .2 .2 .8 .4...
    .2 .2 .4 .4 .4 .4 .8 .4 .2 .2 .6 .2 .4 .2 .2 .8];
% while 1
%     makeTune(n, T, Fs);
%     pause(sum(T)); 
% end
% sound(cos(exp(0:1e-3:2)*2*pi*10)); % R2D2 sound
s = makeTune(n, T, Fs);

% remove last sample to make length(s) even
s(end) = []; 

% plot the fft of s
S = abs(fft(s))/length(s);
S = S(1:length(s)/2 +1);
S(2:end-1) = 2*S(2:end-1);
f = Fs*(0:length(s)/2)/length(s);

figure(); plot(f, S);
xlabel('f [Hz]'); ylabel('|S(f)|');
xlim([0, 800]);
ax = gca;
ax.XTick = [2^(N.C/12)*f0 2^(N.D/12)*f0 2^(N.E/12)*f0 ...
    2^(N.F/12)*f0 2^(N.G/12)*f0 2^(N.A/12)*f0 2^(N.B/12)*f0 2^(N.C5/12)*f0];
ax.XTickLabel = {'261.63','293.66','329.67','349.22','391.99','440','493.88','523.25'};
ax.XTickLabelRotation = 90;

% add noise
load('noise.mat');
sn = s + noise.';
% remove last sample to make length(s) even
sn(end) = []; 

sound(sn);
pause(sum(T));

% plot the fft of s + noise
Sn = abs(fft(sn))/length(sn);
Sn = Sn(1:length(sn)/2 +1);
Sn(2:end-1) = 2*Sn(2:end-1);
f = Fs*(0:length(sn)/2)/length(sn);

figure(); plot(f, Sn);
xlabel('f [Hz]'); ylabel('|S_{n}(f)|');

% filter the noise 
% butterworth filter of order 3, with cuttoff frequency
[b, a] = butter(40, 1500/(Fs/2));  

sf = filter(b, a, sn);
sound(sf);
% plot the fft of s + noise
Sf = abs(fft(sf))/length(sf);
Sf = Sf(1:length(sf)/2 +1);
Sf(2:end-1) = 2*Sf(2:end-1);
f = Fs*(0:length(sf)/2)/length(sf);
figure(); plot(f, Sf);
xlabel('f [Hz]'); ylabel('|S_{f}(f)|');