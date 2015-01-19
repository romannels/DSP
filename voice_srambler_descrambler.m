%-------------------------------------------------------------------------%
%---------------------------VOICE SCRAMBLER-------------------------------%
%----------------------------DSP FALL 2014--------------------------------%
%-----------------------------Roman Nelson--------------------------------%
%-------------------------------------------------------------------------%

clear all;
clc;

Fs = 8e3; % sampling freq
s = wavread('sample.wav');
fo = 2e3; % arbitrary foldover
wo = 2 * pi * fo / Fs; % arbitrary foldover (periodic)
N = 100

sound(s(1:20000),Fs);
pause;
% remove negative frequencies
hilbertaroni = firpm(N,[.01 .99],[1 1],'Hilbert'); % Hilbert filter
sHilbert = filter(hilbertaroni, 1, s); % use Hilbert filter on signal
sHilbertPrime = sHilbert * 1i;
sPrime =  [zeros(1, N/2 + 1) s((N/2 + 2): length(s))']; % delay real signal
sPlus = sPrime' + sHilbertPrime; % add real and transformed signals

 [SPlus, w] = freqz(sPlus, 1, 1024, 'whole'); % get spectrum
 plot(w/pi - 1, 20 * log10(fftshift(SPlus)));
 pause
 periodogram(SPlus);
 pause
 % look at periodograms of signal between various steps
 title('Negative Frequencies Removed (attempt)', 'FontSize', 24, 'FontWeight', 'bold');
xlabel('Frequency (Hz)', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('Magnitude (dB)', 'FontSize', 18, 'FontWeight', 'bold')
saveas(gcf,'Neg_removed_attempt.jpg');
pause
[Hilbertaroni, w] = freqz(hilbertaroni, 1,1024, 'whole');
 plot(w/pi - 1, 20 * log10(fftshift(SPlus)));
 title('Hilbert Transform Freq Response (length = 100)', 'FontSize', 24, 'FontWeight', 'bold');
xlabel('Frequency (Hz)', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('Magnitude (dB)', 'FontSize', 18, 'FontWeight', 'bold')
saveas(gcf,'HilbertXform.jpg');
pause
periodogram(sPlus);
 title('Periodogram of S+', 'FontSize', 24, 'FontWeight', 'bold');
xlabel('Frequency (Hz)', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('Magnitude (dB)', 'FontSize', 18, 'FontWeight', 'bold')
saveas(gcf,'hilbertoneperiodo.jpg');
pause

x = zeros(1,length(sPlus));
% soundsc(sPlus(1:10000));
% frequency shift
trippy = ones(1, length(sPlus)); %rand(1, length(sPlus)); %randomize


for n = 1: length(sPlus)
    x(n) = sPlus(n) * exp((-1i) * wo * n * trippy(n));
end

% soundsc(x(1:10000));

xHilbert = filter(hilbertaroni, 1, x);
xHilbertPrime = 1i * xHilbert;
xPrime = [zeros(1, N/2 + 1) x(N/2 + 2 : length(x))];
xPlus =  xPrime + xHilbertPrime;
xPlus = xPlus / 2;

xMinus = x - xHilbertPrime;
xMinus = xMinus / 2;
scramble = zeros(1,length(xPlus));

for n = 1 : length(xPlus)
    scramble(n) = xPlus(n) + (-1)^n * xMinus(n);
end

sound(real(scramble(1:20000)), Fs);
periodogram(scramble);
 title('Periodogram of scrambled signal', 'FontSize', 24, 'FontWeight', 'bold');
xlabel('Frequency (Hz)', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('Magnitude (dB)', 'FontSize', 18, 'FontWeight', 'bold')
saveas(gcf,'scrambleperiodo.jpg');

 pause
 
 % descrambler. This is the exact same thing as above save for one part...
 [Scramble, w] = freqz(scramble, 1, 1024, 'whole');
 
 foPrime = (Fs / 2) - fo;
 woPrime = 2 * pi * foPrime / Fs; % this is the only difference
 % This time we shift what was the upper positive frequencies before down to
 % the negative frequencies so that when we flip it they return to their
 % original position
 
 scrambleHilbert = filter(hilbertaroni, 1, scramble);
 scrambleHilbertPrime = scrambleHilbert * 1i;
 scramblePrime = [zeros(1,N + 1) scramble(N + 2:length(scramble)) ];
 scramblePlus = scramblePrime + scrambleHilbertPrime;
 z = zeros(1,length(scramblePlus));
 
 for n = 1 : length(scramblePlus)
     z(n) = scramblePlus(n) * exp((-1i) *woPrime * (trippy(n)) * n);
 end
 
 zHilbert = filter(hilbertaroni, 1, z);
 zHilbertPrime = 1i * zHilbert;
 zPrime = [zeros(1,N + 1) z(N + 2 : length(z))];
 zPlus = zPrime + zHilbertPrime;
 zPlus = zPlus / 2;
 
 zMinus = z - zHilbertPrime;
 zMinus = zMinus / 2;
 descramble = zeros(1, length(zPlus));
 for n = 1 : length(zPlus)
     descramble(n) = zPlus(n) + (-1)^n * zMinus(n);
 end
 
  sound(real(descramble(1:40000)), Fs);
  periodogram(descramble);
  pause
  
title('Periodogram of Descrambled signal', 'FontSize', 24, 'FontWeight', 'bold');
xlabel('Frequency (Hz)', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('Magnitude (dB)', 'FontSize', 18, 'FontWeight', 'bold')
saveas(gcf,'descrambleperiodo.jpg');
  
periodogram(s);

title('Periodogram of original signal', 'FontSize', 24, 'FontWeight', 'bold');
xlabel('Frequency (Hz)', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('Magnitude (dB)', 'FontSize', 18, 'FontWeight', 'bold')
saveas(gcf,'origperiodo.jpg');
 
 
