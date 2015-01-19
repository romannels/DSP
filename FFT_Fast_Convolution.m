%-------------------------------------------------------------------------%
%--------------------------------FFT/DFT----------------------------------%
%----------------------------DSP FALL 2014--------------------------------%
%-----------------------------Roman Nelson--------------------------------%
%-------------------------------------------------------------------------%

clc
clear all

Fs = 8e3; % sampling frequency (Hz)
n = 0:1/Fs:1; % samples
messItUp = 2.1; % uh oh not periodic in 128 samples
f = (127 + messItUp) *2; % test signal frequency
x = sin(n * 2 * f * pi); % test signal
%plot(x); % debug
padding = 128; % truncation length 
x_prime = x(1:padding); % truncated signal
x_prime_prime = [x_prime, zeros(1,padding*7)]; % zero padded signal
x_prime_naught = x(1:padding*7 + 128); % sample big chunk
save = 0;
% FFT Fun

X = fft(x);
plot((1:length(X)), abs(X)); % original signal
title('Original Signal Frequency Spectrum', 'FontSize', 24, 'FontWeight', 'bold');
xlabel('Frequency (Hz)', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('Magnitude', 'FontSize', 18, 'FontWeight', 'bold')
if save == 1;
saveas(gcf,'originalFreq.jpg');
end
pause


X = fft(x_prime);
plot((1:length(X)), abs(X)); % truncated signal
title('128 Samples of Original Signal with 128 zero padding Frequency Spectrum', 'FontSize', 24, 'FontWeight', 'bold');
xlabel('Frequency (Hz)', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('Magnitude', 'FontSize', 18, 'FontWeight', 'bold')
if save == 1;
saveas(gcf,'sample128pad.jpg');
end
pause


X = fft(x_prime_prime);
plot((1:length(X)), abs(X)); % lots of zero padding
title('128 Samples of Original Signal with 896 zero padding Frequency Spectrum', 'FontSize', 24, 'FontWeight', 'bold');
xlabel('Frequency (Hz)', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('Magnitude', 'FontSize', 18, 'FontWeight', 'bold')
if save == 1;
saveas(gcf,'sample896pad.jpg');
end
pause

X = fft(x_prime_naught);
plot((1:length(X)), abs(X)); % replace zeros with actual data
title('1024 Samples of Original Signal Frequency Spectrum', 'FontSize', 24, 'FontWeight', 'bold');
xlabel('Frequency (Hz)', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('Magnitude', 'FontSize', 18, 'FontWeight', 'bold')
if save == 1;
saveas(gcf,'sample1024.jpg');
end
pause

x_prime_prime = [x_prime zeros(1,7896)]; % sample a ton
X = fft(x_prime_prime);
plot((1:length(X)), abs(X));
title('128 Samples of Original Signal with 7896 zero padding Frequency Spectrum', 'FontSize', 24, 'FontWeight', 'bold');
xlabel('Frequency (Hz)', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('Magnitude', 'FontSize', 18, 'FontWeight', 'bold')
if save == 1;
saveas(gcf,'sample7896pad.jpg');
end
pause

% Fast datasream convolution
signal = rand(1,10000); % random signal
sigLength = length(signal); % not much shorter, but looks nicer
sampleLength = 2000;
buffer = zeros(1, sampleLength); % allocate
%convo = zeros(1, sigLength);
N = 1000; % filter order
Wn = 0.3; % random cutoff
superFancyFilter = fir1(N, Wn); % low pass
convo = zeros(1, sigLength+ 1); % allocate
 
for i = 1: sigLength
    if mod(i, sampleLength) == 0    
%        Buffer = fft([buffer, zeros(1, N - sampleLength + 1)]);
%        Convo = superFancyFilter .* Buffer';
        convo_prime = conv(buffer, superFancyFilter); % cheated a bit here
        for l = 1 : length(convo_prime)

            
            %convo_prime = ifft(Convo);
            if (i - sampleLength + l) > (sigLength + 1)
                continue
            end
            %add current chunk to old ones
        convo(i - sampleLength + l) = convo(i - sampleLength + l)+ convo_prime(l);
        end
        continue
    end
    
    buffer(mod(i,sampleLength)) = signal(i); % build a buffer
                 
end
y = filter(superFancyFilter, 1, signal); % please be close!!!!!


subplot(2,1,1), plot(1:1000, y(1:1000) / max(y) );
title('Fast Convolution', 'FontSize', 24, 'FontWeight', 'bold');
xlabel('Frequency (Hz)', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('Magnitude', 'FontSize', 18, 'FontWeight', 'bold')

subplot(2,1,2), plot(1:1000, convo(1:1000) / max(convo), 'r');
title('MATLAB filter()', 'FontSize', 24, 'FontWeight', 'bold');
xlabel('Frequency (Hz)', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('Magnitude', 'FontSize', 18, 'FontWeight', 'bold')
if save == 1;
saveas(gcf,'comparison.jpg');
end
pause;

subplot(111);
diff = convo - [filter(superFancyFilter, 1, signal), 0];
plot(abs(diff));
title('Difference between "filter()" and fast convo', 'FontSize', 24, 'FontWeight', 'bold');
xlabel('Frequency (Hz)', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('Magnitude', 'FontSize', 18, 'FontWeight', 'bold')
if save == 1;
saveas(gcf,'diff.jpg');
end
