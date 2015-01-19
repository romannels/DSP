% HW 5
% Roman Nelson
% DSP Fall 2014

clear all
FsOriginal = 8e3; %samples/sec
L = 5; % upsample scale factor 
FsNew = L * FsOriginal; %samples/second
fdc = 1.2 * 4000; % filter slightly above 4k
Fsby2 = FsNew / 2; % nyquist rate
Wn = fdc / Fsby2; %normalized filter cutoff
N = 50; %filter order
signal = wavread('voice_samp_8k.wav');
plotLength = 100; % for making plots more readable
plotStart = 10 % offset for plots to get to good stuff
signal = signal'; %transpose (don't know if this is neccessary)
sigSize = size(signal) * L; % determine new signal size
upSample = kron(signal, ones(1,L)); %upsample
%plot(upSample(plotStart:(plotLength - 1)), 'r-.'); % no interpolation

hold on

Bsh = [0.5,0.5]; % zero order hold interpolation(sample and hold) 
holdIt = filter(Bsh, 1, upSample); % run signal through filter
plot(holdIt(plotStart:(plotLength - 1)), 'g--');
%title('sample and hold, time domain');
hold on

Bli = [0.25, 0.5, 0.25]; % linear interpolation
linearInt = filter(Bli, 1, upSample); %run signal through filter
%plot(linearInt(plotStart:(plotLength - 1)), 'm');
%title('linear interpolation time domain');
hold on

Blp = fir1(N,Wn); % low pass filter, low order
lowPassIt = filter(Blp, 1, upSample); %run signal through filter
%soundsc(lowPassIt, FsNew);
plot(lowPassIt(((N / 2) + plotStart):((plotLength - 1) + (N /2))), 'b');
title('low order low pass interpolation time domain');
hold on

N = 40 * N;
Blp2 = fir1(N,Wn); % low pass filter, high order
lowPassIt = filter(Blp2, 1, upSample); %run signal through filter
%soundsc(lowPassIt, FsNew);
plot(lowPassIt(((N / 2) + plotStart):((plotLength - 1) + (N /2))), 'r:');
axis([1 (plotLength-10) -0.09 .12])
xlabel('time in 1/Fs seconds')
ylabel('amplitude')
%title('time domain signal comparison')

grid
%legend('low order', 'high order');
%legend('No Interp.', 'S & H','Lin. interp.', 'Low Order L.P', 'High Order L.P.');



