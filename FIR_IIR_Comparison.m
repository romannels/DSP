%------------------------------------------------------------------------%-
% HW 6                                                                   %-
% Roman Nelson                                                           %-
% DSP Fall 2014                                                          %-
%                                                                        %-
% Program to demonstrate various FIR and IIR filters                     %-
%------------------------------------------------------------------------%-

clear all;
%--------------------------------------------------------------------------
% assignment requirements
myPhoneNumber = [4 0 2 4 3 2 4 6 1 5];
Fs = 35e3; % Sampling freq
stopBandMin_dB = 40;
mystopBandMin_dB = 1.25 * stopBandMin_dB;
rippleMax_dB = 1;
%--------------------------------------------------------------------------
% initialize global bandpass parameters
bandwidth = 0;
omegaCenter = 0;
    for i = 1:3
        bandwidth = (myPhoneNumber(7-i) * 10^(i-1)) + bandwidth;
    end
    
    for i = 1:4
        omegaCenter = (myPhoneNumber(11-i) * 10^(i-1)) + omegaCenter;
    end  
Wn = [(omegaCenter - (bandwidth / 2)), (omegaCenter + (bandwidth / 2))];
stopBandTestFreqLow = Wn(1) / 2;
stopBandTestFreqHigh = Wn(2) / 2;
passBandTestFreq = omegaCenter;
transitionBandTestFreq = 0;
WnNorm = Wn / (Fs / 2); % normalize cutoff freq's

plotUpperLimit = (Wn(2) + 500);
plotLowerLimit = (Wn(1) - 500);
%--------------------------------------------------------------------------
% band limits
fpass = [(Wn(1)+ 50) (Wn(2) - 50)];
rpassUpper = [(rippleMax_dB / 2) (rippleMax_dB / 2)];
rpassLower = [-(rippleMax_dB / 2) -(rippleMax_dB / 2)];
fstopUpper = [(1.05*Wn(2)) (1.05*Wn(2)) plotUpperLimit]; 
rstopUpper = [rippleMax_dB -mystopBandMin_dB -mystopBandMin_dB];
fstopLower = [plotLowerLimit (.95*Wn(1)) (.95*Wn(1))]; 
rstopLower = [-mystopBandMin_dB -mystopBandMin_dB rippleMax_dB];
%--------------------------------------------------------------------------
filterType = 3;
%--------------------------------------------------------------------------
% FIR Bandpass
N_FIR = 900; % Order
A_FIR = 1;
B_FIR = fir1(N_FIR, WnNorm);
[H, w] = freqz(B_FIR, A_FIR, 2^16, Fs);
if(filterType == 1)
    plot(w, 20*log10(abs([H])), fpass, rpassUpper, 'r--');
    
    
end
hold on
plot(fpass, rpassLower, 'g');
hold on
plot(fstopUpper, rstopUpper, 'r');
hold on
%plot(fstopLower, rstopLower, 'r');
axis([plotLowerLimit plotUpperLimit -65 2]);



%plot(w, angle(H))
%axis([4400 4800 -5 5]);
%axis([(omegaCenter - 4e3) (omegaCenter + 4e3) -100 10])
 hold off
% impulse response
Y = [1 , zeros(1, N_FIR)];
P = filter(B_FIR, A_FIR, Y);
plot(P);
[gd, w] = grpdelay(B_FIR, A_FIR);
plot(gd);
%x = 0:1/(4 * Fs):pi/8;
%Y = sin(2 * pi * (2 * omegaCenter)  * x);
%plot(x, Y, 'g');

%Y = filter(B_FIR, A, Y);
%plot(x, Y);

% Chebyshev 1 Bandpass
N_Cheby = 6;
[B_Cheby, A_Cheby] = cheby1(N_Cheby, rippleMax_dB, WnNorm);
[H, w] = freqz(B_Cheby, A_Cheby, 2^16, Fs);
if(filterType == 2)
    plot(w, 20*log10(abs([H])), fpass, rpassUpper, 'r--');
    hold off
    Y = [1 , zeros(1, N_FIR)];
    P = filter(B_Cheby, A_Cheby, Y);
    %plot(w, angle(H))
    %axis([4400 4800 -5 5]);
   plot(P);
   [gd, w] = grpdelay(B_Cheby, A_Cheby);
    plot(gd);
end

% Elliptical Bandpass
N_Ellip = 4;
Rp = rippleMax_dB;
Rs = 50 * Rp;
[B_Ellip, A_Ellip] = ellip(N_Ellip, Rp, Rs, WnNorm);
[H, w] = freqz(B_Ellip, A_Ellip, 2^18, Fs);
if(filterType == 3)
    plot(w, 20*log10(abs([H])), fpass, rpassUpper, 'r--');
        hold off
        
    Y = [1 , zeros(1, N_FIR)];
    P = filter(B_Ellip, A_Ellip, Y);
    %plot(w, angle(H))
    axis([4400 4800 -5 5]);
    plot(P);
   [gd, w] = grpdelay(B_Ellip, A_Ellip);
    plot(gd);
end
