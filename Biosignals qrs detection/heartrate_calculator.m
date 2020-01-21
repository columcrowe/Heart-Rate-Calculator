function bpm = heartrate_calculator(ecgsignal,fs)
% Outputs heartrate (bpm) calculated at each heart beat in the ECG signal
% Requires input sampling rate of signal fs (Hz)
% MATLAB Version 2017b

% ECG signal
unfiltered_ecgsignal = ecgsignal;

% Sampling Frequency (Hz)

Nyquist=fs/2;
N1=length(ecgsignal); %length of ecg

t=[0:N1-1]./fs; %time (1/Hz) = seconds

ecgsignal = sgolayfilt(ecgsignal, 0, 7); % smoothening filter

ecgsignal = ecgsignal - mean (ecgsignal); % cancel DC conponents

% PAN J., and TOMPKINS W. J. (1985): ‘A Real-Time QRS Detection Algorithm”, IEEE Trans.Biomed. Eng., 32, pp. 230-236
% 5-15 Hz

% Low Pass Filtering
[a,b] = butter(6,15/Nyquist,'low');
ecgsignal = filtfilt(a,b,ecgsignal);

% High Pass Filtering
[a,b] = butter(6,5/Nyquist,'high');
ecgsignal = filtfilt(a,b,ecgsignal);

% Derivative Filtering
a = (fs/8)*[-1 -2 0 2 1];
ecgsignal = filtfilt(a,1,ecgsignal);
ecgsignal = ecgsignal / max(abs(ecgsignal)); % normalise

% Signal Squaring
ecgsignal = ecgsignal .^2;

% Moving Window Integration
% (150 ms)
windowsize = round(0.15*fs);
a = (1/windowsize)*[ones(1,windowsize)];
ecgsignal = filtfilt(a,1,ecgsignal);

ecgsignal = ecgsignal / max(abs(ecgsignal)); % normalise

processed_ecgsignal = ecgsignal;

% setting threshold
threshold=mean(processed_ecgsignal);

% finding where signal exceeds threshold value (heartbeats/qrs complexes)
for i = 1:N1
    if processed_ecgsignal(i)>threshold
        peaks(i)=1;
    end
end

% noting where heartbeat segments start and end
% start point locations 0->1
% end point locations 1->0
c=1;
for i=2:length(peaks)
    if peaks(i)-peaks(i-1)==1
        segment_start(c)=i;
        c=c+1;
    end
end

c=1;
i=length(peaks);
if peaks(i)==1
    segment_end(c)=i;
    c=c+1;
end 
while i>=2
    if peaks(i)-peaks(i-1)==-1
        segment_end(c)=i;
        c=c+1;
    end
    i=i-1;
end
segment_end = fliplr(segment_end);

% # Number of heartbeats in the signal sample
beats=length(segment_end);

% finding maximum within segments = position of r peak
for i=1:beats
    [~, qrs_max_locs(i)] = max(unfiltered_ecgsignal(segment_start(i):segment_end(i)));
end
% locations in relation to full signal
for i=1:beats
    r_peak_locs(i) = qrs_max_locs(i)+segment_start(i)-1; % add offset
end

% R-R intervals
for i=2:length(r_peak_locs)
    r_r_interval(i) = t(r_peak_locs(i)) - t(r_peak_locs(i-1));
    i=i+1;
end
r_r_interval=r_r_interval.';

% calculate heartrate (bpm)
for i=2:length(r_r_interval)
    bpm(i) = 1./(r_r_interval(i)./60);
end
bpm(1)=nan;

% plotting
figure
zoom xon;
hold on
plot(unfiltered_ecgsignal);
hold on
plot(r_peak_locs,unfiltered_ecgsignal(r_peak_locs),'rs','MarkerFaceColor','r');
legend('ecg signal','r _ peaks');
xlabel('samples');
ylabel('amplitude(mV)');
for i=1:length(bpm)
    text(r_peak_locs(i),unfiltered_ecgsignal(r_peak_locs(i)),num2str(round(bpm(i))),'FontSize',8);
end

end




