clc;
close all
clear

load('C:\Users\Navish\Dropbox\Electrorotation Data\17h_17m_35s.mat')

X_signal = UntitledVoltage_0.Data; % Contains X-Photomultiplier signal
Y_signal = UntitledVoltage_1.Data; % Contains Y-Photomultiplier signal
rate = UntitledVoltage_2.Property.wf_samples; % number of samples per sec
F_state = UntitledVoltage_2.Data; % Shows whether the field was ON or OFF
F_state = (F_state > 2); % convert to binary (0/1)
Time = (10^-4)*(1:size(X_signal))';

Switch = [diff(F_state)]; % find if the field was switched
OnOff = find(Switch); % Indices of the switching data points
On_pairs = [OnOff(1:2:end) OnOff(2:2:end)];
Off_pairs = [OnOff(2:2:end-2) OnOff(3:2:end)];

On_cnt = (On_pairs(3,2)-On_pairs(3,1));
Off_cnt = (Off_pairs(3,2)-Off_pairs(3,1));
On_dur = On_cnt/rate;
Off_dur = Off_cnt/rate;
OnOff_duration = Off_duration + On_duration;

Off_pairs = vertcat(Off_pairs,[OnOff(end) OnOff(end)+(Off_pairs(3,2)-Off_pairs(3,1))]);
Rec_array = Off_pairs(end,2):Off_cnt:size(X_signal,1);
Rec_pairs = [Rec_array(1:end-1)' Rec_array(2:end)'];

%Frequency calculation

for Index = 1:size(On_pairs,1)

%ON
[b,a] = butter(6,40/(rate/2),'high');
X_chunk = X_signal(On_pairs(Index,1):On_pairs(Index,2));
Y_chunk = Y_signal(On_pairs(Index,1):On_pairs(Index,2));

X_chunk = filter(b,a,X_chunk);

[pxx,f] = periodogram(X_chunk - mean(X_chunk), [],[],rate);
[~,loc] = max(pxx);
FRE_ON_per(Index) = f(loc);

% [pxx2,f2] = pmtm(X_chunk - mean(X_chunk),1.25,[],rate);
% [~,loc] = max(pxx2);
% FRE_ON_pmtm = f2(loc)

% figure;
% subplot(3,1,1)
% plot(f,pxx)
% xlim([0 250])
% subplot(3,1,2)
% plot(f1,pxx1)
% xlim([150 250])
% subplot(3,1,3)
% plot(f2,pxx2)
% xlim([150 250])


%Off
X_chunk = X_signal(Off_pairs(Index,1):Off_pairs(Index,2));
Y_chunk = Y_signal(Off_pairs(Index,1):Off_pairs(Index,2));

[pxx,f] = periodogram(X_chunk - mean(X_chunk), [],[],rate);
[~,loc] = max(pxx);
FRE_Off_per(Index) = f(loc);

% figure;
% subplot(3,1,1)
% plot(f,pxx)
% xlim([0 20])
% 10, 14, 18, 27

% [pxx2,f2] = pmtm(X_chunk - mean(X_chunk),1.25,[],rate);
% [~,loc] = max(pxx2);
% FRE_Off_pmtm = f2(loc)
end
figure;
plot(FRE_ON_per,'--o')
ylim([0 max(FRE_ON_per)+ (max(FRE_ON_per) - min(FRE_ON_per))])
figure;
plot(FRE_Off_per,'--o')
