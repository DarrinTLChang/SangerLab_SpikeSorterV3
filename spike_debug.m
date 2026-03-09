clc 
clear all
close all

%% look at 510 right arm reach data
% load spikes (spikeWaveform)
load 'E:\s510_202203\Mat Data\Analysis_new_sorting_v1_5cluster\Voluntary\period12\SpikeClusters_5cluster_3std_wav\spikeWaveform.mat'

% load bipolardata (data,fs)
load 'E:\s510_202203\Mat Data\Voluntary\micro_bipolar\period12\microVoSTN_L_2_filtered.mat'
bipolar_data = data;

% load monopolardata (data,fs)
load 'E:\s510_202203\Mat Data\Voluntary\micro\period12\microVoSTN_L_2.mat'
monopolar_data = data;
micro_fs=fs;

% load emg (data,fs)
load 'E:\s510_202203\Mat Data\Voluntary\emg\Run_number_360_Plot_and_Store_Rep_1.12.mat'
emg_data = emgDelsys;
emg_fs=sfEmg;

%% plot raw recordings
channel = 34;
data = monopolar_data;
display_segment = [0.001,100]; % in  seconds
display_segment_index=floor(display_segment*micro_fs);
display_length = display_segment(2)-display_segment(1);
t= display_segment_index(1):display_segment_index(2);
t=t/micro_fs;
spiketimes = spikeWaveform(channel).time;
figure,
ax(1) = subplot(3,1,1)
plot(t , data(display_segment_index(1):display_segment_index(2)))
hold on 
for i=1:length(spiketimes) 
    % add spike on top of time data
    spiketimes = spikeWaveform(channel).time;
    spike_indices = spikeWaveform(channel).time{1,i}(:)>display_segment (1) *1000 & ...
    spikeWaveform(channel).time{1,i}(:)<display_segment(2)*1000;
    plot(spiketimes{1,i}(spike_indices)/1000 , data(round(spiketimes{1,i}(spike_indices)*micro_fs/1000)),'*')
    hold on 
end
yline(4e-5)
yline(-4e-5)
xlim(display_segment)
title('raw input')

filterCutoffs = [300 3000];
[filtData,~] =  butterfilter (filterCutoffs,micro_fs,8,data,'bandpass');

% filterCutoffs = [535 545];
% [filtData,~] =  butterfilter (filterCutoffs,micro_fs,4,filtData,'stop');
% filterCutoffs = [655 665];
% [filtData,~] =  butterfilter (filterCutoffs,micro_fs,4,filtData,'stop');
% filterCutoffs = [895 905];
% [filtData,~] =  butterfilter (filterCutoffs,micro_fs,4,filtData,'stop');
% filterCutoffs = [1015 1025];
% [filtData,~] =  butterfilter (filterCutoffs,micro_fs,4,filtData,'stop');
% filterCutoffs = [1135 1145];
% [filtData,~] =  butterfilter (filterCutoffs,micro_fs,4,filtData,'stop');
% filterCutoffs = [1255 1265];
% [filtData,~] =  butterfilter (filterCutoffs,micro_fs,4,filtData,'stop');
% filterCutoffs = [1375 1385];
% [filtData,~] =  butterfilter (filterCutoffs,micro_fs,4,filtData,'stop');
% filterCutoffs = [1495 1505];
% [filtData,~] =  butterfilter (filterCutoffs,micro_fs,4,filtData,'stop');
% centerfrq= 1740;
% filterCutoffs = [centerfrq-5 centerfrq+5];
% [filtData,~] =  butterfilter (filterCutoffs,micro_fs,4,filtData,'stop');
% centerfrq= 1860;
% filterCutoffs = [centerfrq-5 centerfrq+5];
% [filtData,~] =  butterfilter (filterCutoffs,micro_fs,4,filtData,'stop');
% centerfrq= 2100;
% filterCutoffs = [centerfrq-5 centerfrq+5];
% [filtData,~] =  butterfilter (filterCutoffs,micro_fs,4,filtData,'stop');
% centerfrq= 2220;
% filterCutoffs = [centerfrq-5 centerfrq+5];
% [filtData,~] =  butterfilter (filterCutoffs,micro_fs,4,filtData,'stop');
% centerfrq= 2460;
% filterCutoffs = [centerfrq-5 centerfrq+5];
% [filtData,~] =  butterfilter (filterCutoffs,micro_fs,4,filtData,'stop');
% centerfrq= 2580;
% filterCutoffs = [centerfrq-5 centerfrq+5];
% [filtData,~] =  butterfilter (filterCutoffs,micro_fs,4,filtData,'stop');
% centerfrq= 2820;
% filterCutoffs = [centerfrq-5 centerfrq+5];
% [filtData,~] =  butterfilter (filterCutoffs,micro_fs,4,filtData,'stop');
% centerfrq= 3180;
% filterCutoffs = [centerfrq-5 centerfrq+5];
% [filtData,~] =  butterfilter (filterCutoffs,micro_fs,4,filtData,'stop');
% centerfrq= 5000;
% filterCutoffs = [centerfrq-5 centerfrq+5];
% [filtData,~] =  butterfilter (filterCutoffs,micro_fs,4,filtData,'stop');

ax(2) = subplot(3,1,2)
plot(t , filtData(display_segment_index(1):display_segment_index(2)))
hold on 
for i=1:length(spiketimes) 
    % add spike on top of time data
    spiketimes = spikeWaveform(channel).time;
    spike_indices = spikeWaveform(channel).time{1,i}(:)>display_segment (1) *1000 & ...
    spikeWaveform(channel).time{1,i}(:)<display_segment(2)*1000;
    plot(spiketimes{1,i}(spike_indices)/1000 , filtData(round(spiketimes{1,i}(spike_indices)*micro_fs/1000)),'*')
    hold on 
end
yline(4e-5)
yline(-4e-5)
xlim(display_segment)
title('bandpasssed filter input')

filterCutoffs = [500 3500];
[filtData,~] =  butterfilter (500,micro_fs,4,data,'high');
ax(3) = subplot(3,1,3)
plot(t , filtData(display_segment_index(1):display_segment_index(2)))
hold on 
for i=1:length(spiketimes) 
    % add spike on top of time data
    spiketimes = spikeWaveform(channel).time;
    spike_indices = spikeWaveform(channel).time{1,i}(:)>display_segment (1) *1000 & ...
    spikeWaveform(channel).time{1,i}(:)<display_segment(2)*1000;
    plot(spiketimes{1,i}(spike_indices)/1000 , filtData(round(spiketimes{1,i}(spike_indices)*micro_fs/1000)),'*')
    hold on 
end
yline(4e-5)
yline(-4e-5)
xlim(display_segment)
title('highpasssed filter input')

linkaxes(ax,'x')
%% ploting PSDs
figure,
subplot(2,1,1)
[pxx,f] = periodogram(filtData(display_segment_index(1):display_segment_index(2)),rectwin(length(filtData(display_segment_index(1):display_segment_index(2)))),length(filtData(display_segment_index(1):display_segment_index(2))),micro_fs);
plot(f,pxx)

subplot(2,1,2)
segmentLength = floor(10*micro_fs);
[pxx, w] = pwelch(filtData(display_segment_index(1):display_segment_index(2)),floor(micro_fs),floor(micro_fs/2),[],micro_fs);
plot(w,pxx)
%% new nonlinear operatore
sample_data=filtData(display_segment_index(1):display_segment_index(2));
figure,
ax(1)=subplot(4,1,1)
t=1:length(sample_data);
plot(sample_data)
hold on
title("raw data")
[up,lo] = envelope(sample_data,50,"peak");
hold on
plot(t,up,'-',t,lo,'--')
hold off

ax(2)= subplot(4,1,2)
neo = sqrt(abs(sample_data.^2-circshift(sample_data,1).*circshift(sample_data,-1)));
plot(neo)
title("new data")

ax(3)= subplot(4,1,3)
neo_slow = nthroot(abs(sample_data.^8-circshift(sample_data,1).*circshift(sample_data,-1).* ...
    circshift(sample_data,2).*circshift(sample_data,-2).* ...
    circshift(sample_data,3).*circshift(sample_data,-3).* ...
    circshift(sample_data,4).*circshift(sample_data,-4)),8);
plot(neo_slow)
title("new data slower")

ax(4)= subplot(4,1,4)

% filterCutoffs = [1200 3500];
% [filtData2,~] =  butterfilter (1200,micro_fs,4,sample_data,'high');
% neo = sqrt(abs(filtData2.^2-circshift(filtData2,1).*circshift(filtData2,-1)));
% plot(neo)

s = spectrogram(sample_data,25,'yaxis');
plot(linspace(t(1),t(end), length(s)), sum(abs(s(10:26,1:end))))


title("spectrogram")


linkaxes(ax,'x')


%%


function [x,g] =  butterfilter (fc, fs, order, data, type)
    %[x,g] = butterfilter(fc, fs, order, data, type)
    % Zero-phase Butterworth filter
    % Outputs: x is the filtered data
    %          g is the gain of the filter
    % Inputs:  fc = Cut-off frequency (Hz)
    %          fs = Sampling rate (Hz)
    %          order = Filter order
    %          data = data to be filtered
    %          type: 'low', 'high', 'bandpass' or 'stop'
    
    if strcmpi(type,'bandpass')
        [B,A] = butter(order,2*fc/fs);     % Creates two vectors for TS = B/A
    else
        [B,A] = butter(order,2*fc/fs,type);     % Creates two vectors for TS = B/A
    end
    % [sos,g] = tf2sos(B,A);
    x = filtfilt(B,A,data);                 
    [~,g] = tf2sos(B,A);
end


