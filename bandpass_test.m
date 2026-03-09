figure,
a =load('E:\s509_202202\Mat Data\Analysis_withBipolarcommonremoval_300_3000_4STD\Voluntary\period6\Spikes_4std\microGPi1_L_1_spikes.mat')
plot(a.spikes(2,:))
hold on
%%
load('E:\s509_202202\Mat Data\Analysis_new_sorting_v1_4k\Voluntary\period6\Spikes_3std\microGPi1_L_1_spikes_35k.mat')
plot(spikes(1,:))
hold on 
%%
load('E:\s509_202202\Mat Data\Analysis_new_sorting_v1_4k\Voluntary\period6\Spikes_3std\microGPi1_L_1_spikes_4k.mat')
plot(spikes(1,:))
%%
x= spikes(20,:);
figure,
[pxx,w]=periodogram(x,rectwin(length(x)),length(x),24414);
plot(w,pxx)
hold on
filterCutoffs = [6000];
[filtData,~] =  butterfilter (filterCutoffs,24414,2,x,'low');
[pxx,w]=periodogram(filtData,rectwin(length(filtData)),length(filtData),24414);
plot(w,pxx)

figure,
plot(x)
hold on 
plot(filtData)


%%
load('E:\s509_202202\Mat Data\Analysis_new_sorting_v1_4k\Voluntary\period6\Spikes_3std\microGPi1_L_1_spikes_55k.mat')
plot(spikes(1,:))
%%
load('E:\s509_202202\Mat Data\Analysis_new_sorting_v1_4k\Voluntary\period6\Spikes_3std\microGPi1_L_1_spikes_7k.mat')
plot(spikes(1,:))
%%
load('E:\s509_202202\Mat Data\Analysis_new_sorting_v1_4k\Voluntary\period6\Spikes_3std\microGPi1_L_1_spikes_10k.mat')
plot(spikes(1,:))

%%
%%
load('E:\s509_202202\Mat Data\Voluntary\micro_bipolar\period6\microGPi1_L_5_filtered.mat')
figure,
[pxx,w]=pwelch(data,24414,12207,[],24414);
plot(w,10*log10(pxx))
%% %%%
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

function fs_out = checkfs(fs_in)
    if fs_in<1e3
        fs_out = fs_in*1e3;
    else
        fs_out = fs_in;
    end
end