
%%
load('E:\s509_202202\Mat Data\Analysis_new_sorting_v1_4k\Voluntary\period6\Spikes_3std\microGPi1_L_1_spikes_4k.mat')
plot(spikes(1,:))
%%  our new method for feature selection
a= [];
for i=1:size(spikes,1)
    spike= spikes(i,:);
    wave=['sym3']; %best sym3 - bior3.3 interesting
    [c,l] = wavedec(spike,6,wave);
    firstcoefs=[];
    secondcoefs=[];
    currentpoint=0;
    secondcoefHigh=[];
    for i=1:length(l)-1
        firsthalf = round(l(i)/2);
        secondhalf = l(i)-firsthalf;
        firstcoefs=[firstcoefs, currentpoint+1:currentpoint+firsthalf];
        currentpoint = currentpoint+firsthalf;
        secondcoefs= [secondcoefs, currentpoint+1:currentpoint+secondhalf];
        if i>length(l)/2
            secondcoefHigh=[secondcoefHigh, currentpoint+1:currentpoint+secondhalf];
        end
        currentpoint = currentpoint+secondhalf;
    
    end
    % first half
    % m=mean(abs(c(firstcoefs)));
    % s=std(abs(c(firstcoefs)));
    % c(abs(c(firstcoefs))<m)=0;
    firstHalf=c(firstcoefs);
    m=mean(abs(firstHalf));
    s=std(abs(firstHalf));
    firstHalf(abs(firstHalf)<m)=0;
    c(firstcoefs) =firstHalf;
    
    % second half
    secondHalf=c(secondcoefs);
    m=mean(abs(secondHalf));
    s=std(abs(secondHalf));
    secondHalf(abs(secondHalf)<m+s)=0;
    c(secondcoefs) =secondHalf;
    c(secondcoefHigh) =0;
    
%     % m=mean(abs(c));
%     % s=std(abs(c));
%     % c(abs(c)<m)=0;
%     % c(24:end)=0;
%     % c(1:10)=0;
%     sum((c==0))
%     figure,
%     plot(spike)
%     hold on 
%     plot(waverec(c,l,wave));
    
    if isempty(a)
        a=[c~=0];
    else
    a= a+[c~=0];
    end
end
%% our new method for feature selection part2
[b goodCoefsIndeces]=sort(a,'descend');
goodCoefsIndeces = goodCoefsIndeces(1:5);
spike= spikes(2,:);
wave=['sym3']; %best sym3 - bior3.3 interesting
[c,l] = wavedec(spike,6,wave);

goodCoefs = c(goodCoefsIndeces);
c(:)=0;
c(goodCoefsIndeces)=goodCoefs;
sum((c==0))

figure,
plot(spike)
hold on 
plot(waverec(c,l,wave));

%% eriques method

scales = 6;
inputs = 5;
nspk = size(spikes,1);
ls = size(spikes,2);

cc=zeros(nspk,ls);
if exist('wavedec')                             % Looks for Wavelets Toolbox
    for i=1:nspk                                % Wavelet decomposition
        [c,l]=wavedec(spikes(i,:),scales,'haar'); % was haar
        cc(i,1:ls)=c(1:ls);
    end
end

for i=1:ls                                  % KS test for coefficient selection   
    thr_dist = std(cc(:,i)) * 3;
    thr_dist_min = mean(cc(:,i)) - thr_dist;
    thr_dist_max = mean(cc(:,i)) + thr_dist;
    aux = cc(find(cc(:,i)>thr_dist_min & cc(:,i)<thr_dist_max),i);
    
    if length(aux) > 10;
        [ksstat]=test_ks(aux);
        sd(i)=ksstat;
    else
        sd(i)=0;
    end
end
[max_sd ind]=sort(sd);
coeff(1:inputs)=ind(ls:-1:ls-inputs+1);

% spccoefs = coeff(1:inputs);
% tempcoefs = cc(:,spccoefs);
% cc(:,1:end)=0;
% cc(:,spccoefs)=tempcoefs;

inspk=zeros(nspk,53);
for i=1:nspk
    for j=1:inputs
        inspk(i,coeff(j))=cc(i,coeff(j));
    end
end

sum((inspk(20,:)==0))
figure,
plot(spikes(2,:))
hold on 
plot(waverec(inspk(2,:),l,'haar'));

%%
wname = ['haar'];
[LoD,HiD,LoR,HiR] = wfilters(wname); 
subplot(2,2,1)
stem(LoD)
title('Decomposition Lowpass Filter')
subplot(2,2,2)
stem(HiD)
title('Decomposition Highpass Filter')
subplot(2,2,3)
stem(LoR)
title('Reconstruction Lowpass Filter')
subplot(2,2,4)
stem(HiR)
title('Reconstruction Highpass Filter')
xlabel(['The four filters for ',wname])
%%
[phi,psi,xval] = wavefun('sym1');
plot(xval,psi)
title('Wavelet')
%%
spike= spikes(4,:);
cwt(spike)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [KSmax] = test_ks(x)
% 
% Calculates the CDF (expcdf)
%[y_expcdf,x_expcdf]=cdfcalc(x);

yCDF = [];
xCDF = [];
x = x(~isnan(x));
n = length(x);
x = sort(x(:));
% Get cumulative sums
yCDF = (1:n)' / n;
% Remove duplicates; only need final one with total count
notdup = ([diff(x(:)); 1] > 0);
x_expcdf = x(notdup);
y_expcdf = [0; yCDF(notdup)];

%
% The theoretical CDF (theocdf) is assumed to be normal  
% with unknown mean and sigma

zScores  =  (x_expcdf - mean(x))./std(x);

%theocdf  =  normcdf(zScores , 0 , 1);
mu = 0; 
sigma = 1; 
theocdf = 0.5 * erfc(-(zScores-mu)./(sqrt(2)*sigma));


%
% Compute the Maximum distance: max|S(x) - theocdf(x)|.
%

delta1    =  y_expcdf(1:end-1) - theocdf;   % Vertical difference at jumps approaching from the LEFT.
delta2    =  y_expcdf(2:end)   - theocdf;   % Vertical difference at jumps approaching from the RIGHT.
deltacdf  =  abs([delta1 ; delta2]);

KSmax =  max(deltacdf);
end