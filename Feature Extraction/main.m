%Loading the data
data=readtable('./bidmc_csv/bidmc_05_Signals.csv');
ppg_raw=data(:,3);
ppg_raw=ppg_raw.PLETH;
time=data(:,1);
t=time.Time_s_;

%Path to the file containing true values
path_truth='./bidmc_csv/bidmc_05_Numerics.csv'

%Function to calculate ground truth
[RR_windows,del_rows]=comp_ground_truth(path_truth)
%Preprocessing of raw PPG data
[y_riiv, y_riav, y_rifv,t1]=preprocessing(ppg_raw,t);
%Generating sliding windows
[temp_riiv,temp_riav,temp_rifv]=sliding_windows(y_riiv, y_riav, y_rifv,t1,del_rows);
%Computing AR density
[RR_computed,y,u,o]=AR_density(temp_riiv,temp_riav,temp_rifv)
%Plotting the predicted against ground truth and calculating the MAE
[MAE,diff]=evaluation(RR_windows,RR_computed);
%Saving the MAEs
%writematrix(diff,'05.csv') 

function [y_riiv, y_riav, y_rifv,t1]=preprocessing(ppg_raw,t);

% data=readtable(path);
% ppg_raw=data(:,3);
% ppg_raw=ppg_raw.GreenCount;
% time=data(:,1);
% t=time.Time;

%Converting Datetime to time
% reftime=t(1,1);
% t_temp=etime(datevec(t), repmat(datevec(reftime),numel(t),1));
% t=t_temp;


%Removing DC component
cutoff_freq=0.004166;
norm_cutoff_freq=cutoff_freq*(125/2);
ppg_raw=highpass(ppg_raw,norm_cutoff_freq,125);
% ppg_raw=lowpass(ppg_raw,0.333,125)
ppg_raw=ppg_raw(1:length(ppg_raw),:);
t=t(1:length(t),:);
% ppg_raw=lowpass(ppg_raw,15,125)
%Finding peaks
[pks,locs]=findpeaks(ppg_raw);
ts=t(locs);
y=ppg_raw;
Fs=125;
annot=ts;
ys=y;
% %Calculating IMS
[peaks,onsets,clipp] = adaptPulseSegment(y,Fs,annot);
% Plotting IMS results
% figure(1);
% subplot();
% plot(t,y);
% hold on
% plot(t(peaks),y(peaks));
% hold on
% plot(t(clipp),y(clipp),'r--');
% hold on
% plot(t(onsets),y(onsets));
% title('Raw PPG with its peaks');
tp=t(peaks);
y_riiv_raw=[]
y_riav_raw=[]
y_rifv_raw=[]
y_riiv=[]
y_riav=[]
y_rifv=[]
for i=1:length(peaks);
    if i==length(peaks);
        break
    else
    %Extracting RIIV
    y_riiv_raw(i,1)=y(peaks(i));
    %Extracting RIAV
    if peaks(1)<onsets(1)
        y_riav_raw(i,1)=y(peaks(i))-y(onsets(i));
    else
        y_riav_raw(i,1)=y(peaks(i+1))-y(onsets(i));
    end
        %Extracting RIFV
    y_rifv_raw(i,1)=t(peaks(i+1))-t(peaks(i));
%         y_rifv_raw(i,1)=y(peaks(i+1))-y(onsets(i));
    end
end
%Resampling the signals at Fs=8Hz using linear interpolation
Fs_r=8;
tp=tp(1:length(y_riav_raw),:);
[y_riiv,t1]=resample(y_riiv_raw,tp,Fs_r,'linear');
[y_riav,t2]=resample(y_riav_raw,tp,Fs_r,'linear');
[y_rifv,t3]=resample(y_rifv_raw,tp,Fs_r,'linear');
%Normalizing the data
y_riiv=(normalize(y_riiv)-mean(y_riiv))/std(y_riiv);
y_riav=(normalize(y_riav)-mean(y_riav))/std(y_riav);
y_rifv=(normalize(y_rifv)-mean(y_rifv))/std(y_rifv);
%Plotting the extracted features
% figure(2);
% subplot(3,1,1);
% plot(t1,y_riiv);
% title('RIIV');
% subplot(3,1,2);
% plot(t2,y_riav);
% title('RIAV');
% subplot(3,1,3);
% plot(t3,y_rifv);
% title('RIFV');
%To check if resampled and pre-resampled overlap
% figure(2)
% subplot(3,1,1);
% plot(t1,y_riiv);
% hold on 
% plot(tp,y_riiv_raw)
% title('RIIV');
% subplot(3,1,2);
% plot(t2,y_riav);
% hold on 
% plot(tp,y_riav_raw)
% title('RIAV');
% subplot(3,1,3);
% plot(t3,y_rifv);
% hold on 
% plot(tp,y_rifv_raw)
% title('RIFV')

function [peaks,onsets,clipp] = adaptPulseSegment(y,Fs,annot)
%ADAPTPULSESEGMENT perform adaptive pulse segmentation and artifact detection 
%in ppg signals
%   [peaks,onsets,artif] = adaptPulseSegment(y,annot)
%
% Inputs:
%       y      vector, ppg signal [Lx1] or [1xL], in which L = length(signal)
%       Fs      scalar, sampling frequency in Hz
%       annot   vector, timestamps (in samples) with location of the peaks
%
% Outputs:
%       peaks   vector, locations of peaks (in samples)
%       onsets  vector, locations of onsets of the beats (in samples)
%       artif   vector, locations of peaks classified as artefacts
% 
% References:
%       Karlen et al., Adaptive Pulse Segmentation and Artifact Detection in 
%       Photoplethysmography for Mobile Applications, 34th Internat. Conf. IEEE-EMBS 2012
%       
% Written by Marco A. Pimentel

doOptimise = 1;
doPlot = 0;
% if nargin < 3
%     % no annotations are provided, therefore, no optimisation will take
%     % place
%     doOptimise = 0; 
% end


% The algorithm in the paper is applied to signals sampled at 125 Hz...
% We do not resample our signal
%ys = resample(y,125,Fs);
%Fs = 125;
% if Fs ~= 300
%     ys = resample(y,300,Fs);
%     Fs = 300;
% else
ys = y;
% end

% The paper is not clear about the selection of the range of search for m
% We define the range of m to be between [0.005 - 0.100] secs (5ms to 100ms)
% We define "m" in terms of samples
opt.bounds = 0.005:0.005:0.100;
opt.m = unique(ceil(opt.bounds*Fs));

opt.perf = zeros(length(opt.m),4); % store results of performance for each m

if doOptimise
    % Perform optimisation
    for i = 1 : length(opt.m)
        % Determine peaks and beat onsets
        [peaks,artifs,clipp,linez,linezSig] = pulseSegment(ys,Fs,opt.m(i));
        % Calculate performance of the peak detection
%         opt.perf(i,:) = evalPerf(annot,linez(:,2));
    end
    
else
    % Do not perform optimization; fix m
    opt.m = 10;
    [peaks,artifs,clipp,linez,linezSig] = pulseSegment(ys,Fs,opt.m);
end

if doPlot
    colData = {'g','y','r'};
    figure; 
    h(1) = subplot(211);
    plot(ys); hold on;
    for i = 1 : size(linez,1)
        %if linezSig(i) > -1
        plot(linez(i,:),ys(linez(i,:)),'-x','Color',colData{linezSig(i)+2});
        %end
    end
    
    h(2) = subplot(212);
    plot(ys,'g'); hold on;
    for i = 1 : size(peaks,1)
        plot(peaks(i,:),ys(peaks(i,:)),'-xr');
    end
    if ~isempty(artifs)
    for i = 1 : size(artifs,1)
        plot(artifs(i,:),ys(artifs(i,:)),'--^b');
    end
    end
    if ~isempty(clipp)
    for i = 1 : size(clipp,1)
        plot(clipp(i,:),ys(clipp(i,:)),'-om');
    end
    end
    linkaxes(h,'x');
    
end

% Correct for the downsmapling performed during the peak detection
onsets = peaks(:,1);
peaks  = peaks(:,2);
for i = 1 : size(peaks,1)
    [~,ind]  = min(ys(max([1 onsets(i)-opt.m]):min([length(ys) onsets(i)+opt.m])));
    onsets(i) = max([1 onsets(i)-opt.m]) + ind(1) - 1;
    [~,ind]  = max(ys(max([1 peaks(i)-opt.m]):min([length(ys) peaks(i)+opt.m])));
    peaks(i) = max([1 peaks(i)-opt.m]) + median(ind) - 1;
end

% Correct minimum value of onset of the beat
for i = 2 : length(onsets)
    [~,ind]   = min(ys(peaks(i-1):peaks(i)));
    onsets(i) = peaks(i-1) + ind - 1;
end
end

function [peaks,artifs,clipp,linez,linezSig] = pulseSegment(ys,Fs,m)
% split signal in different segments
nseg = floor(length(ys)/m);   % number of segments
% intialize loop variables
seg = 1;    % segment counter
z = 1;      % line counter 
segInLine = 1;  % line controler
linez = zeros(nseg,2); linez(1,:) = [1,m];
% slope of segment/line
a = zeros(nseg,1); a(1) = slope(ys,linez(1,:));
% "classify" line segment according to the slope
linezSig = zeros(nseg,1); linezSig(1) = sign(a(1));
% Start loop over segments
z = z + 1; seg = seg + 1;
for i = 2 : nseg    % loop over segments
    linez(z,:) = [(seg-1)*m+1 seg*m];
    try
        a(z) = slope(ys,linez(z,:));
    catch
        a = 1;
    end
    linezSig(z) = sign(a(z));
    if sign(a(z)) == sign(a(z-1))
        linez(z-1,:) = [(seg-1-segInLine)*m+1 seg*m];
        seg = seg + 1;
        segInLine = segInLine + 1;
    else
        z = z + 1;
        seg = seg + 1;
        segInLine = 1;
    end
end

% remove extra spaces created in output variables
linezSig(sum(linez,2)==0,:) = [];
linez(sum(linez,2)==0,:) = [];
% Apply adaptive threshold algorithm
% For this algorithm to work, we need to first find a valide line segment 
% in order to intialize the thresholds! In order to this, we define a flag
% to control the intialization in the main loop
FOUND_L1 = 0;

% The algorithm includes the definition of 4 adaptation parameters
% We define the following adaptation parameters
% a =     | a_fast_low    a_fast_high |
%         | a_slow_low    a_slow_high |
% 
a = [0.5 1.6; ...
     0.6 2.0];
 
% Define fixed thresholds described in the paper
ThT  = 0.03 * Fs;    % Duration of the line
ThIB = 0.24 * Fs;    % Interbeat invertal (240 ms) 

% Define parameters used in the main loop
alpha = zeros(size(linez,1),1);
for i = 1 : size(linez,1)
    alpha(i) = slope(ys,linez(i,:));   % slopes of line segments
end
theta = diff(ys(linez),[],2);
durat = diff(linez,[],2);       % duration of line segments (in samples)

% remove lines that do not have the necessary duration
linez(durat<ThT,:) = [];
theta(durat<ThT,:) = [];
alpha(durat<ThT,:) = [];
horiz = horizontalLine(ys,linez,Fs);

FLAG = 0;
artifs = []; clipp = [];
% Select window for detect firs peaks!
wind = theta(theta>0);
try 
    wind = wind(1:10);
catch
    wind = wind;
end
ThAlow  = prctile(wind,95)*0.6;
ThAhigh = prctile(wind,95)*1.8;
peaks = [];
for z = 1 : size(linez,1)-1   % loop over line segments
    if FOUND_L1
        if alpha(z) > 0 && ... % slope must be positive
                alpha(z-1) ~= 0 && ...  % peaks before or after clipping are artefactual
                alpha(z+1) ~= 0
            if theta(z) >= ThAlow && theta(z) <= ThAhigh && ...
                    linez(z,2) >= peaks(end,2) + ThIB
                ThAlow  = (ThAlow + theta(z)*a(2,1))/2;
                ThAhigh = theta(z) * a(2,2);
                FLAG = 0;
                currTheta = [currTheta; theta(z)];
                peaks = [peaks    ; linez(z,:)];
            else
                if FLAG > 0
                    ThAlow  = (ThAlow + min(currTheta(max([1 end-4]):end))*a(1,1))/2;
                    ThAhigh = max(currTheta(max([1 end-4]):end)) * a(1,2);
                    %ThAlow  = (ThAlow + theta(z)*a(1,1))/2;
                    %ThAhigh = theta(z) * a(1,2);
                end
                FLAG = FLAG + 1;
                artifs = [artifs; linez(z,:)];
            end
        elseif theta(z) > 0 && ... 
                ((theta(z-1) ~= 0 || horiz(z-1) ~= 0) && ...
                (theta(z+1) ~= 0 || horiz(z+1) ~= 0))
            if theta(z) >= ThAlow && theta(z) <= ThAhigh && ...
                    linez(z,2) >= peaks(end,2) + ThIB
                ThAlow  = (ThAlow + theta(z)*a(2,1))/2;
                ThAhigh = theta(z) * a(2,2);
                FLAG = 0;
                currTheta = [currTheta; theta(z)];
                peaks = [peaks; linez(z,:)];
            else
                if FLAG > 0
                    %ThAlow  = (ThAlow + currTheta*a(1,1))/2;
                    %ThAhigh = currTheta * a(1,2);
                    ThAlow  = (ThAlow + min(currTheta(max([1 end-4]):end))*a(1,1))/2;
                    ThAhigh = max(currTheta(max([1 end-4]):end)) * a(1,2);
                    %ThAlow  = (ThAlow + theta(z)*a(1,1))/2;
                    %ThAhigh = theta(z) * a(1,2);
                end
                FLAG = FLAG + 1;
                artifs = [artifs; linez(z,:)];
            end
        elseif theta(z) == 0 && horiz(z) == 0
            artifs  = [artifs; linez(z,:)];
            clipp   = [clipp; linez(z,:)];
        end 
    else
        if alpha(z) > 0 && durat(z) >= ThT && ...
                theta(z) >= ThAlow && theta(z) <= ThAhigh 
            FOUND_L1 = 1;
            ThAlow  = theta(z)*0.5;
            ThAhigh = theta(z)*2.0;
            peaks = linez(z,:);    % loaction of onsets and peaks
            currTheta = theta(z);
        end
    end
end
end

function out = slope(ys,interv)
start = interv(1); stop = interv(2);
out = sum(diff(ys([start:stop])))/(stop-start);
%out = median(gradient(ys(start:stop)));
end

function out = horizontalLine(ys,linez,Fs)
% Get horizontal lines from signal given linez
out = zeros(size(linez,1),1);
for i = 1 : size(linez,1)
    out(i) = median(abs(diff(ys(linez(i,1):linez(i,2)))));
    % check duration of the peaks
    if out(i) == 0 && diff(linez(i,:)) <= 0.200*Fs
        out(i) = 0.1;
    end
end

end
function yth = flatline(y,th,hyst)
  yth = false(size(y)); % change #1 - `false` vs `zeros`
  c1 = y > (th - hyst); % change #2 - precompute c1, c2
  c2 = y > (th + hyst);

  for k = 2:numel(y)
      if yth(k - 1) 
          yth(k) = c1(k);
      else
          yth(k) = c2(k);
      end 
  end
end
end


function [temp_riiv,temp_riav,temp_rifv]=sliding_windows(y_riiv, y_riav, y_rifv,t1,del_rows);
% Computing AR
% Splitting the data into 32 seconds window with an 8 second time difference
j=0;
row=0;
col=0;
temp_riiv=0;
temp_riav=0;
temp_rifv=0;
temp1=y_riiv;
temp2=y_riav;
temp3=y_rifv;
[m,n]=size(temp1);
for i = 1:63:length(t1);
    row=row+1;
    col=0;
    for j=i:i+255;
        col=col+1;
        if j>m;
            break
            break
        else
            temp_riiv(row,col)=temp1(j);
            temp_riav(row,col)=temp2(j);
            temp_rifv(row,col)=temp3(j);
        end
    end
end
% temp_riiv=highpass(temp_riiv,0.004166,4);
% temp_riav=highpass(temp_riav,0.004166,4);
% temp_rifv=highpass(temp_rifv,0.004166,4);
temp_riiv=removerows(temp_riiv,del_rows);
temp_riav=removerows(temp_riav,del_rows);
temp_rifv=removerows(temp_rifv,del_rows);

end


function [RR_computed,y,u,o]=AR_density(temp_riiv,temp_riav,temp_rifv);
riiv_sliding=temp_riiv;
riav_sliding=temp_riav;
rifv_sliding=temp_rifv;

RR_riiv=[];
RR_riav=[];
RR_rifv=[];
RR_computed=[];
row=1;
[q,w]=size(riav_sliding);

for i=1:q;
    AR_riiv=AR_density(riiv_sliding(i,:),3);
    RR_riiv=[RR_riiv AR_riiv];
    AR_riav=AR_density(riav_sliding(i,:),4);
    RR_riav=[RR_riav AR_riav];
    AR_rifv=AR_density(rifv_sliding(i,:),5);
    RR_rifv=[RR_rifv AR_rifv];
end
RR_riiv=transpose(RR_riiv);
RR_riav=transpose(RR_riav);
RR_rifv=transpose(RR_rifv);

if RR_riiv(1,1)==0;
    RR_riiv(1,1)=max(RR_riiv);
end
if RR_riav(1,1)==0
    RR_riav(1,1)=max(RR_riav)
end
if RR_rifv(1,1)==0;
    RR_rifv(1,1)=max(RR_rifv);
end
o=[];
for i=1:length(RR_riav);
    if RR_riav(i,:)==0;
%         o(i,1)=max(RR_riav);
        o(i,1)=o(i-1,:);
%         o(i,1)=mean(o)
    else
        o(i,1)=RR_riav(i,1);
    end

end
u=[]
for i=1:length(RR_riiv);
    if RR_riiv(i,:)==0;
%         u(i,1)=max(RR_riiv);
        u(i,1)=u(i-1,:);
%         u(i,1)=mean(u)
    else
        u(i,1)=RR_riiv(i,1);
    end

end
y=[]
for i=1:length(RR_rifv);
    if RR_rifv(i,:)==0;
%         y(i,1)=max(RR_rifv);
        y(i,1)=y(i-1,:);
%         y(i,1)=mean(y)
    else
        y(i,1)=RR_rifv(i,1);
    end

end
RR_computed(:,1)=(o+y+u)/3;

function RR_final=AR_density(sig,fig_no);
    RR_final=[];
    for i=2:2;
        A=ar(sig,i,'burg');
        m=A.A;
%         m=resample(m,);
        [H,F] = freqz(1,m,[],8);
%         [H,F]=freqz(100,.004166,m);
%         RR(:,1)=F*60
        amp=20*log10(abs(H));
        amp=amp(15:512,:);
        F=F(15:512,:);
%         RR=F*60;
%         figure(fig_no);
%         plot(RR,amp);
%         xlabel('RR (bpm)');
%         ylabel('PSD (dB/Hz)');
%         hold on ;
%         [e r]=size(amp);
%         amp=amp(50:e,r);
        [val, index]=max(amp);
        max_freq=F(index);
        temp=max_freq*60;
        RR_final=[RR_final temp];
    end
end
end

function [RR_windows,del_rows]=comp_ground_truth(path_truth)
data=readtable(path_truth);
RR_true=data(:,4);
RR_true=RR_true.RESP;
RR_windows=[];
[m n] =size(RR_true);
row1=1
q=0
for i =1:8:length(RR_true)
    temp=[];
    temp1=0;
    col=0;
    row=1;
    for j=i:i+32;
        col=col+1;
        if j>m
            break
            break
        else
            temp(row,1)= RR_true(j);
            row=row+1;
        end

    end
    temp1=mean(temp);
    RR_windows(row1,1)= temp1;
    row1=row1+1;
    q=q+1;
end

[del_rows, columns] = find(isnan(RR_windows));
RR_windows=removerows(RR_windows,del_rows);

end

function[MAE,diff]=evaluation(RR_windows,RR_computed)
%Comparing the computed results with the ground truth
%Input : RR_windows (ground truth), RR_computed(calculated)
figure(6);
plot(RR_windows,'o--');
hold on 
plot(RR_computed);
ylabel('RR (bpm)')
xlabel('Window samples of 32s')

%In this section we will calculate the MAE and p-error values 

%MAE
final=[]
[v p]=size(RR_windows);
RR_computed=RR_computed(1:v,:);
n=length(RR_computed);
diff=abs(RR_computed-RR_windows);
MAE=(1/n)*sum(diff)
% final=[RR_windows RR_computed]

end
