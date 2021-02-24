%Loading the datasets
data=readtable('bidmc_csv/bidmc_04_Signals.csv');
data_true=readtable('bidmc_csv/bidmc_04_Numerics');
ppg_raw=data(:,3);
ppg_raw=ppg_raw.PLETH;
time=data(:,1);
t=time.Time_s_;

%Converting Datetime to time
% reftime=t(1,1);
% t_temp=etime(datevec(t), repmat(datevec(reftime),numel(t),1));
% t=t_temp;

%Removing DC component
%w=2piF
rate=64
f1=0.1
f2=18
order=4
[B,A] = butter(order,[f1/(rate/2) f2/(rate/2)]);
data_filt=filter(B,A,ppg_raw); 

%Normalizing the data
data_normalized=(normalize(data_filt)-mean(data_filt))/std(data_filt);

% Splitting the data into 32 seconds window with an 8 second time difference
j=0;
row=0;
col=0;
temp=0;
data_windows=[]
[m,n]=size(data_normalized);
for i = 1:1001:length(t);
    row=row+1;
    col=0;
    for j=i:i+4001;
        col=col+1;
        if j>m;
            data_windows(row,col)=temp;
            
        else
            data_windows(row,col)=data_normalized(j);
            temp=data_normalized(j);
        end
    end
end


%For testing csv file
%Loading the dataset
RR_true=data_true(:,4);
RR_true=RR_true.RESP;
RR_windows=[];
[m n] =size(RR_true);
row1=1
q=0
% Splitting the data into 32 seconds window with an 8 second time difference
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

%Making sure both files have the same length
[a b]=size(data_windows);
[c d]=size(RR_windows);
if a>c
    data_windows=data_windows(1:c,:);
else
    RR_windows=RR_windows(1:a,:);
end

%Saving the csv files
writematrix(RR_windows,'bidmc4_true.csv') 
writematrix(data_windows,'bidmc4.csv') 