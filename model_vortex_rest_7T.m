clear all;
path2=[ '../../Turbulence/Basics'];
addpath(genpath(path2));
path3=[ '../../Turbulence/SC_longrange'];
addpath(genpath(path3));
path4=[ '../../TNonequilibrium'];
addpath(genpath(path3));

load (['hcp7t_rfMRI_REST1_PA_schaefer1000.mat']);

load ('schaefer1000to100.mat');
partitiondef=schaefer1000to100;
NP=100;

NSUB=182;
N=1000;
Tmax=900;

Isubdiag = find(tril(ones(NP),-1));

load sc_schaefer_1000.mat;
C=sc_schaefer;
C=C-diag(diag(C));
C1000=C/max(max(C));

for i=1:NP
    partition{i}=find(partitiondef==i);
end

Ctot=C1000;
for i=1:NP
    for j=1:NP
        Caux=Ctot(partition{i},partition{j});
        CtotNP(i,j)=mean(Caux(:));
        if i==j
            CtotNP(i,j)=0;
        end
    end
end

maxC=0.2;
C100=CtotNP/max(max(CtotNP))*maxC;

%%
% Parameters of the data
TR=1;  % Repetition Time (seconds)
% Bandpass filter settings
fnq=1/(2*TR);                 % Nyquist frequency
flp = 0.008;                    % lowpass frequency of filter (Hz)
fhi = 0.08;                    % highpass
Wn=[flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
k=2;                          % 2nd order butterworth filter
[bfilt,afilt]=butter(k,Wn);   % construct the filter

%
FCemp=zeros(NSUB,NP,NP);
for sub=1:NSUB
    sub
    clear signal_filt Phases;
    ts=subject{sub}.schaeferts;
    ts=ts(:,1:Tmax);
    for seed=1:N
        ts(seed,:)=detrend(ts(seed,:)-nanmean(ts(seed,:)));
        if sum(isnan(ts(seed,:)))==0
            signal_filt(seed,:)=filtfilt(bfilt,afilt,ts(seed,:));
            Xanalytic = hilbert(demean(signal_filt(seed,:)));
            Phases(seed,:) = angle(Xanalytic);
        end
    end
    ts=signal_filt(:,50:end-50);
    Phases=Phases(:,50:end-50);

    % Define Partition

    for i=1:NP
        boldpart(i,:)=nanmean(ts(partition{i},:));
        Vorticity(i,:)=abs(sum(complex(cos(Phases(partition{i},:)),sin(Phases(partition{i},:))))/length(partition{i}));
    end

    FCemp(sub,:,:)=corrcoef(Vorticity');
    FCempbold(sub,:,:)=corrcoef(boldpart');
    

    ts=boldpart;
    Tm=size(ts,2);
    [Ns, Tmaxred]=size(ts);
    TT=Tmaxred;
    Ts = TT*TR;
    freq = (0:TT/2-1)/Ts;
    for seed=1:NP
        pw = abs(fft(ts(seed,:)));
        PowSpect(:,seed,sub) = pw(1:floor(TT/2)).^2/(TT/TR);
    end
end

% HOPF parameters Vortex space

Power_Areas=squeeze(mean(PowSpect,3));
for seed=1:NP
    Power_Areas(:,seed)=gaussfilt(freq,Power_Areas(:,seed)',0.01);
end

[maxpowdata,index]=max(Power_Areas);
f_diff = freq(index);

for i=1:NP
    if f_diff(i)<0.02
        f_diff(i)=mean(f_diff);
    end
    if f_diff(i)>0.03
        f_diff(i)=mean(f_diff);
    end
end
%%%
sigma=0.01;
epsFC=0.004;

for sub=1:NSUB
    sub
    FCempsub=squeeze(FCemp(sub,:,:));
    Cnew=C100;
    olderror=100000;
    clear errorFC corrFC;
    for iter=1:5000
        % Linear Hopf FC
        [FCsim]=hopf_int(Cnew,f_diff,sigma);

        errorFC(iter)=nanmean(nanmean((FCempsub-FCsim).^2));
        aux=corrcoef(FCsim(Isubdiag),FCempsub(Isubdiag),'Rows','complete');
        corrFC(iter)=aux(2);

        if mod(iter,10)<0.1
            errornow=nanmean(nanmean((FCempsub-FCsim).^2));
            if  (olderror-errornow)/errornow<0.01
                break;
            end
            if  olderror<errornow
                break;
            end
            olderror=errornow;
        end

        for i=1:NP  %% learning
            for j=1:NP
                if (C100(i,j)>0 || abs(j-i)==NP/2)
                    Cnew(i,j)=Cnew(i,j)+epsFC*(FCempsub(i,j)-FCsim(i,j));
                    if Cnew(i,j)<0
                        Cnew(i,j)=0;
                    end
                end
            end
        end
        Cnew = Cnew/max(max(Cnew))*maxC;
    end
    mean(corrFC(iter-10:iter))
    Ceff(sub,:,:)=Cnew;
end

% save results_Cbad_7T.mat Ceff FCemp FCempbold;
save results_model_vortex_rest_7T.mat Ceff FCemp FCempbold;