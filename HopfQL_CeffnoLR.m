clear all;
path2=[ '../../Nonequilibrium/'];
addpath(genpath(path2));
path3=[ '../../Tenet/TENET/'];
addpath(genpath(path3));
path4=[ '../../TaskResrvoir/DataHCP100ordered'];
addpath(genpath(path4));
path5=[ '../../Turbulence/Basics'];
addpath(genpath(path5));
path6=[ '../../Turbulence/SC_longrange'];
addpath(genpath(path6));
path7=[ '../../VortexModel/PartitionModel'];
addpath(genpath(path7));

load (['hcp_REST1_LR_schaefer1000.mat']);

N=1000;
indexN=1:N;
sigma=0.01;
NSUB=971;
Tau=3;
Isubdiag = find(tril(ones(N),-1));

%% FCemp
% Parameters of the data
TR=0.72;  % Repetition Time (seconds)
% Bandpass filter settings
fnq=1/(2*TR);                 % Nyquist frequency
flp = 0.008;                    % lowpass frequency of filter (Hz)
fhi = 0.08;                    % highpass
Wn=[flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
k=2;                          % 2nd order butterworth filter
[bfilt,afilt]=butter(k,Wn);   % construct the filter

FC=zeros(NSUB,N,N);
COVtau=zeros(NSUB,N,N);
for nsub=1:NSUB
    nsub
    ts=subject{nsub}.schaeferts;  % fMRI
    ts=ts(indexN,:);
    clear signal_filt Phases;
    for seed=1:N
        ts(seed,:)=detrend(ts(seed,:)-nanmean(ts(seed,:)));
        if sum(isnan(ts(seed,:)))==0
            signal_filt(seed,:)=filtfilt(bfilt,afilt,ts(seed,:));
        end
    end
    ts2=signal_filt(:,100:end-100);
    FCemp=corrcoef(ts2');
    FC(nsub,:,:)=FCemp;
    COVemp=cov(ts2');
    % COV(tau)
    tst=ts2';
    for i=1:N
        for j=1:N
            sigratio(i,j)=1/sqrt(COVemp(i,i))/sqrt(COVemp(j,j));
            [clag lags] = xcov(tst(:,i),tst(:,j),Tau);
            indx=find(lags==Tau);
            COVtauemp(i,j)=clag(indx)/size(tst,1);
        end
    end
    COVtauemp=COVtauemp.*sigratio;
    COVtau(nsub,:,:)=COVtauemp;
end

FCemp=squeeze(mean(FC));
COVtauemp=squeeze(mean(COVtau));

%% Compute Crest

load hpcdata1003_f_diff_fce.mat;
f_diff=f_diff';

sig = 0.01;
dt=0.1*TR/2;
dsig = sqrt(dt)*sig;

AA=-0.02;

a=AA*ones(N,1);
a=repmat(a,1,2);
wo = f_diff*(2*pi);
omega = repmat(wo,1,2);
omega(:,1) = -omega(:,1);

load sc_schaefer_1000.mat;
C=sc_schaefer;
C = C/max(max(C));
C=C-diag(diag(C));

load COG_schaefer1000.mat;
dd=cog;

for i=1:N
    for j=1:N
        distancia(i,j)=sqrt(sum((dd(i,:)-dd(j,:)).^2));
        if distancia(i,j)>80
            dmask(i,j)=0;
        else
            dmask(i,j)=1;
        end
    end
end

C=dmask.*C;

epsFC=0.0004;
epsFCtau=0.0001;

maxC=0.2;

C=C*maxC;

Tmax=1200
nn=1;
Cnew=C;
nng=1;

olderror=100000;
for iter=1:160
    iter
    % Linear Hopf FC

    [FCsim,COVsim,COVsimtotal,A]=hopf_int(Cnew,f_diff,sigma);
    COVtausim=expm((Tau*TR)*A)*COVsimtotal;
    COVtausim=COVtausim(1:N,1:N);
    for i=1:N
        for j=1:N
            sigratiosim(i,j)=1/sqrt(COVsim(i,i))/sqrt(COVsim(j,j));
        end
    end
    COVtausim=COVtausim.*sigratiosim;

    errorFC(iter)=nanmean(nanmean((FCemp-FCsim).^2));
    aux=corrcoef(FCsim(Isubdiag),FCemp(Isubdiag),'Rows','complete');
    corrFC(iter)=aux(2)
    errorCOVtau(iter)=nanmean(nanmean((COVtauemp-COVtausim).^2));

    % if mod(iter,10)<0.1
    %     errornow=nanmean(nanmean((FCemp-FCsim).^2))+nanmean(nanmean((COVtauemp-COVtausim).^2));
    %     if  (olderror-errornow)/errornow<0.0001
    %         break;
    %     end
    %     if  olderror<errornow
    %         break;
    %     end
    %     olderror=errornow;
    % end

    for i=1:N  %% learning
        for j=1:N
            if (C(i,j)>0 || abs(j-i)==N/2)
                if ~isnan(FCemp(i,j))
                    Cnew(i,j)=Cnew(i,j)+epsFC*(FCemp(i,j)-FCsim(i,j)) ...
                        +epsFCtau*(COVtauemp(i,j)-COVtausim(i,j));
                end
                if Cnew(i,j)<0
                    Cnew(i,j)=0;
                end
            end
        end
    end
    Cnew = Cnew/max(max(Cnew))*maxC;

    if iter>140 
        CeffnoLR(nng,:,:)=Cnew;
        nng=nng+1;
    end
end

save results_HopfQL_CeffnoLR.mat CeffnoLR;