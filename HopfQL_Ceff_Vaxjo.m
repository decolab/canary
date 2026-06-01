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

load ('schaefer1000to100.mat');
partitiondef=schaefer1000to100;
NP=100;

for i=1:NP
    partition{i}=find(partitiondef==i);
end

N=1000;
indexN=1:N;
sigma=0.01;
NSUB=971;
Tau=3;
Isubdiag = find(tril(ones(N),-1));
IsubdiagNP = find(tril(ones(NP),-1));

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

clear FC;
vortexcenteremp=[];
for nsub=1:NSUB
    nsub
    ts=subject{nsub}.schaeferts;  % fMRI
    ts=ts(indexN,:);
    clear signal_filt Phases;
    for seed=1:N
        ts(seed,:)=detrend(ts(seed,:)-nanmean(ts(seed,:)));
        if sum(isnan(ts(seed,:)))==0
            signal_filt(seed,:)=filtfilt(bfilt,afilt,ts(seed,:));
            Xanalytic = hilbert(demean(signal_filt(seed,:)));
            Phases(seed,:) = angle(Xanalytic);
        end
    end
    Phasesemp=Phases(:,100:end-100);
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

    %%
    for i=1:NP
        Vorticity(i,:)=abs(nansum(complex(cos(Phasesemp(partition{i},:)),sin(Phasesemp(partition{i},:))))/length(partition{i}));
    end

    Tmax2=size(Vorticity,2);
    TurbulenceEmpirical(nsub)=std(Vorticity(:));

    %% Probability in Vortex Space
    for i=1:NP
        for j=1:Tmax2
            if Vorticity(i,j)<0.2
                binvortex(i,j)=0;
            end
            if Vorticity(i,j)>=0.2 && Vorticity(i,j)<0.4
                binvortex(i,j)=1;
            end
            if Vorticity(i,j)>=0.4 && Vorticity(i,j)<0.6
                binvortex(i,j)=2;
            end
            if Vorticity(i,j)>=0.6 && Vorticity(i,j)<0.8
                binvortex(i,j)=3;
            end
            if Vorticity(i,j)>=0.8
                binvortex(i,j)=4;
            end
        end
    end

    npa=1;
    for i=1:NP
        for j=1:NP
            if i~=j
                bin1=binvortex(i,:);
                bin2=binvortex(j,:);
                P1_0=length(find(bin1==0))/length(bin1);
                P1_1=length(find(bin1==1))/length(bin1);
                P1_2=length(find(bin1==2))/length(bin1);
                P1_3=length(find(bin1==3))/length(bin1);
                P1_4=length(find(bin1==4))/length(bin1);
                P2_0=length(find(bin2==0))/length(bin1);
                P2_1=length(find(bin2==1))/length(bin1);
                P2_2=length(find(bin2==2))/length(bin1);
                P2_3=length(find(bin2==3))/length(bin1);
                P2_4=length(find(bin2==4))/length(bin1);
                P12_00=length(intersect(find(bin1==0),find(bin2==0)))/length(find(bin2==0));
                P12_01=length(intersect(find(bin1==0),find(bin2==1)))/length(find(bin2==1));
                P12_02=length(intersect(find(bin1==0),find(bin2==2)))/length(find(bin2==2));
                P12_03=length(intersect(find(bin1==0),find(bin2==3)))/length(find(bin2==3));
                P12_04=length(intersect(find(bin1==0),find(bin2==4)))/length(find(bin2==4));
                P12_10=length(intersect(find(bin1==1),find(bin2==0)))/length(find(bin2==0));
                P12_11=length(intersect(find(bin1==1),find(bin2==1)))/length(find(bin2==1));
                P12_12=length(intersect(find(bin1==1),find(bin2==2)))/length(find(bin2==2));
                P12_13=length(intersect(find(bin1==1),find(bin2==3)))/length(find(bin2==3));
                P12_14=length(intersect(find(bin1==1),find(bin2==4)))/length(find(bin2==4));
                P12_20=length(intersect(find(bin1==2),find(bin2==0)))/length(find(bin2==0));
                P12_21=length(intersect(find(bin1==2),find(bin2==1)))/length(find(bin2==1));
                P12_22=length(intersect(find(bin1==2),find(bin2==2)))/length(find(bin2==2));
                P12_23=length(intersect(find(bin1==2),find(bin2==3)))/length(find(bin2==3));
                P12_24=length(intersect(find(bin1==2),find(bin2==4)))/length(find(bin2==4));
                P12_30=length(intersect(find(bin1==3),find(bin2==0)))/length(find(bin2==0));
                P12_31=length(intersect(find(bin1==3),find(bin2==1)))/length(find(bin2==1));
                P12_32=length(intersect(find(bin1==3),find(bin2==2)))/length(find(bin2==2));
                P12_33=length(intersect(find(bin1==3),find(bin2==3)))/length(find(bin2==3));
                P12_34=length(intersect(find(bin1==3),find(bin2==4)))/length(find(bin2==4));
                P12_40=length(intersect(find(bin1==4),find(bin2==0)))/length(find(bin2==0));
                P12_41=length(intersect(find(bin1==4),find(bin2==1)))/length(find(bin2==1));
                P12_42=length(intersect(find(bin1==4),find(bin2==2)))/length(find(bin2==2));
                P12_43=length(intersect(find(bin1==4),find(bin2==3)))/length(find(bin2==3));
                P12_44=length(intersect(find(bin1==4),find(bin2==4)))/length(find(bin2==4));

                deltapair(npa)=abs(P1_0-P2_0*P12_00-P2_1*P12_01-P2_2*P12_02-P2_3*P12_03-P2_4*P12_04) ...
                    /2/sqrt(P2_0*P12_00*P2_1*P12_01*P2_2*P12_02*P2_3*P12_03*P2_4*P12_04) ...
                    +abs(P1_1-P2_0*P12_10-P2_1*P12_11-P2_2*P12_12-P2_3*P12_13-P2_4*P12_14) ...
                    /2/sqrt(P2_0*P12_10*P2_1*P12_11*P2_2*P12_12*P2_3*P12_13*P2_4*P12_14) ...
                    +abs(P1_2-P2_0*P12_20-P2_1*P12_21-P2_2*P12_22-P2_3*P12_23-P2_4*P12_24) ...
                    /2/sqrt(P2_0*P12_20*P2_1*P12_21*P2_2*P12_22*P2_3*P12_23*P2_4*P12_24) ...
                    +abs(P1_3-P2_0*P12_30-P2_1*P12_31-P2_2*P12_32-P2_3*P12_33-P2_4*P12_34) ...
                    /2/sqrt(P2_0*P12_30*P2_1*P12_31*P2_2*P12_32*P2_3*P12_33*P2_4*P12_34) ...
                    +abs(P1_4-P2_0*P12_40-P2_1*P12_41-P2_2*P12_42-P2_3*P12_43-P2_4*P12_44) ...
                    /2/sqrt(P2_0*P12_40*P2_1*P12_41*P2_2*P12_42*P2_3*P12_43*P2_4*P12_44);
                npa=npa+1;
            end
        end
    end
    deltapair(deltapair==Inf)=[];
    interferencevaxjo_emp(nsub)=nanmean(deltapair);
    %%
    if nsub<101
        vortexcenteremp=[vortexcenteremp Vorticity];
    end
end

FCemp=squeeze(mean(FC));
COVtauemp=squeeze(mean(COVtau));
[auxidx Center]=kmeans(vortexcenteremp',15,'Replicates',10,'MaxIter',10000,'Display','off');

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

epsFC=0.0004;
epsFCtau=0.0001;

maxC=0.2;

C=C*maxC;

Tmax=1200
nn=1;
Cnew=C;
nnb=1;
nng=1;

olderror=100000;
for iter=1:5000
    iter
    %% Meta
    wC=Cnew;
    sumC = repmat(sum(wC,2),1,2);
    xs=zeros(Tmax,N);
    z = 0.1*ones(N,2);
    nn=0;
    % discard first 2000 time steps
    for t=0:dt:2000
        suma = wC*z - sumC.*z; % sum(Cij*xi) - sum(Cij)*xj
        zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
        z = z + dt*(a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma) + dsig*randn(N,2);
    end
    % actual modeling (x=BOLD signal (Interpretation), y some other oscillation)
    for t=0:dt:((Tmax-1)*TR)
        suma = wC*z - sumC.*z; % sum(Cij*xi) - sum(Cij)*xj
        zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
        z = z + dt*(a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma) + dsig*randn(N,2);
        if abs(mod(t,TR))<0.01
            nn=nn+1;
            xs(nn,:)=z(:,1)';
        end
    end
    ts=xs';
    clear signal_filt Phases;
    for seed=1:N
        ts(seed,:)=detrend(ts(seed,:)-mean(ts(seed,:)));
        signal_filt(seed,:) =filtfilt(bfilt,afilt,ts(seed,:));
        Xanalytic = hilbert(demean(signal_filt(seed,:)));
        Phases(seed,:) = angle(Xanalytic);
    end
    timeseries=signal_filt(:,100:end-100);
    Phasessim=Phases(:,100:end-100);

    for i=1:NP
        Vorticitysim(i,:)=abs(nansum(complex(cos(Phasessim(partition{i},:)),sin(Phasessim(partition{i},:))))/length(partition{i}));
    end
    Tmax2=size(Vorticitysim,2);
    Turbulence(iter)=std(Vorticitysim(:));

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

    %% Energy
    Theta0=COVsimtotal;
    A=-A;
    D = 0.5*(sigma^2)*eye(2*N);
    iD=inv(D);
    Theta=(Theta0+Theta0')/2;
    DTA=D*inv(Theta)-A;
    [P L]=eig(Theta);
    AP=A*P;
    TT=P'*Theta*P;
    XTX=AP*TT*AP';

    for node=1:2*N
        EntroFlow2(node)=iD(node,node)*XTX(node,node)-A(node,node);
    end
    Energy(iter)=mean(EntroFlow2);

    %%

    errorFC(iter)=nanmean(nanmean((FCemp-FCsim).^2));
    aux=corrcoef(FCsim(Isubdiag),FCemp(Isubdiag),'Rows','complete');
    corrFC(iter)=aux(2)
    errorCOVtau(iter)=nanmean(nanmean((COVtauemp-COVtausim).^2));

    %% Spectral gap and entangelment
    [V D]=eig(Cnew);
    D=diag(real(D));
    Dall=D;
    [dmax]=sort(D,'descend');
    spacings = abs(dmax);
    gapThreshold = mean(spacings) + std(spacings);
    gapIdx = find(spacings > gapThreshold); 
    sumgap=0;
    if ~isempty(gapIdx)
        for g = 1:numel(gapIdx)
            sumgap=sumgap+spacings(gapIdx(g));
        end
    end
    spectral_gap(iter)=sumgap;


    %% Coherence and entropy

    NR=15;
    IsubdiagNR = find(tril(ones(NR),-1));

    clear Probm;
    for t=1:Tmax2
        for i=1:NR
            dd(i)=norm(Center(i,:)'-Vorticitysim(:,t));
        end
        [aux id]=min(dd);
        IDX(t)=id;
    end

    for i=1:NR
        Probm(i)=length(find(IDX==i));
    end

    Probm=(1./sqrt(sum((Probm.^2)))).*Probm;

    density=abs(Probm'*Probm);
    coherence_mean(iter)=mean(density(IsubdiagNR));

    Probm=Probm.^2;
    Probm(Probm==0)=[];
    EntropyVN(iter)=-sum(Probm.*log(Probm));

    %%

    %% Vaxjol
    for i=1:NP
        for j=1:Tmax2
            if Vorticitysim(i,j)<0.2
                binvortex(i,j)=0;
            end
            if Vorticitysim(i,j)>=0.2 && Vorticitysim(i,j)<0.4
                binvortex(i,j)=1;
            end
            if Vorticitysim(i,j)>=0.4 && Vorticitysim(i,j)<0.6
                binvortex(i,j)=2;
            end
            if Vorticitysim(i,j)>=0.6 && Vorticitysim(i,j)<0.8
                binvortex(i,j)=3;
            end
            if Vorticitysim(i,j)>=0.8
                binvortex(i,j)=4;
            end
        end
    end

    npa=1;
    for i=1:NP
        for j=1:NP
            if i~=j
                bin1=binvortex(i,:);
                bin2=binvortex(j,:);

                P1_0=length(find(bin1==0))/length(bin1);
                P1_1=length(find(bin1==1))/length(bin1);
                P1_2=length(find(bin1==2))/length(bin1);
                P1_3=length(find(bin1==3))/length(bin1);
                P1_4=length(find(bin1==4))/length(bin1);
                P2_0=length(find(bin2==0))/length(bin1);
                P2_1=length(find(bin2==1))/length(bin1);
                P2_2=length(find(bin2==2))/length(bin1);
                P2_3=length(find(bin2==3))/length(bin1);
                P2_4=length(find(bin2==4))/length(bin1);
                P12_00=length(intersect(find(bin1==0),find(bin2==0)))/length(find(bin2==0));
                P12_01=length(intersect(find(bin1==0),find(bin2==1)))/length(find(bin2==1));
                P12_02=length(intersect(find(bin1==0),find(bin2==2)))/length(find(bin2==2));
                P12_03=length(intersect(find(bin1==0),find(bin2==3)))/length(find(bin2==3));
                P12_04=length(intersect(find(bin1==0),find(bin2==4)))/length(find(bin2==4));
                P12_10=length(intersect(find(bin1==1),find(bin2==0)))/length(find(bin2==0));
                P12_11=length(intersect(find(bin1==1),find(bin2==1)))/length(find(bin2==1));
                P12_12=length(intersect(find(bin1==1),find(bin2==2)))/length(find(bin2==2));
                P12_13=length(intersect(find(bin1==1),find(bin2==3)))/length(find(bin2==3));
                P12_14=length(intersect(find(bin1==1),find(bin2==4)))/length(find(bin2==4));
                P12_20=length(intersect(find(bin1==2),find(bin2==0)))/length(find(bin2==0));
                P12_21=length(intersect(find(bin1==2),find(bin2==1)))/length(find(bin2==1));
                P12_22=length(intersect(find(bin1==2),find(bin2==2)))/length(find(bin2==2));
                P12_23=length(intersect(find(bin1==2),find(bin2==3)))/length(find(bin2==3));
                P12_24=length(intersect(find(bin1==2),find(bin2==4)))/length(find(bin2==4));
                P12_30=length(intersect(find(bin1==3),find(bin2==0)))/length(find(bin2==0));
                P12_31=length(intersect(find(bin1==3),find(bin2==1)))/length(find(bin2==1));
                P12_32=length(intersect(find(bin1==3),find(bin2==2)))/length(find(bin2==2));
                P12_33=length(intersect(find(bin1==3),find(bin2==3)))/length(find(bin2==3));
                P12_34=length(intersect(find(bin1==3),find(bin2==4)))/length(find(bin2==4));
                P12_40=length(intersect(find(bin1==4),find(bin2==0)))/length(find(bin2==0));
                P12_41=length(intersect(find(bin1==4),find(bin2==1)))/length(find(bin2==1));
                P12_42=length(intersect(find(bin1==4),find(bin2==2)))/length(find(bin2==2));
                P12_43=length(intersect(find(bin1==4),find(bin2==3)))/length(find(bin2==3));
                P12_44=length(intersect(find(bin1==4),find(bin2==4)))/length(find(bin2==4));

                deltapair(npa)=abs(P1_0-P2_0*P12_00-P2_1*P12_01-P2_2*P12_02-P2_3*P12_03-P2_4*P12_04) ...
                    /2/sqrt(P2_0*P12_00*P2_1*P12_01*P2_2*P12_02*P2_3*P12_03*P2_4*P12_04) ...
                    +abs(P1_1-P2_0*P12_10-P2_1*P12_11-P2_2*P12_12-P2_3*P12_13-P2_4*P12_14) ...
                    /2/sqrt(P2_0*P12_10*P2_1*P12_11*P2_2*P12_12*P2_3*P12_13*P2_4*P12_14) ...
                    +abs(P1_2-P2_0*P12_20-P2_1*P12_21-P2_2*P12_22-P2_3*P12_23-P2_4*P12_24) ...
                    /2/sqrt(P2_0*P12_20*P2_1*P12_21*P2_2*P12_22*P2_3*P12_23*P2_4*P12_24) ...
                    +abs(P1_3-P2_0*P12_30-P2_1*P12_31-P2_2*P12_32-P2_3*P12_33-P2_4*P12_34) ...
                    /2/sqrt(P2_0*P12_30*P2_1*P12_31*P2_2*P12_32*P2_3*P12_33*P2_4*P12_34) ...
                    +abs(P1_4-P2_0*P12_40-P2_1*P12_41-P2_2*P12_42-P2_3*P12_43-P2_4*P12_44) ...
                    /2/sqrt(P2_0*P12_40*P2_1*P12_41*P2_2*P12_42*P2_3*P12_43*P2_4*P12_44);

                Prob=[P12_00 P12_10 P12_20 P12_30 P12_40];
                Prob(Prob==0)=[];
                Entropy0=-sum(Prob.*log10(Prob));
                Prob=[P12_01 P12_11 P12_21 P12_31 P12_41];
                Prob(Prob==0)=[];
                Entropy1=-sum(Prob.*log10(Prob));
                Prob=[P12_02 P12_12 P12_22 P12_32 P12_42];
                Prob(Prob==0)=[];
                Entropy2=-sum(Prob.*log10(Prob));
                Prob=[P12_03 P12_13 P12_23 P12_33 P12_43];
                Prob(Prob==0)=[];
                Entropy3=-sum(Prob.*log10(Prob));
                Prob=[P12_04 P12_14 P12_24 P12_34 P12_44];
                Prob(Prob==0)=[];
                Entropy4=-sum(Prob.*log10(Prob));

                Entropy(npa)=(Entropy0*P2_0+Entropy1*P2_1+Entropy2*P2_2+Entropy3*P2_3+Entropy4*P2_4);
                npa=npa+1;
            end
        end
    end
    deltapair(deltapair==Inf)=[];
    interferencevaxjo(iter)=nanmean(deltapair);
    uncertainty(iter)=nanmean(Entropy);

    %%

    if mod(iter,10)<0.1
        errornow=nanmean(nanmean((FCemp-FCsim).^2))+nanmean(nanmean((COVtauemp-COVtausim).^2));
        if  (olderror-errornow)/errornow<0.0001
            break;
        end
        if  olderror<errornow
            break;
        end
        olderror=errornow;
    end

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
    if iter<21
        Cbad(nnb,:,:)=Cnew;
        nnb=nnb+1;
    end
    if iter>140 
        Cgood(nng,:,:)=Cnew;
        nng=nng+1;
    end
end
Crest=Cnew;

figure(1)
subplot(4,2,1)
plot(errorFC)
subplot(4,2,2)
plot(corrFC)
subplot(4,2,3)
plot(interferencevaxjo)
subplot(4,2,4)
plot(uncertainty)
subplot(4,2,5)
plot(coherence_mean)
subplot(4,2,6)
plot(EntropyVN)
subplot(4,2,7)
plot(spectral_gap)
subplot(4,2,8)
plot(Turbulence);
hold on;
plot(mean(TurbulenceEmpirical)*ones(1,length(Turbulence)));

figure(2)
plot(Energy)

%%
for nsub=1:NSUB
    nsub
    nsel=randperm(20);
    C=squeeze(Cbad(nsel(1),:,:));
    wC=C;
    sumC = repmat(sum(wC,2),1,2);
    xs=zeros(Tmax,N);
    z = 0.1*ones(N,2);
    nn=0;
    % discard first 2000 time steps
    for t=0:dt:2000
        suma = wC*z - sumC.*z; % sum(Cij*xi) - sum(Cij)*xj
        zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
        z = z + dt*(a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma) + dsig*randn(N,2);
    end
    % actual modeling (x=BOLD signal (Interpretation), y some other oscillation)
    for t=0:dt:((Tmax-1)*TR)
        suma = wC*z - sumC.*z; % sum(Cij*xi) - sum(Cij)*xj
        zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
        z = z + dt*(a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma) + dsig*randn(N,2);
        if abs(mod(t,TR))<0.01
            nn=nn+1;
            xs(nn,:)=z(:,1)';
        end
    end
    ts=xs';
    clear signal_filt Phases;
    for seed=1:N
        ts(seed,:)=detrend(ts(seed,:)-mean(ts(seed,:)));
        signal_filt(seed,:) =filtfilt(bfilt,afilt,ts(seed,:));
        Xanalytic = hilbert(demean(signal_filt(seed,:)));
        Phases(seed,:) = angle(Xanalytic);
    end
    Phasessim=Phases(:,100:end-100);

    for i=1:NP
        Vorticitysim(i,:)=abs(nansum(complex(cos(Phasessim(partition{i},:)),sin(Phasessim(partition{i},:))))/length(partition{i}));
    end
    Tmax2=size(Vorticitysim,2);

    %% Energy
    [FCsim,COVsim,COVsimtotal,A]=hopf_int(C,f_diff,sigma);

    Theta0=COVsimtotal;
    A=-A;
    D = 0.5*(sigma^2)*eye(2*N);
    iD=inv(D);
    Theta=(Theta0+Theta0')/2;
    DTA=D*inv(Theta)-A;
    [P L]=eig(Theta);
    AP=A*P;
    TT=P'*Theta*P;
    XTX=AP*TT*AP';

    for node=1:2*N
        EntroFlow2(node)=iD(node,node)*XTX(node,node)-A(node,node);
    end
    Energy_bad(nsub)=mean(EntroFlow2);

    %% Spectral gap and entangelment
    [V D]=eig(C);
    D=diag(real(D));
    Dall=D;
    [dmax]=sort(D,'descend');
    spacings = abs(dmax);
    gapThreshold = mean(spacings) + std(spacings);
    gapIdx = find(spacings > gapThreshold); 
    sumgap=0;
    if ~isempty(gapIdx)
        for g = 1:numel(gapIdx)
            sumgap=sumgap+spacings(gapIdx(g));
        end
    end
    spectral_gap_bad(nsub)=sumgap;


    %% Coherence and entropy

    NR=15;
    IsubdiagNR = find(tril(ones(NR),-1));

    clear Probm;
    for t=1:Tmax2
        for i=1:NR
            dd(i)=norm(Center(i,:)'-Vorticitysim(:,t));
        end
        [aux id]=min(dd);
        IDX(t)=id;
    end

    for i=1:NR
        Probm(i)=length(find(IDX==i));
    end

    Probm=(1./sqrt(sum((Probm.^2)))).*Probm;

    density=abs(Probm'*Probm);
    coherence_mean_bad(nsub)=mean(density(IsubdiagNR));

    Probm=Probm.^2;
    Probm(Probm==0)=[];
    EntropyVN_bad(nsub)=-sum(Probm.*log(Probm));

    %%
    %% Vaxjol
    for i=1:NP
        for j=1:Tmax2
            if Vorticitysim(i,j)<0.2
                binvortex(i,j)=0;
            end
            if Vorticitysim(i,j)>=0.2 && Vorticitysim(i,j)<0.4
                binvortex(i,j)=1;
            end
            if Vorticitysim(i,j)>=0.4 && Vorticitysim(i,j)<0.6
                binvortex(i,j)=2;
            end
            if Vorticitysim(i,j)>=0.6 && Vorticitysim(i,j)<0.8
                binvortex(i,j)=3;
            end
            if Vorticitysim(i,j)>=0.8
                binvortex(i,j)=4;
            end
        end
    end

    npa=1;
    for i=1:NP
        for j=1:NP
            if i~=j
                bin1=binvortex(i,:);
                bin2=binvortex(j,:);

                P1_0=length(find(bin1==0))/length(bin1);
                P1_1=length(find(bin1==1))/length(bin1);
                P1_2=length(find(bin1==2))/length(bin1);
                P1_3=length(find(bin1==3))/length(bin1);
                P1_4=length(find(bin1==4))/length(bin1);
                P2_0=length(find(bin2==0))/length(bin1);
                P2_1=length(find(bin2==1))/length(bin1);
                P2_2=length(find(bin2==2))/length(bin1);
                P2_3=length(find(bin2==3))/length(bin1);
                P2_4=length(find(bin2==4))/length(bin1);
                P12_00=length(intersect(find(bin1==0),find(bin2==0)))/length(find(bin2==0));
                P12_01=length(intersect(find(bin1==0),find(bin2==1)))/length(find(bin2==1));
                P12_02=length(intersect(find(bin1==0),find(bin2==2)))/length(find(bin2==2));
                P12_03=length(intersect(find(bin1==0),find(bin2==3)))/length(find(bin2==3));
                P12_04=length(intersect(find(bin1==0),find(bin2==4)))/length(find(bin2==4));
                P12_10=length(intersect(find(bin1==1),find(bin2==0)))/length(find(bin2==0));
                P12_11=length(intersect(find(bin1==1),find(bin2==1)))/length(find(bin2==1));
                P12_12=length(intersect(find(bin1==1),find(bin2==2)))/length(find(bin2==2));
                P12_13=length(intersect(find(bin1==1),find(bin2==3)))/length(find(bin2==3));
                P12_14=length(intersect(find(bin1==1),find(bin2==4)))/length(find(bin2==4));
                P12_20=length(intersect(find(bin1==2),find(bin2==0)))/length(find(bin2==0));
                P12_21=length(intersect(find(bin1==2),find(bin2==1)))/length(find(bin2==1));
                P12_22=length(intersect(find(bin1==2),find(bin2==2)))/length(find(bin2==2));
                P12_23=length(intersect(find(bin1==2),find(bin2==3)))/length(find(bin2==3));
                P12_24=length(intersect(find(bin1==2),find(bin2==4)))/length(find(bin2==4));
                P12_30=length(intersect(find(bin1==3),find(bin2==0)))/length(find(bin2==0));
                P12_31=length(intersect(find(bin1==3),find(bin2==1)))/length(find(bin2==1));
                P12_32=length(intersect(find(bin1==3),find(bin2==2)))/length(find(bin2==2));
                P12_33=length(intersect(find(bin1==3),find(bin2==3)))/length(find(bin2==3));
                P12_34=length(intersect(find(bin1==3),find(bin2==4)))/length(find(bin2==4));
                P12_40=length(intersect(find(bin1==4),find(bin2==0)))/length(find(bin2==0));
                P12_41=length(intersect(find(bin1==4),find(bin2==1)))/length(find(bin2==1));
                P12_42=length(intersect(find(bin1==4),find(bin2==2)))/length(find(bin2==2));
                P12_43=length(intersect(find(bin1==4),find(bin2==3)))/length(find(bin2==3));
                P12_44=length(intersect(find(bin1==4),find(bin2==4)))/length(find(bin2==4));

                deltapair(npa)=abs(P1_0-P2_0*P12_00-P2_1*P12_01-P2_2*P12_02-P2_3*P12_03-P2_4*P12_04) ...
                    /2/sqrt(P2_0*P12_00*P2_1*P12_01*P2_2*P12_02*P2_3*P12_03*P2_4*P12_04) ...
                    +abs(P1_1-P2_0*P12_10-P2_1*P12_11-P2_2*P12_12-P2_3*P12_13-P2_4*P12_14) ...
                    /2/sqrt(P2_0*P12_10*P2_1*P12_11*P2_2*P12_12*P2_3*P12_13*P2_4*P12_14) ...
                    +abs(P1_2-P2_0*P12_20-P2_1*P12_21-P2_2*P12_22-P2_3*P12_23-P2_4*P12_24) ...
                    /2/sqrt(P2_0*P12_20*P2_1*P12_21*P2_2*P12_22*P2_3*P12_23*P2_4*P12_24) ...
                    +abs(P1_3-P2_0*P12_30-P2_1*P12_31-P2_2*P12_32-P2_3*P12_33-P2_4*P12_34) ...
                    /2/sqrt(P2_0*P12_30*P2_1*P12_31*P2_2*P12_32*P2_3*P12_33*P2_4*P12_34) ...
                    +abs(P1_4-P2_0*P12_40-P2_1*P12_41-P2_2*P12_42-P2_3*P12_43-P2_4*P12_44) ...
                    /2/sqrt(P2_0*P12_40*P2_1*P12_41*P2_2*P12_42*P2_3*P12_43*P2_4*P12_44);

                Prob=[P12_00 P12_10 P12_20 P12_30 P12_40];
                Prob(Prob==0)=[];
                Entropy0=-sum(Prob.*log10(Prob));
                Prob=[P12_01 P12_11 P12_21 P12_31 P12_41];
                Prob(Prob==0)=[];
                Entropy1=-sum(Prob.*log10(Prob));
                Prob=[P12_02 P12_12 P12_22 P12_32 P12_42];
                Prob(Prob==0)=[];
                Entropy2=-sum(Prob.*log10(Prob));
                Prob=[P12_03 P12_13 P12_23 P12_33 P12_43];
                Prob(Prob==0)=[];
                Entropy3=-sum(Prob.*log10(Prob));
                Prob=[P12_04 P12_14 P12_24 P12_34 P12_44];
                Prob(Prob==0)=[];
                Entropy4=-sum(Prob.*log10(Prob));

                Entropy(npa)=(Entropy0*P2_0+Entropy1*P2_1+Entropy2*P2_2+Entropy3*P2_3+Entropy4*P2_4);
                npa=npa+1;
            end
        end
    end
    deltapair(deltapair==Inf)=[];
    interferencevaxjo_bad(nsub)=nanmean(deltapair);
    uncertainty_bad(nsub)=nanmean(Entropy);
end

%%

for nsub=1:NSUB
    nsub
    nsel=randperm(20);
    C=squeeze(Cgood(nsel(1),:,:));
    wC=C;
    sumC = repmat(sum(wC,2),1,2);
    xs=zeros(Tmax,N);
    z = 0.1*ones(N,2);
    nn=0;
    % discard first 2000 time steps
    for t=0:dt:2000
        suma = wC*z - sumC.*z; % sum(Cij*xi) - sum(Cij)*xj
        zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
        z = z + dt*(a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma) + dsig*randn(N,2);
    end
    % actual modeling (x=BOLD signal (Interpretation), y some other oscillation)
    for t=0:dt:((Tmax-1)*TR)
        suma = wC*z - sumC.*z; % sum(Cij*xi) - sum(Cij)*xj
        zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
        z = z + dt*(a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma) + dsig*randn(N,2);
        if abs(mod(t,TR))<0.01
            nn=nn+1;
            xs(nn,:)=z(:,1)';
        end
    end
    ts=xs';
    clear signal_filt Phases;
    for seed=1:N
        ts(seed,:)=detrend(ts(seed,:)-mean(ts(seed,:)));
        signal_filt(seed,:) =filtfilt(bfilt,afilt,ts(seed,:));
        Xanalytic = hilbert(demean(signal_filt(seed,:)));
        Phases(seed,:) = angle(Xanalytic);
    end
    Phasessim=Phases(:,100:end-100);

    for i=1:NP
        Vorticitysim(i,:)=abs(nansum(complex(cos(Phasessim(partition{i},:)),sin(Phasessim(partition{i},:))))/length(partition{i}));
    end
    Tmax2=size(Vorticitysim,2);

        %% Energy
    [FCsim,COVsim,COVsimtotal,A]=hopf_int(C,f_diff,sigma);

    Theta0=COVsimtotal;
    A=-A;
    D = 0.5*(sigma^2)*eye(2*N);
    iD=inv(D);
    Theta=(Theta0+Theta0')/2;
    DTA=D*inv(Theta)-A;
    [P L]=eig(Theta);
    AP=A*P;
    TT=P'*Theta*P;
    XTX=AP*TT*AP';

    for node=1:2*N
        EntroFlow2(node)=iD(node,node)*XTX(node,node)-A(node,node);
    end
    Energy_good(nsub)=mean(EntroFlow2);

    %% Spectral gap and entangelment
    [V D]=eig(C);
    D=diag(real(D));
    Dall=D;
    [dmax]=sort(D,'descend');
    spacings = abs(dmax);
    gapThreshold = mean(spacings) + std(spacings);
    gapIdx = find(spacings > gapThreshold); 
    sumgap=0;
    if ~isempty(gapIdx)
        for g = 1:numel(gapIdx)
            sumgap=sumgap+spacings(gapIdx(g));
        end
    end
    spectral_gap_good(nsub)=sumgap;


    %% Coherence and entropy

    NR=15;
    IsubdiagNR = find(tril(ones(NR),-1));

    clear Probm;
    for t=1:Tmax2
        for i=1:NR
            dd(i)=norm(Center(i,:)'-Vorticitysim(:,t));
        end
        [aux id]=min(dd);
        IDX(t)=id;
    end

    for i=1:NR
        Probm(i)=length(find(IDX==i));
    end

    Probm=(1./sqrt(sum((Probm.^2)))).*Probm;

    density=abs(Probm'*Probm);
    coherence_mean_good(nsub)=mean(density(IsubdiagNR));

    Probm=Probm.^2;
    Probm(Probm==0)=[];
    EntropyVN_good(nsub)=-sum(Probm.*log(Probm));

    %%
    %% Vaxjol
    for i=1:NP
        for j=1:Tmax2
            if Vorticitysim(i,j)<0.2
                binvortex(i,j)=0;
            end
            if Vorticitysim(i,j)>=0.2 && Vorticitysim(i,j)<0.4
                binvortex(i,j)=1;
            end
            if Vorticitysim(i,j)>=0.4 && Vorticitysim(i,j)<0.6
                binvortex(i,j)=2;
            end
            if Vorticitysim(i,j)>=0.6 && Vorticitysim(i,j)<0.8
                binvortex(i,j)=3;
            end
            if Vorticitysim(i,j)>=0.8
                binvortex(i,j)=4;
            end
        end
    end

    npa=1;
    for i=1:NP
        for j=1:NP
            if i~=j
                bin1=binvortex(i,:);
                bin2=binvortex(j,:);

                P1_0=length(find(bin1==0))/length(bin1);
                P1_1=length(find(bin1==1))/length(bin1);
                P1_2=length(find(bin1==2))/length(bin1);
                P1_3=length(find(bin1==3))/length(bin1);
                P1_4=length(find(bin1==4))/length(bin1);
                P2_0=length(find(bin2==0))/length(bin1);
                P2_1=length(find(bin2==1))/length(bin1);
                P2_2=length(find(bin2==2))/length(bin1);
                P2_3=length(find(bin2==3))/length(bin1);
                P2_4=length(find(bin2==4))/length(bin1);
                P12_00=length(intersect(find(bin1==0),find(bin2==0)))/length(find(bin2==0));
                P12_01=length(intersect(find(bin1==0),find(bin2==1)))/length(find(bin2==1));
                P12_02=length(intersect(find(bin1==0),find(bin2==2)))/length(find(bin2==2));
                P12_03=length(intersect(find(bin1==0),find(bin2==3)))/length(find(bin2==3));
                P12_04=length(intersect(find(bin1==0),find(bin2==4)))/length(find(bin2==4));
                P12_10=length(intersect(find(bin1==1),find(bin2==0)))/length(find(bin2==0));
                P12_11=length(intersect(find(bin1==1),find(bin2==1)))/length(find(bin2==1));
                P12_12=length(intersect(find(bin1==1),find(bin2==2)))/length(find(bin2==2));
                P12_13=length(intersect(find(bin1==1),find(bin2==3)))/length(find(bin2==3));
                P12_14=length(intersect(find(bin1==1),find(bin2==4)))/length(find(bin2==4));
                P12_20=length(intersect(find(bin1==2),find(bin2==0)))/length(find(bin2==0));
                P12_21=length(intersect(find(bin1==2),find(bin2==1)))/length(find(bin2==1));
                P12_22=length(intersect(find(bin1==2),find(bin2==2)))/length(find(bin2==2));
                P12_23=length(intersect(find(bin1==2),find(bin2==3)))/length(find(bin2==3));
                P12_24=length(intersect(find(bin1==2),find(bin2==4)))/length(find(bin2==4));
                P12_30=length(intersect(find(bin1==3),find(bin2==0)))/length(find(bin2==0));
                P12_31=length(intersect(find(bin1==3),find(bin2==1)))/length(find(bin2==1));
                P12_32=length(intersect(find(bin1==3),find(bin2==2)))/length(find(bin2==2));
                P12_33=length(intersect(find(bin1==3),find(bin2==3)))/length(find(bin2==3));
                P12_34=length(intersect(find(bin1==3),find(bin2==4)))/length(find(bin2==4));
                P12_40=length(intersect(find(bin1==4),find(bin2==0)))/length(find(bin2==0));
                P12_41=length(intersect(find(bin1==4),find(bin2==1)))/length(find(bin2==1));
                P12_42=length(intersect(find(bin1==4),find(bin2==2)))/length(find(bin2==2));
                P12_43=length(intersect(find(bin1==4),find(bin2==3)))/length(find(bin2==3));
                P12_44=length(intersect(find(bin1==4),find(bin2==4)))/length(find(bin2==4));

                deltapair(npa)=abs(P1_0-P2_0*P12_00-P2_1*P12_01-P2_2*P12_02-P2_3*P12_03-P2_4*P12_04) ...
                    /2/sqrt(P2_0*P12_00*P2_1*P12_01*P2_2*P12_02*P2_3*P12_03*P2_4*P12_04) ...
                    +abs(P1_1-P2_0*P12_10-P2_1*P12_11-P2_2*P12_12-P2_3*P12_13-P2_4*P12_14) ...
                    /2/sqrt(P2_0*P12_10*P2_1*P12_11*P2_2*P12_12*P2_3*P12_13*P2_4*P12_14) ...
                    +abs(P1_2-P2_0*P12_20-P2_1*P12_21-P2_2*P12_22-P2_3*P12_23-P2_4*P12_24) ...
                    /2/sqrt(P2_0*P12_20*P2_1*P12_21*P2_2*P12_22*P2_3*P12_23*P2_4*P12_24) ...
                    +abs(P1_3-P2_0*P12_30-P2_1*P12_31-P2_2*P12_32-P2_3*P12_33-P2_4*P12_34) ...
                    /2/sqrt(P2_0*P12_30*P2_1*P12_31*P2_2*P12_32*P2_3*P12_33*P2_4*P12_34) ...
                    +abs(P1_4-P2_0*P12_40-P2_1*P12_41-P2_2*P12_42-P2_3*P12_43-P2_4*P12_44) ...
                    /2/sqrt(P2_0*P12_40*P2_1*P12_41*P2_2*P12_42*P2_3*P12_43*P2_4*P12_44);

                Prob=[P12_00 P12_10 P12_20 P12_30 P12_40];
                Prob(Prob==0)=[];
                Entropy0=-sum(Prob.*log10(Prob));
                Prob=[P12_01 P12_11 P12_21 P12_31 P12_41];
                Prob(Prob==0)=[];
                Entropy1=-sum(Prob.*log10(Prob));
                Prob=[P12_02 P12_12 P12_22 P12_32 P12_42];
                Prob(Prob==0)=[];
                Entropy2=-sum(Prob.*log10(Prob));
                Prob=[P12_03 P12_13 P12_23 P12_33 P12_43];
                Prob(Prob==0)=[];
                Entropy3=-sum(Prob.*log10(Prob));
                Prob=[P12_04 P12_14 P12_24 P12_34 P12_44];
                Prob(Prob==0)=[];
                Entropy4=-sum(Prob.*log10(Prob));

                Entropy(npa)=(Entropy0*P2_0+Entropy1*P2_1+Entropy2*P2_2+Entropy3*P2_3+Entropy4*P2_4);
                npa=npa+1;
            end
        end
    end
    deltapair(deltapair==Inf)=[];
    interferencevaxjo_good(nsub)=nanmean(deltapair);
    uncertainty_good(nsub)=nanmean(Entropy);
end


%% GBCs for rendering

for nsub=1:100
    nsub
    C=Crest;
    wC=C;
    sumC = repmat(sum(wC,2),1,2);
    xs=zeros(Tmax,N);
    z = 0.1*ones(N,2);
    nn=0;
    % discard first 2000 time steps
    for t=0:dt:2000
        suma = wC*z - sumC.*z; % sum(Cij*xi) - sum(Cij)*xj
        zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
        z = z + dt*(a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma) + dsig*randn(N,2);
    end
    % actual modeling (x=BOLD signal (Interpretation), y some other oscillation)
    for t=0:dt:((Tmax-1)*TR)
        suma = wC*z - sumC.*z; % sum(Cij*xi) - sum(Cij)*xj
        zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
        z = z + dt*(a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma) + dsig*randn(N,2);
        if abs(mod(t,TR))<0.01
            nn=nn+1;
            xs(nn,:)=z(:,1)';
        end
    end
    ts=xs';
    clear signal_filt Phases;
    for seed=1:N
        ts(seed,:)=detrend(ts(seed,:)-mean(ts(seed,:)));
        signal_filt(seed,:) =filtfilt(bfilt,afilt,ts(seed,:));
        Xanalytic = hilbert(demean(signal_filt(seed,:)));
        Phases(seed,:) = angle(Xanalytic);
    end
    Phasessim=Phases(:,100:end-100);

    for i=1:NP
        Vorticitysim(i,:)=abs(nansum(complex(cos(Phasessim(partition{i},:)),sin(Phasessim(partition{i},:))))/length(partition{i}));
    end
    Tmax2=size(Vorticitysim,2);

    %%
    %% Vaxjol
    for i=1:NP
        for j=1:Tmax2
            if Vorticitysim(i,j)<0.2
                binvortex(i,j)=0;
            end
            if Vorticitysim(i,j)>=0.2 && Vorticitysim(i,j)<0.4
                binvortex(i,j)=1;
            end
            if Vorticitysim(i,j)>=0.4 && Vorticitysim(i,j)<0.6
                binvortex(i,j)=2;
            end
            if Vorticitysim(i,j)>=0.6 && Vorticitysim(i,j)<0.8
                binvortex(i,j)=3;
            end
            if Vorticitysim(i,j)>=0.8
                binvortex(i,j)=4;
            end
        end
    end

    npa=1;
    for i=1:NP
        for j=1:NP
            if i~=j
                bin1=binvortex(i,:);
                bin2=binvortex(j,:);

                P1_0=length(find(bin1==0))/length(bin1);
                P1_1=length(find(bin1==1))/length(bin1);
                P1_2=length(find(bin1==2))/length(bin1);
                P1_3=length(find(bin1==3))/length(bin1);
                P1_4=length(find(bin1==4))/length(bin1);
                P2_0=length(find(bin2==0))/length(bin1);
                P2_1=length(find(bin2==1))/length(bin1);
                P2_2=length(find(bin2==2))/length(bin1);
                P2_3=length(find(bin2==3))/length(bin1);
                P2_4=length(find(bin2==4))/length(bin1);
                P12_00=length(intersect(find(bin1==0),find(bin2==0)))/length(find(bin2==0));
                P12_01=length(intersect(find(bin1==0),find(bin2==1)))/length(find(bin2==1));
                P12_02=length(intersect(find(bin1==0),find(bin2==2)))/length(find(bin2==2));
                P12_03=length(intersect(find(bin1==0),find(bin2==3)))/length(find(bin2==3));
                P12_04=length(intersect(find(bin1==0),find(bin2==4)))/length(find(bin2==4));
                P12_10=length(intersect(find(bin1==1),find(bin2==0)))/length(find(bin2==0));
                P12_11=length(intersect(find(bin1==1),find(bin2==1)))/length(find(bin2==1));
                P12_12=length(intersect(find(bin1==1),find(bin2==2)))/length(find(bin2==2));
                P12_13=length(intersect(find(bin1==1),find(bin2==3)))/length(find(bin2==3));
                P12_14=length(intersect(find(bin1==1),find(bin2==4)))/length(find(bin2==4));
                P12_20=length(intersect(find(bin1==2),find(bin2==0)))/length(find(bin2==0));
                P12_21=length(intersect(find(bin1==2),find(bin2==1)))/length(find(bin2==1));
                P12_22=length(intersect(find(bin1==2),find(bin2==2)))/length(find(bin2==2));
                P12_23=length(intersect(find(bin1==2),find(bin2==3)))/length(find(bin2==3));
                P12_24=length(intersect(find(bin1==2),find(bin2==4)))/length(find(bin2==4));
                P12_30=length(intersect(find(bin1==3),find(bin2==0)))/length(find(bin2==0));
                P12_31=length(intersect(find(bin1==3),find(bin2==1)))/length(find(bin2==1));
                P12_32=length(intersect(find(bin1==3),find(bin2==2)))/length(find(bin2==2));
                P12_33=length(intersect(find(bin1==3),find(bin2==3)))/length(find(bin2==3));
                P12_34=length(intersect(find(bin1==3),find(bin2==4)))/length(find(bin2==4));
                P12_40=length(intersect(find(bin1==4),find(bin2==0)))/length(find(bin2==0));
                P12_41=length(intersect(find(bin1==4),find(bin2==1)))/length(find(bin2==1));
                P12_42=length(intersect(find(bin1==4),find(bin2==2)))/length(find(bin2==2));
                P12_43=length(intersect(find(bin1==4),find(bin2==3)))/length(find(bin2==3));
                P12_44=length(intersect(find(bin1==4),find(bin2==4)))/length(find(bin2==4));

                deltapair(npa)=abs(P1_0-P2_0*P12_00-P2_1*P12_01-P2_2*P12_02-P2_3*P12_03-P2_4*P12_04) ...
                    /2/sqrt(P2_0*P12_00*P2_1*P12_01*P2_2*P12_02*P2_3*P12_03*P2_4*P12_04) ...
                    +abs(P1_1-P2_0*P12_10-P2_1*P12_11-P2_2*P12_12-P2_3*P12_13-P2_4*P12_14) ...
                    /2/sqrt(P2_0*P12_10*P2_1*P12_11*P2_2*P12_12*P2_3*P12_13*P2_4*P12_14) ...
                    +abs(P1_2-P2_0*P12_20-P2_1*P12_21-P2_2*P12_22-P2_3*P12_23-P2_4*P12_24) ...
                    /2/sqrt(P2_0*P12_20*P2_1*P12_21*P2_2*P12_22*P2_3*P12_23*P2_4*P12_24) ...
                    +abs(P1_3-P2_0*P12_30-P2_1*P12_31-P2_2*P12_32-P2_3*P12_33-P2_4*P12_34) ...
                    /2/sqrt(P2_0*P12_30*P2_1*P12_31*P2_2*P12_32*P2_3*P12_33*P2_4*P12_34) ...
                    +abs(P1_4-P2_0*P12_40-P2_1*P12_41-P2_2*P12_42-P2_3*P12_43-P2_4*P12_44) ...
                    /2/sqrt(P2_0*P12_40*P2_1*P12_41*P2_2*P12_42*P2_3*P12_43*P2_4*P12_44);

                Prob=[P12_00 P12_10 P12_20 P12_30 P12_40];
                Prob(Prob==0)=[];
                Entropy0=-sum(Prob.*log10(Prob));
                Prob=[P12_01 P12_11 P12_21 P12_31 P12_41];
                Prob(Prob==0)=[];
                Entropy1=-sum(Prob.*log10(Prob));
                Prob=[P12_02 P12_12 P12_22 P12_32 P12_42];
                Prob(Prob==0)=[];
                Entropy2=-sum(Prob.*log10(Prob));
                Prob=[P12_03 P12_13 P12_23 P12_33 P12_43];
                Prob(Prob==0)=[];
                Entropy3=-sum(Prob.*log10(Prob));
                Prob=[P12_04 P12_14 P12_24 P12_34 P12_44];
                Prob(Prob==0)=[];
                Entropy4=-sum(Prob.*log10(Prob));

                Entropy(npa)=(Entropy0*P2_0+Entropy1*P2_1+Entropy2*P2_2+Entropy3*P2_3+Entropy4*P2_4);
                if deltapair(npa)==Inf
                    Matrix_Vaxjo(nsub,i,j)=NaN;
                else
                    Matrix_Vaxjo(nsub,i,j)=deltapair(npa);
                end
                if Entropy(npa)==Inf
                    Matrix_Uncertainty(nsub,i,j)=NaN;
                else
                    Matrix_Uncertainty(nsub,i,j)=Entropy(npa);
                end
                npa=npa+1;
            end
        end
    end
end
Matrix_Vaxjo=squeeze(nanmean(Matrix_Vaxjo));
GBC_Vaxjo=nanmean(Matrix_Vaxjo);
Matrix_Uncertainty=squeeze(nanmean(Matrix_Uncertainty));
GBC_Uncertainty=nanmean(Matrix_Uncertainty);

%%


figure(3)
boxplot([((interferencevaxjo_bad-mean(interferencevaxjo_emp)).^2)' ((interferencevaxjo_good-mean(interferencevaxjo_emp)).^2)']);
a=(interferencevaxjo_bad-mean(interferencevaxjo_emp)).^2;
b=(interferencevaxjo_good-mean(interferencevaxjo_emp)).^2;
pp=ranksum(a,b)

figure(4)
boxplot([interferencevaxjo_bad' interferencevaxjo_good']);
a=interferencevaxjo_bad;
b=interferencevaxjo_good;
pp=ranksum(a,b)

load results_EmpiricalUncertainty.mat;

gg1=(uncertainty_bad-mean(UncertaintyEmpirical)).^2;
[gg1 idx]=rmoutliers(gg1,"ThresHold",0.9);
gg2=(uncertainty_good-mean(UncertaintyEmpirical)).^2;
[gg2 idx]=rmoutliers(gg2,"ThresHold",0.9);
MIN=min(length(gg1),length(gg2));
gg1=gg1(1:MIN);
gg2=gg2(1:MIN);

figure(5)
boxplot([gg1' gg2']);
a=gg1;
b=gg2;
pp=ranksum(a,b)

gg1=uncertainty_bad;
[gg1 idx]=rmoutliers(gg1,"ThresHold",0.9);
gg2=uncertainty_good;
[gg2 idx]=rmoutliers(gg2,"ThresHold",0.9);
MIN=min(length(gg1),length(gg2));
gg1=gg1(1:MIN);
gg2=gg2(1:MIN);
figure(6)
boxplot([gg1' gg2']);
a=gg1;
b=gg2;
pp=ranksum(a,b)

figure(7)
boxplot([coherence_mean_bad' coherence_mean_good']);
a=coherence_mean_bad;
b=coherence_mean_good;
pp=ranksum(a,b)

figure(8)
boxplot([EntropyVN_bad' EntropyVN_good']);
a=EntropyVN_bad;
b=EntropyVN_good;
pp=ranksum(a,b)

figure(9)
boxplot([spectral_gap_bad' spectral_gap_good'])
a=spectral_gap_bad;
b=spectral_gap_good;
pp=ranksum(a,b)

figure(10)
boxplot([Energy_bad' Energy_good'])
a=Energy_bad;
b=Energy_good;
pp=ranksum(a,b)

save results_HopfQL_Ceff_Vaxjo.mat errorFC corrFC ...
    interferencevaxjo interferencevaxjo_emp interferencevaxjo_good interferencevaxjo_bad Crest ...
    uncertainty uncertainty_good uncertainty_bad ...
    TurbulenceEmpirical Turbulence ...
    coherence_mean EntropyVN EntropyVN_good EntropyVN_bad coherence_mean_bad coherence_mean_good ...
    spectral_gap spectral_gap_good spectral_gap_bad ...
    Energy Energy_bad Energy_good;

save results_GBCs_Vaxjo_Uncertainty.mat GBC_Uncertainty GBC_Vaxjo Matrix_Uncertainty Matrix_Vaxjo;
