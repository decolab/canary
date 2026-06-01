clear all;
path2=[ '../../Turbulence/Basics'];
addpath(genpath(path2));
path3=[ '../../Turbulence/SC_longrange'];
addpath(genpath(path3));
path4=[ '../../Nonequilibrium'];
addpath(genpath(path4));
path5=[ '../../BANDAdata'];
addpath(genpath(path5));
path6=[ '../../BANDAdata/sFC_Schaefer2018_1000parcels_7Networks_order_REST_AP'];
addpath(genpath(path6));

load ('schaefer1000to100.mat');
partitiondef=schaefer1000to100;
NP=100;
N=1000;

Isubdiag = find(tril(ones(N),-1));

%% Loading the SC, anatomy of HCP in schaefer 1000

for i=1:NP
    partition{i}=find(partitiondef==i);
end

%%
% Parameters of the data
TR=0.8;  % Repetition Time (seconds)
% Bandpass filter settings
fnq=1/(2*TR);                 % Nyquist frequency
flp = 0.008;                    % lowpass frequency of filter (Hz)
fhi = 0.08;                    % highpass
Wn=[flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
k=2;                          % 2nd order butterworth filter
[bfilt,afilt]=butter(k,Wn);   % construct the filter

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


%% Loading IDs and timeseries
load ("Original_groups_ID.mat")
load ("subjects_BANDA_char.mat")
load ("sFC_Schaefer2018_1000Parcels_7Networks_order_REST1_AP.mat")
%% Extracting parcitipants ids from fmri file order
% subjects is 202x11 char -> convert to clean string array (202x1)
subjNames_full = strtrim(string(cellstr(subjects)));
subjNames_base = extractBefore(subjNames_full, "_");
%% Group ID
groupIDs = struct();
groupIDs.CNT = strtrim(string(CNT_ID(:)));
groupIDs.ANX = strtrim(string(ANX_ID(:)));
groupIDs.DEP = strtrim(string(DEP_ID(:)));
%% Get indices in subject{} for each group
groupIdx = struct();
fn = fieldnames(groupIDs);
for i = 1:numel(fn)
    g = fn{i};
    ids = groupIDs.(g);
    [tf, loc] = ismember(ids, subjNames_base);   % loc are indices into subject{}
    groupIdx.(g) = loc(tf);                % valid indices
end
%% Example: loop and read timeseries per group
for i = 1:numel(fn)
    g = fn{i};
    clear vaxjo2 Interferencevaxjo_submat;
    vaxjo2NULLall=[];
    idxs = groupIdx.(g);
    for nsub = 1:numel(idxs)
        nsub
        s = subject{idxs(nsub)};          % 1x1 struct
        ts = s.schaeferts;            % 100 x T (from your screenshot)
        clear signal_filt Phases;
        for seed=1:N
            ts(seed,:)=detrend(ts(seed,:)-nanmean(ts(seed,:)));
            if sum(isnan(ts(seed,:)))==0
                signal_filt(seed,:)=filtfilt(bfilt,afilt,ts(seed,:));
                Xanalytic = hilbert(demean(signal_filt(seed,:)));
                Phases(seed,:) = angle(Xanalytic);
            end
        end

        ts=signal_filt(:,20:end-20);
        Phases=Phases(:,20:end-20);

        for i=1:NP
            boldpart(i,:)=nanmean(ts(partition{i},:));
            Vorticity(i,:)=abs(sum(complex(cos(Phases(partition{i},:)),sin(Phases(partition{i},:))))/length(partition{i}));
        end

        FC=corrcoef(boldpart');

        FCemp.(g)=FC;

        Tmax2=size(Vorticity,2);
        NullVorticity=Vorticity(randperm(NP),randperm(Tmax2));

        %%%%
        %% Vaxjol
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

                %% surro
                if NullVorticity(i,j)<0.2
                    binvortexNULL(i,j)=0;
                end
                if NullVorticity(i,j)>=0.2 && NullVorticity(i,j)<0.4
                    binvortexNULL(i,j)=1;
                end
                if NullVorticity(i,j)>=0.4 && NullVorticity(i,j)<0.6
                    binvortexNULL(i,j)=2;
                end
                if NullVorticity(i,j)>=0.6 && NullVorticity(i,j)<0.8
                    binvortexNULL(i,j)=3;
                end
                if NullVorticity(i,j)>=0.8
                    binvortexNULL(i,j)=4;
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

                    %% Vaxjo

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

                    vaxjo2(nsub,i,j)=deltapair(npa);

                    %% Surrogates
                    bin1=binvortexNULL(i,:);
                    bin2=binvortexNULL(j,:);

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

                    %% Vaxjo

                    vaxjo2NULL(npa)=abs(P1_0-P2_0*P12_00-P2_1*P12_01-P2_2*P12_02-P2_3*P12_03-P2_4*P12_04) ...
                        /2/sqrt(P2_0*P12_00*P2_1*P12_01*P2_2*P12_02*P2_3*P12_03*P2_4*P12_04) ...
                        +abs(P1_1-P2_0*P12_10-P2_1*P12_11-P2_2*P12_12-P2_3*P12_13-P2_4*P12_14) ...
                        /2/sqrt(P2_0*P12_10*P2_1*P12_11*P2_2*P12_12*P2_3*P12_13*P2_4*P12_14) ...
                        +abs(P1_2-P2_0*P12_20-P2_1*P12_21-P2_2*P12_22-P2_3*P12_23-P2_4*P12_24) ...
                        /2/sqrt(P2_0*P12_20*P2_1*P12_21*P2_2*P12_22*P2_3*P12_23*P2_4*P12_24) ...
                        +abs(P1_3-P2_0*P12_30-P2_1*P12_31-P2_2*P12_32-P2_3*P12_33-P2_4*P12_34) ...
                        /2/sqrt(P2_0*P12_30*P2_1*P12_31*P2_2*P12_32*P2_3*P12_33*P2_4*P12_34) ...
                        +abs(P1_4-P2_0*P12_40-P2_1*P12_41-P2_2*P12_42-P2_3*P12_43-P2_4*P12_44) ...
                        /2/sqrt(P2_0*P12_40*P2_1*P12_41*P2_2*P12_42*P2_3*P12_43*P2_4*P12_44);

                    %%

                    npa=npa+1;

                end
            end
        end

        vaxjo2NULL(isinf(vaxjo2NULL))=[];
        vaxjo2NULL(isnan(vaxjo2NULL))=[];
        vaxjo2NULLall=[vaxjo2NULLall vaxjo2NULL];
        FCemp_submat(nsub,:,:)=FC;
    end

    mNULL=nanmean(vaxjo2NULLall);
    sNULL=nanstd(vaxjo2NULLall);

    for nsub=1:numel(idxs)
        for i=1:NP
            for j=1:NP
                if i~=j
                    Interferencevaxjo_submat(nsub,i,j)=(vaxjo2(nsub,i,j)-mNULL)/sNULL;
                    if isinf(Interferencevaxjo_submat(nsub,i,j))
                        Interferencevaxjo_submat(nsub,i,j)=randn*0.001;
                    end
                    if isnan(Interferencevaxjo_submat(nsub,i,j))
                        Interferencevaxjo_submat(nsub,i,j)=randn*0.001;
                    end
                end
            end
        end
    end

    Interferencevaxjo_sub.(g)=Interferencevaxjo_submat;
    FCemp_sub.(g)=FCemp_submat;
end

save results_BANDA_Vaxjo_emp.mat FCemp FCemp_sub Interferencevaxjo_sub;