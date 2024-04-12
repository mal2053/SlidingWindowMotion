%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up directories
ddir = '/Users/martinlindquist/Desktop/SlidingWindowMotion/imagefiles';
subjdir = filenames(fullfile(ddir,'KKI2009*'));
motiondir = '/Users/martinlindquist/Desktop/SlidingWindowMotion/motionfiles';
addpath('/Users/martinlindquist/Desktop/SlidingWindowMotion/DC_toolbox/')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% NOTE: You will need the CanlabCore tools to use this script: https://github.com/canlab


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% List of subjects
sub = ['01'; '02'; '03'; '04'; '05'; '06'; '07'; '09'; '10'; '12'; '13'; '14'; '15'; '16'; '18'; '23'; '28'; '30'; '32'; '39'];
subjdir(8) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Settings
TR = 2;                  % TR in sec
num_regions = 268;       % Number of regions
T = 210;                 % Number of time points
win_length = 30;         % Window length
wT = T - win_length + 1; % Number of time points - window length
time=1:210;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup filter matrix
TR = 2;     W = zeros(T,T+win_length-1);
for i=1:T, W(i,i:i+win_length-1) = 1; end
W = W(1:T,win_length:T+win_length-1);
W = win_length.*W./(sum(W,2)*ones(1,T));


% Setup noise smoother
win_length2 = 10;
Wn = zeros(T,T+win_length2-1);
for i=1:T, Wn(i,i:i+win_length2-1) = 1; end
Wn = Wn(1:T,win_length2:T+win_length2-1);
Wn = win_length2.*Wn./(sum(Wn,2)*ones(1,T));


Fitall = zeros(20,6);
NNall = zeros(wT,20);
ResNNall = zeros(20,7);
ResNNcorr = zeros(20,6);
ResNNsubjects = [];
ResfXYall = zeros(20,7);
ResfXYcorr = zeros(20,6);
ResfXYsubjects = [];
ResNSall = zeros(20,7);
ResNScorr = zeros(20,6);
ResNSsubjects = [];


for subject =1:20

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Load motion parametrs and create projection matrix for motion regression 
    X = load(fullfile(motiondir,strcat('rp_aKKI2009-',sub(subject,:),'-fMRI_ext24.txt')));   
    PX = X*inv(X'*X)*X';
    
    X2 = X(:,1:7);
    PX2 = X2*inv(X2'*X2)*X2';               % Project on space spanned by nuisance
    iPX2 = eye(size(PX2)) - PX2;            % Project on space orthogonal to the one spanned by nuisance
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute nuisance norms (filtered/non-filtered)
    NN = zeros(wT,1);
    for i=1:wT
        tmp = X(i:i+(win_length-1),2:7);
        tmp2 = tmp - repmat(mean(tmp),win_length,1);
        NN(i) = sqrt(sum(sum(tmp2.^2)));
    end
    
    
    NN0 = zeros(T,1);
    for i=1:T
        NN0(i) = sqrt(sum(X(i,2:7).^2));
    end
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Create windowing matrix orthogonal to space spanned by the motion
    W2 = (eye(T) - NN0*inv(NN0'*NN0)*NN0')*W;


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Read data and parcellate
    obj = fmri_data(fullfile(subjdir{subject},strcat('wraKKI2009-',sub(subject,:),'-fMRI.nii')));
    atlas_obj = load_atlas('shen');
    regions = extract_roi_averages(obj,atlas_obj);
    
    
    % True dynamic connectivity 
    rho = .5*sin(2*pi*time/150);
        
    reps = 1000;            % Number of repetitions
    

    MotionCorr = zeros(reps,6,2);
    ResNN = zeros(reps,7,6);
    ResX = zeros(reps,7,6);
    ResY = zeros(reps,7,6);
    ResfXY = zeros(reps,7,6);
    ResNS = zeros(reps,7,6);
    Fit = zeros(reps, 6);
    Nuis = zeros(reps, 3);
      
    for r=1:reps 
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Generate data
        
        % Extract two random motion time courses
        [a, b] = sort(rand(268,1));
        eps_x = PX*regions(b(1)).dat;
        eps_x = zscore(Wn*eps_x);
        eps_y = PX*regions(b(2)).dat;
        eps_y = zscore(Wn*eps_y);

        % Generate two simulated time courses with motion
        x0 = zeros(T,1);
        y0 = zeros(T,1); 

        for k=2:T
            tmp = mvnrnd([0; 0], 0.25*[[1 rho(k)]; [rho(k) 1]]);

            x0(k) = tmp(1);
            y0(k) = tmp(2);

        end

        x = x0 + eps_x;
        y = y0 + eps_y;     
        xy = x.*y;
    
        % Calculate motion corrected signal and motion signal and filtered
        % motion signal
        mcx = iPX2*x;
        motion_x = PX2*x;
        mcy = iPX2*y;
        motion_y = PX2*y;

        mcxy = iPX2*xy;
        motion_xy = motion_x.*motion_y;
    
        Nuis(r,1) = corr(motion_x, motion_y);
        Nuis(r,2) = corr(motion_x, mcy);
        Nuis(r,3) = corr(mcx, motion_y);
    
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Windowed motion
        fmotion_x = zeros(wT,1);
        fmotion_y = zeros(wT,1);
        fmotion_xy = zeros(wT,1);
        frho = zeros(wT,1);
        
        for k=1:wT
            fmotion_x(k) = sum(motion_x(k:k+(win_length-1)));
            fmotion_y(k) = sum(motion_y(k:k+(win_length-1)));
            fmotion_xy(k) = sum(motion_xy(k:k+(win_length-1)));
            frho(k) = sum(rho(k:k+(win_length-1)));
        end
    
    
        MotionCorr(r,:,1) = Matrix2Vector(corr([NN0 motion_x motion_y motion_xy]), tril(ones(4),-1))';
        MotionCorr(r,:,2) = Matrix2Vector(corr([NN fmotion_x fmotion_y fmotion_xy]), tril(ones(4),-1))'; 
    
        C = zeros(wT,6);
        Cov = zeros(wT,6);
    
        D1 = zeros(wT,6);
        D2 = zeros(wT,6);
        D3 = zeros(wT,6);
        D4 = zeros(wT,6);
        D5 = zeros(wT,6);
        D6 = zeros(wT,6);
    
        for t=1:wT
           
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Set-up for block regression
            % 6 df case
            XC = X(t:t+(win_length-1),2:7);     
            XC = XC - repmat(mean(XC), win_length, 1);
            H = XC*pinv(XC);                    % Hat matrix     
            xInt = x(t:t+(win_length-1));
            yInt = y(t:t+(win_length-1));
            mxInt = H*xInt;
            myInt = H*yInt;
            mxyInt = H*(xInt.*yInt);
 

            % 12 df case
            XC2 = X(t:t+(win_length-1),2:13);     
            XC2 = XC2 - repmat(mean(XC2), win_length, 1);
            H2 = XC2*pinv(XC2);                    % Hat matrix
            xInt = x(t:t+(win_length-1));
            yInt = y(t:t+(win_length-1));
            mxInt2 = H2*xInt;
            myInt2 = H2*yInt;
    
    
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % I. Pipeline 1: Original signal
                        
            xInt = xInt - mean(xInt);
            yInt = yInt - mean(yInt);
            xsq = sum(xInt.^2);
            ysq = sum(yInt.^2);
            D1(t,1) = sum(xInt);
            D1(t,2) = sum(yInt);
            D1(t,3) = sum(xInt.*yInt);
            D1(t,4) = xsq;
            D1(t,5) = ysq;
            Cov(t,1) = (sum(xInt.*yInt) - sum(xInt)*sum(yInt)); 
            C(t,1) = (sum(xInt.*yInt) - sum(xInt)*sum(yInt))/sqrt(xsq*ysq);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
            % II. Pipeline 2: Motion-corrected signal
    
            xmcInt = mcx(t:t+(win_length-1));
            xmcInt = xmcInt - mean(xmcInt);
            ymcInt = mcy(t:t+(win_length-1));
            ymcInt = ymcInt - mean(ymcInt);
            xsq = sum(xmcInt.^2);
            ysq = sum(ymcInt.^2);    
            D2(t,1) = sum(xmcInt);
            D2(t,2) = sum(ymcInt);
            D2(t,3) = sum(xmcInt.*ymcInt);
            D2(t,4) = xsq;
            D2(t,5) = ysq;
            Cov(t,2) = (sum(xmcInt.*ymcInt) - sum(xmcInt)*sum(ymcInt)); 
            C(t,2) = (sum(xmcInt.*ymcInt) - sum(xmcInt)*sum(ymcInt))/sqrt(xsq*ysq);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
            % III. Pipeline 3: Motion signal
    
            xmInt = motion_x(t:t+(win_length-1));
            xmInt = xmInt - mean(xmInt);
            ymInt = motion_y(t:t+(win_length-1));
            ymInt = ymInt - mean(ymInt);
            xsq = sum(xmInt.^2);
            ysq = sum(ymInt.^2);    
            D3(t,1) = sum(xmInt);
            D3(t,2) = sum(ymInt);
            D3(t,3) = sum(xmInt.*ymInt);
            D3(t,4) = xsq;
            D3(t,5) = ysq;
            Cov(t,3) = (sum(xmInt.*ymInt) - sum(xmInt)*sum(ymInt)); 
            C(t,3) = corr(xmInt, ymInt);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % IV. Pipeline 4: Block regression 
       
            mcxInt = xInt - mxInt;
            mcyInt = yInt - myInt;      
            D4(t,1) = sum(mcxInt);
            D4(t,2) = sum(mcyInt);
            xsq = sum(mcxInt.^2);
            ysq = sum(mcyInt.^2);   
            D4(t,4) = xsq;
            D4(t,5) = ysq;
            mcxyInt = mcxInt.*mcyInt;
            mcxyInt = mcxyInt - mean(mcxyInt);
            D4(t,3) = sum(mcxyInt);
            Cov(t,4)= sum(mcxyInt) - sum(mcxInt)*sum(mcyInt);
            C(t,4)= corr(mcxInt,mcyInt);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % VI. Pipeline 5: Block Partial Correlation
     
            mcxInt2 = xInt - mxInt2;
            mcyInt2 = yInt - myInt2;              
            D6(t,1) = sum(mcxInt2);
            D6(t,2) = sum(mcyInt2);
            xsq = sum(mcxInt2.^2);
            ysq = sum(mcyInt2.^2);   
            D6(t,4) = xsq;
            D6(t,5) = ysq;
            mcxyInt2 = mcxInt2.*mcyInt2;
            mcxyInt2 = mcxyInt2 - mean(mcxyInt2);
            D6(t,3) = sum(mcxyInt2);
            Cov(t,6)= sum(mcxyInt2) - sum(mcxInt2)*sum(mcyInt2);
%            C(t,6)= corr(mcxInt2,mcyInt2);
            C(t,6) = partialcorr(mcxInt,mcyInt, X(t:t+(win_length-1),2:13) );
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % V. Pipeline 6: Motion Corrected Filter 
    
        f2mcx = W2*mcx;
        f2mcx =  f2mcx(win_length:end);
        D5(:,1) =  f2mcx;
        f2mcy = W2*mcy;
        f2mcy =  f2mcy(win_length:end);
        D5(:,2) =  f2mcy;
        f2mcxy = W2*mcxy;
        f2mcxy = f2mcxy(win_length:end);
        D5(:,3) = f2mcxy;
        f2mcxsq = W2*(mcx.*mcx);
        f2mcxsq =  f2mcxsq(win_length:end);
        D5(:,4) = f2mcxsq - f2mcx.*f2mcx;
        f2mcysq = W2*(mcx.*mcx);
        f2mcysq =  f2mcysq(win_length:end);
        D5(:,5) = f2mcysq - f2mcy.*f2mcy;
        f2mcxy = W2*(mcx.*mcxy);
        f2mcxy =  f2mcxy(win_length:end);
        Cov(:,5)= f2mcxy - f2mcx.*f2mcy;
        f2corr = zeros(T,1);
        for t=1:T
            D = [mcx mcy].* repmat(W2(t,:)',1,2);    
            tmp = corr(D);
            f2corr(t) = tmp(2,1); 
        end
        C(:,5) = f2corr(win_length:end); 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Check correlation between estimated and true TVFC
        Fit(r,1) = corr(C(:,1),frho);
        Fit(r,2) = corr(C(:,2),frho);
        Fit(r,3) = corr(C(:,3),frho);
        Fit(r,4) = corr(C(:,4),frho);
        Fit(r,5) = corr(C(:,5),frho);
        Fit(r,6) = corr(C(:,6),frho);
    
    
        % Check correlation between components of TVFC and nuisance signal

        ResNN(r,:,1) = [corr(C(:,1),NN) corr(Cov(:,1),NN) corr(D1(:,1),NN) corr(D1(:,2),NN) corr(D1(:,3),NN) corr(D1(:,4),NN) corr(D1(:,5),NN)];
        ResNN(r,:,2) = [corr(C(:,2),NN) corr(Cov(:,2),NN) corr(D2(:,1),NN) corr(D2(:,2),NN) corr(D2(:,3),NN) corr(D2(:,4),NN) corr(D2(:,5),NN)];
        ResNN(r,:,3) = [corr(C(:,3),NN) corr(Cov(:,3),NN) corr(D3(:,1),NN) corr(D3(:,2),NN) corr(D3(:,3),NN) corr(D3(:,4),NN) corr(D3(:,5),NN)];
        ResNN(r,:,4) = [corr(C(:,4),NN) corr(Cov(:,4),NN) corr(D4(:,1),NN) corr(D4(:,2),NN) corr(D4(:,3),NN) corr(D4(:,4),NN) corr(D4(:,5),NN)];
        ResNN(r,:,5) = [corr(C(:,5),NN) corr(Cov(:,5),NN) corr(D5(:,1),NN) corr(D5(:,2),NN) corr(D5(:,3),NN) corr(D5(:,4),NN) corr(D5(:,5),NN)];
        ResNN(r,:,6) = [corr(C(:,6),NN) corr(Cov(:,6),NN) corr(D6(:,1),NN) corr(D6(:,2),NN) corr(D6(:,3),NN) corr(D6(:,4),NN) corr(D6(:,5),NN)];
    
        ResX(r,:,1) = [corr(C(:,1),fmotion_x) corr(Cov(:,1),fmotion_x) corr(D1(:,1),fmotion_x) corr(D1(:,2),fmotion_x) corr(D1(:,3),fmotion_x) corr(D1(:,4),fmotion_x) corr(D1(:,5),fmotion_x)];
        ResX(r,:,2) = [corr(C(:,2),fmotion_x) corr(Cov(:,2),fmotion_x) corr(D2(:,1),fmotion_x) corr(D2(:,2),fmotion_x) corr(D2(:,3),fmotion_x) corr(D2(:,4),fmotion_x) corr(D2(:,5),fmotion_x)];
        ResX(r,:,3) = [corr(C(:,3),fmotion_x) corr(Cov(:,3),fmotion_x) corr(D3(:,1),fmotion_x) corr(D3(:,2),fmotion_x) corr(D3(:,3),fmotion_x) corr(D3(:,4),fmotion_x) corr(D3(:,5),fmotion_x)];
        ResX(r,:,4) = [corr(C(:,4),fmotion_x) corr(Cov(:,4),fmotion_x) corr(D4(:,1),fmotion_x) corr(D4(:,2),fmotion_x) corr(D4(:,3),fmotion_x) corr(D4(:,4),fmotion_x) corr(D4(:,5),fmotion_x)];
        ResX(r,:,5) = [corr(C(:,5),fmotion_x) corr(Cov(:,5),fmotion_x) corr(D5(:,1),fmotion_x) corr(D5(:,2),fmotion_x) corr(D5(:,3),fmotion_x) corr(D5(:,4),fmotion_x) corr(D5(:,5),fmotion_x)];
        ResX(r,:,6) = [corr(C(:,6),fmotion_x) corr(Cov(:,6),fmotion_x) corr(D6(:,1),fmotion_x) corr(D6(:,2),fmotion_x) corr(D6(:,3),fmotion_x) corr(D6(:,4),fmotion_x) corr(D6(:,5),fmotion_x)];
    
        ResY(r,:,1) = [corr(C(:,1),fmotion_y) corr(Cov(:,1),fmotion_y) corr(D1(:,1),fmotion_y) corr(D1(:,2),fmotion_y) corr(D1(:,3),fmotion_y) corr(D1(:,4),fmotion_y) corr(D1(:,5),fmotion_y)];
        ResY(r,:,2) = [corr(C(:,2),fmotion_y) corr(Cov(:,2),fmotion_y) corr(D2(:,1),fmotion_y) corr(D2(:,2),fmotion_y) corr(D2(:,3),fmotion_y) corr(D2(:,4),fmotion_y) corr(D2(:,5),fmotion_y)];
        ResY(r,:,3) = [corr(C(:,3),fmotion_y) corr(Cov(:,3),fmotion_y) corr(D3(:,1),fmotion_y) corr(D3(:,2),fmotion_y) corr(D3(:,3),fmotion_y) corr(D3(:,4),fmotion_y) corr(D3(:,5),fmotion_y)];
        ResY(r,:,4) = [corr(C(:,4),fmotion_y) corr(Cov(:,4),fmotion_y) corr(D4(:,1),fmotion_y) corr(D4(:,2),fmotion_y) corr(D4(:,3),fmotion_y) corr(D4(:,4),fmotion_y) corr(D4(:,5),fmotion_y)];
        ResY(r,:,5) = [corr(C(:,5),fmotion_y) corr(Cov(:,5),fmotion_y) corr(D5(:,1),fmotion_y) corr(D5(:,2),fmotion_y) corr(D5(:,3),fmotion_y) corr(D5(:,4),fmotion_y) corr(D5(:,5),fmotion_y)];
        ResY(r,:,6) = [corr(C(:,6),fmotion_y) corr(Cov(:,6),fmotion_y) corr(D6(:,1),fmotion_y) corr(D6(:,2),fmotion_y) corr(D6(:,3),fmotion_y) corr(D6(:,4),fmotion_y) corr(D6(:,5),fmotion_y)];
    
        ResfXY(r,:,1) = [corr(C(:,1),fmotion_xy) corr(Cov(:,1),fmotion_xy) corr(D1(:,1),fmotion_xy) corr(D1(:,2),fmotion_xy) corr(D1(:,3),fmotion_xy) corr(D1(:,4),fmotion_xy) corr(D1(:,5),fmotion_xy)];
        ResfXY(r,:,2) = [corr(C(:,2),fmotion_xy) corr(Cov(:,2),fmotion_xy) corr(D2(:,1),fmotion_xy) corr(D2(:,2),fmotion_xy) corr(D2(:,3),fmotion_xy) corr(D2(:,4),fmotion_xy) corr(D2(:,5),fmotion_xy)];
        ResfXY(r,:,3) = [corr(C(:,3),fmotion_xy) corr(Cov(:,3),fmotion_xy) corr(D3(:,1),fmotion_xy) corr(D3(:,2),fmotion_xy) corr(D3(:,3),fmotion_xy) corr(D3(:,4),fmotion_xy) corr(D3(:,5),fmotion_xy)];
        ResfXY(r,:,4) = [corr(C(:,4),fmotion_xy) corr(Cov(:,4),fmotion_xy) corr(D4(:,1),fmotion_xy) corr(D4(:,2),fmotion_xy) corr(D4(:,3),fmotion_xy) corr(D4(:,4),fmotion_xy) corr(D4(:,5),fmotion_xy)];
        ResfXY(r,:,5) = [corr(C(:,5),fmotion_xy) corr(Cov(:,5),fmotion_xy) corr(D5(:,1),fmotion_xy) corr(D5(:,2),fmotion_xy) corr(D5(:,3),fmotion_xy) corr(D5(:,4),fmotion_xy) corr(D5(:,5),fmotion_xy)];
        ResfXY(r,:,6) = [corr(C(:,6),fmotion_xy) corr(Cov(:,6),fmotion_xy) corr(D6(:,1),fmotion_xy) corr(D6(:,2),fmotion_xy) corr(D6(:,3),fmotion_xy) corr(D6(:,4),fmotion_xy) corr(D6(:,5),fmotion_xy)];
    
        ResNS(r,:,1) = [corr(C(:,1),C(:,3)) corr(Cov(:,1),C(:,3)) corr(D1(:,1),C(:,3)) corr(D1(:,2),C(:,3)) corr(D1(:,3),C(:,3)) corr(D1(:,4),C(:,3)) corr(D1(:,5),C(:,3))];
        ResNS(r,:,2) = [corr(C(:,2),C(:,3)) corr(Cov(:,2),C(:,3)) corr(D2(:,1),C(:,3)) corr(D2(:,2),C(:,3)) corr(D2(:,3),C(:,3)) corr(D2(:,4),C(:,3)) corr(D2(:,5),C(:,3))];
        ResNS(r,:,3) = [corr(C(:,3),C(:,3)) corr(Cov(:,3),C(:,3)) corr(D3(:,1),C(:,3)) corr(D3(:,2),C(:,3)) corr(D3(:,3),C(:,3)) corr(D3(:,4),C(:,3)) corr(D3(:,5),C(:,3))];
        ResNS(r,:,4) = [corr(C(:,4),C(:,3)) corr(Cov(:,4),C(:,3)) corr(D4(:,1),C(:,3)) corr(D4(:,2),C(:,3)) corr(D4(:,3),C(:,3)) corr(D4(:,4),C(:,3)) corr(D4(:,5),C(:,3))];
        ResNS(r,:,5) = [corr(C(:,5),C(:,3)) corr(Cov(:,5),C(:,3)) corr(D5(:,1),C(:,3)) corr(D5(:,2),C(:,3)) corr(D5(:,3),C(:,3)) corr(D5(:,4),C(:,3)) corr(D5(:,5),C(:,3))];
        ResNS(r,:,6) = [corr(C(:,6),C(:,3)) corr(Cov(:,6),C(:,3)) corr(D6(:,1),C(:,3)) corr(D6(:,2),C(:,3)) corr(D6(:,3),C(:,3)) corr(D6(:,4),C(:,3)) corr(D6(:,5),C(:,3))];
 
        
    end

    ResNNsubjects{subject} = ResNN;
    ResNNall(subject,:) = mean(squeeze(ResNN(:,:,1)));
    ResNNcorr(subject,:) = mean(squeeze(ResNN(:,1,:)));
    Fitall(subject,:) = mean(squeeze(Fit));
    NNall(:,subject) = NN;

   
    ResfXYsubjects{subject} = ResfXY;
    ResfXYall(subject,:) = mean(squeeze(ResfXY(:,:,1)));
    ResfXYcorr(subject,:) = mean(squeeze(ResfXY(:,1,:)));

    ResNSsubjects{subject} = ResNS;
    ResNSall(subject,:) = mean(squeeze(ResNS(:,:,1)));
    ResNScorr(subject,:) = mean(squeeze(ResNS(:,1,:)));
    
    
    disp(subject)

end

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure: Simulation set-up and results

figure
subplot 131
plot(rho,'LineWidth',2)
title('True time-varrying connectivity \rho(t)')
axis([0 211 -.6 .6])
xlabel('Time (TR)');
ylabel('Correlation');
set(gca,'FontSize',16);

subplot 132
plot(NNall,'LineWidth',1.5)
title('Nuisance norms')
xlabel('Time (TR)');
axis([0 182 0 3.2])
set(gca,'FontSize',16);

subplot 133
%Fitall(:,3) = [];
violin(Fitall)
xticks([1:6])
xticklabels({'P1','P2', 'P3', 'P4', 'P5','P6'})
title('Correlation with \rho(t)')
ylabel('Correlation');
xlabel('Pipeline');
axis([0 7 -1 1])
set(gca,'FontSize',16);


%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure: Simulation results

figure
subplot 231
ResNN = ResNNsubjects{subject};
dat = [ResNN(:,1,1) ResNN(:,2,1) ResNN(:,3,1) ResNN(:,4,1) ResNN(:,5,1) ResNN(:,6,1) ResNN(:,7,1)];
grp = [ones(1000,1); 2*ones(1000,1); 3*ones(1000,1); 4*ones(1000,1); 5*ones(1000,1); 6*ones(1000,1); 7*ones(1000,1)];
violin(dat)
xticks([1:7])
xticklabels({'Corr','Cov', 'x*w', 'y*w','(xy)*w','(x^2)*w','(y^2)*w'})
axis([0 8 -1.04 1.04])
title('Correlation with nuisance norm (Subject 2)')
ylabel('Correlation');
set(gca,'FontSize',16);

subplot 232
violin(ResNNall)
xticks([1:7])
xticklabels({'Corr','Cov', 'x*w', 'y*w','(xy)*w','(x^2)*w','(y^2)*w'})
axis([0 8 -1.04 1.04])
title('Correlation with nuisance norm (All subjects)')
ylabel('Correlation');
set(gca,'FontSize',16);

subplot 233
violin(ResNNcorr)
xticks([1:6])
xticklabels({'P1','P2', 'P3', 'P4', 'P5','P6'})
title('Correlation with nuisance norm')
ylabel('Correlation');
xlabel('Pipeline');
axis([0 7 -1.04 1.04])
set(gca,'FontSize',16);

subplot 234
ResNS = ResNSsubjects{subject};
dat = [ResNS(:,1,1) ResNS(:,2,1) ResNS(:,3,1) ResNS(:,4,1) ResNS(:,5,1) ResNS(:,6,1) ResNS(:,7,1)];
grp = [ones(1000,1); 2*ones(1000,1); 3*ones(1000,1); 4*ones(1000,1); 5*ones(1000,1); 6*ones(1000,1); 7*ones(1000,1)];
violin(dat)
xticks([1:7])
xticklabels({'Corr','Cov', 'x*w', 'y*w','(xy)*w','(x^2)*w','(y^2)*w'})
axis([0 8 -1.04 1.04])
title('Correlation with nuisance correlation (Subject 2)')
ylabel('Correlation');
set(gca,'FontSize',16);

subplot 235
violin(ResNSall)
xticks([1:7])
xticklabels({'Corr','Cov', 'x*w', 'y*w','(xy)*w','(x^2)*w','(y^2)*w'})
axis([0 8 -1.04 1.04])
title('Correlation with nuisance correlation (All subjects)')
ylabel('Correlation');
set(gca,'FontSize',16);

subplot 236
ResNScorr(:,3) = ResNScorr(:,3) + normrnd(0,0.00000001,20,1);
xticks([1:6])
violin(ResNScorr)
xticks([1:6])
xticklabels({'P1','P2', 'P3', 'P4', 'P5','P6'})
title('Correlation with nuisance correlation')
ylabel('Correlation');
xlabel('Pipeline');
axis([0 7 -1.04 1.04])
set(gca,'FontSize',16);


