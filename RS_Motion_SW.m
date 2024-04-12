%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up directories
ddir = '/Users/martinlindquist/Desktop/SlidingWindowMotion/imagefiles';
subjdir = filenames(fullfile(ddir,'KKI2009*'));
motiondir = '/Users/martinlindquist/Desktop/SlidingWindowMotion/motionfiles';
addpath('DC_toolbox/')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% NOTE: You will need the CanlabCore tools to use this script: https://github.com/canlab

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Chose which subjects to run

Session = 1;                % Set  Session to either 1 or 2.

if Session == 1
    sub = ['01'; '02'; '03'; '04'; '05'; '06'; '07'; '09'; '10'; '12'; '13'; '14'; '15'; '16'; '18'; '23'; '28'; '30'; '32'; '39'];
elseif Session == 2
    sub = ['25'; '37'; '22'; '11'; '31'; '20'; '34'; '42'; '21'; '19'; '24'; '17'; '26'; '35'; '38'; '27'; '40'; '33'; '36'; '41'];
else
    error('Wrong session number')
end
subjdir(8) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Settings
TR = 2;                  % TR in sec
num_regions = 268;       % Number of regions
num_connections = 35778; % Number of connections
T = 210;                 % Number of time points
win_length = 30;         % Window length
wT = T - win_length + 1; % Number of time points - window length
num_sub = length(sub);   % Number of subjects
radius = 50;             % Radius for FWD calculation
type = 1;                % Investigate correlation with nuisance norm (1) or motion correlation (2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% List of Pipelines
%
% Pipeline 1: No motion regression
% Pipeline 2: Motion-corrected 
% Pipeline 3: Motion
% Pipeline 4: Windowed motion-corrected (6 dof)
% Pipeline 5: Motion-corrected window
% Pipeline 6: Windowed motion-corrected (12 dof)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Additional set-up for saving results
Q = tril(ones(num_regions, num_regions),-1);

cm0 = zeros(num_regions, num_sub);
cm1 = zeros(num_regions, num_sub);


% Correlations between motion time courses and processed timecourses
corrP1 = zeros(num_regions, num_sub,3);
corrP2 = zeros(num_regions, num_sub,3);
corrP3 = zeros(num_regions, num_sub,3);
corrP4 = zeros(num_regions, num_sub,3);
corrP5 = zeros(num_regions, num_sub,3);
corrP6 = zeros(num_regions, num_sub,3);

% Correlations between nuisance norm and dynamic covariance matrix
corrCovxyP1 = zeros(num_regions, num_regions, num_sub);
corrCovxyP2 = zeros(num_regions, num_regions, num_sub);
corrCovxyP3 = zeros(num_regions, num_regions, num_sub);
corrCovxyP4 = zeros(num_regions, num_regions, num_sub);
corrCovxyP5 = zeros(num_regions, num_regions, num_sub);
corrCovxyP6 = zeros(num_regions, num_regions, num_sub);

% Correlations between nuisance norm and windowed signal
corrMUXP1 = zeros(num_regions, num_sub);
corrMUXP2 = zeros(num_regions, num_sub);
corrMUXP3 = zeros(num_regions, num_sub);
corrMUXP4 = zeros(num_regions, num_sub);
corrMUXP5 = zeros(num_regions, num_sub);
corrMUXP6 = zeros(num_regions, num_sub);

% Correlations between nuisance norm and windowed standard deviation
corrSDXP1 = zeros(num_regions, num_sub);
corrSDXP2 = zeros(num_regions, num_sub);
corrSDXP3 = zeros(num_regions, num_sub);
corrSDXP4 = zeros(num_regions, num_sub);
corrSDXP5 = zeros(num_regions, num_sub);
corrSDXP6 = zeros(num_regions, num_sub);

% Correlations between nuisance norm and windowd E(X)E(Y) product
corrEXEYP1 = zeros(num_regions, num_regions, num_sub);
corrEXEYP2 = zeros(num_regions, num_regions, num_sub);
corrEXEYP3 = zeros(num_regions, num_regions, num_sub);
corrEXEYP4 = zeros(num_regions, num_regions, num_sub);
corrEXEYP5 = zeros(num_regions, num_regions, num_sub);
corrEXEYP6 = zeros(num_regions, num_regions, num_sub);

% Correlations between nuisance norm and windowd E(XY) product
corrEXYP1 = zeros(num_regions, num_regions, num_sub);
corrEXYP2 = zeros(num_regions, num_regions, num_sub);
corrEXYP3 = zeros(num_regions, num_regions, num_sub);
corrEXYP4 = zeros(num_regions, num_regions, num_sub);
corrEXYP5 = zeros(num_regions, num_regions, num_sub);
corrEXYP6 = zeros(num_regions, num_regions, num_sub);

% Correlations between nuisance norm and dynamic correlation matrix
corrRxyP1 = zeros(35778, num_sub);
corrRxyP2 = zeros(35778, num_sub);
corrRxyP3 = zeros(35778, num_sub);
corrRxyP4 = zeros(35778, num_sub);
corrRxyP5 = zeros(35778, num_sub);
corrRxyP6 = zeros(35778, num_sub);

% Dynamic correlation matrix
CP1 = zeros(35778, wT, num_sub);
CP2 = zeros(35778, wT, num_sub);
CP3 = zeros(35778, wT, num_sub);
CP4 = zeros(35778, wT, num_sub);
CP5 = zeros(35778, wT, num_sub);
CP6 = zeros(35778, wT, num_sub);

% Compute proportion cbvariation due to motion
Comp1 = zeros(num_regions, num_regions, num_sub);
Comp2 = zeros(num_regions, num_regions, num_sub);

% Framewsie displacement
FDMat = zeros(T, num_sub);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create windowing matrix 
W = zeros(T,T+win_length-1);
for i=1:T, W(i,i:i+win_length-1) = 1; end
W = W(1:T,win_length:T+win_length-1);
W = win_length.*W./(sum(W,2)*ones(1,T));


NNall = zeros(wT,num_sub);
MTCall = zeros(T,num_sub); 
MTCconvall = zeros(T,num_sub); 
 
for s=1:num_sub

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Load motion parametrs and create projection matrix for motion regression 
    X = load(fullfile(motiondir,strcat('rp_aKKI2009-',sub(s,:),'-fMRI_ext24.txt')));%
    X = X(:,1:13);
    PX = X*inv(X'*X)*X';
    iPX = eye(size(PX)) - PX;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Read data and parcellate
    obj = fmri_data(fullfile(subjdir{s},strcat('wraKKI2009-',sub(s,:),'-fMRI.nii')));
    atlas_obj = load_atlas('shen');
    r = extract_roi_averages(obj,atlas_obj);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute frame-wise displacement
    ts = X(:,2:7);
    tmp = ts(:,4:6); 
    tmp = (2*radius*pi/360)*tmp;
    ts(:,4:6) = tmp;
    dts = diff(ts);
    fwd = sum(abs(dts),2);
    
    if (length(fwd) == (T-1))
        FDMat(2:T,s) = fwd;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute nuisance norm
    NN = zeros(wT,1);
    for i=1:wT
        tmp = X(i:i+(win_length-1),2:7);
        tmp2 = tmp - repmat(mean(tmp),win_length,1);
        NN(i) = sqrt(sum(sum(tmp2.^2)));
    end
    NNall(:,s) = NN;
  
    NN0 = zeros(T,1);
    for i=1:T
        NN0(i) = sqrt(sum(X(i,2:7).^2));
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Time course for Pipelines 1-3

    TCnoMC = zeros(T, num_regions);     % Pipeline 1
    TC = zeros(T, num_regions);         % Pipeline 2
    TCM = zeros(T, num_regions);        % Pipeline 3
    
    for i=1:num_regions
        y = r(i).dat;
        yhat = PX*y;
        TCnoMC(:,i) = y;      % Time course without motion correction
        TCM(:,i) = yhat;      % Motion time course                                                                                        e 
        TC(:,i) = y-yhat;     % Motion-corrected time course
    end

    % Compute windowed motion time course
    MTC = mean(TCM- repmat(mean(TCM),210,1),2);
    MTCconv = W*MTC;

    MTCall(:,s) = MTC;
    MTCconvall(:,s) = MTCconv;


    % Correlation between motion time course and motion-corrected time course
    for i=1:num_regions, cm1(i,s) = corr(TCM(:,i), TC(:,i)); end
    
    % Correlation between motion time course and time course without motion correction
    for i=1:num_regions, cm0(i,s) = corr(TCM(:,i), TCnoMC(:,i)); end
  
    TCnoMCconv  = W*TCnoMC;                 % Windowed original time course; Pipeline 1
    TCconv  = W*TC;                         % Windowed motion-corrected time course; Pipeline 2
    TCMconv = W*TCM;                        % Windowed motion time course; Pipeline 3

       
    % Correlation between motion time courses and windowed non-motion-corrected time course
    for i=1:num_regions, corrP1(i,s,1) = corr(TCM(:,i), TCnoMCconv(:,i)); end
    for i=1:num_regions, corrP1(i,s,2) = corr(TCMconv(:,i), TCnoMCconv(:,i)); end
    for i=1:num_regions, corrP1(i,s,3) = corr(NN, TCnoMCconv(win_length:end,i)); end

    % Correlation between motion time courses and windowed motion-corrected time course
    for i=1:num_regions, corrP2(i,s,1) = corr(TCM(:,i), TCconv(:,i)); end
    for i=1:num_regions, corrP2(i,s,2) = corr(TCMconv(:,i), TCconv(:,i)); end
    for i=1:num_regions, corrP2(i,s,3) = corr(NN, TCconv(win_length:end,i)); end

   % Correlation between motion time courses and windowed motion time course
    for i=1:num_regions, corrP3(i,s,1) = corr(TCM(:,i), TCMconv(:,i)); end
    for i=1:num_regions, corrP3(i,s,2) = corr(TCMconv(:,i), TCMconv(:,i)); end
    for i=1:num_regions, corrP3(i,s,3) = corr(NN, TCMconv(win_length:end,i)); end
     


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Block regression (pipeline 4 and 6)
    %
    TCMCconv = zeros(wT, num_regions);
    TCMCconv2 = zeros(wT, num_regions);
    MC = zeros(num_regions, num_regions, wT);
    MC2 = zeros(num_regions, num_regions, wT);
    TCtmp = zeros(num_regions, win_length);
    TCtmp2 = zeros(num_regions, win_length);
    
    for t=1:wT

        XC = X(t:t+(win_length-1),2:7);    
        XC = XC - repmat(mean(XC), win_length, 1);
        H = XC*pinv(XC);                    % Hat matrix

        XC2 = X(t:t+(win_length-1),2:13);
        XC2 = XC2 - repmat(mean(XC2),win_length,1);
        H2 = XC2*pinv(XC2);

        for j=1:num_regions       
            y = TC(t:t+(win_length-1),j);
            y = y - mean(y);

            yhat = H*y;
            TCMCconv(t,j) = sum(y-yhat);        % Pipeline 4
            TCtmp(j,:) = y-yhat;
      
            yhat2 = H2*y;
            TCMCconv2(t,j) = sum(y-yhat2);       % Pipeline 6
            TCtmp2(j,:) = y-yhat2;
        end
        
        MC(:,:,t) = corr(TCtmp');
        MC2(:,:,t) = corr(TCtmp2');
        
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Create windowing matrix orthogonal to space spanned by the motion
    % Pipeline 5
    W2 =(eye(T) - PX)*W;
    W2 = W2/win_length;

    TCMCconvf  = W2*TC;                 % Windowed signal using motion-corrected window; Pipeline 5
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Investigate different terms of the TVFC

    Cov_P1 = zeros(num_regions,num_regions,wT);
    Cov_P2 = zeros(num_regions,num_regions,wT);
    Cov_P3 = zeros(num_regions,num_regions,wT);
    Cov_P4 = zeros(num_regions,num_regions,wT);
    Cov_P6 = zeros(num_regions,num_regions,wT);

    MUX_P1 = zeros(num_regions,wT);             
    MUX_P2 = zeros(num_regions,wT); 
    MUX_P3 = zeros(num_regions,wT); 
    MUX_P4 = zeros(num_regions,wT); 
    MUX_P5 = zeros(num_regions,wT); 
    MUX_P6 = zeros(num_regions,wT); 

    SDX_P1 = zeros(num_regions,wT);  
    SDX_P2 = zeros(num_regions,wT); 
    SDX_P3 = zeros(num_regions,wT); 
    SDX_P4 = zeros(num_regions,wT); 
    SDX_P5 = zeros(num_regions,wT); 
    SDX_P6 = zeros(num_regions,wT); 

    EXY_P1 = zeros(num_regions,num_regions,wT);
    EXY_P2 = zeros(num_regions,num_regions,wT);
    EXY_P3 = zeros(num_regions,num_regions,wT);
    EXY_P4 = zeros(num_regions,num_regions,wT);
    EXY_P5 = zeros(num_regions,num_regions,wT);
    EXY_P6 = zeros(num_regions,num_regions,wT);

    EXEY_P1 = zeros(num_regions,num_regions,wT);
    EXEY_P2 = zeros(num_regions,num_regions,wT);
    EXEY_P3 = zeros(num_regions,num_regions,wT);
    EXEY_P4 = zeros(num_regions,num_regions,wT);
    EXEY_P5 = zeros(num_regions,num_regions,wT);
    EXEY_P6 = zeros(num_regions,num_regions,wT);

    for t=1:wT

        XC = X(t:t+(win_length-1),2:7);
        XC = XC - repmat(mean(XC),win_length,1);
        H = XC*pinv(XC);

        XC2 = X(t:t+(win_length-1),2:13);
        XC2 = XC2 - repmat(mean(XC2),win_length,1);
        H2 = XC2*pinv(XC2);


        for i=1:num_regions
                        
            xP1 = TCnoMC(t:t+(win_length-1),i);
            xP1 = xP1 - mean(xP1);
                      
            xP2 = TC(t:t+(win_length-1),i);
            xP2 = xP2 - mean(xP2);

            tmpx = H*xP1;
            xP4 = xP1 - tmpx;

            xP3 = TCM(t:t+(win_length-1),i);
            xP3 = xP3 - mean(xP3);

            tmpx2 = H2*xP1;
            xP6 = xP1 - tmpx2;


            MUX_P1(i,t) = sum(xP1); 
            MUX_P2(i,t) = sum(xP2);
            MUX_P3(i,t) = sum(xP3);
            MUX_P4(i,t) = sum(xP4);
            MUX_P6(i,t) = sum(xP6);


            SDX_P1(i,t) = std(xP1); 
            SDX_P2(i,t) = std(xP2);
            SDX_P3(i,t) = std(xP3);
            SDX_P4(i,t) = std(xP4);
            SDX_P6(i,t) = std(xP6);

            for j=1:num_regions
                          
                yP1 = TCnoMC(t:t+(win_length-1),j);
                yP1 = yP1 - mean(yP1);

                yP2 = TC(t:t+(win_length-1),j);
                yP2 = yP2 - mean(yP2);
                
                tmpy = H*yP1;
                yP4 = yP1 - tmpy;

                yP3 = TCM(t:t+(win_length-1),j);
                yP3 = yP3 - mean(yP3);

                xy = xP4.*yP4;
                xyc = xy;

              
                tmpy2 = H2*yP1;
                yP6 = yP1 - tmpy2;

                Cov_P1(i,j,t) = sum(xP1.*yP1) - sum(xP1)*sum(yP1); 
                Cov_P2(i,j,t) = sum(xP2.*yP2) - sum(xP2)*sum(yP2);
                Cov_P3(i,j,t) = sum(xP3.*yP3) - sum(xP3)*sum(yP3);
                Cov_P4(i,j,t) = sum(xP4.*yP4) - sum(xP4)*sum(yP4);
                Cov_P6(i,j,t) = sum(xP6.*yP6) - sum(xP6)*sum(yP6);


                EXY_P1(i,j,t) = sum(xP1.*yP1); 
                EXY_P2(i,j,t) = sum(xP2.*yP2);
                EXY_P3(i,j,t) = sum(xP3.*yP3);
                EXY_P4(i,j,t) = sum(xP4.*yP4);
                EXY_P6(i,j,t) = sum(xP6.*yP6);
                
                
                EXEY_P1(i,j,t) = sum(xP1)*sum(yP1); 
                EXEY_P2(i,j,t) = sum(xP2)*sum(yP2);
                EXEY_P3(i,j,t) = sum(xP3)*sum(yP3);
                EXEY_P4(i,j,t) = sum(xP4)*sum(yP4);
                EXEY_P6(i,j,t) = sum(xP6)*sum(yP6);
                
            end
        end
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Pipeline 5
    
    Cov_P5 = zeros(num_regions,num_regions,T);
    for i=1:num_regions
        x = TC(:,i);

        tmp1 = W2*x;
        MUX_P5(i,:) = tmp1(win_length:end);
        tmp2 = W2*(x.*x)/(win_length-1);
        SDX_P5(i,:) = tmp2(win_length:end); 
  
        for j=1:num_regions
            y = TC(:,j);
            Cov_P5(i,j,:) = W2*(x.*y) - (W2*x).*(W2*y);
                           
            tmp1 = W2*(x.*y);
            EXY_P5(i,j,:) = tmp1(win_length:end);
            tmp2 = (W2*x).*(W2*y);
            EXEY_P5(i,j,:) = tmp2(win_length:end); 

        end
    end

    Cov_P5 = Cov_P5(:,:,win_length:end);
 
    % Compute dynamic connectivity with motion corrected window
    MCF = zeros(num_regions,num_regions,T);
    for i=1:T
        D = TC.* repmat(W2(i,:)',1,num_regions);    
        MCF(:,:,i)=corr(D);
    end
    MCF = MCF(:,:,win_length:end); 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute correlations needed to compute proportion motion
    for i=1:num_regions
        for j=1:num_regions    
            Comp1(i,j,s) = corr(TCM(:,i), TCconv(:,j));
            Comp2(i,j,s) = corr(TC(:,i), TCconv(:,j));
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Correlation between windowed time course with motion correction within each
    % window and motion time courses (6 regressors)
    for i=1:num_regions, corrP4(i,s,1) = corr(TCMCconv(:,i), TCM(win_length:end,i)); end
    for i=1:num_regions, corrP4(i,s,2) = corr(TCMCconv(:,i), TCMconv(win_length:end,i)); end
    for i=1:num_regions, corrP4(i,s,3) = corr(TCMCconv(:,i), NN); end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Correlation between filtered windowed time course and motion time course
    for i=1:num_regions, corrP5(i,s,1) = corr(TCMCconvf(:,i), TCM(:,i)); end
    for i=1:num_regions, corrP5(i,s,2) = corr(TCMCconvf(:,i), TCMconv(:,i)); end
    for i=1:num_regions, corrP5(i,s,3) = corr(TCMCconvf(win_length:end,i), NN); end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Correlation between windowed time course with motion correction within each
    % window and motion time courses (12 regressors)
    for i=1:num_regions, corrP6(i,s,1) = corr(TCMCconv2(:,i), TCM(win_length:end,i)); end
    for i=1:num_regions, corrP6(i,s,2) = corr(TCMCconv2(:,i), TCMconv(win_length:end,i)); end
    for i=1:num_regions, corrP6(i,s,3) = corr(TCMCconv2(:,i), NN); end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute dynamic connectivity with no motion correction
    NMC = sliding_window(TCnoMC, win_length);
    NMC = NMC(:,:,win_length:end); 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute dynamic connectivity with motion-corrected time course
    C = sliding_window(TC,win_length);
    C = C(:,:,win_length:end);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute dynamic connectivity with motion time course
    M = sliding_window(TCM,win_length);
    M = M(:,:,win_length:end);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Correlation between motion time courses and components of TVFC
    for i=1:num_regions

        corrMUXP1(i,s) = corr(NN, squeeze(MUX_P1(i,:)')); 
        corrMUXP2(i,s) = corr(NN, squeeze(MUX_P2(i,:)')); 
        corrMUXP3(i,s) = corr(NN, squeeze(MUX_P3(i,:)')); 
        corrMUXP4(i,s) = corr(NN, squeeze(MUX_P4(i,:)'));
        corrMUXP5(i,s) = corr(NN, squeeze(MUX_P5(i,:)')); 
        corrMUXP6(i,s) = corr(NN, squeeze(MUX_P6(i,:)')); 
                    
        corrSDXP1(i,s) = corr(NN, squeeze(SDX_P1(i,:)')); 
        corrSDXP2(i,s) = corr(NN, squeeze(SDX_P2(i,:)')); 
        corrSDXP3(i,s) = corr(NN, squeeze(SDX_P3(i,:)')); 
        corrSDXP4(i,s) = corr(NN, squeeze(SDX_P4(i,:)'));
        corrSDXP5(i,s) = corr(NN, squeeze(SDX_P5(i,:)')); 
        corrSDXP6(i,s) = corr(NN, squeeze(SDX_P6(i,:)')); 

        for j=1:num_regions

            if (type == 1)

                corrCovxyP1(i,j,s) = corr(NN, squeeze(Cov_P1(i,j,:))); 
                corrCovxyP2(i,j,s) = corr(NN, squeeze(Cov_P2(i,j,:))); 
                corrCovxyP3(i,j,s) = corr(NN, squeeze(Cov_P3(i,j,:))); 
                corrCovxyP4(i,j,s) = corr(NN, squeeze(Cov_P4(i,j,:)));
                corrCovxyP5(i,j,s) = corr(NN, squeeze(Cov_P5(i,j,:))); 
                corrCovxyP6(i,j,s) = corr(NN, squeeze(Cov_P6(i,j,:))); 
    
                corrEXEYP1(i,j,s) = corr(NN, squeeze(EXEY_P1(i,j,:))); 
                corrEXEYP2(i,j,s) = corr(NN, squeeze(EXEY_P2(i,j,:))); 
                corrEXEYP3(i,j,s) = corr(NN, squeeze(EXEY_P3(i,j,:))); 
                corrEXEYP4(i,j,s) = corr(NN, squeeze(EXEY_P4(i,j,:)));
                corrEXEYP5(i,j,s) = corr(NN, squeeze(EXEY_P5(i,j,:))); 
                corrEXEYP6(i,j,s) = corr(NN, squeeze(EXEY_P6(i,j,:))); 
    
                corrEXYP1(i,j,s) = corr(NN, squeeze(EXY_P1(i,j,:))); 
                corrEXYP2(i,j,s) = corr(NN, squeeze(EXY_P2(i,j,:))); 
                corrEXYP3(i,j,s) = corr(NN, squeeze(EXY_P3(i,j,:))); 
                corrEXYP4(i,j,s) = corr(NN, squeeze(EXY_P4(i,j,:)));
                corrEXYP5(i,j,s) = corr(NN, squeeze(EXY_P5(i,j,:))); 
                corrEXYP6(i,j,s) = corr(NN, squeeze(EXY_P6(i,j,:))); 

            elseif (type == 2)
   
                corrCovxyP1(i,j,s) = corr(squeeze(M(i,j,:)), squeeze(Cov_P1(i,j,:))); 
                corrCovxyP2(i,j,s) = corr(squeeze(M(i,j,:)), squeeze(Cov_P2(i,j,:))); 
                corrCovxyP3(i,j,s) = corr(squeeze(M(i,j,:)), squeeze(Cov_P3(i,j,:))); 
                corrCovxyP4(i,j,s) = corr(squeeze(M(i,j,:)), squeeze(Cov_P4(i,j,:)));
                corrCovxyP5(i,j,s) = corr(squeeze(M(i,j,:)), squeeze(Cov_P5(i,j,:))); 
                corrCovxyP6(i,j,s) = corr(squeeze(M(i,j,:)), squeeze(Cov_P6(i,j,:))); 
    
                corrEXEYP1(i,j,s) = corr(squeeze(M(i,j,:)), squeeze(EXEY_P1(i,j,:))); 
                corrEXEYP2(i,j,s) = corr(squeeze(M(i,j,:)), squeeze(EXEY_P2(i,j,:))); 
                corrEXEYP3(i,j,s) = corr(squeeze(M(i,j,:)), squeeze(EXEY_P3(i,j,:))); 
                corrEXEYP4(i,j,s) = corr(squeeze(M(i,j,:)), squeeze(EXEY_P4(i,j,:)));
                corrEXEYP5(i,j,s) = corr(squeeze(M(i,j,:)), squeeze(EXEY_P5(i,j,:))); 
                corrEXEYP6(i,j,s) = corr(squeeze(M(i,j,:)), squeeze(EXEY_P6(i,j,:))); 
    
                corrEXYP1(i,j,s) = corr(squeeze(M(i,j,:)), squeeze(EXY_P1(i,j,:))); 
                corrEXYP2(i,j,s) = corr(squeeze(M(i,j,:)), squeeze(EXY_P2(i,j,:))); 
                corrEXYP3(i,j,s) = corr(squeeze(M(i,j,:)), squeeze(EXY_P3(i,j,:))); 
                corrEXYP4(i,j,s) = corr(squeeze(M(i,j,:)), squeeze(EXY_P4(i,j,:)));
                corrEXYP5(i,j,s) = corr(squeeze(M(i,j,:)), squeeze(EXY_P5(i,j,:))); 
                corrEXYP6(i,j,s) = corr(squeeze(M(i,j,:)), squeeze(EXY_P6(i,j,:))); 

            end

        end
    end

    for i=1:wT, tmp = squeeze(NMC(:,:,i)); CP1(:,i,s) = tmp(Q == 1); end
    for i=1:wT, tmp = squeeze(C(:,:,i)); CP2(:,i,s) = tmp(Q == 1); end
    for i=1:wT, tmp = squeeze(M(:,:,i)); CP3(:,i,s) = tmp(Q == 1); end
    for i=1:wT, tmp = squeeze(MC(:,:,i)); CP4(:,i,s) = tmp(Q == 1); end
    for i=1:wT, tmp = squeeze(MCF(:,:,i)); CP5(:,i,s) = tmp(Q == 1); end
    for i=1:wT, tmp = squeeze(MC2(:,:,i)); CP6(:,i,s) = tmp(Q == 1); end
   
    for i=1:35778
        if (type == 1)

            corrRxyP1(i,s) = corr(NN, squeeze(CP1(i,:,s))'); 
            corrRxyP2(i,s) = corr(NN, squeeze(CP2(i,:,s))'); 
            corrRxyP3(i,s) = corr(NN, squeeze(CP3(i,:,s))'); 
            corrRxyP4(i,s) = corr(NN, squeeze(CP4(i,:,s))'); 
            corrRxyP5(i,s) = corr(NN, squeeze(CP5(i,:,s))'); 
            corrRxyP6(i,s) = corr(NN, squeeze(CP6(i,:,s))'); 

        elseif (type == 2)

            corrRxyP1(i,s) = corr(squeeze(CP3(i,:,s))', squeeze(CP1(i,:,s))'); 
            corrRxyP2(i,s) = corr(squeeze(CP3(i,:,s))', squeeze(CP2(i,:,s))'); 
            corrRxyP3(i,s) = corr(squeeze(CP3(i,:,s))', squeeze(CP3(i,:,s))'); 
            corrRxyP4(i,s) = corr(squeeze(CP3(i,:,s))', squeeze(CP4(i,:,s))'); 
            corrRxyP5(i,s) = corr(squeeze(CP3(i,:,s))', squeeze(CP5(i,:,s))'); 
            corrRxyP6(i,s) = corr(squeeze(CP3(i,:,s))', squeeze(CP6(i,:,s))'); 
            
        end
    end



 disp(s)
end



%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine clusters for each Pipeline


load(which('Network.mat'))
[a,b] = sort(Network);
Q = tril(ones(num_regions,num_regions),-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pipeline 1

% Concatenate data
CP1all = [];
for s=1:num_sub, CP1all = [CP1all squeeze(CP1(:,:,s))]; end

% Perform k-means
[IDXP1, CCP1] = kmeans(CP1all', 2, 'Replicates', 100);

% Organize clusters
Clust1_P1 = Vector2Matrix(CCP1(1,:), Q);
Clust2_P1 = Vector2Matrix(CCP1(2,:), Q);


%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pipeline 2

% Concatenate data
CP2all = [];
for s=1:num_sub, CP2all = [CP2all squeeze(CP2(:,:,s))]; end

% Perform k-means
[IDXP2, CCP2] = kmeans(CP2all', 2, 'Replicates', 100);

% Organize clusters
Clust1_P2 = Vector2Matrix(CCP2(1,:), Q);
Clust2_P2 = Vector2Matrix(CCP2(2,:), Q);


%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pipeline 3

% Concatenate data
CP3all = [];
for s=1:num_sub, CP3all = [CP3all squeeze(CP3(:,:,s))]; end

% Perform k-means
[IDXP3, CCP3] = kmeans(CP3all', 2, 'Replicates', 100);

% Organize clusters
Clust1_P3 = Vector2Matrix(CCP3(1,:), Q);
Clust2_P3 = Vector2Matrix(CCP3(2,:), Q);


%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pipeline 4

% Concatenate data
CP4all = [];
for s=1:num_sub, CP4all = [CP4all squeeze(CP4(:,:,s))]; end

% Perform k-means
[IDXP4, CCP4] = kmeans(CP4all', 2, 'Replicates', 100);

% Organize clusters
Clust1_P4 = Vector2Matrix(CCP4(1,:), Q);
Clust2_P4 = Vector2Matrix(CCP4(2,:), Q);


%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pipeline 5

% Concatenate data
CP5all = [];
for s=1:num_sub, CP5all = [CP5all squeeze(CP5(:,:,s))]; end

% Perform k-means
[IDXP5, CCP5] = kmeans(CP5all', 2, 'Replicates', 100);

% Organize clusters
Clust1_P5 = Vector2Matrix(CCP5(1,:), Q);
Clust2_P5 = Vector2Matrix(CCP5(2,:), Q);


%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pipeline 6

% Concatenate data
CP6all = [];
for s=1:num_sub, CP6all = [CP6all squeeze(CP6(:,:,s))]; end

% Perform k-means
[IDXP6, CCP6] = kmeans(CP6all', 2, 'Replicates', 100);

% Organize clusters
Clust1_P6 = Vector2Matrix(CCP6(1,:), Q);
Clust2_P6 = Vector2Matrix(CCP6(2,:), Q);



%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure: Multidimensional scaling

D = [CCP1; CCP2; CCP3; CCP4; CCP5; CCP6];
DD = corr(D');
Y = cmdscale(DD);

figure
plot(Y(1,1),Y(1,2),'.b', 'MarkerSize', 30)
hold
plot(Y(2,1),Y(2,2),'s', 'MarkerSize', 10,'MarkerFaceColor','blue')
plot(Y(3,1),Y(3,2),'.r', 'MarkerSize', 30)
plot(Y(4,1),Y(4,2),'s', 'MarkerSize', 10,'MarkerFaceColor','red')
plot(Y(5,1),Y(5,2),'.g', 'MarkerSize', 30)
plot(Y(6,1),Y(6,2),'s', 'MarkerSize', 10,'MarkerFaceColor','green')
plot(Y(7,1),Y(7,2),'.c', 'MarkerSize', 30)
plot(Y(8,1),Y(8,2),'s', 'MarkerSize', 10, 'MarkerFaceColor','cyan')
plot(Y(9,1),Y(9,2),'.y', 'MarkerSize', 30)
plot(Y(10,1),Y(10,2),'s', 'MarkerSize', 10, 'MarkerFaceColor','yellow')
plot(Y(11,1),Y(11,2),'.k', 'MarkerSize', 30)
plot(Y(12,1),Y(12,2),'s', 'MarkerSize', 10, 'MarkerFaceColor','black')

xlabel('Dimension 1');
ylabel('Dimension 2');
title('Session 1')
set(gca,'FontSize',18);
axis([-.5 .5 -.5 .5])

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure: Brain States

figure
subplot(6,2,1)
imagesc(Clust1_P1(b,b), [-1 1])
set(gca,'xtick',[])
set(gca,'ytick',[])
%title('State 1','FontSize',16);

subplot(6,2,2)
imagesc(Clust2_P1(b,b), [-1 1])
set(gca,'xtick',[])
set(gca,'ytick',[])
%title('State 2','FontSize',16);

subplot(6,2,3)
imagesc(Clust1_P2(b,b), [-1 1])
set(gca,'xtick',[])
set(gca,'ytick',[])
%title('State 1','FontSize',16);

subplot(6,2,4)
imagesc(Clust2_P2(b,b), [-1 1])
set(gca,'xtick',[])
set(gca,'ytick',[])
%title('State 2','FontSize',16);

subplot(6,2,5)
imagesc(Clust1_P3(b,b), [-1 1])
set(gca,'xtick',[])
set(gca,'ytick',[])
%title('Nuisance State 1','FontSize',16);

subplot(6,2,6)
imagesc(Clust2_P3(b,b), [-1 1])
set(gca,'xtick',[])
set(gca,'ytick',[])
%title('Nuisance State 2','FontSize',16);

subplot(6,2,7)
imagesc(Clust1_P4(b,b), [-1 1])
set(gca,'xtick',[])
set(gca,'ytick',[])
%title('State 1','FontSize',16);

subplot(6,2,8)
imagesc(Clust2_P4(b,b), [-1 1])
set(gca,'xtick',[])
set(gca,'ytick',[])
%title('State 2','FontSize',16);

subplot(6,2,9)
imagesc(Clust1_P5(b,b), [-1 1])
set(gca,'xtick',[])
set(gca,'ytick',[])
%title('State 1','FontSize',16);

subplot(6,2,10)
imagesc(Clust2_P5(b,b), [-1 1])
set(gca,'xtick',[])
set(gca,'ytick',[])
%title('State 2','FontSize',16);

subplot(6,2,11)
imagesc(Clust1_P6(b,b), [-1 1])
set(gca,'xtick',[])
set(gca,'ytick',[])
%title('State 1','FontSize',16);

subplot(6,2,12)
imagesc(Clust2_P6(b,b), [-1 1])
set(gca,'xtick',[])
set(gca,'ytick',[])
%title('State 2','FontSize',16);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure: Transition plots

figure
subplot 611
imagesc(reshape(IDXP1,wT,num_sub)')
xlabel('Time (TR)');
ylabel('Subject');
set(gca,'FontSize',16);

subplot 612
imagesc(reshape(IDXP2,wT,num_sub)')
xlabel('Time (TR)');
ylabel('Subject');
set(gca,'FontSize',16);

subplot 613
imagesc(reshape(IDXP3,wT,num_sub)')
xlabel('Time (TR)');
ylabel('Subject');
set(gca,'FontSize',16);

subplot 614
imagesc(reshape(IDXP4,wT,num_sub)')
xlabel('Time (TR)');
ylabel('Subject');
set(gca,'FontSize',16);

subplot 615
imagesc(reshape(IDXP5,wT,num_sub)')
xlabel('Time (TR)');
ylabel('Subject');
set(gca,'FontSize',16);

subplot 616
imagesc(reshape(IDXP6,wT,num_sub)')
xlabel('Time (TR)');
ylabel('Subject');
set(gca,'FontSize',16);


%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modeling transitions and transition times

% Find transition points
[numTrans1, Dwell1] = FindTransition(IDXP1, wT, num_sub);
[numTrans2, Dwell2] = FindTransition(IDXP2, wT, num_sub);
[numTrans3, Dwell3] = FindTransition(IDXP3, wT, num_sub);
[numTrans4, Dwell4] = FindTransition(IDXP4, wT, num_sub);
[numTrans5, Dwell5] = FindTransition(IDXP5, wT, num_sub);
[numTrans6, Dwell6] = FindTransition(IDXP6, wT, num_sub);

% Model number of transitions
[lambdahat1,lambdaci1] = poissfit(numTrans1);
[lambdahat2,lambdaci2] = poissfit(numTrans2);
[lambdahat3,lambdaci3] = poissfit(numTrans3);
[lambdahat4,lambdaci4] = poissfit(numTrans4);
[lambdahat5,lambdaci5] = poissfit(numTrans5);
[lambdahat6,lambdaci6] = poissfit(numTrans6);

% Model time to transition
[muhat1,muhatCI1] = expfit([Dwell1{:}]);
[muhat2,muhatCI2] = expfit([Dwell2{:}]);
[muhat3,muhatCI3] = expfit([Dwell3{:}]);
[muhat4,muhatCI4] = expfit([Dwell4{:}]);
[muhat5,muhatCI5] = expfit([Dwell5{:}]);
[muhat6,muhatCI6] = expfit([Dwell6{:}]);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure: Modeling transitions and transition times

t=0:1:20;
figure
subplot 121
plot(poisspdf(t, lambdahat1),'LineWidth',2)
hold
plot(poisspdf(t, lambdahat2),'LineWidth',2)
plot(poisspdf(t, lambdahat3),'LineWidth',2)
plot(poisspdf(t, lambdahat4),'LineWidth',2)
plot(poisspdf(t, lambdahat5),'LineWidth',2)
plot(poisspdf(t, lambdahat6),'LineWidth',2)
legend('Pipeline 1', 'Pipeline 2', 'Pipeline 3','Pipeline 4','Pipeline 5','Pipeline 6')
title('Number of Transitions')
xlabel('Time (TR)');
set(gca,'FontSize',18);

t=0:1:150;
subplot 122
plot(exppdf(t, muhat1),'LineWidth',2)
hold
plot(exppdf(t, muhat2),'LineWidth',2)
plot(exppdf(t, muhat3),'LineWidth',2)
plot(exppdf(t, muhat4),'LineWidth',2)
plot(exppdf(t, muhat5),'LineWidth',2)
plot(exppdf(t, muhat6),'LineWidth',2)
legend('Pipeline 1', 'Pipeline 2', 'Pipeline 3','Pipeline 4','Pipeline 5','Pipeline 6')
title('Transition Time')
xlabel('Time (TR)');
set(gca,'FontSize',18);


%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure: Correlations between motion time courses and processed timecourses

figure

subplot 311
plot((1:num_regions),squeeze(corrP1(b,:,1)),'.b', 'MarkerSize', .1)
hold
plot(((num_regions+1):2*num_regions),squeeze(corrP2(b,:,1)),'.r', 'MarkerSize', .1)
plot(((2*num_regions+1):3*num_regions),squeeze(corrP3(b,:,1)),'.g', 'MarkerSize', .1)
plot(((3*num_regions+1):4*num_regions),squeeze(corrP4(b,:,1)),'.c', 'MarkerSize', .1)
plot(((4*num_regions+1):5*num_regions),squeeze(corrP5(b,:,1)),'.y', 'MarkerSize', .1)
plot(((5*num_regions+1):6*num_regions),squeeze(corrP6(b,:,1)),'.black', 'MarkerSize', .1)
axis([1 6*num_regions -1.02 1.02])

xticks([50 150 250 num_regions+[50 150 250 ] 2*num_regions+[50 150 250 ] 3*num_regions+[50 150 250 ] 4*num_regions+[50 150 250] 5*num_regions+[50 150 250 ]])
xticklabels([50 150 250 50 150 250 50 150 250 50 150 250 50 150 250 50 150 250])
xlabel('Region');
title('Correlation between motion and processed timecourses')
set(gca,'FontSize',16);


subplot 312
plot((1:num_regions),squeeze(corrP1(b,:,2)),'.b', 'MarkerSize', .1)
hold
plot(((num_regions+1):2*num_regions),squeeze(corrP2(b,:,2)),'.r', 'MarkerSize', .1)
plot(((2*num_regions+1):3*num_regions),squeeze(corrP3(b,:,2)),'.g', 'MarkerSize', .1)
plot(((3*num_regions+1):4*num_regions),squeeze(corrP4(b,:,2)),'.c', 'MarkerSize', .1)
plot(((4*num_regions+1):5*num_regions),squeeze(corrP5(b,:,2)),'.y', 'MarkerSize', .1)
plot(((5*num_regions+1):6*num_regions),squeeze(corrP6(b,:,2)),'.black', 'MarkerSize', .1)
axis([1 6*num_regions -1.02 1.02])

xticks([50 150 250 num_regions+[50 150 250 ] 2*num_regions+[50 150 250 ] 3*num_regions+[50 150 250 ] 4*num_regions+[50 150 250 ] 5*num_regions+[50 150 250 ]])
xticklabels([50 150 250 50 150 250 50 150 250 50 150 250 50 150 250 50 150 250 ])
xlabel('Region');
title('Correlation between windowed motion and processed timecourses')
set(gca,'FontSize',16);


subplot 313
plot((1:num_regions),squeeze(corrP1(b,:,3)),'.b', 'MarkerSize', .1)
hold
plot(((num_regions+1):2*num_regions),squeeze(corrP2(b,:,3)),'.r', 'MarkerSize', .1)
plot(((2*num_regions+1):3*num_regions),squeeze(corrP3(b,:,3)),'.g', 'MarkerSize', .1)
plot(((3*num_regions+1):4*num_regions),squeeze(corrP4(b,:,3)),'.c', 'MarkerSize', .1)
plot(((4*num_regions+1):5*num_regions),squeeze(corrP5(b,:,3)),'.y', 'MarkerSize', .1)
plot(((5*num_regions+1):6*num_regions),squeeze(corrP6(b,:,3)),'.black', 'MarkerSize', .1)
axis([1 6*num_regions -1.02 1.02])

xticks([50 150 250 num_regions+[50 150 250 ] 2*num_regions+[50 150 250 ] 3*num_regions+[50 150 250 ] 4*num_regions+[50 150 250 ] 5*num_regions+[50 150 250 ]])
xticklabels([50 150 250 50 150 250 50 150 250 50 150 250 50 150 250 50 150 250 ])
xlabel('Region');
title('Correlation between nuisance norm and processed timecourses')
set(gca,'FontSize',16);





%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure: Correlations between nuisance norm and dynamic correlation matrix
% (Manhattan Plots)


figure

subplot 411
plot((1:num_connections),corrRxyP1,'.b', 'MarkerSize', .1)
hold
plot(((num_connections+1):2*num_connections),corrRxyP2,'.r', 'MarkerSize', .1)
plot(((2*num_connections+1):3*num_connections),corrRxyP3,'.g', 'MarkerSize', .1)
plot(((3*num_connections+1):4*num_connections),corrRxyP4,'.c', 'MarkerSize', .1)
plot(((4*num_connections+1):5*num_connections),corrRxyP5,'.y', 'MarkerSize', .1)
plot(((5*num_connections+1):6*num_connections),corrRxyP6,'.black', 'MarkerSize', .1)
axis([1 6*num_connections -1 1])

xticks([15000 30000 num_connections+[15000 30000 ] 2*num_connections+[15000 30000 ] 3*num_connections+[15000 30000 ] 4*num_connections+[15000 30000 ] 5*num_connections+[15000 30000 ]])
xticklabels([15000 30000 15000 30000 15000 30000 15000 30000 15000 30000 15000 30000])
xlabel('Region');
title('Correlation between nuisance norm and time-varying correlation')
set(gca,'FontSize',16);


subplot 412
XYP1 = zeros(num_connections, num_sub);
XYP2 = zeros(num_connections, num_sub);
XYP3 = zeros(num_connections, num_sub);
XYP4 = zeros(num_connections, num_sub);
XYP5 = zeros(num_connections, num_sub);
XYP6 = zeros(num_connections, num_sub);

for i=1:num_sub,
    tmp = corrEXYP1(:,:,i);
    XYP1(:,i) = tmp(Q==1);
    tmp = corrEXYP2(:,:,i);
    XYP2(:,i) = tmp(Q==1); 
    tmp = corrEXYP3(:,:,i);
    XYP3(:,i) = tmp(Q==1);
    tmp = corrEXYP4(:,:,i);
    XYP4(:,i) = tmp(Q==1);     
    tmp = corrEXYP5(:,:,i);
    XYP5(:,i) = tmp(Q==1);     
    tmp = corrEXYP6(:,:,i);
    XYP6(:,i) = tmp(Q==1);     
end

plot((1:num_connections),XYP1,'.b', 'MarkerSize', .1)
hold
plot(((num_connections+1):2*num_connections),XYP2,'.r', 'MarkerSize', .1)
plot(((2*num_connections+1):3*num_connections),XYP3,'.g', 'MarkerSize', .1)
plot(((3*num_connections+1):4*num_connections),XYP4,'.c', 'MarkerSize', .1)
plot(((4*num_connections+1):5*num_connections),XYP5,'.y', 'MarkerSize', .1)
plot(((5*num_connections+1):6*num_connections),XYP6,'.black', 'MarkerSize', .1)
axis([1 6*num_connections -1 1])

xticks([15000 30000 num_connections+[15000 30000 ] 2*num_connections+[15000 30000 ] 3*num_connections+[15000 30000 ] 4*num_connections+[15000 30000 ] 5*num_connections+[15000 30000 ]])
xticklabels([15000 30000 15000 30000 15000 30000 15000 30000 15000 30000 15000 30000])
xlabel('Region');
title('Correlation between nuisance norm and E(XY)')
set(gca,'FontSize',16);


subplot 413
pP1 = zeros(num_connections, num_sub);
pP2 = zeros(num_connections, num_sub);
pP3 = zeros(num_connections, num_sub);
pP4 = zeros(num_connections, num_sub);
pP5 = zeros(num_connections, num_sub);
pP6 = zeros(num_connections, num_sub);

for i=1:num_sub,
    tmp = corrEXEYP1(:,:,i);
    pP1(:,i) = tmp(Q==1);
    tmp = corrEXEYP2(:,:,i);
    pP2(:,i) = tmp(Q==1); 
    tmp = corrEXEYP3(:,:,i);
    pP3(:,i) = tmp(Q==1);
    tmp = corrEXEYP4(:,:,i);
    pP4(:,i) = tmp(Q==1);     
    tmp = corrEXEYP5(:,:,i);
    pP5(:,i) = tmp(Q==1);    
    tmp = corrEXEYP6(:,:,i);
    pP6(:,i) = tmp(Q==1);    
end

plot((1:num_connections),pP1,'.b', 'MarkerSize', 1)
hold
plot(((num_connections+1):2*num_connections),pP2,'.r', 'MarkerSize', 1)
plot(((2*num_connections+1):3*num_connections),pP3,'.g', 'MarkerSize', 1)
plot(((3*num_connections+1):4*num_connections),pP4,'.c', 'MarkerSize', 1)
plot(((4*num_connections+1):5*num_connections),pP5,'.y', 'MarkerSize', 1)
plot(((5*num_connections+1):6*num_connections),pP6,'.black', 'MarkerSize', 1)
axis([1 6*num_connections -1 1])

xticks([15000 30000 num_connections+[15000 30000 ] 2*num_connections+[15000 30000 ] 3*num_connections+[15000 30000 ] 4*num_connections+[15000 30000 ] 5*num_connections+[15000 30000 ]])
xticklabels([15000 30000 15000 30000 15000 30000 15000 30000 15000 30000 15000 30000])
xlabel('Region');
title('Correlation between nuisance norm and E(X)E(Y)')
set(gca,'FontSize',16);


subplot 427

plot((1:num_regions),corrMUXP1,'.b', 'MarkerSize', 1)
hold
plot(((num_regions+1):2*num_regions),corrMUXP2,'.r', 'MarkerSize', 1)
plot(((2*num_regions+1):3*num_regions),corrMUXP3,'.g', 'MarkerSize', 1)
plot(((3*num_regions+1):4*num_regions),corrMUXP4,'.c', 'MarkerSize', 1)
plot(((4*num_regions+1):5*num_regions),corrMUXP5,'.y', 'MarkerSize', 1)
plot(((5*num_regions+1):6*num_regions),corrMUXP6,'.black', 'MarkerSize', 1)
axis([1 6*num_regions -1 1])

xticks([50 150 250 num_regions+[50 150 250 ] 2*num_regions+[50 150 250 ] 3*num_regions+[50 150 250 ] 4*num_regions+[50 150 250 ] 5*num_regions+[50 150 250 ]])
xticklabels([50 150 250 50 150 250 50 150 250 50 150 250 50 150 250 50 150 250])
xlabel('Region');
title('Correlation between nuisance norm and EX')
set(gca,'FontSize',16);


subplot 428
plot((1:num_regions),corrSDXP1,'.b', 'MarkerSize', 1)
hold
plot(((num_regions+1):2*num_regions),corrSDXP2,'.r', 'MarkerSize', 1)
plot(((2*num_regions+1):3*num_regions),corrSDXP3,'.g', 'MarkerSize', 1)
plot(((3*num_regions+1):4*num_regions),corrSDXP4,'.c', 'MarkerSize', 1)
plot(((4*num_regions+1):5*num_regions),corrSDXP5,'.y', 'MarkerSize', 1)
plot(((5*num_regions+1):6*num_regions),corrSDXP6,'.black', 'MarkerSize', 1)
axis([1 6*num_regions -1 1])

xticks([50 150 250 num_regions+[50 150 250 ] 2*num_regions+[50 150 250 ] 3*num_regions+[50 150 250 ] 4*num_regions+[50 150 250 ] 5*num_regions+[50 150 250 ]])
xticklabels([50 150 250 50 150 250 50 150 250 50 150 250 50 150 250 50 150 250 ])
xlabel('Region');
title('Correlation between nuisance norm and SD(X)')
set(gca,'FontSize',16);




%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure: Proportion motion

corrM = zeros(num_regions, num_regions);
corrT = zeros(num_regions, num_regions);
for i=1:num_regions
    for j=1:num_regions 
        tmp = mean((squeeze(Comp1(i,j,:)))); 
        corrM(i,j) = tmp; 
        tmp = mean((squeeze(Comp2(i,j,:)))); corrT(i,j) = tmp; 
    end
end

Prop = (corrM.^2)./((corrM.^2) + (corrT.^2));

Q1 = zeros(8,8);
for i=1:8, for j=1:8, Q1(i,j) = mean(mean(corrM(Network == i, Network == j))); end; end

Q2 = zeros(8,8);
for i=1:8, for j=1:8, Q2(i,j) = mean(mean(corrT(Network == i, Network == j))); end; end

PropNet = (Q1.^2)./((Q1.^2) + (Q2.^2));

figure
subplot 121
imagesc(Prop,[0 1])
set(gca,'xtick',[])
set(gca,'ytick',[])

subplot 122
imagesc(PropNet,[0 1])
set(gca,'xtick',[])
set(gca,'ytick',[])

colorbar
set(gca,'FontSize',16);

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure: Correlation between nuisance norm and time-varying correlation
% componets with Violin plots


D1 = zeros(num_connections*20,6);
D1(:,1) = reshape(corrRxyP1,num_connections*20,1);
D1(:,2) = reshape(corrRxyP2,num_connections*20,1);
D1(:,3) = reshape(corrRxyP3,num_connections*20,1) + normrnd(0,0.000001, num_connections*20, 1);
D1(:,4) = reshape(corrRxyP4,num_connections*20,1);
D1(:,5) = reshape(corrRxyP5,num_connections*20,1);
D1(:,6) = reshape(corrRxyP6,num_connections*20,1);

figure
violin(D1,'Labels',{'P1','P2', 'P3', 'P4', 'P5','P6'})
title('Correlation between nuisance norm and time-varying correlation')
ylabel('Correlation');
xlabel('Pipeline');
axis([0 7 -1.04 1.04])
set(gca,'FontSize',16);



XYP1 = zeros(num_connections, num_sub);
XYP2 = zeros(num_connections, num_sub);
XYP3 = zeros(num_connections, num_sub);
XYP4 = zeros(num_connections, num_sub);
XYP5 = zeros(num_connections, num_sub);
XYP6 = zeros(num_connections, num_sub);

for i=1:num_sub,
    tmp = corrEXYP1(:,:,i);
    XYP1(:,i) = tmp(Q==1);
    tmp = corrEXYP2(:,:,i);
    XYP2(:,i) = tmp(Q==1); 
    tmp = corrEXYP3(:,:,i);
    XYP3(:,i) = tmp(Q==1);
    tmp = corrEXYP4(:,:,i);
    XYP4(:,i) = tmp(Q==1);     
    tmp = corrEXYP5(:,:,i);
    XYP5(:,i) = tmp(Q==1);     
    tmp = corrEXYP6(:,:,i);
    XYP6(:,i) = tmp(Q==1);     
end


D2 = zeros(num_connections*20,6);
D2(:,1) = reshape(XYP1,num_connections*20,1);
D2(:,2) = reshape(XYP2,num_connections*20,1);
D2(:,3) = reshape(XYP3,num_connections*20,1);
D2(:,4) = reshape(XYP4,num_connections*20,1);
D2(:,5) = reshape(XYP5,num_connections*20,1);
D2(:,6) = reshape(XYP6,num_connections*20,1);


figure
violin(D2,'Labels',{'P1','P2', 'P3', 'P4', 'P5','P6'})
title('Correlation between nuisance norm and (xy)*w')
ylabel('Correlation');
xlabel('Pipeline');
axis([0 7 -1.04 1.04])
set(gca,'FontSize',16);


pP1 = zeros(num_connections, num_sub);
pP2 = zeros(num_connections, num_sub);
pP3 = zeros(num_connections, num_sub);
pP4 = zeros(num_connections, num_sub);
pP5 = zeros(num_connections, num_sub);
pP6 = zeros(num_connections, num_sub);

for i=1:num_sub,
    tmp = corrEXEYP1(:,:,i);
    pP1(:,i) = tmp(Q==1);
    tmp = corrEXEYP2(:,:,i);
    pP2(:,i) = tmp(Q==1); 
    tmp = corrEXEYP3(:,:,i);
    pP3(:,i) = tmp(Q==1);
    tmp = corrEXEYP4(:,:,i);
    pP4(:,i) = tmp(Q==1);     
    tmp = corrEXEYP5(:,:,i);
    pP5(:,i) = tmp(Q==1);    
    tmp = corrEXEYP6(:,:,i);
    pP6(:,i) = tmp(Q==1);    
end

D3 = zeros(num_connections*20,6);
D3(:,1) = reshape(pP1,num_connections*20,1);
D3(:,2) = reshape(pP2,num_connections*20,1);
D3(:,3) = reshape(pP3,num_connections*20,1);
D3(:,4) = reshape(pP4,num_connections*20,1);
D3(:,5) = reshape(pP5,num_connections*20,1);
D3(:,6) = reshape(pP6,num_connections*20,1);


figure
violin(D3,'Labels',{'P1','P2', 'P3', 'P4', 'P5','P6'})
title('Correlation between nuisance norm and (x*w)E(y*w)')
ylabel('Correlation');
xlabel('Pipeline');
axis([0 7 -1.04 1.04])
set(gca,'FontSize',16);


D4 = zeros(num_regions*20,6);
D4(:,1) = reshape(corrMUXP1,num_regions*20,1);
D4(:,2) = reshape(corrMUXP2,num_regions*20,1);
D4(:,3) = reshape(corrMUXP3,num_regions*20,1);
D4(:,4) = reshape(corrMUXP4,num_regions*20,1);
D4(:,5) = reshape(corrMUXP5,num_regions*20,1);
D4(:,6) = reshape(corrMUXP6,num_regions*20,1);


figure
violin(D4,'Labels',{'P1','P2', 'P3', 'P4', 'P5','P6'})
title('Correlation between nuisance norm and x*w')
xlabel('Pipeline');
axis([0 7 -1.04 1.04])
set(gca,'FontSize',16);


D5 = zeros(num_regions*20,6);
D5(:,1) = reshape(corrSDXP1,num_regions*20,1);
D5(:,2) = reshape(corrSDXP2,num_regions*20,1);
D5(:,3) = reshape(corrSDXP3,num_regions*20,1);
D5(:,4) = reshape(corrSDXP4,num_regions*20,1);
D5(:,5) = reshape(corrSDXP5,num_regions*20,1);
D5(:,6) = reshape(corrSDXP6,num_regions*20,1);

figure
violin(D5,'Labels',{'P1','P2', 'P3', 'P4', 'P5','P6'})
title('Correlation between nuisance norm and (x^2)*w')
xlabel('Pipeline');
axis([0 7 -1.04 1.04])
set(gca,'FontSize',16);
