function [dPsi, M, SO_time, hLatent, RLatent] = cuttingPlaneFull_Latent_f_t1(k, X, hLatent, RLatent, ...
    PsiLatent, numSign, W, Ypos, Yneg, batchSize, Template, Nrotate, SAMPLES, ClassScores)
% [dPsi, M, SO_time] = cuttingPlaneFull(k, X, W, Yp, Yn, batchSize, SAMPLES, ClassScores)
%
%   k           = k parameter for the SO
%   X           = d*n data matrix
%   W           = d*d PSD metric
%   Yp          = cell-array of relevant results for each point
%   Yn          = cell-array of irrelevant results for each point
%   batchSize   = number of points to use in the constraint batch
%   SAMPLES     = indices of valid points to include in the batch
%   ClassScores = structure for synthetic constraints
%
%   dPsi        = dPsi vector for this batch
%   M           = mean loss on this batch
%   SO_time     = time spent in separation oracle

    global SO PSI PSIALL DISTANCE CPGRADIENT;

    [d,n,m] = size(getLatent(numSign, X, hLatent));
    L = zeros(size(W));
    parfor i = 1:size(W,3)
        [vecs,vals] = eig(0.5 * (W(:,:,i) + W(:,:,i)'));
        L(:,:,i) = real(abs(vals)).^0.5 * vecs';
    end  
    clear 'W'
%     dis = (D1(:,2:end) - Dwgt);max(abs(dis(:)))
    Xtemplate = cell2mat(X(Template));
    XTest = (X(setdiff(1:n, Template)));
    Xtest = (cell2mat(XTest));
    clear 'X';clear 'W';clear 'XTest'

    M       = 0;
    S       = zeros(n);
    
%     S1      = zeros(n);
    
    
    dIndex  = sub2ind([n n], 1:n, 1:n);
    nbase = length(Template) / Nrotate;
    SO_time = 0;
    ntrain = length(Template);
    Range = zeros(Nrotate, 2);plus = zeros(Nrotate, 1);
    for jj = 1:Nrotate 
        Range(jj,1) = nbase*(jj-1) + 1;
        Range(jj,2) = nbase*jj;
        plus(jj) = (jj - 1)*nbase;      
    end
                   
    Tnum = numSign(setdiff((1:n)', Template));
    
    
    AH_index = cell(length(Tnum), 1)';
    AC_index = cell(length(Tnum), 1)';
    for i = 1:length(Tnum)
        AH_index{i} = [1:Tnum(i)]';
        AC_index{i} = (i+ntrain)*ones(Tnum(i), 1);
    end    
    
    if isempty(ClassScores)
        TS  = zeros(batchSize, n);

        
        ntest = size(Xtest, 2);
        D1 = setDistanceFullMKL_f([Xtemplate Xtest], L, ...
                ntrain + (1:ntest), 1:ntrain);
        TSAMPLES = [1:ntest*Nrotate];
            yik = zeros(nbase, length(TSAMPLES));
            lik = zeros(1, length(TSAMPLES));
            sik = zeros(1, length(TSAMPLES));
            
            xx = cell2mat(AC_index);
            C_index = repmat(xx(:),  [Nrotate, 1]);
            R_index = repmat([1:Nrotate], [ntest, 1]);R_index = R_index(:);
            xx = cell2mat(AH_index);
            H_index = repmat(xx(:),  [Nrotate, 1]);
            
            B_index = [0:Nrotate-1]*ntest;
            
            SO_start        = tic();
                 
            for j = TSAMPLES
                Ci = C_index(j);
                if isempty(Ypos{Ci})
                    continue;
                end
                if isempty(Yneg)
                    Ynegative = setdiff((1:n)', [Ci ; Ypos{Ci}]);
                else
                    Ynegative = Yneg{Ci};
                end
                D    = D1([1:nbase]+plus(R_index(j)), j-B_index(R_index(j)));
                [yik(:,j), lik(j), sik(j)]   =   SO(1, D, Ypos{Ci}, Ynegative, k);
            end
            SO_time         = SO_time + toc(SO_start);
            
             parfor j = SAMPLES'
                if isempty(Yneg)
                    Ynegative = setdiff((1:n)', [j ; Ypos{j}]);
                else
                    Ynegative = Yneg{j};
                end
                index = find(C_index == j);
                [si, ord] = max(sik(index));
                hLatent(j) = H_index(index(ord(1)));
                RLatent(j) = R_index(index(ord(1)));

                yi = yik(:,index(ord(1)));li = lik(index(ord(1)));
                M               = M + li /batchSize;
                T1 = (RLatent(j)-1) * nbase;
                TS(j-ntrain,:)         = metricPsiPOALL(j, yi'+T1, n, Ypos{j}+...
                    T1, Ynegative+T1);   
             end
             
        % Reconstruct the S matrix from TS
        S(SAMPLES,:)    = TS;
        S(:,SAMPLES)    = S(:,SAMPLES) + TS';
        S(dIndex)       = S(dIndex) - sum(TS, 1);
        
%         %%%for debug
%         S1(SAMPLES,:)   = TS1;
%         S1(:,SAMPLES)   = S1(:,SAMPLES) + TS1';
%         S1(dIndex)      = S1(dIndex) - sum(TS1, 1);
%         %%%end debug
        
    else

        % Do it class-wise for efficiency
        batchSize = 0;
        for j = 1:length(ClassScores.classes)
            c       = ClassScores.classes(j);
            points  = find(ClassScores.Y == c);

            Yneg    = find(ClassScores.Yneg{j});
            yp      = ClassScores.Ypos{j};
            
            if length(points) <= 1
                continue;
            end

            batchSize = batchSize + length(points);
            TS      = zeros(length(points), n);
            parfor x = 1:length(points)
                i           = points(x);
                yl          = yp;
                yl(i)       = 0;
                Ypos        = find(yl);
                SO_start    = tic();
                    [yi, li]    = SO(i, D, Ypos, Yneg, k);
                SO_time     = SO_time + toc(SO_start);

                M           = M + li;
                TS(x,:)     = PSI(i, yi', n, Ypos, Yneg);
            end

            S(points,:) = S(points,:) + TS;
            S(:,points) = S(:,points) + TS';
            S(dIndex)   = S(dIndex) - sum(TS, 1);
        end
        M   = M / batchSize;
    end
    
%     %%%When hLatent = 1
%     %%%for debug
%     dPsi1   = CPGRADIENT(getLatent(numSign, X, hLatent), S1, batchSize);
%     
%     dPsi    = CPGRADIENT(getLatent(numSign, X, hLatent), S, batchSize);
%     dPsi2   = PsiLatent - dPsi;
%     dis = dPsi1 - dPsi2;
%     max(abs(dis(:)))
%     pause
%     %%%end debug
%     
    dPsi    = PsiLatent - CPGRADIENT(getLatent(numSign, mat2cell([Xtemplate, Xtest], ...
        size(Xtemplate,1), numSign), hLatent), S, batchSize);
end
