function [dPsi, M, SO_time, hLatent, RLatent, X] = cuttingPlaneFull_latent_Large_f1(k, X, hLatent, RLatent, ...
    PsiLatent, numSign, CIndex, W, datastr, Ypos, Yneg, batchSize, Template, Nrotate, SAMPLES, ClassScores)
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
L = zeros(size(W));
    for i = 1:size(W,3)
        [vecs,vals] = eig(0.5 * (W(:,:,i) + W(:,:,i)'));
        L(:,:,i) = real(abs(vals)).^0.5 * vecs';
    end
%     clear 'W'
    
    [d,n,m] = size(X);
    dPsi = PsiLatent;
    M       = 0;
    BatchRound = size(CIndex, 1);
    
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

        smid = 0;
        for R = 1:BatchRound
            index = [CIndex(R, 1):CIndex(R, 2)];
            
            load(fullfile(datastr,[num2str(R) '.mat']), 'Xcell'); 
            Xtest = (cell2mat(Xcell))';
            clear 'Xcell'
            ntest = size(Xtest, 2);
            
            D1 = setDistanceFullMKL_f([X(:, Template) Xtest], L, ...
                ntrain + (1:ntest), 1:ntrain);
            
            TSAMPLES = [1:ntest*Nrotate];
            yik = zeros(nbase, length(TSAMPLES));
            lik = zeros(1, length(TSAMPLES));
            sik = zeros(1, length(TSAMPLES));
            
            xx = cell2mat(AC_index(index));
            C_index = repmat(xx(:),  [Nrotate, 1]);
            R_index = repmat([1:Nrotate], [ntest, 1]);R_index = R_index(:);
            xx = cell2mat(AH_index(index));
            H_index = repmat(xx(:),  [Nrotate, 1]);
            B_index = [0:Nrotate-1]*ntest;
            SO_start        = tic();
            
            
            MRR  = ((1:(nbase)).^-1)';
            Ypos1 = cell2mat(Ypos(C_index));R =1 /(nbase-1);
            D1 = permute(reshape(D1, [nbase, Nrotate, ntest]), [1 3 2]);
            D1 = -reshape(D1, [nbase, length(TSAMPLES)]);
            parfor j = TSAMPLES
                D = D1(:,j);
                ind1 = Ypos1(j);         
                Dpos = D(ind1);
                D(ind1) = inf;
                [Dneg,ord] = sort(D, 'descend');ord = ord(2:end);
                Dneg(1,:) = Dpos';
                Td = Dpos - Dneg;
                std = cumsum(Td)*R; 
                I = -2*std-MRR;
                [I1, Index] = max(I);
                ord1 = [0;ord];
                ord1(Index) =ind1;
                ord1(1:Index-1) = ord(1:Index-1);
                yik(:,j) = ord1;
                lik(j) = (1 - MRR(Index))';
                sik(j) = std(end,:) + I1+1;
            end
                
%             for j = TSAMPLES
%                 Ci = C_index(j);
%                 Ynegative = Yneg{Ci};
%                 D    = D1([1:nbase]+plus(R_index(j)), j-B_index(R_index(j)));
%                 [yik(:,j), lik(j), sik(j)]   =   SO(1, D, Ypos{Ci}, Ynegative, k);
%             end
            SO_time         = SO_time + toc(SO_start);
            
            TSAMPLES = SAMPLES(index);
            
            S       = zeros(ntrain+length(TSAMPLES));
            n1 = ntrain+length(TSAMPLES);
            TS  = zeros(length(TSAMPLES), n1);
            ntrainS = ntrain+smid;
            Tsum = [0, cumsum(numSign(TSAMPLES'))];

            Cs = 0;
            for j = TSAMPLES'
                Ynegative = Yneg{j};
                index = find(C_index == j);
                [si, ord] = max(sik(index));
                hLatent(j) = H_index(index(ord(1)));
                RLatent(j) = R_index(index(ord(1)));

                yi = yik(:,index(ord(1)));li = lik(index(ord(1)));
                M               = M + li /batchSize;
                X(:, j) = Xtest(:, hLatent(j)+Tsum(j-ntrainS));
                TS(j-ntrainS,:)         = PSIALL(j, yi'+plus(RLatent(j)), n1, Ypos{j}+...
                    plus(RLatent(j)), Ynegative+plus(RLatent(j)));   
                
                Cs = Cs + si - li;
                
            end

            
            % Reconstruct the S matrix from TS
            S(TSAMPLES-smid,:)    = TS;
            S(:,TSAMPLES-smid)    = S(:,TSAMPLES-smid) + TS';
            dIndex1  = sub2ind([n1 n1], 1:n1, 1:n1);
            S(dIndex1)       = S(dIndex1) - sum(TS, 1);
            smid = smid + length(TSAMPLES);
            dPsi    = dPsi - CPGRADIENT(X(:, [1:ntrain, TSAMPLES']), S, batchSize);

            
            
            if abs(sum(sum(W .* CPGRADIENT(X(:, [1:ntrain, TSAMPLES']), S, batchSize))) - Cs / batchSize) > 1e-6
                t = 1;
                pause
            end
        end
        
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
            for x = 1:length(points)
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

end


