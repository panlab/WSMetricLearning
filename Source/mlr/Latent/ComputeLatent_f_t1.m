function [Psi, hLatent, Rlatent, tord] = ComputeLatent_f_t1(k, X, hLatent, numSign, W, Ypos, Yneg, ...
    batchSize, Template, Nrotate, SAMPLES, ClassScores)
%%% Psi(x_i,y_i)
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

    global PSIALL CPGRADIENT DISTANCE;
    L = zeros(size(W));
    for i = 1:size(W,3)
        [vecs,vals] = eig(0.5 * (W(:,:,i) + W(:,:,i)'));
        L(:,:,i) = real(abs(vals)).^0.5 * vecs';
    end
    tord = [];
    clear 'W'

    [d,n,m] = size(getLatent(numSign, X, hLatent));
    Rlatent = zeros([length(hLatent), 1]);

    Xtemplate = cell2mat(X(Template));
    XTest = (X(setdiff(1:n, Template)));
    clear 'X'

    S       = zeros(n);
    dIndex  = sub2ind([n n], 1:n, 1:n);
    nbase = length(Template) / Nrotate;
    ntrain = length(Template);
    v_ps = zeros(Nrotate, 1);
    for jj = 1:Nrotate 
        v_ps(jj) = (jj - 1)*nbase;
    end
%     rangeS = [0,cumsum(numSign(1:end-1))];
%     rangeE = [cumsum(numSign)];
    
    Tnum = numSign(setdiff((1:n)', Template));
    
    if isempty(ClassScores)
        TS  = zeros(batchSize, n);
        Xtest = (cell2mat(XTest));
        
        D1 = mat2cell(setDistanceFullMKL_f([Xtemplate Xtest], L, ...
            ntrain + (1:size(Xtest, 2)), 1:ntrain), ntrain,Tnum);

        
        parfor j = SAMPLES'
                if isempty(Ypos{j})
                    continue;
                end
                if isempty(Yneg)
                    % Construct a negative set 
                    Ynegative = setdiff((1:n)', [j ; Ypos{j}]);
                else
                    Ynegative = Yneg{j};
                end
                D    = D1{j-ntrain};
                [YposT, YnegativeT] = MultiY(Ypos{j}, Ynegative, nbase, Nrotate);
                ScorePos        = - D(YposT,:);
                ScoreNeg        = - D(YnegativeT,:);
                
                ScoreNeg = reshape(ScoreNeg, [size(ScoreNeg,1) / Nrotate, Nrotate, size(ScoreNeg, 2)]);
                C1 = ScorePos - reshape(sum(ScoreNeg, 1) / size(ScoreNeg, 1), size(ScorePos));                
                [C, ord] = sort(C1(:));
                [Rlatent(j), hLatent(j)] = ind2sub(size(C1), ord(end));
                T1 =  (Rlatent(j)- 1)*nbase;
                TS(j-ntrain,:)         = metricPsiPOALL(j, [Ypos{j}; Ynegative'] + T1,  ...
                    n, Ypos{j} + T1, Ynegative + T1);
        end
        % Reconstruct the S matrix from TS
        S(SAMPLES,:)    = TS;
        S(:,SAMPLES)    = S(:,SAMPLES) + TS';
        S(dIndex)       = S(dIndex) - sum(TS, 1);
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
                TS(x,:)     = PSIALL(i, n, Ypos, Yneg);
            end

            S(points,:) = S(points,:) + TS;
            S(:,points) = S(:,points) + TS';
            S(dIndex)   = S(dIndex) - sum(TS, 1);
        end
    end
    Psi    = CPGRADIENT(getLatent(numSign, mat2cell([Xtemplate, Xtest], ...
        size(Xtemplate,1), numSign), hLatent), S, batchSize);
end