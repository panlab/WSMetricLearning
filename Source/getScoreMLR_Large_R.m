function [AllScore, ranktime, tr_size, Metric, PCACOEFF] = getScoreMLR_Large(...
    TrainFun, feattype, ts_fold_idx, ts_idx, tr_idx, setting, fdatabase, ...
    Mstr, Samplevoted, NotRatio, cindex, normmerge, knnpara, cpara)   
tr_size = [];

% % % % setting.Flatent = 1;
% % % % % setting.Mclassfy = 0;
% % % % % setting.testSOutput = 1;
% % % % % 
% ts_idx = ts_idx(1:300);
% ts_fold_idx = ts_fold_idx(1:300);
% Samplevoted = Samplevoted(1:300);
% NotRatio = NotRatio(1:300,:);
% setting.ts_idx_conf = setting.ts_idx_conf(1:300);

% setting.testSOutput = 0;
% setting.Flatent = 1;

Mstr = setting.Mstr;
addpath(genpath('mlr/'));
PCACOEFF = [];


% 4GB file limit
bytelimit = setting.GBlimit*2^30;


tr_label = zeros(length(tr_idx), 1);ts_label = zeros(length(ts_idx), 1);
clabel = unique(fdatabase{1}.label);
if setting.Mclassfy
    Lscore = 0;
else
    Lscore = Inf * ones(length(ts_idx), length(tr_idx));
end
labelmap = setting.labelmap;
ranktime = 0;
    try
        load(fullfile(setting.Modelresult,['ctr_fea-', Mstr]), 'ctr_fea', 'tr_label');
    catch
        [ctr_fea, tr_label] = GetFeastr(fdatabase, feattype, tr_idx, setting, labelmap);
        save(fullfile(setting.Modelresult,['ctr_fea-', Mstr]), 'ctr_fea', 'tr_label');
    end
    setting.fsize = [];
    if mod(setting.latent, 5)
    jj = 1;fsize = [];
    for i = 1:length(feattype)
        bookfeattype = length(fdatabase{i}.path);
        for j = 1:length(bookfeattype)
            fpath = fdatabase{i}.path{j}{tr_idx(jj)};load(fpath, 'fea', 'label');
            fsize = [fsize; size(fea{1}{1})];
        end
    end
    setting.fsize = [fsize(1,1:2), sum(fsize(:,3))];
    end
    
    Ntrain = size(ctr_fea,1);
    if ~setting.isKNNlatent
        KNNlatent = Ntrain - 1;
    else
        KNNlatent = setting.KNNlatent(1);
    end
    if KNNlatent < 1
        if KNNlatent > 0
            KNNlatent = ceil(abs(KNNlatent) * length(labelmap));
        else
            setting.rankKNN = 1;
            KNNlatent = abs(KNNlatent);
        end
    end
    

    setting.tr_label = tr_label;
    tr_imname = GetTrname(tr_idx, fdatabase);
    bitsize1 = size(ctr_fea, 1) * size(ctr_fea, 1) * 8;
    bitsize = numel(ctr_fea) * 8; 
    if setting.Mclassfy
        try
            if setting.PCA > 0
                load(fullfile(setting.Modelresult,['Model-', Mstr]), 'Metric', 'PCACOEFF');
            else
                load(fullfile(setting.Modelresult,['Model-', Mstr]), 'Metric');
            end
        catch
            Metric = eye(prod(setting.fsize));
            
            setting.datastr = fullfile(setting.latentresult,'XCdata');
            if ~exist(setting.datastr)
                mkdir(setting.datastr)
            end
            if setting.PCA > 0
                setting.datastr = fullfile(setting.latentresult,'XCdata_PCA');
                if ~exist(setting.datastr)
                    mkdir(setting.datastr)
                end
            end
            pcastr = '';
            if isempty(setting.DRstr)
    if setting.PCAenergy ~= -0.9
        pcastr = [pcastr '_E' num2str(setting.PCA)];
    end
            else
    if setting.PCAenergy ~= -0.9
        pcastr = [pcastr '_E' setting.DRstr];
    end
            end
            if setting.LocalPCA
    pcastr = [pcastr, '_L'];
            end
            if ~isempty(setting.RRatiostr)
    pcastr = [setting.RRatiostr, pcastr];
            end 
            
            th = tic;

            [dataX, dataY, numSign, CIndex, latentinfo, ts_imname] = GetSampleLatent(1, ts_fold_idx, ts_idx, tr_idx, Samplevoted, NotRatio, ...
                ctr_fea, Metric, feattype, fdatabase, setting, Lscore, bytelimit, bitsize, bitsize1);           
            

            th = toc(th);
            fprintf('Time for Get latent Samples : %d \n', th)
            Ntrain = size(ctr_fea, 1);           
            numD = size(CIndex, 1);
            
            th = tic;
            if setting.PCA > 0 %%%Computing traininG PCA
                fname = [setting.TFstr, '_', setting.strfea, setting.featNormstr, pcastr,'_PCACOEFF'];
                setting1.NormFea = setting.NormFea;setting1.fname = fname;
                setting1.featsize = setting.featsize;setting1.LocalPCA = setting.LocalPCA;
            
                [PCACOEFF{1}, PCACOEFF{2}, dataX] = GetPCAfea(setting1, [ctr_fea, tr_label], ...
                    setting.PCA, [cell2mat(dataX), zeros(size(cell2mat(dataX), 1), 1)]);
                setting.PCACOEFF = PCACOEFF;
            end
            
            if setting.PCA < 0
                TrainPCA = (ctr_fea - repmat(setting.PCACOEFF{2}, ...
                    [size(ctr_fea,1),1])) * setting.PCACOEFF{1};
                dataX = [TrainPCA;cell2mat(dataX)];
            end
            
            if setting.PCA > 0
                bitsize = size(ctr_fea,1)*size(dataX, 2) * 8; 
                CIndex = GetPCADATA(setting, numD, numSign, CIndex, bytelimit, bitsize, bitsize1);
            end
            
            dataX = mat2cell(dataX, [ones(1, size(dataX,1))], size(dataX, 2));
            
            th = toc(th);
            fprintf('Time for PCA : %d \n', th)
            
            th = tic;
            setting.numSign = [ones(1, size(ctr_fea,1)), numSign];
            Tlatentinfo = latentinfo;
            
            latentfile = fullfile(setting.latentresult,['latent', Mstr(1:end-4), setting.svotestr], 'Train');
            
            dataX([1:Ntrain]) = ChangeRotateOrd(dataX([1:Ntrain]), length(setting.rotate));
            th = tic;
            dataX1= dataX;
            
            if setting.selfpaced(1) == 1
                load(fullfile(setting.ModelresultO,['Model-', Mstr]), 'Metric');
                Winit = Metric;
            else
                Winit= [];
            end
                
            if size(CIndex, 1) == 1
                dataX(Ntrain+1:end) = GetCellX(setting.datastr, CIndex);
                dataX1 = dataX;
            else
            end

            [Metric, Xi, D] = TrainFun(setting.selfpaced, Winit, dataX1, [cell(Ntrain, 2); dataY], ...
                    setting.datastr, setting.numSign, CIndex, Tlatentinfo, latentfile, setting.rotate, setting.rawsize, ...
                    KNNlatent, setting.SnippetRatio, ts_imname, tr_imname, setting.Asubname, setting.C,  setting.LOSS,  setting.k,  setting.REG);
%                

            
            ranktime = toc(th);
            fprintf('Time for metric learning : %d \n', th)
            if setting.PCA > 0
                save(fullfile(setting.Modelresult,['Model-', Mstr]), 'Metric', 'PCACOEFF');
            else
                save(fullfile(setting.Modelresult,['Model-', Mstr]), 'Metric');
            end
            PlotFIG(D.ObjFun, fullfile(setting.Modelresult,['Iteration-', Mstr(1:end-4)]), 1);

            rmdir((fullfile(setting.latentresult,'XCdata')), 's');
            if setting.PCA > 0
                rmdir((fullfile(setting.latentresult,'XCdata_PCA')), 's');
            end
            
        end
    end
    if ~setting.Mclassfy   %%%%testing
        try
            if ~isempty(TrainFun)  %%%if use metric
                if setting.PCA > 0
                    load(fullfile(setting.Modelresult,['Model-', Mstr]), 'Metric', 'PCACOEFF');
                else
                    load(fullfile(setting.Modelresult,['Model-', Mstr]), 'Metric');
                    PCACOEFF = setting.PCACOEFF;
                end
            else
                Metric = [];
                PCACOEFF = setting.PCACOEFF;
            end
            load(fullfile(setting.latentresult,['latent', Mstr(1:end-4), setting.svotestr, '.mat']), 'LX', 'LY', 'LL', 'LR', 'Lscore', 'ranktime');
        catch
            setting.PCACOEFF = PCACOEFF;
            if setting.PCA
                ctr_fea = (ctr_fea - repmat(setting.PCACOEFF{2}, [size(ctr_fea,1),1])) * setting.PCACOEFF{1};
            end
            th = tic;
            
            if ~isempty(TrainFun)
                try
                    load(fullfile(setting.Modelresult,['Model-', Mstr]), 'Metric');
                catch
                    fprintf('Error Load for Model\n');
                    pause;
                end
            else
                Metric = [];
            end
            
            latentfile = fullfile(setting.latentresult,['latent', Mstr(1:end-4), setting.svotestr], 'Test');
            
            if setting.Flatent
                sizet = size(ctr_fea);Nrotate = length(setting.rotate);
                ctr_fea = permute(reshape(ctr_fea, [Nrotate, size(ctr_fea, 1) / Nrotate,...
                    size(ctr_fea, 2)]), [2 1 3]);
                ctr_fea = reshape(ctr_fea, sizet);
            end
            
            [X, Y, numSign, CIndex, latentinfo, ts_imname, LX, LY, LL, LR, Lscore, Lmargin, SignBased] = GetSampleLatent(0, ts_fold_idx, ts_idx, tr_idx,...
                Samplevoted, NotRatio, ctr_fea, Metric, feattype, fdatabase, setting, Lscore, bytelimit, bitsize, bitsize1);
            ranktime = toc(th);
            
            save(fullfile(setting.latentresult,['latent', Mstr(1:end-4), setting.svotestr, '.mat']), 'LX', 'LY', 'LL', 'LR', 'Lscore', 'Lmargin', 'ranktime');      
            ffile = fullfile(setting.latentresult,['latent', Mstr(1:end-4), setting.svotestr], 'Test');
            
%             if SignBased
%                 ShowLatent(ffile, latentinfo, ts_imname, tr_imname, setting.Asubname, setting.rotate, setting.rawsize, 2);
%             else
%                 ShowLatent(ffile, latentinfo, ts_imname, tr_imname, setting.Asubname, setting.rotate, setting.rawsize, 3);
%             end
            
        end
    end
AllScore = Lscore;

function feat = Combine(feat, fea, latent)
for i = 1:length(feat)
    for j = 1:length(fea{i})
        if latent == 1
    feat{i}{j}(:,:, end+1:end+size(fea{i}{j}, 3)) = fea{i}{j};
        end
        if latent == 2
    feat{i}{j}(end+1:end+length(fea{i}{j}),:) = fea{i}{j};
        end
    end
end

function Map = ComputeOrgDistance(ts_idx, ctr_fea, feattype, fdatabase, setting)
Map = cell(length(ts_idx), length(setting.displace));
for jj = 1:length(ts_idx),  %%%compute for each sample
    fprintf('Compute latent position and rotation %d / %d\n', jj, length(ts_idx))
    for i = 1:length(feattype)
        bookfeattype = length(fdatabase{i}.path);
        for j = 1:length(bookfeattype)
            fpath = fdatabase{i}.path{j}{ts_idx(jj)};
            load(fpath, 'fea', 'label');
            if i == 1 && j== 1
                feat = setting.normweigh(i)*fea{1}; 
            else
                feat = Combine(feat, setting.normweigh(i)*fea{1}, setting.latent);
            end
            if setting.latent == 2
                if i == 1 && j== 1
                    feat = setting.normweigh(i)*fea{2};
                else
                    feat = Combine(feat, setting.normweigh(i)*fea{2}, setting.latent);
                end
            end
        end
    end
    for hh = 1:length(setting.displace) %%%for each position, if no dispard , then 0
        Map{jj,hh} = ComputeDistance(feat{hh}, ctr_fea, setting.fsize);             
    end
end


function [X, Y, numSign, CIndex, latentinfo, ts_imname, LX, LY, LL, LR, Lscore, Lmargin, SignBased] = GetSampleLatent(islatent, ts_fold_idx, ts_idx, tr_idx, ...
    Samplevoted, NotRatio, ctr_fea, W, feattype, fdatabase, setting,Lscore, bytelimit, bitsize, bitsize1)
if (setting.bestS || setting.KNNlatentTrain) && ~setting.latent
    [X, Y, numSign, CIndex, latentinfo, ts_imname, LX, LY, LL, LR, Lscore, Lmargin, SignBased] = GetSampleLatent_bestS(islatent, ts_fold_idx, ts_idx, tr_idx, ...
        Samplevoted, NotRatio, ctr_fea, W, feattype, fdatabase, setting,Lscore, bytelimit, bitsize, bitsize1);
    return;
end


numSign = [];
latentinfo = [];
ts_imname = {};
CIndex = [];
if islatent
    LX =0;LY = 0;LL = 0;LR = 0;
else
    
LX = zeros(length(ts_idx), length(tr_idx));
LY = zeros(length(ts_idx), length(tr_idx));
LL = zeros(length(ts_idx), length(tr_idx));
LR = zeros(length(ts_idx), length(tr_idx));
end
X = [];Y = [];
Range = setting.ts_idx_conf;
if isfield(setting, 'ts_fold_idx')
    ts_Fold_label = setting.ts_fold_idx(Range);
else
    ts_Fold_label = ts_fold_idx;
end
labelmap = setting.labelmap;
SignBased = 1;
if setting.bestS
    if islatent || (~islatent && setting.testSOutput)
        [idd, b, c] = unique(ts_Fold_label);
        IDMap = c;
        MaxFlod =(max(unique(c))); 
        numSign = zeros(1, MaxFlod);
    else
        MaxFlod = length(ts_idx);
        IDMap = [1:length(ts_idx)];
        numSign = zeros(1, MaxFlod);
        SignBased = 0;
    end
else
    MaxFlod = length(ts_idx);
    IDMap = [1:length(ts_idx)];
    numSign = zeros(1, MaxFlod);
    SignBased = 0;
end
if islatent
%     Samplevoted(:) = 1;   %%%let the 
    NotRatio(:) = 1;
    X = cell(length(ts_idx), 1);
    Y = cell(length(ts_idx), 2);
    latentinfo = cell(length(ts_idx), 1);
    if setting.bestS
        X = cell(MaxFlod, 1);
        Y = cell(MaxFlod, 2); 
        latentinfo = cell(MaxFlod, 1);
    end
end
ts_imname = cell(1, length(ts_idx));
Lmargin = zeros(1, length(ts_idx));
iused = 0;
round = 0;istart = 1;
orglength = sqrt(bitsize1 / 8)';
orglength1 = orglength;
usedbit = bitsize;
usedbit1 = orglength*orglength*8;
for jj = 1:length(ts_idx),
    if Samplevoted(jj) == 0
        fprintf('Not voted snippets %d\n', jj)
        continue;
    end
    iused  = iused+1;
    aindex = find(NotRatio(jj,:) ~= 0);
    fprintf('Compute latent position and rotation %d / %d\n', jj, length(ts_idx))
    for i = 1:length(feattype)
        bookfeattype = length(fdatabase{i}.path);
        for j = 1:length(bookfeattype)
            fpath = fdatabase{i}.path{j}{ts_idx(jj)};
            load(fpath, 'fea', 'label');
            ts_imname{iused} = fdatabase{i}.imgpath{ts_idx(jj)};
            if setting.latent == 2
                if i == 1 && j== 1
% %                     feat = setting.normweigh(i)*fea{2};
                    feat = fea{2};
                    for jjj = 1:length(feat)
                        for kkk = 1:length(feat{jjj})
                            feat{jjj}{kkk} = setting.normweigh(i)*feat{jjj}{kkk};
                        end
                    end
                else
%                     feat = Combine(feat, setting.normweigh(i)*fea{2});
                    for jjj = 1:length(fea{2})
                        for kkk = 1:length(fea{2}{jjj})
                            fea{2}{jjj}{kkk} = setting.normweigh(i)*fea{2}{jjj}{kkk};
                        end
                    end
                    feat = Combine(feat, fea{2}, setting.latent);
                end
            else
             tmp = fea{1};
            if i == 1 && j== 1
                feat = setting.normweigh(i)*tmp; 
            else
                feat = Combine(feat, setting.normweigh(i)*tmp, setting.latent);
            end
            end
        end
    end
    aindex1 = find(NotRatio(jj,setting.ConsiderID) ~= 0);
    nrotate =  length(setting.rotate);
    if islatent   %%trining
        if setting.latent == 1
            for li = 1:length(setting.displace) %%%for each position, if no dispard , then 0
            xdim = size(feat{li}, 1) - setting.fsize(1);
            ydim = size(feat{li}, 2) - setting.fsize(2);
            for xi = 1:xdim+1
                for yi = 1:ydim+1
                    latfeat = feat{li}(xi:xi+setting.fsize(1)-1, yi:yi+setting.fsize(2)-1,:);
                    numSign(iused) = numSign(iused) + 1;
                    X{iused}(numSign(iused),:) =  reshape(latfeat, [1, numel(latfeat)]); 
                    index = setting.MapSub2ind{li}(xi, yi);
                    latentinfo{iused}(numSign(iused),:) = [setting.MapRange(index, :),jj,ts_Fold_label(jj)];
                end
            end
            end
        else
            for li = 1:length(feat) %%%for each position, if no dispard , then 0
                for xi = 1:length(feat{li})
                    latfeat = feat{li}{xi};
                    numSign(iused) = numSign(iused) + 1;
                    X{iused}(numSign(iused),:) =  reshape(latfeat, [1, numel(latfeat)]); 
                    index = setting.MapSub2ind{li}(setting.Map{li}(xi,1),setting.Map{li}(xi,2));
                    latentinfo{iused}(numSign(iused),:) = [setting.MapRange{li}(index, :),jj,ts_Fold_label(jj)];
                end
            end     
        end
        
        idx = find(setting.tr_label == labelmap(label));
        Y{iused, 1} = ceil(idx(end)/nrotate);
        Y{iused, 2} = setdiff([1:length(setting.tr_label)/nrotate], Y{iused, 1});

        if setting.PCA < 0
            X{(iused)} = (X{(iused)} - repmat(setting.PCACOEFF{2}, ...
                [size(X{(iused)},1),1])) * setting.PCACOEFF{1};
        end
        
        Xtmp = (X{(iused)});
        
        
        thislength = size(Xtmp, 1);
        bitthis = numel(Xtmp)*8;
        bitthis1 = 2*(thislength*orglength)*8+(2*orglength1+1)*8;
        EndRound = 0;
        if (usedbit + bitthis + usedbit1 + bitthis1 > bytelimit)
            EndRound = 1;
            round = round + 1;iend = iused-1;
            ffile = fullfile(setting.latentresult,'XCdata', [num2str(round), '.mat']);
            Xcell = X(istart:iend);
            save(ffile, 'Xcell', '-v7.3');clear 'Xcell'
            CIndex = [CIndex; [istart, iend]];
            X(istart:iend) = cell(iend-istart+1,1);
            istart = iend+1;
            usedbit = bitsize;
            usedbit1 = orglength*orglength*8;
            
        end
        usedbit = usedbit + bitthis;
        usedbit1 = usedbit1 + bitthis1;
        orglength1 = orglength1 + 1;
        XX{iused} = Xtmp(1,:); %%%%%for the first one
    else
        score = zeros(length(feat), length(aindex));
        xpos = zeros(length(feat), length(aindex));
        ypos = zeros(length(feat), length(aindex));
        rpos = zeros(length(feat), length(aindex));
        Index = intersect([1:length(labelmap)], aindex1);
        
        if setting.Flatent
            Tindex = repmat(Index(:), [1, nrotate]) + repmat([0:nrotate-1]*length(labelmap), [length(Index), 1]); 
            Tindex = Tindex(:);
        else
            TindexY = repmat(Index, [nrotate, 1]);
            TindexX = repmat([1:nrotate]', [1, length(Index)]);
            Tindex = sub2ind([nrotate, length(labelmap)], TindexX, TindexY);
            Tindex = Tindex(:);
        end
        
        ComputeID = setting.ConsiderID(aindex1);
        
        margin = zeros(1, length(setting.displace));
        for hh = 1:length(setting.displace) %%%for each position, if no dispard , then 0
            if setting.Flatent
                [rscore, rxpos, rypos, rrpos, margin(hh)] = ...
                    bestmatchingW_MRR(feat{hh}, ctr_fea(Tindex, :), setting.Map{hh}, ...
                    setting.PCACOEFF, setting.fsize, W, nrotate);
            else
                [rscore, rxpos, rypos, rrpos] = ...
                    bestmatchingW(feat{hh}, ctr_fea(Tindex, :), setting.Map{hh}, ...
                    setting.PCACOEFF, setting.fsize, W, nrotate);
            end
            
            tmp = Inf * ones(1, length(labelmap));
            tmp(ComputeID) = rscore; score(hh,:) = tmp(aindex); 
            tmp = Inf * ones(1, length(labelmap));
            tmp(ComputeID) = rxpos; xpos(hh,:) = tmp(aindex); 
            tmp = Inf * ones(1, length(labelmap));
            tmp(ComputeID) = rypos; ypos(hh,:) = tmp(aindex); 
            tmp = Inf * ones(1, length(labelmap));
            tmp(ComputeID) = rrpos; rpos(hh,:) = tmp(aindex);  
        end
        score(length(setting.displace)+1:end,:) = Inf;
        
        if setting.Flatent
            [Lmargin(jj), ord] = max(margin);
            index = sub2ind(size(score), ord*ones(size(score,2), 1), [1:size(score,2)]');
            score = score(ord, :);
            LL(jj,aindex) = ord;Ltmp = ord.*ones(length(index),1);
        else
            [score, ord] = min(score, [], 1);
            index = sub2ind(size(xpos), ord(:), [1:length(ord)]');
            LL(jj,aindex) = ord;Ltmp = ord;
        end

        LX(jj,aindex) = xpos(index);Xtmp = xpos(index);
        LY(jj,aindex) = ypos(index);Ytmp = ypos(index);
        
        LR(jj,aindex) = rpos(index);Rtmp = rpos(index);
        
        Lscore(jj,aindex) = score;
        [ss, ord] = min(score);
        sindex = aindex(ord(1));
        
        index = setting.MapSub2ind{Ltmp(ord(1))}(Xtmp(ord(1)), Ytmp(ord(1)));
        latentinfo(jj,:) = [setting.MapRange{Ltmp(ord(1))}(index, :),jj,...
            ts_Fold_label(jj), Rtmp(ord(1)), sindex];
    end
end

if islatent
    numSign = numSign(1:iused);
    latentinfo = latentinfo(1:iused);
    Y =Y(1:iused,:);
    ts_imname = ts_imname(1:iused);
    
    round = round + 1;iend = iused;
    ffile = fullfile(setting.latentresult,'XCdata', [num2str(round), '.mat']);
    Xcell = X(istart:iend);
    save(ffile, 'Xcell', '-v7.3')
    CIndex = [CIndex; [istart, iend]];
    clear 'X'
    X = XX';
end
    
nbase = size(LL, 2);
if SignBased && ~islatent
    for tt = 1:MaxFlod
        index = find(IDMap == tt);
        LXT = LX(index, :);
        LYT = LY(index, :);
        LRT = LR(index, :);
        LLT = LL(index, :);
        LscoreT = Lscore(index, :);
        
        if setting.Flatent
            [marginB, ord] = max(Lmargin(index));
            ind = sub2ind(size(LscoreT), ord * ones(1,nbase), [1:nbase]);
            ss = LscoreT(ord,:);
        else
            [ss, ord] = min(LscoreT, [], 1);
            ind = sub2ind(size(LscoreT), ord, [1:nbase]);
        end
        
        nn = size(LXT, 1);
        
        LX(index,:) = repmat(reshape(LXT(ind), [1, nbase]), [nn, 1]);
        LY(index,:) = repmat(reshape(LYT(ind), [1, nbase]), [nn, 1]);
        LL(index,:) = repmat(reshape(LLT(ind), [1, nbase]), [nn, 1]);
        LR(index,:) = repmat(reshape(LRT(ind), [1, nbase]), [nn, 1]);
        Lscore(index,:) = repmat(reshape(LscoreT(ind), [1, nbase]), [nn, 1]);
        
        [ss, ord1] = min(ss); %%%%
        if setting.Flatent
            row = ord;
        else
            row = ord(ord1);
        end
        xi = LXT(row, ord1);
        yi = LYT(row, ord1);
        li = LLT(row, ord1);
        ri = LRT(row, ord1);
        tindex = setting.MapSub2ind{li}(xi, yi);
        jjj = ts_Fold_label(index);
        latentinfo(index,:) = [repmat([setting.MapRange{li}(tindex, :), index(row)], [nn,1]), jjj(:), repmat([ri, ord1], [nn, 1])];
    end
end


function CIndex = GetPCADATA(setting, numD, numSign, CIndex, bytelimit, bitsize, bitsize1)
round = 0;
XcellT = {};CIndexT = [];
istart = 1;

usedbit = bitsize;
% usedbit1 = 0;
orglength = sqrt(bitsize1 / 8)';
usedbit1 = orglength*orglength*8;
orglength1 = orglength;


for jj = 1:numD
    tsize = numSign(CIndex(jj, 1):CIndex(jj, 2));
    ffile1 = fullfile(setting.latentresult,'XCdata', [num2str(jj), '.mat']);
    load(ffile1, 'Xcell')
    Xtmp = cell2mat(Xcell);
    Xtmp = (Xtmp - repmat(setting.PCACOEFF{2}, [size(Xtmp,1),1])) * setting.PCACOEFF{1};

    EndRound = 0;
    Xcell = mat2cell(Xtmp, tsize, size(Xtmp, 2));
    XXcell = Xcell;
    
    thislength = size(Xtmp, 1);
    bitthis = numel(Xtmp)*8;
    bitthis1 = 2*(thislength*orglength)*8 + (2*orglength1+1)*8;   
    
%     if (usedbit + bitthis > bytelimit) || (usedbit1 + bitthis1 > bytelimit)
    if (usedbit + bitthis + usedbit1 + bitthis1 > bytelimit)
        
        EndRound = 1;
        round = round + 1;
        ffile = fullfile(setting.latentresult,'XCdata_PCA', [num2str(round), '.mat']);
        iend = jj - 1;
        Xcell = XcellT';save(ffile, 'Xcell', '-v7.3')
%         usedbit = 0;
        XcellT = {};
        CIndexT(end+1,:) = [CIndex(istart,1), CIndex(iend,2)];
        istart = iend+1;
        
        usedbit = bitsize;
        usedbit1 = orglength*orglength*8;
        
    end
%     usedbit = usedbit + bitthis;
        usedbit = usedbit + bitthis;
        usedbit1 = usedbit1 + bitthis1;
        orglength1 = orglength1 + length(Xcell);
        

    XcellT(end+1:end+length(XXcell)) = XXcell; 
    
end

if (~EndRound) || (EndRound && jj == numD)
    round = round + 1;iend = jj;
    ffile = fullfile(setting.latentresult,'XCdata_PCA', [num2str(round), '.mat']);
    Xcell = XcellT';save(ffile, 'Xcell', '-v7.3')
    CIndexT(end+1,:) = [CIndex(istart,1), CIndex(iend,2)];
end
CIndex = CIndexT;


function Map = ComputeDistance(feat, ctr_fea, fsize)
[m, n] = size(ctr_fea);
xdim = size(feat, 1) - fsize(1);
ydim = size(feat, 2) - fsize(2);
Map = cell(xdim+1, ydim+1);
for i = 1:xdim+1
    for j = 1:ydim+1
        Map{i,j} = cell(1, m);
        feattmp = feat(i:i+fsize(1)-1,j:j+fsize(2)-1,:);
        feattmp = feattmp(:);
        for t = 1:m
            dis = feattmp' - ctr_fea(t,:);
            Map{i,j}{t} = dis'*dis;      
        end
    end
end

function [ctr_fea, tr_label] = GetFeastr(fdatabase, feattype, tr_idx, setting, labelmap)
ctr_fea = []; 
for i = 1:length(feattype)
            bookfeattype = length(fdatabase{i}.path);
            for j = 1:length(bookfeattype)
                ctr_fea_t = []; tr_label = [];
                for jj = 1:length(tr_idx),
                    if ~mod(jj, 5),
                        fprintf('.');
                    end
                    if ~mod(jj, 100),
                        fprintf(' %d images processed\n', jj);
                    end
                    fpath = fdatabase{i}.path{j}{tr_idx(jj)};
                    load(fpath, 'fea', 'label');
                    if ~ismember(labelmap(label), setting.Consider)
                        continue;
                    end
                    if mod(setting.latent, 5)
                        for tt = 1:length(setting.rotate)   %%%if no rotate then for one
                            tmp = fea{1}{tt}(:);
                            if setting.latent == 2
                                tmp = fea{2}{tt}(:);
                            end
                            ctr_fea_t = [ctr_fea_t; setting.normweigh(i)*tmp'];
                            tr_label = [tr_label; labelmap(label)];
                        end
                    else
                        tmp = fea;
                        ctr_fea_t = [ctr_fea_t; setting.normweigh(i)*tmp'];
                        tr_label = [tr_label; labelmap(label)];
                    end
                    
                end
                ctr_fea = [ctr_fea,ctr_fea_t]; 
             end
        end
        clear 'ctr_fea_t'

function tr_imname = GetTrname(tr_idx, fdatabase)
tr_imname = cell(length(tr_idx), 1);
i = 1;
for jj = 1:length(tr_idx),
    tr_imname{jj} = fdatabase{i}.imgpath{tr_idx(jj)};
end

function dataX = GetCellX(datastr, CIndex)
BatchRound = size(CIndex, 1);
dataX = cell(CIndex(end, 2), 1);
for R = 1:BatchRound
    index = [CIndex(R, 1):CIndex(R, 2)];
    load(fullfile(datastr,[num2str(R) '.mat']), 'Xcell'); 
    dataX(index) = Xcell;
end
            
