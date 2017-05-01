function Ymatched = GetMatched(Ypos, SAMPLES, NT)
Ymatched = cell2mat((Ypos(SAMPLES))');
tnum = cell2mat(cellfun(@size, Ypos(SAMPLES), 'UniformOutput', false));
        tnum = tnum(:, 2);
        cnum = mat2cell([1:length(tnum)]', ones(1, length(tnum)), 1);
        ctnum = mat2cell(tnum, ones(1, length(tnum)), 1);
        TT = cell2mat((cellfun(@(x,y) y*ones(1, x), ctnum, cnum, 'UniformOutput', false))');
        Ymatched = sub2ind([length(SAMPLES), NT], TT(:), Ymatched(:));
        Ymatched = sort(Ymatched);