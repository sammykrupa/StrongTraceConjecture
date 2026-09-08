function [pivotCols, pivotColumns] = symRREFPivotCols(R)
    % symRREFPivotCols returns pivot column indices and columns
    % for a symbolic matrix R that is already in RREF.

    R = sym(R);
    [m, n] = size(R);

    pivotCols = [];

    for i = 1:m
        for j = 1:n
            entry = simplify(R(i,j));

            % Entry is treated as nonzero unless it is provably zero
            if ~isAlways(entry == 0, 'Unknown', 'false')
                pivotCols(end+1) = j;
                break
            end
        end
    end

    % In case of redundant detections, preserve order
    pivotCols = unique(pivotCols, 'stable');

    % The actual pivot columns of R
    pivotColumns = R(:, pivotCols);
end