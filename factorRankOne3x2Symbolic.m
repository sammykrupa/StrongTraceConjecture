function [aCell, bCell] = factorRankOne3x2Symbolic(C)
% factorRankOne3x2Symbolic
%
% Factors each symbolic rank-one 3x2 matrix C(:,:,i) as
%
%       C(:,:,i) = aCell{i} * bCell{i}.'
%
% where
%
%       aCell{i} is 3x1
%       bCell{i} is 2x1
%
% Input:
%       C : 3x2xN symbolic array
%
% Output:
%       aCell : 1xN cell array, each entry is 3x1 symbolic vector
%       bCell : 1xN cell array, each entry is 2x1 symbolic vector
%
% Notes:
%       The factorization of a rank-one matrix is not unique.
%       If C = a*b.', then also C = (lambda*a)*(b/lambda).'
%       for any nonzero scalar lambda.
%
%       This code chooses a pivot entry C(p,q) that is not identically zero.
%       It then sets
%
%           a = C(:,q)
%           b = C(p,:).' / C(p,q)
%
%       For a rank-one matrix, this guarantees
%
%           C = a*b.'
%
%       as long as the chosen pivot is nonzero.

    % Make sure input is symbolic.
    C = sym(C);

    % Check that the first two dimensions are 3x2.
    if size(C,1) ~= 3 || size(C,2) ~= 2
        error('Input C must have size 3x2xN.');
    end

    % Number of matrices.
    N = size(C,3);

    % Preallocate output cell arrays.
    aCell = cell(1,N);
    bCell = cell(1,N);

    % Process each 3x2 matrix C(:,:,i).
    for i = 1:N

        % Extract the i-th 3x2 symbolic matrix.
        M = simplify(C(:,:,i));

        % Find a symbolic pivot entry M(p,q) that is not identically zero.
        [p, q] = findSymbolicPivot(M);

        % If no nonzero pivot was found, then M is the zero matrix.
        % The zero matrix can be written as zeros(3,1)*zeros(2,1).'
        if isempty(p)
            aCell{i} = sym(zeros(3,1));
            bCell{i} = sym(zeros(2,1));
            continue;
        end

        % Store the pivot value.
        pivot = M(p,q);

        % For a rank-one matrix, choosing
        %
        %       a = M(:,q)
        %       b = M(p,:).' / M(p,q)
        %
        % gives M = a*b.'.
        a = simplify(M(:,q));
        b = simplify(M(p,:).' / pivot);

        % Store the result.
        aCell{i} = a;
        bCell{i} = b;

        % Optional verification.
        % This checks whether a*b.' - M simplifies to the zero matrix.
        D = simplify(a*b.' - M);

        verified = true;
        for k = 1:numel(D)
            if ~isIdenticallyZero(D(k))
                verified = false;
                break;
            end
        end

        if ~verified
            warning(['Matrix C(:,:,%d) was not verified as rank one, ', ...
                     'or simplification was inconclusive.'], i);
        end
    end
end


function [p, q] = findSymbolicPivot(M)
% findSymbolicPivot
%
% Finds the first entry of M that is not identically zero.
%
% Output:
%       p, q : row and column of pivot
%
% If M is identically zero, returns p = [] and q = [].

    p = [];
    q = [];

    for col = 1:size(M,2)
        for row = 1:size(M,1)

            if ~isIdenticallyZero(M(row,col))
                p = row;
                q = col;
                return;
            end

        end
    end
end


function tf = isIdenticallyZero(expr)
% isIdenticallyZero
%
% Returns true if expr is symbolically equal to zero.
% Returns false otherwise.
%
% For symbolic expressions containing parameters, this checks whether the
% expression is identically zero, not whether it can become zero for some
% parameter values.

    expr = simplify(expr);

    % Direct symbolic zero test.
    if isequal(expr, sym(0))
        tf = true;
        return;
    end

    % Additional symbolic test.
    % If MATLAB cannot prove the expression is always zero, isAlways returns
    % false in most versions.
    try
        tf = isAlways(expr == 0);
    catch
        tf = false;
    end
end