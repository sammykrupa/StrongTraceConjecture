function out = tn_tangent_basis_symbolic(P, aCell, bCell, kappa)
%TN_TANGENT_BASIS_SYMBOLIC Symbolic basis for the tangent space to the T_N manifold.
%
%   out = tn_tangent_basis_symbolic(P, aCell, bCell, kappa)
%
% Given a T_N parametrization in (R^{3x2})^N,
%
%   X_i = P + C_1 + ... + C_{i-1} + kappa_i * C_i,
%   C_i = a_i * b_i.',
%   sum_i C_i = 0,
%
% this function constructs a symbolic basis of the tangent space to the
% local T_N-manifold at the point (X_1,...,X_N).
%
% -------------------------------------------------------------------------
% INPUTS
% -------------------------------------------------------------------------
%   P     : 3x2 symbolic or numeric matrix.
%
%   aCell : 1xN or Nx1 cell array. Each entry aCell{i} must be a 3x1 vector.
%
%   bCell : 1xN or Nx1 cell array. Each entry bCell{i} must be a 2x1 vector.
%           (Row vectors are also accepted; they are converted to columns.)
%
%   kappa : 1xN or Nx1 symbolic/numeric vector of the parameters kappa_i.
%
% -------------------------------------------------------------------------
% OUTPUT
% -------------------------------------------------------------------------
%   out is a struct with the following fields:
%
%     out.basis       - cell array containing an actual basis of the image of
%                       the linearized T_N parametrization.
%                       out.basis{m} is a 1xN cell array, and out.basis{m}{i}
%                       is the 3x2 matrix H_i of the m-th basis vector.
%
%     out.basisMatrix - symbolic matrix whose columns are vec(out.basis{m}),
%                       obtained by stacking H_1(:), H_2(:), ..., H_N(:).
%
%     out.labels      - descriptive labels for the returned basis vectors.
%
%     out.candidateBasis       - the unreduced family of 6 + N + dim ker(A)
%                                tangent directions coming from the parameter
%                                variations. At a regular point this family
%                                already has size 5N and is linearly independent.
%
%     out.candidateBasisMatrix - the corresponding symbolic matrix whose
%                                columns are the candidate directions.
%
%     out.candidateLabels      - labels for the candidate directions.
%
%     out.basisPivots - pivot columns selected from out.candidateBasisMatrix.
%
%     out.basisRank   - number of returned basis vectors.
%
%     out.expectedDim - expected manifold dimension 5N at a regular point.
%
%     out.X           - 1xN cell array containing the configuration X_i.
%
%     out.C           - 1xN cell array containing C_i = a_i * b_i.'.
%
%     out.E           - N-by-4 cell array.  E{i,alpha} is a chosen basis of
%                       T_{C_i}(rank-one manifold), namely
%                           E{i,1} = e1 \otimes b_i,
%                           E{i,2} = e2 \otimes b_i,
%                           E{i,3} = e3 \otimes b_i,
%                           E{i,4} = a_i \otimes b_i^perp,
%                       with b_i^perp = [-b_i(2); b_i(1)].
%
%     out.A           - 6-by-(4N) symbolic matrix for the linearized closure
%                       map  (dC_1,...,dC_N) -> sum_i dC_i.
%
%     out.R           - reduced row echelon form of out.A.
%
%     out.pivots      - pivot columns of out.A.
%
%     out.freeCols    - free columns of out.A.
%
%     out.kernelCoeff - cell array of coefficient vectors spanning ker(out.A).
%                       Each coefficient vector is expressed in the basis
%                       E{i,alpha}.
%
%     out.colMap      - (4N)-by-2 array.  Row c stores [i, alpha], meaning
%                       column c of out.A corresponds to E{i,alpha}.
%
% -------------------------------------------------------------------------
% MATHEMATICAL BACKGROUND
% -------------------------------------------------------------------------
% An admissible first variation of the parameters is given by
%
%   dP in R^{3x2},
%   dkappa_i in R,
%   dC_i in T_{C_i}(rank-one manifold),
%   sum_i dC_i = 0.
%
% The induced first variation of X_i is
%
%   dX_i = dP + sum_{j < i} dC_j + kappa_i * dC_i + dkappa_i * C_i.
%
% This code builds a basis of the tangent space by combining:
%
%   (1) the 6 translation directions coming from dP,
%   (2) the N directions coming from dkappa_i,
%   (3) a basis of the kernel of the linearized closure map sum_i dC_i = 0.
%
% At a regular point one has rank(A) = 6, so dim ker(A) = 4N - 6, and the
% total tangent-space dimension is
%
%   6 + N + (4N - 6) = 5N.
%
% -------------------------------------------------------------------------
% IMPORTANT NOTE
% -------------------------------------------------------------------------
% This routine assumes you are at a REGULAR T_N point, meaning the matrix A
% below has full row rank 6.  If rank(A) < 6, the usual local 5N-dimensional
% manifold chart is singular at this point, and the routine stops.
%
% Because pivot selection matters, this routine is best used at a concrete
% point (numeric or symbolic constants).  If your input contains abstract
% symbolic parameters, you may need to impose assumptions with ASSUME/ASSUMEALSO
% so that the symbolic row-reduction can choose pivots correctly.
%
% -------------------------------------------------------------------------
% EXAMPLE
% -------------------------------------------------------------------------
%   syms k1 k2 k3 k4 real
%   P = sym(zeros(3,2));
%   a = { [1;0;0], [0;1;0], [-1;0;0], [0;-1;0] };
%   b = { [1;1], [1;-1], [1;1], [1;-1] };
%   out = tn_tangent_basis_symbolic(P, a, b, [k1 k2 k3 k4]);
%
%   % Access the 7th tangent basis vector, which is an N-tuple of 3x2 matrices:
%   H = out.basis{7};
%   H1 = H{1};
%   H2 = H{2};
%
%   % The full basis as columns of one symbolic matrix:
%   B = out.basisMatrix;
%
%   % If you also want the unreduced candidate family before column-selection:
%   Bcand = out.candidateBasisMatrix;
%
% -------------------------------------------------------------------------

    %------------------------------------------------------------------
    % Basic input checks and symbolic conversion
    %------------------------------------------------------------------
    if nargin ~= 4
        error('tn_tangent_basis_symbolic requires exactly 4 inputs: P, aCell, bCell, kappa.');
    end

    if ~iscell(aCell) || ~iscell(bCell)
        error('aCell and bCell must both be cell arrays.');
    end

    N = numel(aCell);
    if N < 1
        error('aCell must contain at least one entry.');
    end

    if numel(bCell) ~= N
        error('aCell and bCell must have the same length N.');
    end

    if numel(kappa) ~= N
        error('kappa must have length N = numel(aCell).');
    end

    P = sym(P);
    if ~isequal(size(P), [3, 2])
        error('P must be a 3x2 matrix.');
    end

    kappa = sym(kappa(:).');  % store kappa as a row vector for convenience

    % Standard basis of R^3.
    e = cell(1,3);
    e{1} = sym([1;0;0]);
    e{2} = sym([0;1;0]);
    e{3} = sym([0;0;1]);

    %------------------------------------------------------------------
    % Build C_i, the local bases E{i,alpha}, and the configuration X_i
    %------------------------------------------------------------------
    C = cell(1,N);
    E = cell(N,4);
    X = cell(1,N);

    % We also keep the original a_i and b_i (converted to symbolic columns),
    % since they are sometimes useful for debugging or post-processing.
    aSym = cell(1,N);
    bSym = cell(1,N);

    sumC = sym(zeros(3,2));

    for i = 1:N
        ai = sym(aCell{i});
        bi = sym(bCell{i});

        ai = ai(:);
        bi = bi(:);

        if numel(ai) ~= 3
            error('aCell{%d} must contain a vector with 3 entries.', i);
        end
        if numel(bi) ~= 2
            error('bCell{%d} must contain a vector with 2 entries.', i);
        end

        aSym{i} = ai;
        bSym{i} = bi;

        % Rank-one increment C_i = a_i \otimes b_i = a_i * b_i.'
        C{i} = simplify(ai * bi.');
        sumC = simplify(sumC + C{i});

        % A convenient 4-element basis of T_{C_i}(rank-one manifold):
        %   E{i,1} = e1 \otimes b_i,
        %   E{i,2} = e2 \otimes b_i,
        %   E{i,3} = e3 \otimes b_i,
        %   E{i,4} = a_i \otimes b_i^perp.
        bPerp = [-bi(2); bi(1)];

        E{i,1} = simplify(e{1} * bi.');
        E{i,2} = simplify(e{2} * bi.');
        E{i,3} = simplify(e{3} * bi.');
        E{i,4} = simplify(ai * bPerp.');
    end

    % Optional sanity check of the closure sum_i C_i = 0.
    % If the input is symbolic with unspecified parameters, isAlways may not be
    % able to decide; in that case we simply skip the warning.
    try
        if ~all(isAlways(sumC(:) == 0))
            warning(['The supplied rank-one increments do not appear to satisfy ', ...
                     'sum_i C_i = 0 identically. Check the input parameterization.']);
        end
    catch
        % Do nothing; undecidable symbolic conditions are common.
    end

    % Reconstruct the T_N configuration X_i from the parameterization.
    runningBase = P;
    for i = 1:N
        X{i} = simplify(runningBase + kappa(i) * C{i});
        runningBase = simplify(runningBase + C{i});
    end

    %------------------------------------------------------------------
    % Build the 6 x (4N) symbolic matrix A encoding the linearized closure
    % constraint   sum_i dC_i = 0.
    %
    % Each admissible dC_i is expanded as
    %   dC_i = sum_{alpha=1}^4 t_{i,alpha} E{i,alpha}.
    %
    % Then the constraint sum_i dC_i = 0 becomes A * t = 0,
    % where t is the column vector of coefficients t_{i,alpha}.
    %------------------------------------------------------------------
    A = sym(zeros(6, 4*N));
    colMap = zeros(4*N, 2);   % row c = [i, alpha]

    col = 0;
    for i = 1:N
        for alpha = 1:4
            col = col + 1;
            A(:, col) = E{i,alpha}(:);   % MATLAB stacks columns, i.e. vec(M) = M(:)
            colMap(col,:) = [i, alpha];
        end
    end

    % Compute the reduced row echelon form and the pivot structure.
    R = rref(A);
    [pivots, ~] = symRREFPivotCols(R)
    rankA = numel(pivots);

    if rankA < 6
        error(['The closure matrix A has rank %d < 6. ', ...
               'This point is not regular for the standard 5N-dimensional ', ...
               'local T_N manifold chart.'], rankA);
    end

    freeCols = setdiff(1:(4*N), pivots, 'stable');

    %------------------------------------------------------------------
    % Build an explicit basis of ker(A) using the RREF data.
    %
    % For each free column f, we set the free variable x_f = 1, the other
    % free variables = 0, and solve for the pivot variables.
    %------------------------------------------------------------------
    kernelCoeff = cell(1, numel(freeCols));

    for m = 1:numel(freeCols)
        f = freeCols(m);

        % z is the coefficient vector in R^{4N} describing a shape variation.
        z = sym(zeros(4*N, 1));
        z(f) = sym(1);

        % In RREF form, the pivot variables satisfy x_pivot = -R(:,free)*x_free.
        % The pivot rows are exactly the first rankA rows of R.
        z(pivots) = -R(1:rankA, f);

        kernelCoeff{m} = simplify(z);
    end

    %------------------------------------------------------------------
    % Assemble a CANDIDATE generating family for the tangent space.
    %
    % Block 1: 6 directions coming from dP.
    % Block 2: N directions coming from dkappa_i.
    % Block 3: (4N - 6) directions coming from the kernel of A.
    %
    % In the generic regular case this family is already a basis of size 5N.
    % To be fully robust, we will later vectorize these directions and perform
    % one more symbolic column-reduction, returning only independent columns.
    %------------------------------------------------------------------
    numCandidate = 6 + N + numel(kernelCoeff);   % equals 5N when rankA = 6
    candidateBasis = cell(1, numCandidate);
    candidateLabels = cell(1, numCandidate);

    %------------------------------
    % Block 1: translation directions (dP)
    %------------------------------
    idx = 0;
    for q = 1:2
        for p = 1:3
            idx = idx + 1;

            dirTuple = zeroTuple(N);
            Epq = sym(zeros(3,2));
            Epq(p,q) = sym(1);

            for i = 1:N
                dirTuple{i} = Epq;
            end

            candidateBasis{idx} = dirTuple;
            candidateLabels{idx} = sprintf('dP(%d,%d)', p, q);
        end
    end

    %------------------------------
    % Block 2: kappa directions (dkappa_i)
    %------------------------------
    for i = 1:N
        idx = idx + 1;

        dirTuple = zeroTuple(N);
        dirTuple{i} = C{i};   % since dX_i = dkappa_i * C_i when all other variations vanish

        candidateBasis{idx} = dirTuple;
        candidateLabels{idx} = sprintf('dkappa_%d', i);
    end

    %------------------------------
    % Block 3: shape directions from ker(A)
    %------------------------------
    for m = 1:numel(kernelCoeff)
        idx = idx + 1;
        z = kernelCoeff{m};

        % Recover the corresponding dC_i from the coefficient vector z.
        dC = zeroTuple(N);
        for c = 1:(4*N)
            coeff = z(c);
            if ~isAlwaysZero(coeff)
                i = colMap(c, 1);
                alpha = colMap(c, 2);
                dC{i} = simplify(dC{i} + coeff * E{i,alpha});
            end
        end

        % Push the variation forward to the configuration space:
        %   dX_i = sum_{j < i} dC_j + kappa_i * dC_i
        % since here dP = 0 and dkappa_i = 0.
        dirTuple = zeroTuple(N);
        cumulative = sym(zeros(3,2));

        for i = 1:N
            dirTuple{i} = simplify(cumulative + kappa(i) * dC{i});
            cumulative = simplify(cumulative + dC{i});
        end

        candidateBasis{idx} = dirTuple;
        candidateLabels{idx} = sprintf('shape_%d_freeCol_%d', m, freeCols(m));
    end

    %------------------------------------------------------------------
    % Pack the candidate family into one big symbolic matrix.
    % Each column is vec(H_1,...,H_N).
    %------------------------------------------------------------------
    candidateBasisMatrix = sym(zeros(6*N, numCandidate));
    for m = 1:numCandidate
        candidateBasisMatrix(:,m) = tupleToVector(candidateBasis{m});
    end

    %------------------------------------------------------------------
    % Select an actual basis by symbolic column reduction.
    %
    % If the point is a regular point of the T_N manifold, then the rank here
    % should be exactly 5N. In that case every candidate column survives.
    % If the rank drops, we still return a basis of the image of the linearized
    % parametrization, together with a warning.
    %------------------------------------------------------------------
    RB = rref(candidateBasisMatrix);  %#ok<ASGLU>
    [basisPivots, ~]=symRREFPivotCols(RB)
    basisRank = numel(basisPivots);

    expectedDim = 5*N;

    if basisRank < numCandidate
        warning(['The pushed-forward candidate directions are not all independent. ', ...
                 'Returning a basis with %d columns instead of %d.'], ...
                 basisRank, numCandidate);
    end

    if basisRank ~= expectedDim
        warning(['The returned tangent-space basis has dimension %d, while the regular ', ...
                 'T_N manifold dimension would be %d. Check whether the point is singular ', ...
                 'or whether additional symbolic assumptions are needed.'], ...
                 basisRank, expectedDim);
    end

    basis = candidateBasis(basisPivots);
    labels = candidateLabels(basisPivots);
    basisMatrix = candidateBasisMatrix(:, basisPivots);

    %------------------------------------------------------------------
    % Return everything in a struct.
    %------------------------------------------------------------------
    out = struct();
    out.basis       = basis;
    out.basisMatrix = basisMatrix;
    out.labels      = labels;
    out.candidateBasis = candidateBasis;
    out.candidateBasisMatrix = candidateBasisMatrix;
    out.candidateLabels = candidateLabels;
    out.basisPivots = basisPivots;
    out.basisRank   = basisRank;
    out.expectedDim = expectedDim;
    out.X           = X;
    out.C           = C;
    out.E           = E;
    out.A           = A;
    out.R           = R;
    out.pivots      = pivots;
    out.freeCols    = freeCols;
    out.kernelCoeff = kernelCoeff;
    out.colMap      = colMap;
    out.a           = aSym;
    out.b           = bSym;
end

%=========================================================================%
% Helper function: create a 1xN cell array of symbolic 3x2 zero matrices.
%=========================================================================%
function tuple = zeroTuple(N)
    tuple = cell(1, N);
    for i = 1:N
        tuple{i} = sym(zeros(3,2));
    end
end

%=========================================================================%
% Helper function: stack a tuple (H_1,...,H_N) into one column vector.
% MATLAB vectorizes a matrix column-by-column via H(:).
%=========================================================================%
function v = tupleToVector(tuple)
    N = numel(tuple);
    v = sym(zeros(6*N, 1));
    for i = 1:N
        rows = (6*(i-1)+1):(6*i);
        v(rows) = tuple{i}(:);
    end
end

%=========================================================================%
% Helper function: robust symbolic zero test.
% For symbolic expressions, isAlways(expr == 0) is preferable when possible.
%=========================================================================%
function tf = isAlwaysZero(expr)
    try
        tf = isAlways(expr == 0);
    catch
        % Fall back to a direct comparison if symbolic decision fails.
        tf = isequal(expr, sym(0));
    end
end
