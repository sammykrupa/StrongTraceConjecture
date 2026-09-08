%% export_solution_latex.m
% This script generates the appendix with the values of the T_N parameters 
% It writes publication-ready LaTeX to solution_values.tex
% Run this script after kappa and sol have been created in the workspace.
% 

assert(isa(kappa, 'sym'), ...
    'kappa must be a symbolic array.');

assert(isfield(sol, 'a') && isa(sol.a, 'sym'), ...
    'sol.a must be a symbolic array.');

assert(isfield(sol, 'n') && isa(sol.n, 'sym'), ...
    'sol.n must be a symbolic array.');

assert(isfield(sol, 'P') && isa(sol.P, 'sym'), ...
    'sol.P must be a symbolic array.');

assert(isequal(size(kappa), [6, 1]), ...
    'kappa must have size 6-by-1.');

assert(isequal(size(sol.a), [3, 1, 6]), ...
    'sol.a must have size 3-by-1-by-6.');

assert(isequal(size(sol.n), [1, 2, 6]), ...
    'sol.n must have size 1-by-2-by-6.');

assert(isequal(size(sol.P), [3, 2]), ...
    'sol.P must have size 3-by-2.');

filename = 'solution_values.tex';

[fid, message] = fopen(filename, 'w');

if fid == -1
    error('Could not open %s: %s', filename, message);
end

try
    fprintf(fid, '%s\n', ...
        '% This file was generated automatically by MATLAB.');
    fprintf(fid, '%s\n', ...
        '% It requires the amsmath package.');
    fprintf(fid, '\n');

    % Write P.
    fprintf(fid, '%s\n', '\begin{equation*}');
    fprintf(fid, '%s\n', ...
        ['P = ', symbolicPmatrix(sol.P), '.']);
    fprintf(fid, '%s\n', '\end{equation*}');
    fprintf(fid, '\n');

    % Write kappa_i, a_i, and n_i for i = 1,...,6.
    for i = 1:6
        ai = sol.a(:, :, i);    % 3-by-1 column vector
        ni = sol.n(:, :, i);    % 1-by-2 row vector

        fprintf(fid, '%s\n', '\begin{equation*}');
        fprintf(fid, '%s\n', '\begin{aligned}');

        fprintf(fid, '%s\n', ...
            ['\kappa_{', num2str(i), '} &= ', ...
             char(latex(kappa(i))), ',\\']);

        fprintf(fid, '%s\n', ...
            ['a_{', num2str(i), '} &= ', ...
             symbolicPmatrix(ai), ',\\']);

        fprintf(fid, '%s\n', ...
            ['n_{', num2str(i), '} &= ', ...
             symbolicPmatrix(ni), '.']);

        fprintf(fid, '%s\n', '\end{aligned}');
        fprintf(fid, '%s\n', '\end{equation*}');
        fprintf(fid, '\n');
    end

    fclose(fid);

catch ME
    fclose(fid);
    rethrow(ME);
end

fprintf('LaTeX values written to %s\n', filename);


%% Local function: convert a symbolic matrix to a LaTeX pmatrix
function output = symbolicPmatrix(M)

    M = sym(M);
    [numberOfRows, numberOfColumns] = size(M);

    latexRows = cell(numberOfRows, 1);

    for row = 1:numberOfRows
        latexEntries = cell(1, numberOfColumns);

        for column = 1:numberOfColumns
            latexEntries{column} = char(latex(M(row, column)));
        end

        latexRows{row} = strjoin(latexEntries, ' & ');
    end

    matrixBody = strjoin(latexRows, ' \\ ');

    output = [
        '\begin{pmatrix}', ...
        matrixBody, ...
        '\end{pmatrix}'
    ];
end