function result = checkIncreasingShape(a, b, c, d)
%CHECKINCREASINGSHAPE Classify admissible increasing convex/concave data.
%
%   result = checkIncreasingShape(a,b,c,d)
%
% Given endpoint and integral data
%
%       f(0) = a,
%       f(c) = b,
%       integral_0^c f(x) dx = d,
%
% this function checks whether there exists:
%
%   1. a strictly increasing, strictly convex function f;
%   2. a strictly increasing, strictly concave function f;
%   3. neither.
%
% The function accepts exact numerical values or symbolic expressions.
%
% The relevant necessary and sufficient conditions are:
%
%   Increasing and strictly convex:
%
%       c > 0,
%       b > a,
%       c*a < d < c*(a+b)/2.
%
%   Increasing and strictly concave:
%
%       c > 0,
%       b > a,
%       c*(a+b)/2 < d < c*b.
%
% If MATLAB cannot decide one of the symbolic inequalities from the
% current assumptions, the classification is reported as "Undetermined".
%
% OUTPUT:
%   result.classification     Main classification
%   result.midpointIntegral   c*(a+b)/2
%   result.convexCondition    Complete symbolic convexity condition
%   result.concaveCondition   Complete symbolic concavity condition
%   result.convexStatus       "true", "false", or "undetermined"
%   result.concaveStatus      "true", "false", or "undetermined"

    % Convert numerical inputs to symbolic objects. This allows exact
    % symbolic comparisons and preserves symbolic inputs.
    a = sym(a);
    b = sym(b);
    c = sym(c);
    d = sym(d);

    % Integral of the chord joining (0,a) and (c,b).
    midpointIntegral = c*(a + b)/2;

    % Necessary and sufficient condition for the existence of a strictly
    % increasing, strictly convex function.
    convexCondition = ...
        (c > 0)                & ...
        (b > a)                & ...
        (d > c*a)              & ...
        (d < midpointIntegral);

    % Necessary and sufficient condition for the existence of a strictly
    % increasing, strictly concave function.
    concaveCondition = ...
        (c > 0)                & ...
        (b > a)                & ...
        (d > midpointIntegral) & ...
        (d < c*b);

    % Determine whether each condition is provably true, provably false,
    % or undecidable from the current assumptions.
    cState       = symbolicTruthState(c > 0);
    convexState  = symbolicTruthState(convexCondition);
    concaveState = symbolicTruthState(concaveCondition);

    % Classify the data.
    if cState == 0

        classification = "Invalid input: c must be strictly positive";

    elseif convexState == 1

        classification = ...
            "Strictly increasing and strictly convex";

    elseif concaveState == 1

        classification = ...
            "Strictly increasing and strictly concave";

    elseif convexState == 0 && concaveState == 0

        classification = "Neither";

    else

        classification = ...
            "Undetermined from the current symbolic assumptions";

    end

    % Store all useful information in the output structure.
    result = struct();

    result.classification   = classification;
    result.midpointIntegral = midpointIntegral;

    result.convexCondition  = convexCondition;
    result.concaveCondition = concaveCondition;

    result.convexStatus  = truthStateLabel(convexState);
    result.concaveStatus = truthStateLabel(concaveState);

    % Display the classification.
    fprintf('\nClassification:\n  %s\n', char(classification));

    fprintf('\nRelevant integral thresholds:\n');
    fprintf('  Lower endpoint integral: ');
    disp(c*a);

    fprintf('  Chord integral:          ');
    disp(midpointIntegral);

    fprintf('  Upper endpoint integral: ');
    disp(c*b);

    fprintf('Convex condition:\n  ');
    disp(convexCondition);

    fprintf('Concave condition:\n  ');
    disp(concaveCondition);
end


function state = symbolicTruthState(condition)
%SYMBOLICTRUTHSTATE Three-valued test for a symbolic condition.
%
%   state = 1   means the condition is provably true.
%   state = 0   means the condition is provably false.
%   state = -1  means MATLAB cannot decide from the assumptions.

    if isAlways(condition, 'Unknown', 'false')

        state = 1;

    elseif isAlways(~condition, 'Unknown', 'false')

        state = 0;

    else

        state = -1;

    end
end


function label = truthStateLabel(state)
%TRUTHSTATELABEL Convert the numerical truth state to a readable string.

    switch state
        case 1
            label = "true";

        case 0
            label = "false";

        otherwise
            label = "undetermined";
    end
end