function [sizes, forwardError, inverseError] = accuracy_comparison_Bunck2()
%ACCURACY_COMPARISON_BUNCK2  Compare Bunck and Bunck2 as in Figure 4 of
% Webb and Maierhofer, "Stable Hermite transforms via the Golub-Welsch algorithm".
% Requires Chebfun's hermpts and Python 3 on PATH (standard library only).
% Both error arrays contain relative matrix 2-norms, with Bunck in column 1
% and Bunck2 in column 2. All sizes and both transform directions are saved.

root = fileparts(mfilename('fullpath'));
addpath(fullfile(root, 'modules'));
work = tempname;
mkdir(work);
cleanup = onCleanup(@() rmdir(work, 's'));
sizes = round(2.^(2:0.5:12))';
methods = {@initialise_Hermite_transform_Bunck, ...
           @initialise_Hermite_transform_Bunck2};
forwardError = zeros(numel(sizes), 2);
inverseError = zeros(numel(sizes), 2);

for i = 1:numel(sizes)
    N = sizes(i);
    [dRef, QRef] = referenceFactors(N, root, work);
    % Like the paper's MATLAB experiment, assemble both matrices from
    % reference factors rounded to double. norm(A,2) uses the full matrix.
    T = dRef .* QRef.';
    Tinv = QRef .* (1 ./ dRef).';
    clear QRef
    normT = norm(T, 2);
    normTinv = norm(Tinv, 2);
    for m = 1:numel(methods)
        [d, Q] = methods{m}(N);
        forwardError(i,m) = norm(d .* Q.' - T, 2) / normT;
        inverseError(i,m) = norm(Q .* (1 ./ d).' - Tinv, 2) / normTinv;
    end
    fprintf('N = %4d   T: %.3e %.3e   T inverse: %.3e %.3e\n', ...
            N, forwardError(i,:), inverseError(i,:));
end

methodNames = {'Bunck', 'Bunck2'};
save(fullfile(root, 'data', 'accuracy_Bunck2.mat'), ...
     'sizes', 'methodNames', 'forwardError', 'inverseError');

fig = figure('Color', 'white', 'Position', [100 100 1030 440]);
layout = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
errors = {forwardError, inverseError};
titles = {'Forward $T$ (coefficients to values)', ...
          'Inverse $T^{-1}$ (values to coefficients)'};
allErrors = [forwardError(:); inverseError(:)];
limits = 10.^[floor(log10(min(allErrors))), ceil(log10(max(allErrors)))];
for panel = 1:2
    ax = nexttile(layout);
    bunck2 = loglog(ax, sizes, errors{panel}(:,2), '-o', 'Color', '#00785D', ...
                   'LineWidth', 2.5, 'MarkerSize', 5, 'MarkerFaceColor', 'white');
    hold(ax, 'on');
    % Keep both curves visible where their errors agree.
    bunck = loglog(ax, sizes, errors{panel}(:,1), '-.', 'Color', '#edb120', ...
                  'LineWidth', 1.8);
    grid(ax, 'on');
    xlim(ax, sizes([1 end]));
    ylim(ax, limits);
    set(ax, 'FontSize', 13);
    xlabel(ax, '$N$', 'Interpreter', 'latex', 'FontSize', 17);
    title(ax, titles{panel}, 'Interpreter', 'latex', 'FontSize', 16);
    if panel == 1
        ylabel(ax, 'Relative spectral-norm error', 'Interpreter', 'latex', 'FontSize', 16);
        legend(ax, [bunck bunck2], methodNames, 'Interpreter', 'latex', ...
               'Location', 'northwest', 'Box', 'off');
    end
end
output = fullfile(root, 'images', 'accuracy_Bunck2');
exportgraphics(fig, [output '.pdf'], 'ContentType', 'vector');
exportgraphics(fig, [output '.png'], 'Resolution', 160);
end

function [d, Q] = referenceFactors(N, root, work)
% Refine the nodes and evaluate ordinary monic Hermite polynomials at
% 64 decimal digits. Only the final nonnegative-node factors are exported.
nodesFile = fullfile(work, 'nodes.txt');
factorsFile = fullfile(work, 'factors.bin');
logFile = fullfile(work, 'reference.log');
writematrix(hermpts(N), nodesFile);
python = 'python3';
if ispc
    python = 'python';
end
helper = fullfile(root, 'dev', 'hermite_accuracy_reference.py');
command = java.lang.ProcessBuilder({python, helper, nodesFile, factorsFile});
command.redirectErrorStream(true);
command.redirectOutput(java.io.File(logFile));
process = command.start();
if process.waitFor() ~= 0
    error('Reference generation failed for N=%d. %s', N, fileread(logFile));
end

file = fopen(factorsFile, 'r', 'ieee-le');
closeFile = onCleanup(@() fclose(file));
[records, count] = fread(file, [N+1, ceil(N/2)], '*double');
assert(count == (N+1)*ceil(N/2), 'Incomplete reference factors for N=%d.', N);
d = records(1,:).';
Q = records(2:end,:);
% Reflect the positive nodes. The middle column for odd N occurs once.
mirror = numel(d):-1:(1+mod(N,2));
parity = (-1).^(0:N-1).';
d = [d(mirror); d];
Q = [parity .* Q(:,mirror), Q];
end
