function [sizes, seconds] = timing_comparison_Bunck2()
%TIMING_COMPARISON_BUNCK2  Reproduce the four-method assembly timing figure.
% Requires Chebfun's hermpts on the MATLAB path. Columns of seconds are
% Direct, Bunck, Bunck2, and Golub-Welsch, in legend order.

root = fileparts(mfilename('fullpath'));
addpath(fullfile(root, 'modules'), fullfile(root, 'modules', 'quadrature'));
methods = {@initialise_Hermite_transform_unstable, ...
           @initialise_Hermite_transform_Bunck, ...
           @initialise_Hermite_transform_Bunck2, ...
           @initialise_Hermite_transform_Golub_Welsch};
sizes = round(2.^(2:0.5:12))';
repeats = 20;
seconds = zeros(numel(sizes), numel(methods));

% Warm the small, recurrence, and asymptotic node-generation paths.
for m = 1:numel(methods)
    for n = [4, 32, 256]
        [~, ~] = methods{m}(n);
    end
end

for i = 1:numel(sizes)
    % Rotate which initializer runs first as N increases.
    for m = circshift(1:numel(methods), i-1)
        assemble = methods{m};
        started = tic;
        for trial = 1:repeats
            [~, ~] = assemble(sizes(i));
        end
        seconds(i,m) = toc(started) / repeats;
    end
    fprintf('Timed N = %d\n', sizes(i));
end

fig = figure('Color', 'white', 'Position', [100 100 860 555]);
% Draw Direct over Golub-Welsch so both are visible where they coincide.
h(4) = loglog(sizes, seconds(:,4), '-', 'Color', 'blue', 'LineWidth', 3);
hold on
h(1) = loglog(sizes, seconds(:,1), '--', 'Color', '#c74cb9', 'LineWidth', 3);
h(2) = loglog(sizes, seconds(:,2), '-.', 'Color', '#edb120', 'LineWidth', 3);
h(3) = loglog(sizes, seconds(:,3), '-o', 'Color', '#00785D', 'LineWidth', 2.5, ...
              'MarkerSize', 4, 'MarkerFaceColor', 'white');
guide = 1e-6 * sizes.^2;
h(5) = loglog(sizes, guide, '-.', 'Color', 'red', 'LineWidth', 1);
grid on
set(gca, 'FontSize', 16);
xlim([sizes(1), sizes(end)]);
times = [seconds(:); guide];
decades = floor(log10(min(times))):ceil(log10(max(times)));
ylim(10.^decades([1 end]));
yticks(10.^decades);
xlabel('$N$', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('Assembly time (s)', 'Interpreter', 'latex', 'FontSize', 22);
legend(h, 'Direct', 'Bunck', 'Bunck2', 'Golub--Welsch', ...
       '$\mathcal{O}(N^2)$', 'Interpreter', 'latex', ...
       'FontSize', 16, 'Location', 'northwest', 'Box', 'off');

output = fullfile(root, 'images', 'assembly_time_Bunck2');
exportgraphics(fig, [output '.pdf'], 'ContentType', 'vector');
exportgraphics(fig, [output '.png'], 'Resolution', 160);
end
