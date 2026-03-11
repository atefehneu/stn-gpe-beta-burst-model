%% Four Signal Lempel-Ziv Complexity Comparison - Nature Communications Quality
% Professional publication-ready figures with optimized colors and styling

%% Step 1: Load four signals
fprintf('=== LOADING FOUR SIGNALS ===\n');

signal1 = nondelaypaperd4params.VarName2(:);    % Signal 1: nondelaypaperd4params
signal2 = identicaldelays.VarName2(:);          % Signal 2: identicaldelays
signal3 = nondbscsgparams.VarName2(:);          % Signal 3: nondbscsgparams
signal4 = csg.VarName2(:);                     % Signal 4: csg

signals = {signal1, signal2, signal3, signal4};
signal_names = {'NonDelayPaperD4', 'IdenticalDelays', 'NonDBSCSG', 'csg'};
signal_labels = {'Non-Delay Paper D4', 'Identical Delays', 'Non-DBS CSG', 'csg'};

% Display signal info
for i = 1:4
    fprintf('Signal %d (%s) - Length: %d, Range: [%.3f, %.3f]\n', ...
        i, signal_names{i}, length(signals{i}), min(signals{i}), max(signals{i}));
end

%% Step 2: Extract burst durations from all signals
fprintf('\n=== BURST DETECTION FOR ALL SIGNALS ===\n');
burst_durations = cell(4, 1);
burst_info = cell(4, 1);

for i = 1:4
    [burst_durations{i}, burst_info{i}] = extract_burst_durations_enhanced(signals{i}, signal_names{i});
    
    if length(burst_durations{i}) < 10
        warning('Signal %d (%s): Only %d bursts detected. Results may be unreliable.', i, signal_names{i}, length(burst_durations{i}));
    end
end

% Check minimum requirements
min_bursts = min(cellfun(@length, burst_durations));
if min_bursts < 5
    error('Insufficient bursts detected in one or more signals (minimum 5 required, got %d)', min_bursts);
end

%% Step 3: LZC Analysis for all signals
fprintf('\n=== LZC COMPUTATION FOR ALL SIGNALS ===\n');
bin_range = 5:80;
lzc_data = cell(4, 1);

for i = 1:4
    [lzc_data{i}] = compute_lzc_series(burst_durations{i}, bin_range, signal_names{i});
end

%% Step 4: Comparative Analysis Summary
fprintf('\n=== BURST CHARACTERISTICS COMPARISON (ALL SIGNALS) ===\n');
fprintf('%-15s\t%s\t%s\t%s\t%s\n', 'Property', signal_names{:});
fprintf('%s\n', repmat('-', 1, 80));

% Burst counts
counts = cellfun(@length, burst_durations);
fprintf('%-15s\t%d\t\t%d\t\t%d\t\t%d\n', 'Burst Count', counts);

% Mean durations
means = cellfun(@mean, burst_durations);
fprintf('%-15s\t%.1f\t\t%.1f\t\t%.1f\t\t%.1f\n', 'Mean Dur (ms)', means);

% CVs
cvs = cellfun(@(x) std(x)/mean(x), burst_durations);
fprintf('%-15s\t%.3f\t\t%.3f\t\t%.3f\t\t%.3f\n', 'Duration CV', cvs);

%% Step 5: Create Nature Communications Quality Visualization
create_publication_quality_plots(burst_durations, lzc_data, bin_range, signal_labels);

%% Enhanced Burst Detection Function
function [durations, info] = extract_burst_durations_enhanced(signal, signal_name)
    fprintf('Analyzing %s...\n', signal_name);

    ind_thres = find(signal > -50);
    if isempty(ind_thres)
        durations = [];
        info = struct('count',0,'mean_dur',NaN,'std_dur',NaN,'cv',NaN,'range',[NaN NaN]);
        return;
    end

    gapSamples = 600;                 % 60 ms at 0.1 ms sampling
    gapIdx   = find(diff(ind_thres) > gapSamples);
    startPos = [1; gapIdx+1];
    endPos   = [gapIdx; numel(ind_thres)];

    bs = ind_thres(startPos);
    be = ind_thres(endPos);

    durations = 0.1 * (be - bs + 1);  % ms

    % Match your main pipeline: discard very short events
    durations = durations(durations >= 100);

    durations = durations(durations > 0);

    info.count = numel(durations);
    info.mean_dur = mean(durations);
    info.std_dur = std(durations);
    info.cv = info.std_dur / info.mean_dur;
    info.range = [min(durations), max(durations)];

    fprintf('  Detected bursts: %d\n', info.count);
    if info.count > 0
        fprintf('  Duration range: [%.1f, %.1f] ms\n', info.range(1), info.range(2));
        fprintf('  Mean duration: %.1f ± %.1f ms\n', info.mean_dur, info.std_dur);
        fprintf('  CV: %.3f\n', info.cv);
    end
end

%% LZC Computation Function
function lzc_data = compute_lzc_series(burst_durations, bin_range, signal_name)
    fprintf('Computing LZC for %s...\n', signal_name);
    
    lzc_data.bins = bin_range;
    lzc_data.raw = zeros(size(bin_range));
    lzc_data.normalized = zeros(size(bin_range));
    
    for i = 1:length(bin_range)
        discrete_seq = discretize_durations_local(burst_durations, bin_range(i));
        [raw_lzc, norm_lzc] = lempel_ziv_complexity(discrete_seq);
        lzc_data.raw(i) = raw_lzc;
        lzc_data.normalized(i) = norm_lzc;
    end
    
    fprintf('  LZC range: %.3f - %.3f (normalized)\n', ...
        min(lzc_data.normalized), max(lzc_data.normalized));
end

%% NATURE COMMUNICATIONS QUALITY VISUALIZATION
function create_publication_quality_plots(burst_durations, lzc_data, bins, signal_labels)
    
    % PROFESSIONAL NEUTRAL PALETTE (Nature Communications style)
    % Using soft, neutral colors for elegant scientific publications
    colors = [
        0.50, 0.50, 0.50;    % Gray (Signal 1)
        0.70, 0.80, 1.00;    % Pastel Blue (Signal 2) 
        1.00, 0.70, 0.70;    % Pastel Red (Signal 3)
        0.90, 0.90, 0.70;    % Pastel Beige (Signal 4)
    ];
    
    % Lighter versions for fills
    colors_light = colors + 0.3 * (1 - colors);  % 30% lighter
    
    % PUBLICATION SETTINGS for 2x3 layout with professional spacing
    fig_width_cm = 25;    % Optimized width for elegant proportions
    fig_height_cm = 16;   % Better height ratio for professional appearance
    dpi = 300;            % Publication quality DPI
    
    % Convert to pixels for MATLAB
    fig_width_px = fig_width_cm * dpi / 2.54;
    fig_height_px = fig_height_cm * dpi / 2.54;
    
    % Create main figure with elegant professional spacing
    figure('Position', [100, 100, fig_width_px, fig_height_px], ...
           'Color', 'white', 'PaperPositionMode', 'auto', ...
           'PaperUnits', 'centimeters', 'PaperSize', [fig_width_cm, fig_height_cm]);
    
    % PROFESSIONAL SUBPLOT LAYOUT with generous spacing (like your examples)
    % Left margin, bottom margin, subplot width, subplot height, horizontal gap, vertical gap
    left_margin = 0.08;
    bottom_margin = 0.12;  
    subplot_width = 0.26;   % Generous width for each subplot
    subplot_height = 0.35;  % Good height for readability
    h_gap = 0.06;          % Professional horizontal spacing
    v_gap = 0.08;          % Professional vertical spacing
    
    % Calculate positions for 2x3 layout with adjacent heatmaps
    positions = zeros(6, 4);  % [x, y, width, height] for each subplot
    
    % Top row: A, B, C
    positions(1,:) = [left_margin, bottom_margin + subplot_height + v_gap, subplot_width, subplot_height];  % A
    positions(2,:) = [left_margin + subplot_width + h_gap, bottom_margin + subplot_height + v_gap, subplot_width, subplot_height];  % B  
    positions(3,:) = [left_margin + 2*(subplot_width + h_gap), bottom_margin + subplot_height + v_gap, subplot_width, subplot_height];  % C (Heatmap 1)
    
    % Bottom row: D, E, F
    positions(4,:) = [left_margin, bottom_margin, subplot_width, subplot_height];  % D
    positions(5,:) = [left_margin + subplot_width + h_gap, bottom_margin, subplot_width, subplot_height];  % E
    positions(6,:) = [left_margin + 2*(subplot_width + h_gap), bottom_margin, subplot_width, subplot_height];  % F (Heatmap 2)
    
    %% Panel A: Burst Duration Distributions (Clean overlaid histograms)
    subplot('Position', positions(1,:));
    create_duration_histograms(burst_durations, signal_labels, colors, colors_light);
    title('A', 'FontName', 'Arial', 'FontSize', 16, 'FontWeight', 'bold', ...
          'Units', 'normalized', 'Position', [-0.12, 1.08, 0]);
    
    %% Panel B: Duration Sequences (Clean overlaid sequences)  
    subplot('Position', positions(2,:));
    create_duration_sequences(burst_durations, signal_labels, colors);
    title('B', 'FontName', 'Arial', 'FontSize', 16, 'FontWeight', 'bold', ...
          'Units', 'normalized', 'Position', [-0.12, 1.08, 0]);
    
    %% Panel C: Heatmap 1 (Adjacent to F for comparison)
    subplot('Position', positions(3,:));
    create_placeholder_panel('C', 'Heatmap 1');
    title('C', 'FontName', 'Arial', 'FontSize', 16, 'FontWeight', 'bold', ...
          'Units', 'normalized', 'Position', [-0.12, 1.08, 0]);
    
    %% Panel D: Lempel-Ziv Complexity Analysis
    subplot('Position', positions(4,:));
    create_lzc_analysis(lzc_data, bins, signal_labels, colors);
    title('D', 'FontName', 'Arial', 'FontSize', 16, 'FontWeight', 'bold', ...
          'Units', 'normalized', 'Position', [-0.12, 1.08, 0]);
    
    %% Panel E: Coefficient of Variation Comparison
    subplot('Position', positions(5,:));
    create_cv_comparison(burst_durations, signal_labels, colors);
    title('E', 'FontName', 'Arial', 'FontSize', 16, 'FontWeight', 'bold', ...
          'Units', 'normalized', 'Position', [-0.12, 1.08, 0]);
    
    %% Panel F: Heatmap 2 (Adjacent to C for easy comparison)
    subplot('Position', positions(6,:));
    create_placeholder_panel('F', 'Heatmap 2');
    title('F', 'FontName', 'Arial', 'FontSize', 16, 'FontWeight', 'bold', ...
          'Units', 'normalized', 'Position', [-0.12, 1.08, 0]);
    
    % Adjust overall layout for publication
    set(gcf, 'Renderer', 'painters');  % Vector graphics for crisp text
    
    fprintf('\n=== NATURE COMMUNICATIONS QUALITY FIGURE COMPLETE ===\n');
    fprintf('Generated elegant publication layout with professional spacing:\n');
    fprintf('Top row:    (A) Overlaid histograms | (B) Overlaid sequences | (C) Heatmap 1\n');
    fprintf('Bottom row: (D) LZC complexity     | (E) CV comparison       | (F) Heatmap 2\n');
    fprintf('Note: Heatmaps C & F are vertically adjacent for easy comparison\n');
    fprintf('Fixed: Clean, readable overlaid plots instead of cramped mini-grids\n');
    fprintf('Specifications: %g×%g cm, %d DPI, elegant proportions, neutral palette\n', ...
            fig_width_cm, fig_height_cm, dpi);
end

%% Panel A: Clean Overlaid Histograms (READABLE!)
function create_duration_histograms(burst_durations, signal_labels, colors, colors_light)
    hold on;
    
    % Fix: force each to row vector before concatenating
    all_data = cellfun(@(x) x(:)', burst_durations, 'UniformOutput', false);
    all_data = [all_data{:}];
    max_dur = max(all_data);
    n_bins = 25;
    
    for i = 1:4
        histogram(burst_durations{i}(:), n_bins, ...
                 'FaceAlpha', 0.6, 'EdgeAlpha', 0.8, ...
                 'FaceColor', colors(i,:), ...
                 'EdgeColor', colors(i,:) * 0.6, ...
                 'LineWidth', 1.5, ...
                 'Normalization', 'probability', ...
                 'DisplayName', signal_labels{i});
    end
    
    xlabel('Burst Duration (ms)', 'FontName', 'Arial', 'FontSize', 14);
    ylabel('Probability', 'FontName', 'Arial', 'FontSize', 14);
    legend('Location', 'northeast', 'FontName', 'Arial', 'FontSize', 12, ...
           'Box', 'off', 'Interpreter', 'none');
    set(gca, 'FontName', 'Arial', 'FontSize', 12, 'LineWidth', 1.2, ...
             'TickDir', 'out', 'Box', 'on');
    grid on; set(gca, 'GridAlpha', 0.3);
    xlim([0, max_dur * 1.05]);
end

%% Panel B: Clean Overlaid Duration Sequences (READABLE!)
function create_duration_sequences(burst_durations, signal_labels, colors)
    hold on;
    
    % Plot parameters
    markers = {'o', 's', '^', 'd'};
    line_styles = {'-', '--', '-.', ':'};
    n_plot = min(50, min(cellfun(@length, burst_durations)));  % Show first 50 for clarity
    
    % Plot with clear, readable styling
    for i = 1:4
        n_points = min(length(burst_durations{i}), n_plot);
        x_vals = 1:n_points;
        
        plot(x_vals, burst_durations{i}(1:n_points), ...
             'Color', colors(i,:), 'Marker', markers{i}, ...
             'MarkerSize', 6, 'LineStyle', line_styles{i}, ...
             'LineWidth', 2.0, 'MarkerFaceColor', colors(i,:), ...
             'MarkerEdgeColor', colors(i,:) * 0.6, ...
             'DisplayName', signal_labels{i});
    end
    
    % Large, readable labels
    xlabel('Burst Index', 'FontName', 'Arial', 'FontSize', 14, 'FontWeight', 'normal');
    ylabel('Duration (ms)', 'FontName', 'Arial', 'FontSize', 14, 'FontWeight', 'normal');
    
    % Clean legend
    legend('Location', 'best', 'FontName', 'Arial', 'FontSize', 12, ...
           'Box', 'off', 'Interpreter', 'none');
    
    % Professional styling
    set(gca, 'FontName', 'Arial', 'FontSize', 12, 'LineWidth', 1.2, ...
             'TickDir', 'out', 'TickLength', [0.01, 0.01], 'Box', 'on');
    grid on;
    set(gca, 'GridAlpha', 0.3, 'GridColor', [0.6, 0.6, 0.6]);
    
    % Add annotation if truncated
    if max(cellfun(@length, burst_durations)) > n_plot
        text(0.02, 0.98, sprintf('First %d bursts shown', n_plot), ...
             'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 11, ...
             'BackgroundColor', 'white', 'EdgeColor', [0.4, 0.4, 0.4], ...
             'Margin', 3);
    end
end

%% Elegant Placeholder for Future Heatmap Panels
function create_placeholder_panel(panel_letter, description)
    % Create professional placeholder with elegant styling (like your examples)
    hold on;
    
    % Create subtle background
    rectangle('Position', [0.05, 0.05, 0.9, 0.9], ...
              'FaceColor', [0.98, 0.98, 0.98], ...
              'EdgeColor', [0.85, 0.85, 0.85], ...
              'LineWidth', 1.5, 'LineStyle', '-');
    
    % Add elegant text
    text(0.5, 0.6, description, ...
         'Units', 'normalized', 'HorizontalAlignment', 'center', ...
         'VerticalAlignment', 'middle', 'FontName', 'Arial', ...
         'FontSize', 13, 'FontWeight', 'normal', 'Color', [0.4, 0.4, 0.4]);
    
    text(0.5, 0.4, '[Reserved for comparison]', ...
         'Units', 'normalized', 'HorizontalAlignment', 'center', ...
         'VerticalAlignment', 'middle', 'FontName', 'Arial', ...
         'FontSize', 11, 'FontWeight', 'normal', 'Color', [0.6, 0.6, 0.6], ...
         'FontAngle', 'italic');
    
    % Professional axes
    set(gca, 'XTick', [], 'YTick', [], 'XLim', [0, 1], 'YLim', [0, 1], ...
             'Box', 'on', 'LineWidth', 1.0);
end

%% Panel D: Elegant LZC Analysis
function create_lzc_analysis(lzc_data, bins, signal_labels, colors)
    hold on;
    
    markers = {'o', 's', '^', 'd'};
    
    % Plot with elegant styling
    for i = 1:4
        plot(bins, lzc_data{i}.normalized, ...
             'Color', colors(i,:), 'Marker', markers{i}, ...
             'MarkerSize', 6, 'LineWidth', 2.5, ...
             'MarkerFaceColor', colors(i,:), ...
             'MarkerEdgeColor', colors(i,:) * 0.6, ...
             'DisplayName', signal_labels{i});
    end
    
    xlabel('Number of Bins', 'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'normal');
    ylabel('Normalized LZC', 'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'normal');
    
    legend('Location', 'best', 'FontName', 'Arial', 'FontSize', 10, ...
           'Box', 'off', 'Interpreter', 'none');
    
    set(gca, 'FontName', 'Arial', 'FontSize', 11, 'LineWidth', 1.2, ...
             'TickDir', 'out', 'TickLength', [0.01, 0.01], 'Box', 'on');
    grid on;
    set(gca, 'GridAlpha', 0.2, 'GridColor', [0.7, 0.7, 0.7]);
    
    % Set elegant limits
    ylim([0, max(cellfun(@(x) max(x.normalized), lzc_data)) * 1.1]);
end

%% Panel E: Elegant CV Comparison
function create_cv_comparison(burst_durations, signal_labels, colors)
    cvs = cellfun(@(x) std(x)/mean(x), burst_durations);
    
    % Create elegant bar plot
    b = bar(cvs, 'FaceAlpha', 0.8, 'EdgeColor', 'black', 'LineWidth', 1.2);
    
    % Set individual colors
    for i = 1:4
        b.FaceColor = 'flat';
        b.CData(i,:) = colors(i,:);
    end
    
    % Add elegant value labels on bars
    for i = 1:4
        text(i, cvs(i) + max(cvs)*0.03, sprintf('%.3f', cvs(i)), ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
             'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 10, ...
             'Color', colors(i,:) * 0.7);
    end
    
    ylabel('Coefficient of Variation', 'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'normal');
    
    % Elegant abbreviated labels for better fit
    signal_labels_short = {'Non-Delay D4', 'Identical Del.', 'Non-DBS CSG', 'csg'};
    set(gca, 'XTick', 1:4, 'XTickLabel', signal_labels_short, ...
             'FontName', 'Arial', 'FontSize', 10, 'FontWeight', 'normal');
    xtickangle(45);
    
    set(gca, 'FontName', 'Arial', 'FontSize', 11, 'LineWidth', 1.2, ...
             'TickDir', 'out', 'TickLength', [0.01, 0.01], 'Box', 'on');
    grid on;
    set(gca, 'GridAlpha', 0.2, 'GridColor', [0.7, 0.7, 0.7]);
    
    % Set elegant limits
    ylim([0, max(cvs) * 1.2]);
end

%% Helper Functions
function discrete_seq = discretize_durations_local(durations, num_bins)
    durations = durations(:);
    edges = quantile(durations, linspace(0, 1, num_bins + 1));
    edges(1) = edges(1) - eps;
    discrete_seq = discretize(durations, edges);
    discrete_seq(isnan(discrete_seq)) = 1;
    discrete_seq = discrete_seq(:)';
end

%%
function [lzc, lzc_normalized] = lempel_ziv_complexity(sequence, normalize)
    if nargin < 2
        normalize = true;
    end
    
    sequence = round(sequence(:)');
    n = length(sequence);
    
    if n == 0
        lzc = 0;
        lzc_normalized = 0;
        return;
    end
    
    % --- your raw LZC parser ---
    complexity = 1;
    i = 1;
    vocab = {};
    
    while i <= n
        match_length = 0;
        for j = i:n
            test_substring = sequence(i:j);
            test_string = num2str(test_substring);
            if any(strcmp(vocab, test_string))
                match_length = j - i + 1;
            else
                break;
            end
        end
        
        if i + match_length <= n
            new_pattern = sequence(i:(i + match_length));
            vocab{end+1} = num2str(new_pattern);
            complexity = complexity + 1;
            i = i + match_length + 1;
        else
            if match_length == 0
                vocab{end+1} = num2str(sequence(i));
                complexity = complexity + 1;
            end
            break;
        end
    end
    
    lzc = complexity;
    
    if normalize
        alphabet_size = length(unique(sequence));
        if alphabet_size > 1 && n > 1
            lzc_normalized = lzc * (log(n) / log(alphabet_size)) / n;
        else
            lzc_normalized = 0;
        end
    else
        lzc_normalized = lzc;
    end
end

%% EXPORT FUNCTION FOR PUBLICATION
% Add this function call at the end to save publication-ready files
function export_publication_figures()
    % Export as high-quality PDF (vector graphics)
    print(gcf, 'LZC_Analysis_Figure', '-dpdf', '-r300', '-painters');
    
    % Export as high-quality PNG (raster backup)
    print(gcf, 'LZC_Analysis_Figure', '-dpng', '-r300');
    
    % Export as EPS for some journals
    print(gcf, 'LZC_Analysis_Figure', '-depsc2', '-r300', '-painters');
    
    fprintf('Publication files exported:\n');
    fprintf('- LZC_Analysis_Figure.pdf (vector, preferred)\n');
    fprintf('- LZC_Analysis_Figure.png (raster backup)\n');
    fprintf('- LZC_Analysis_Figure.eps (alternative vector)\n');
end