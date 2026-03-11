%%Clinical vs Computational Comparison Figure (v3: 6 OFF+ON datasets, Hedges g per panel)

%% 1. Prepare Clinical Data
clinical_OFF = [
    11.38, 197, 487, 632, 1.01, 768;
    12.50, 187, 457, 746, 0.98, 794;
    11.59, 202, 463, 694, 0.92, 829;
    11.39, 191, 443, 720, 1.05, 709;
    13.7,  178, 496, 712, 0.97, 776;
    7.69,  185, 531, 724, 1.04, 749;
    5.26,  196, 489, 804, 1.01, 778;
    13.64, 203, 485, 680, 0.88, 875;
    8.43,  187, 455, 664, 1.1,  692;
    18.64, 216, 459, 648, 0.78, 945;
    16.92, 213, 495, 649, 0.87, 894;
    12.86, 187, 509, 744, 0.93, 847
];

clinical_ON = [
    20.31, 202, 462, 676, 0.86, 866;
    5.13,  205, 479, NaN, 1.05, 735;
    3.75,  210, 445, NaN, 1.06, 716;
    8.45,  206, 449, 728, 0.94, 830;
    3.85,  196, 442, 612, 1.04, 759;
    4.94,  175, 502, NaN, 1.08, 733;
    11.84, 190, 497, 732, 1.01, 752;
    15.15, 203, 483, 832, 0.88, 869;
    3.33,  192, 427, NaN, 1.2,  624;
    11.76, 184, 526, 776, 0.90, 864;
    4.92,  198, 488, NaN, 1.07, 732;
    0.00,  167, NaN, NaN, 1.44, 533
];

%% 2. Define All Computational Datasets
% Replace off_ds2…off_ds5 and on_ds2…on_ds5 with your actual MATLAB table variables.
off_datasets = {neunon1, neunon21, neunon31, neunon41, neunon51, neunon61};
on_datasets  = {neudbs1, neudbs28, neudbs32, neudbs413, neudbs516, neudbs64};

num_datasets = numel(off_datasets);   % = 6
num_neurons  = 10;                    % VarName2 … VarName11
varCols      = 2:11;

%% 3. Extract Voltage Matrices
fprintf('Extracting computational data...\n');

OFF_volt_cell = cell(1, num_datasets);
ON_volt_cell  = cell(1, num_datasets);

for d = 1:num_datasets
    tbl_off = off_datasets{d};
    tbl_on  = on_datasets{d};
    cols_off = zeros(height(tbl_off), num_neurons);
    cols_on  = zeros(height(tbl_on),  num_neurons);
    for v = 1:num_neurons
        vname = sprintf('VarName%d', varCols(v));
        cols_off(:,v) = tbl_off.(vname);
        cols_on(:,v)  = tbl_on.(vname);
    end
    OFF_volt_cell{d} = cols_off;
    ON_volt_cell{d}  = cols_on;
end

%% 4. Analyse Each Neuron → Per-Dataset Mean Metrics
% comp_OFF_all / comp_ON_all : [num_datasets × 6]  (mean across neurons per dataset)

comp_OFF_all = zeros(num_datasets, 6);
comp_ON_all  = zeros(num_datasets, 6);

for d = 1:num_datasets
    neuron_off = zeros(num_neurons, 6);
    neuron_on  = zeros(num_neurons, 6);
    for n = 1:num_neurons
        m = analyze_bursts(OFF_volt_cell{d}(:,n));
        neuron_off(n,:) = [m.ls_ratio, m.mean_100_400, m.mean_400_600, ...
                           m.mean_600_900, m.burst_rate, m.inter_burst_interval];
        m = analyze_bursts(ON_volt_cell{d}(:,n));
        neuron_on(n,:)  = [m.ls_ratio, m.mean_100_400, m.mean_400_600, ...
                           m.mean_600_900, m.burst_rate, m.inter_burst_interval];
    end
    comp_OFF_all(d,:) = nanmean(neuron_off, 1);
    comp_ON_all(d,:)  = nanmean(neuron_on,  1);
end

fprintf('Data extraction complete!\n\n');

% Alias for downstream code (5 paired dataset-mean rows)
comp_OFF = comp_OFF_all;   % [5×6]
comp_ON  = comp_ON_all;    % [5×6]

%% 5. Effect-Size Analysis (paired Hedges' g) — used by all panels
rng(1);
nBoot = 5000;
metric_short = {'L/S','100-400','400-600','600-900','Rate','IBI'};

dzClin  = nan(1,6); ciDzClin  = nan(6,2); nClin  = nan(1,6);
dzModel = nan(1,6); ciDzModel = nan(6,2); nModel = nan(1,6);
gClin   = nan(1,6); ciClin    = nan(6,2);
gModel  = nan(1,6); ciModel   = nan(6,2);

for k = 1:6
    [dzClin(k),  gClin(k),  ciDzClin(k,:),  ciClin(k,:),  nClin(k)]  = ...
        paired_dz_gz_boot(clinical_OFF(:,k), clinical_ON(:,k), nBoot, 0.05);
    [dzModel(k), gModel(k), ciDzModel(k,:), ciModel(k,:), nModel(k)] = ...
        paired_dz_gz_boot(comp_OFF(:,k),     comp_ON(:,k),     nBoot, 0.05);
end

fprintf('\n=== EFFECT SIZES (paired Hedges g; OFF-ON) ===\n');
for k = 1:6
    fprintf('%-9s  Clinical: g=%.2f [%.2f, %.2f] (n=%d)   Model: g=%.2f [%.2f, %.2f] (n=%d)\n', ...
        metric_short{k}, ...
        gClin(k),  ciClin(k,1),  ciClin(k,2),  nClin(k), ...
        gModel(k), ciModel(k,1), ciModel(k,2), nModel(k));
end

ok = ~isnan(gClin) & ~isnan(gModel);
pearson_r  = corr(gClin(ok)', gModel(ok)', 'Type','Pearson');
spearman_r = corr(gClin(ok)', gModel(ok)', 'Type','Spearman');
rmse_es    = sqrt(mean((gModel(ok) - gClin(ok)).^2));
fprintf('\nAcross metrics (M=%d): Pearson r=%.3f, Spearman rho=%.3f, RMSE=%.3f\n', ...
    sum(ok), pearson_r, spearman_r, rmse_es);

%% 6. Colors
color_clinical_off = [0.85, 0.85, 0.85];
color_clinical_on  = [0.80, 0.90, 1.00];
color_comp_off     = [1.00, 0.70, 0.70];
color_comp_on      = [1.00, 0.85, 0.85];

slope_colors = {
    [0.8, 0.9, 1.0];   % 100-400 ms
    [1.0, 0.8, 0.8];   % 400-600 ms
    [0.9, 0.9, 0.8]    % 600-900 ms
};

metric_names = {'L/S ratio', 'Burst rate (bursts/s)', 'IBI (ms)'};

%% 7. Figure  (3-row × 4-col layout)
% Row 1 : (A) L/S Clinical  (B) L/S Model  (C) Rate Clinical  (D) Rate Model
% Row 2 : (E) IBI Clinical  (F) IBI Model  (G) Duration slope [spans cols 3-4]
% Row 3 : (H) Calibration   [spans cols 1-2]   (I) Effect sizes [spans cols 3-4]

fig = figure('Position', [50, 50, 1800, 1200], 'Color', 'w');

bar_indices  = [1, 5, 6];   % columns: L/S, Rate, IBI
panel_labels = {'A','B','C','D','E','F'};
label_idx    = 1;

%% ---- Rows 1–2 : Box plots (A–F) with Hedges' g annotation ----
% metric_idx 1→L/S (col 1), 2→Rate (col 5), 3→IBI (col 6)

for metric_idx = 1:3
    col_i = bar_indices(metric_idx);   % column index into clinical/comp matrices

    % Determine subplot positions
    if metric_idx <= 2
        sp_clin = (metric_idx-1)*2 + 1;   % row 1: subplots 1,3
        sp_mod  = (metric_idx-1)*2 + 2;   % row 1: subplots 2,4
    else
        sp_clin = 5;                        % row 2: subplot 5
        sp_mod  = 6;                        % row 2: subplot 6
    end

    % ---- Clinical subplot ----
    subplot(3, 4, sp_clin);
    clin_off_all = clinical_OFF(~isnan(clinical_OFF(:,col_i)), col_i);
    clin_on_all  = clinical_ON(~isnan(clinical_ON(:,col_i)),   col_i);
    hold on;
    drawBoxPlot(1, clin_off_all, color_clinical_off);
    drawBoxPlot(2, clin_on_all,  color_clinical_on);

    % Significance bracket
    valid_p = ~isnan(clinical_OFF(:,col_i)) & ~isnan(clinical_ON(:,col_i));
    [p_val, y_top] = sigBracket(clinical_OFF(valid_p,col_i), clinical_ON(valid_p,col_i), ...
                                [clin_off_all; clin_on_all]);

    % Hedges' g annotation (bottom-right corner)
    addHedgesG(gClin(col_i), ciClin(col_i,:));

    xlim([0.5, 2.5]);
    set(gca,'XTick',[1,2],'XTickLabel',{'OFF','ON'},'FontSize',11,...
        'Box','off','TickDir','out','LineWidth',1);
    ylabel([metric_names{metric_idx} ' (Clinical)'],'FontWeight','bold','FontSize',11);
    grid on; set(gca,'GridAlpha',0.3);
    title(['(' panel_labels{label_idx} ') Clinical'],'FontSize',13,'FontWeight','bold');
    label_idx = label_idx + 1;

    % ---- Model subplot ----
    subplot(3, 4, sp_mod);
    comp_off_all = comp_OFF(~isnan(comp_OFF(:,col_i)), col_i);
    comp_on_all  = comp_ON(~isnan(comp_ON(:,col_i)),   col_i);
    hold on;
    drawBoxPlot(1, comp_off_all, color_comp_off);
    drawBoxPlot(2, comp_on_all,  color_comp_on);

    % Significance bracket
    valid_p = ~isnan(comp_OFF(:,col_i)) & ~isnan(comp_ON(:,col_i));
    [p_val, y_top] = sigBracket(comp_OFF(valid_p,col_i), comp_ON(valid_p,col_i), ...
                                [comp_off_all; comp_on_all]);

    % Hedges' g annotation
    addHedgesG(gModel(col_i), ciModel(col_i,:));

    xlim([0.5, 2.5]);
    set(gca,'XTick',[1,2],'XTickLabel',{'OFF','ON'},'FontSize',11,...
        'Box','off','TickDir','out','LineWidth',1);
    ylabel([metric_names{metric_idx} ' (Model mean)'],'FontWeight','bold','FontSize',11);
    grid on; set(gca,'GridAlpha',0.3);
    title(['(' panel_labels{label_idx} ') Model'],'FontSize',13,'FontWeight','bold');
    label_idx = label_idx + 1;
end

%% ---- Row 2 cols 3-4 : Duration Slope Chart (G) ----
subplot(3, 4, [7, 8]);
duration_ranges = {'100-400ms','400-600ms','600-900ms'};
hold on;

for range_idx = 1:3
    col_idx = range_idx + 1;   % columns 2,3,4

    % Clinical
    pc = ~isnan(clinical_OFF(:,col_idx)) & ~isnan(clinical_ON(:,col_idx));
    off_c = clinical_OFF(pc, col_idx); on_c = clinical_ON(pc, col_idx);
    for j = 1:length(off_c)
        plot([1,2],[off_c(j),on_c(j)],'-','Color',[0.6,0.6,0.6,0.3],'LineWidth',1);
    end
    mu_oc = mean(off_c); mu_nc = mean(on_c);
    plot([1,2],[mu_oc,mu_nc],'-o','Color',slope_colors{range_idx},'LineWidth',3,...
        'MarkerSize',10,'MarkerFaceColor',slope_colors{range_idx},'MarkerEdgeColor','k');
    text(0.85,mu_oc,sprintf('%.0f',mu_oc),'HorizontalAlignment','right','FontSize',9,'FontWeight','bold');
    text(2.15,mu_nc,sprintf('%.0f',mu_nc),'HorizontalAlignment','left', 'FontSize',9,'FontWeight','bold');

    % Model
    pm = ~isnan(comp_OFF(:,col_idx)) & ~isnan(comp_ON(:,col_idx));
    off_m = comp_OFF(pm, col_idx); on_m = comp_ON(pm, col_idx);
    for j = 1:length(off_m)
        plot([3.5,4.5],[off_m(j),on_m(j)],'-','Color',[1.0,0.7,0.7,0.3],'LineWidth',1);
    end
    mu_om = mean(off_m); mu_nm = mean(on_m);
    plot([3.5,4.5],[mu_om,mu_nm],'-o','Color',slope_colors{range_idx},'LineWidth',3,...
        'MarkerSize',10,'MarkerFaceColor',slope_colors{range_idx},'MarkerEdgeColor','k');
    text(3.35,mu_om,sprintf('%.0f',mu_om),'HorizontalAlignment','right','FontSize',9,'FontWeight','bold');
    text(4.65,mu_nm,sprintf('%.0f',mu_nm),'HorizontalAlignment','left', 'FontSize',9,'FontWeight','bold');
end

all_dur = [clinical_OFF(:,2:4); clinical_ON(:,2:4); comp_OFF(:,2:4); comp_ON(:,2:4)];
all_dur = all_dur(~isnan(all_dur));
ylim([min(all_dur)-50, max(all_dur)+50]);
xlim([0.5, 5]);
set(gca,'XTick',[1.5,4],'XTickLabel',{'Clinical','Model'},'FontSize',11,...
    'Box','off','TickDir','out','LineWidth',1);
ylabel('Duration (ms)','FontSize',12,'FontWeight','bold');
title('(G) Duration Range Transitions','FontSize',13,'FontWeight','bold');
grid on; set(gca,'GridAlpha',0.3);
lh = [plot(NaN,NaN,'-o','Color',slope_colors{1},'LineWidth',2.5,'MarkerSize',8,'MarkerFaceColor',slope_colors{1}), ...
      plot(NaN,NaN,'-o','Color',slope_colors{2},'LineWidth',2.5,'MarkerSize',8,'MarkerFaceColor',slope_colors{2}), ...
      plot(NaN,NaN,'-o','Color',slope_colors{3},'LineWidth',2.5,'MarkerSize',8,'MarkerFaceColor',slope_colors{3})];
legend(lh, duration_ranges,'Location','best','FontSize',9);

%% ---- Row 3 : Calibration (H) + Effect Sizes (I) ----
subplot(3, 4, [9, 10]);
calibrationAgreementPanel(metric_short, gClin, ciClin, gModel, ciModel);
title(sprintf('(H) Calibration: ES  (Pearson r=%.2f, Spearman \\rho=%.2f, RMSE=%.2f)', ...
    pearson_r, spearman_r, rmse_es),'FontSize',13,'FontWeight','bold');

subplot(3, 4, [11, 12]);
hold on;
x_positions = 1:6;
width = 0.35;
for i = 1:6
    plotBar(x_positions(i) - width/2, width, dzClin(i),  color_clinical_off);
    plotBar(x_positions(i) + width/2, width, dzModel(i), color_comp_off);

    % Value labels
    yc = dzClin(i)  + sign(dzClin(i))*0.08;
    ym = dzModel(i) + sign(dzModel(i))*0.08;
    if dzClin(i) == 0,  yc = 0.08; end
    if dzModel(i) == 0, ym = 0.08; end
    text(x_positions(i)-width/2, yc, sprintf('%.2f', dzClin(i)), ...
        'HorizontalAlignment','center','FontSize',10,'FontWeight','bold');
    text(x_positions(i)+width/2, ym, sprintf('%.2f', dzModel(i)), ...
        'HorizontalAlignment','center','FontSize',10,'FontWeight','bold');
end

yline(0,   'k-',  'LineWidth', 2.5);
yline( 0.5,'k:',  'LineWidth', 1.5);
yline(-0.5,'k:',  'LineWidth', 1.5);
yline( 0.8,'k:',  'LineWidth', 1, 'Alpha', 0.5);
yline(-0.8,'k:',  'LineWidth', 1, 'Alpha', 0.5);
xlim([0.5, 6.5]);
xlabel('Metric','FontWeight','bold','FontSize',13);
ylabel('Effect size (paired Cohen''s d_z)','FontWeight','bold','FontSize',13);
set(gca,'XTick',1:6,'XTickLabel',metric_short,'FontSize',11,...
    'Box','off','TickDir','out','LineWidth',1);
grid on; set(gca,'GridAlpha',0.3);
h1 = patch(NaN,NaN, color_clinical_off, 'EdgeColor','k','LineWidth',2);
h2 = patch(NaN,NaN, color_comp_off,     'EdgeColor','k','LineWidth',2);
legend([h1,h2],{'Clinical','Model'},'Location','northeast','FontSize',11,'Box','off');
title('(I) Effect Sizes: OFF vs ON DBS','FontSize',13,'FontWeight','bold');

fprintf('\n=== FIGURE COMPLETE ===\n');

%% =========================================================================
%% HELPER FUNCTIONS
%% =========================================================================

%% addHedgesG 
function addHedgesG(g_val, ci_val)
% Adds a small text box: "g = X.XX [lo, hi]" in the lower-right of the current axes.
% Must be called while hold is on, BEFORE xlim/ylim are finalized
% (uses 'Units','normalized' so position is axes-relative).
    if isnan(g_val)
        txt = 'g = NaN';
    else
        txt = sprintf('g = %.2f\n[%.2f, %.2f]', g_val, ci_val(1), ci_val(2));
    end
    text(0.97, 0.04, txt, ...
        'Units',              'normalized', ...
        'HorizontalAlignment','right', ...
        'VerticalAlignment',  'bottom', ...
        'FontSize',           9, ...
        'FontWeight',         'bold', ...
        'BackgroundColor',    'white', ...
        'EdgeColor',          [0.6 0.6 0.6], ...
        'Margin',             3, ...
        'LineWidth',          1);
end

%% sigBracket
function [p_val, y_top] = sigBracket(off_paired, on_paired, all_vals)
    p_val = NaN; y_top = max(all_vals(~isnan(all_vals)));
    if numel(off_paired) >= 3
        [~, p_val] = ttest(off_paired, on_paired);
        y_top = max(all_vals(~isnan(all_vals)));
        if p_val < 0.05
            y_sig = y_top * 1.08;
            plot([1,2],[y_sig,y_sig],'k-','LineWidth',2);
            if     p_val < 0.001, sig_txt = '***';
            elseif p_val < 0.01,  sig_txt = '**';
            else,                 sig_txt = '*';
            end
            text(1.5, y_sig*1.02, sig_txt, 'HorizontalAlignment','center', ...
                'FontSize',12,'FontWeight','bold');
        end
        ylim([0, y_top*1.15]);
    end
end

%% plotBar 
function plotBar(x_ctr, w, val, col)
    if isnan(val), return; end
    if val >= 0
        rectangle('Position',[x_ctr-w/2, 0,   w, val],  'FaceColor',col,'EdgeColor','k','LineWidth',2);
    else
        rectangle('Position',[x_ctr-w/2, val, w, abs(val)], 'FaceColor',col,'EdgeColor','k','LineWidth',2);
    end
end

%% drawBoxPlot
function drawBoxPlot(x_pos, data, color)
    data = data(~isnan(data));
    if length(data) < 2, return; end
    q1 = prctile(data,25); q2 = median(data); q3 = prctile(data,75);
    iqr_v = q3 - q1;
    lw = max(min(data), q1 - 1.5*iqr_v);
    uw = min(max(data), q3 + 1.5*iqr_v);
    bw = 0.4;
    rectangle('Position',[x_pos-bw/2, q1, bw, q3-q1],...
        'FaceColor',color,'EdgeColor','k','LineWidth',2);
    plot([x_pos-bw/2, x_pos+bw/2],[q2,q2],'k-','LineWidth',3);
    plot([x_pos,x_pos],[q3,uw],'k-','LineWidth',2);
    plot([x_pos,x_pos],[q1,lw],'k-','LineWidth',2);
    plot([x_pos-bw/4, x_pos+bw/4],[uw,uw],'k-','LineWidth',2);
    plot([x_pos-bw/4, x_pos+bw/4],[lw,lw],'k-','LineWidth',2);
    jit = (rand(length(data),1)-0.5)*0.15;
    scatter(x_pos+jit, data, 30,'k','filled','MarkerFaceAlpha',0.4);
end

%% analyze_bursts — burst detection and metric extraction
function metrics = analyze_bursts(data)
    empty_m = struct('dur_t',[],'total_bursts',0,'long_bursts_count',0,...
        'short_bursts_count',0,'long_proportion',0,'short_proportion',0,...
        'ls_ratio',0,'burst_rate',0,'mean_duration',NaN,'median_duration',NaN,...
        'std_duration',NaN,'mean_100_400',NaN,'mean_400_600',NaN,...
        'mean_600_900',NaN,'inter_burst_interval',NaN,'mean_long_duration',0);

    ind_thres = find(data > -50);
    if isempty(ind_thres), metrics = empty_m; return; end

    gapSamples = 600;
    gapIdx   = find(diff(ind_thres) > gapSamples);
    startPos = [1; gapIdx+1]; endPos = [gapIdx; numel(ind_thres)];
    bs = ind_thres(startPos); be = ind_thres(endPos);
    dur_t = 0.1*(be - bs + 1);

    keep  = dur_t >= 100;
    dur_t = dur_t(keep); bs = bs(keep); be = be(keep);
    if isempty(dur_t), metrics = empty_m; return; end

    rec_s = length(data)*0.1/1000;
    long_b   = dur_t(dur_t >= 600); short_b = dur_t(dur_t < 600);
    d100_400 = dur_t(dur_t>=100 & dur_t<400);
    d400_600 = dur_t(dur_t>=400 & dur_t<600);
    d600_900 = dur_t(dur_t>=600 & dur_t<900);

    metrics.dur_t              = dur_t;
    metrics.total_bursts       = length(dur_t);
    metrics.long_bursts_count  = length(long_b);
    metrics.short_bursts_count = length(short_b);
    metrics.long_proportion    = length(long_b)/length(dur_t)*100;
    metrics.short_proportion   = length(short_b)/length(dur_t)*100;
    if metrics.short_bursts_count > 0
        metrics.ls_ratio = metrics.long_bursts_count/metrics.short_bursts_count*100;
    else
        metrics.ls_ratio = metrics.long_bursts_count*100;
    end
    metrics.burst_rate      = length(dur_t)/rec_s;
    metrics.mean_duration   = mean(dur_t);
    metrics.median_duration = median(dur_t);
    metrics.std_duration    = std(dur_t);
    metrics.mean_100_400 = NaN; metrics.mean_400_600 = NaN; metrics.mean_600_900 = NaN;
    if ~isempty(d100_400), metrics.mean_100_400 = mean(d100_400); end
    if ~isempty(d400_600), metrics.mean_400_600 = mean(d400_600); end
    if ~isempty(d600_900), metrics.mean_600_900 = mean(d600_900); end
    if length(dur_t) > 1
        metrics.inter_burst_interval = mean((bs(2:end)-be(1:end-1))*0.1);
    else
        metrics.inter_burst_interval = NaN;
    end
    metrics.mean_long_duration = 0;
    if ~isempty(long_b), metrics.mean_long_duration = mean(long_b); end
end

%% paired_dz_gz_boot
function [dz, gz, ci_dz, ci_gz, n] = paired_dz_gz_boot(off, on, nBoot, alpha)
    if nargin<3||isempty(nBoot), nBoot=5000; end
    if nargin<4||isempty(alpha), alpha=0.05;  end
    paired = ~isnan(off) & ~isnan(on);
    d = off(paired) - on(paired); n = numel(d);
    dz=NaN; gz=NaN; ci_dz=[NaN NaN]; ci_gz=[NaN NaN];
    if n<3, return; end
    sd = std(d,0);
    if sd<=0||~isfinite(sd), return; end
    dz = mean(d)/sd;
    df = n-1; J = 1-3/(4*df-1); gz = J*dz;
    dzBoot = NaN(nBoot,1);
    for b = 1:nBoot
        idx = randi(n,n,1); db = d(idx); sdb = std(db,0);
        if sdb>0&&isfinite(sdb), dzBoot(b)=mean(db)/sdb; end
    end
    dzBoot = dzBoot(isfinite(dzBoot));
    if isempty(dzBoot), return; end
    ci_dz = prctile(dzBoot, [100*(alpha/2), 100*(1-alpha/2)]);
    ci_gz = J*ci_dz;
end

function calibrationAgreementPanel(metricNames, gClin, ciClin, gModel, ciModel)
    hold on;
    x = gClin(:); y = gModel(:);
    ok = ~isnan(x) & ~isnan(y);
    x=x(ok); y=y(ok); ciX=ciClin(ok,:); ciY=ciModel(ok,:);
    names=metricNames(ok); origIdx=find(ok);
    if isempty(x)
        text(0.5,0.5,'No valid metrics','Units','normalized','HorizontalAlignment','center');
        axis off; return;
    end

    % Axis limits based on point estimates only (not CIs) + small padding
    allPts = [x; y];
    pad = 0.5;
    lim = [min(allPts)-pad, max(allPts)+pad];

    mcol=[0.00,0.45,0.74; 0.85,0.33,0.10; 0.93,0.69,0.13;
          0.49,0.18,0.56; 0.47,0.67,0.19; 0.30,0.75,0.93];

    % Identity line
    plot(lim, lim, 'k--', 'LineWidth', 1.5);

    legH = [];
    for i = 1:numel(x)
        col = mcol(origIdx(i),:);

        % Clip CI whiskers to axis limits for display only
        ciX_plot = max(min(ciX(i,:), lim(2)), lim(1));
        ciY_plot = max(min(ciY(i,:), lim(2)), lim(1));

        plot([ciX_plot(1), ciX_plot(2)], [y(i), y(i)], '-', 'Color', col, 'LineWidth', 2);
        plot([x(i), x(i)], [ciY_plot(1), ciY_plot(2)], '-', 'Color', col, 'LineWidth', 2);

        % Arrow tip if CI exceeds axis
        if ciY(i,1) < lim(1)
            plot(x(i), lim(1)+0.1, 'v', 'Color', col, 'MarkerSize', 7, 'MarkerFaceColor', col);
        end
        if ciY(i,2) > lim(2)
            plot(x(i), lim(2)-0.1, '^', 'Color', col, 'MarkerSize', 7, 'MarkerFaceColor', col);
        end
        if ciX(i,1) < lim(1)
            plot(lim(1)+0.1, y(i), '<', 'Color', col, 'MarkerSize', 7, 'MarkerFaceColor', col);
        end

        h = scatter(x(i), y(i), 120, col, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
        text(x(i), y(i), [' ' names{i}], 'FontSize', 9, 'FontWeight', 'bold', 'VerticalAlignment', 'middle');
        legH(i) = h;
    end

    pr = corr(x, y, 'Type', 'Pearson');
    sr = corr(x, y, 'Type', 'Spearman');
    rm = sqrt(mean((y - x).^2));
    txt = sprintf('Pearson r = %.2f\nSpearman \\rho = %.2f\nRMSE = %.2f', pr, sr, rm);
    text(lim(1)+0.05*range(lim), lim(2)-0.05*range(lim), txt, ...
        'FontSize', 11, 'FontWeight', 'bold', 'VerticalAlignment', 'top', ...
        'BackgroundColor', 'white', 'Margin', 4, 'EdgeColor', [0.2,0.2,0.2], 'LineWidth', 1.2);

    xlabel('Clinical effect size (Hedges g)',  'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Model effect size (Hedges g)',      'FontSize', 12, 'FontWeight', 'bold');
    xlim(lim); ylim(lim); axis square;
    grid on; set(gca, 'GridAlpha', 0.3, 'Box', 'off', 'TickDir', 'out', 'LineWidth', 1);
    legend(legH, names, 'Location', 'southeast', 'FontSize', 9);
end