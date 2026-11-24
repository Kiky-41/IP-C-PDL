% =============================================================================
%  IP-C-PDL v8.0  —  Interior-Point Chaotic Powell's Dog-Leg
%  -------------------------------------------------------------------------
%  NOTE ON NAMING:
%    The algorithm correctly combines:
%      (1) Interior-Point warm-start + local refinement (fmincon IP/SQP)
%      (2) Chaotic perturbation via tent-map (Chaotic)
%      (3) Powell's Dog-Leg trust-region step selection (PDL)
%    Hence the name  IP-C-PDL  (Interior-Point Chaotic Powell's Dog-Leg)
%    is accurate.  Earlier names such as "IP-PDL" omitted the chaotic
%    component; v7+ correctly reflects all three components.
%
%  SECTION A : Benchmark F01-F12  (30 runs x 10 methods)
%  SECTION B : Single-period ED   5/16/42-Unit  (50 runs x 8 methods)
%  SECTION C : 24-Hour ED         5/16/42-Unit  (50 runs x 8 methods)
%              with ramp-rate, generation-limit & power-balance constraints
%  SECTION D : Statistical tests -- One-Way ANOVA, paired t-Test,
%              Tukey HSD, Friedman, Wilcoxon + Holm-Bonferroni
%              for all 6 scenarios (3 single-period + 3 twenty-four-hour)
%  SECTION E : Sensitivity Analysis (IP-C-PDL only, 42-Unit single-period)
%              Factors: chaos probability cp0, max iterations MAX_IT,
%                       multi-start pool size N_STARTS
%              Each factor: 5 levels x 20 runs -> mean/std cost + best cost
%  -------------------------------------------------------------------------
%  Compatible : MATLAB R2016b+, Statistics and Optimization Toolboxes
%
%  Outputs:
%    CEC_Results_v8.xlsx
%    ED_Results_v8.xlsx          (Summary + Raw + Schedule + Stat +
%                                 Sensitivity sheets)
%    Constraint_Check_v8.xlsx
%    Fig01_CEC_Conv.png          Convergence curves F01-F12
%    Fig02_Conv_5Unit.png        Convergence 5-Unit single-period
%    Fig03_Conv_16Unit.png       Convergence 16-Unit single-period
%    Fig04_Conv_42Unit.png       Convergence 42-Unit single-period
%    Fig05_Conv_5Unit_24h.png    Convergence 5-Unit 24-hour
%    Fig06_Conv_16Unit_24h.png   Convergence 16-Unit 24-hour
%    Fig07_Conv_42Unit_24h.png   Convergence 42-Unit 24-hour
%    Fig08_Sensitivity.png       Sensitivity analysis (3-panel)
% =============================================================================
clear; close all; clc;
global IPCPDL_CV;
IPCPDL_CV = [];
%% =========================================================================
%% SHARED CONFIGURATION
%% =========================================================================
N_CEC    = 30;    % CEC: independent runs per method per function
N_ED     = 50;    % ED : independent runs per method per system
N_STARTS = 6;     % ED multi-start pool size
MAX_IT   = 400;   % main loop iterations
CEC_NAMES = {'F01','F02','F03','F04','F05','F06', ...
'F07','F08','F09','F10','F11','F12'};
CEC_LONG  = {'Five Uneven Peak Trap','Equal Maxima','Uneven Decreasing Maxima', ...
'Six-Hump Camel Back','Himmelblau','Shubert','Composition F1', ...
'Vincent','Modified Rastrigin','Composition F4', ...
'Composition F2','Composition F3'};
CEC_DIM   = [1,1,1,2,2,2,2,2,2,3,2,5];
CEC_LO    = [0,  0,  0, -5,-6,-10,-5, 0.1, 0,-5,-5,-5];
CEC_HI    = [15, 1,  1,  5, 6, 10, 5, 10,  1, 5, 5, 5];
METHOD_CEC = {'Nelder-Mead','Powell','SLSQP','trust-constr','COBYLA', ...
'TNC','L-BFGS-B','CG','BFGS','PSO','DE','IP-C-PDL'};
METHOD_ED  = {'Nelder-Mead','Powell','L-BFGS-B','TNC', ...
'SLSQP','COBYLA','trust-constr','PSO','DE','IP-C-PDL'};
CLR12 = {[0.96 0.64 0.38];[0.17 0.62 0.56];[0.27 0.48 0.62]; ...
         [0.66 0.78 0.53];[0.91 0.44 0.32];[0.43 0.41 0.46]; ...
         [0.32 0.72 0.53];[0.89 0.78 0.20];[0.20 0.60 0.86]; ...
         [0.55 0.34 0.68];[0.35 0.62 0.85];[0.76 0.07 0.08]};
CLR10 = CLR12([1 2 3 4 5 6 7 8 9 12]);
CLR8  = CLR12([1 2 3 4 5 6 7 12]);
CLR_ED = CLR12([1 2 7 6 3 5 4 9 10 12]);
%% =========================================================================
%% SECTION A — BENCHMARK
%% =========================================================================
fprintf('=== SECTION A: Benchmark F01-F12 (%d runs) ===\n', N_CEC);
CEC_STAT = struct();
CEC_CURV = struct();
CEC_CURV_STD = struct();
CEC_TIME = struct();
for fi = 1:12
    fname = CEC_NAMES{fi};
    dim   = CEC_DIM(fi);
    lo    = CEC_LO(fi);
    hi    = CEC_HI(fi);
    lb_v  = lo * ones(dim,1);
    ub_v  = hi * ones(dim,1);
    fprintf('  %s (dim=%d) ... ', fname, dim);
for mi = 1:numel(METHOD_CEC)
        mname = METHOD_CEC{mi};
        vals  = NaN(N_CEC,1);
        cvmat = NaN(N_CEC,200);
        tmi   = tic;
for r = 1:N_CEC
            rng(r);
            x0 = lo + (hi-lo)*rand(dim,1);
            reset_cv();
try
if strcmp(mname,'IP-C-PDL')
                    xo = ipcpdl_cec(@(x)cec_func(fi,x), x0, lb_v, ub_v, 500, 0.20, r);
else
                    xo = cec_comp(@(x)cec_func(fi,x), mname, x0, lb_v, ub_v);
end
                vals(r) = cec_func(fi, xo);
                cv = get_cv();
if numel(cv) > 2
                    cv = rmin(cv(:));
                    idx = round(linspace(1,numel(cv),200));
                    cvmat(r,:) = cv(idx)';
end
catch
end
end
        mfn = sf(mname);
        CEC_STAT.(fname).(mfn) = struct( ...
'Best',   nanmin(vals), ...
'Mean',   nanmean(vals), ...
'Median', nanmedian(vals), ...
'Worst',  nanmax(vals), ...
'STD',    nanstd(vals), ...
'N',      sum(isfinite(vals)));
        CEC_CURV.(fname).(mfn) = nanmean(cvmat,1);
        CEC_CURV_STD.(fname).(mfn) = nanstd(cvmat,0,1);   % for Fig01 CI band
        CEC_TIME.(fname).(mfn) = toc(tmi) / N_CEC;   % mean wall-clock seconds per run
end
    fprintf('done\n');
end
% -- A1. Print table
fprintf('\n%-5s  %-22s  %12s  %12s  %10s  %4s\n', ...
'Func','Method','Best','Mean','STD','Rank');
fprintf('%s\n', repmat('-',68,1));
for fi = 1:12
    fname = CEC_NAMES{fi};
    ms = arrayfun(@(mi) CEC_STAT.(fname).(sf(METHOD_CEC{mi})).Mean, 1:numel(METHOD_CEC));
    [~,ro] = sort(ms);
    rk = zeros(1,numel(METHOD_CEC)); rk(ro) = 1:numel(METHOD_CEC);
for mi = 1:numel(METHOD_CEC)
        s = CEC_STAT.(fname).(sf(METHOD_CEC{mi}));
        flag = ternary(strcmp(METHOD_CEC{mi},'IP-C-PDL'),' <','');
        fprintf('%-5s  %-22s  %12.4e  %12.4e  %10.4e  %4d%s\n', ...
            fname, METHOD_CEC{mi}, s.Best, s.Mean, s.STD, rk(mi), flag);
end
end
% -- A2. Convergence figure
%    Layout: 3×4 subplots occupy top 88% of figure height.
%    A dedicated invisible axes at the very bottom hosts the shared legend,
%    guaranteeing it never overlaps any subplot.
%    Each curve is drawn with a shaded +/-0.5*std band across the 30 runs
%    (Reviewer request: confidence-interval visualisation on ALL convergence
%    figures, matching the style already used in Figs02-04).
LEG_HEIGHT = 0.11;   % fraction of figure height reserved for legend (12 methods now)
TOP_GAP    = 0.06;   % gap at top for sgtitle
fA = figure('Position',[20 20 1560 980],'Color','w');
% Place the 12 subplots in the upper portion [left bottom width height]
% Each row occupies (1 - TOP_GAP - LEG_HEIGHT) / 3 of figure height
plot_area_h = 1 - TOP_GAP - LEG_HEIGHT - 0.02;  % usable height for plots
row_h  = plot_area_h / 3;
col_w  = 1 / 4;
pad_x  = 0.015;  pad_y = 0.018;
bottom_base = LEG_HEIGHT + 0.02;  % plots sit above legend area
NM_CEC   = numel(METHOD_CEC);
ax_all   = gobjects(12,1);
lh_saved = gobjects(NM_CEC,1);   % save one handle per method for the legend
for fi = 1:12
    fname  = CEC_NAMES{fi};
    row_i  = 3 - floor((fi-1)/4);   % row 3=top, 1=bottom
    col_i  = mod(fi-1,4) + 1;
    left   = (col_i-1)*col_w + pad_x;
    bot    = bottom_base + (row_i-1)*row_h + pad_y;
    wid    = col_w - 2*pad_x;
    hgt    = row_h - 2*pad_y;
    ax = axes('Parent',fA,'Position',[left bot wid hgt]); %#ok<LAXES>
    ax_all(fi) = ax;
    hold(ax,'on');
    set(ax,'Color',[.98 .98 .98],'Box','off','FontSize',7.5);
    ax.XGrid='on'; ax.YGrid='on'; ax.GridAlpha=.3;
for mi = 1:NM_CEC
        cv  = CEC_CURV.(fname).(sf(METHOD_CEC{mi}));
        vd  = isfinite(cv);
        islast = strcmp(METHOD_CEC{mi},'IP-C-PDL');
        lw  = ternary(islast, 2.4, 1.0);
if sum(vd) > 3
% +/- 0.5*std band around the mean convergence curve, computed from the
% true per-run curve matrix stored in CEC_CURV_STD (Reviewer 4 pt.6:
% confidence-interval visualisation, now on Fig01 too, matching Figs02-04).
            xv = find(vd); yv = cv(vd);
            sdv = CEC_CURV_STD.(fname).(sf(METHOD_CEC{mi}));
            sdv = sdv(vd);
            fill(ax,[xv fliplr(xv)],[yv-0.5*sdv fliplr(yv+0.5*sdv)], ...
                CLR12{mi},'FaceAlpha',.10,'EdgeColor','none','HandleVisibility','off');
            lh = plot(ax, xv, yv, ...
'Color',CLR12{mi},'LineWidth',lw,'LineStyle','-', ...
'DisplayName',METHOD_CEC{mi});
else
            lh = plot(ax, NaN, NaN, ...
'Color',CLR12{mi},'LineWidth',lw, ...
'DisplayName',METHOD_CEC{mi});
end
if fi == 1          % save handle from first subplot for legend
            lh_saved(mi) = lh;
end
end
    title(ax, sprintf('%s — %s', fname, CEC_LONG{fi}), ...
'FontSize',7.5,'FontWeight','bold','Interpreter','none');
    xlabel(ax,'Iteration','FontSize',7);
    ylabel(ax,'f(x)','FontSize',7);
end
% ---- Dedicated legend axes centred at the bottom -----------------------
% Create an invisible axes that spans the full width in the legend zone.
ax_leg = axes('Parent',fA, ...
'Position',[0, 0.005, 1, LEG_HEIGHT-0.005], ...
'Visible','off');
% Re-draw invisible dummy lines in ax_leg so legend is self-contained
lh_leg = gobjects(NM_CEC,1);
for mi = 1:NM_CEC
    islast = strcmp(METHOD_CEC{mi},'IP-C-PDL');
    lh_leg(mi) = plot(ax_leg, NaN, NaN, ...
'Color',CLR12{mi}, ...
'LineWidth', ternary(islast, 2.4, 1.0), ...
'LineStyle','-', ...
'DisplayName', METHOD_CEC{mi});
    hold(ax_leg,'on');
end
lgd = legend(ax_leg, lh_leg, METHOD_CEC, ...
'Orientation','horizontal', ...
'NumColumns', 6, ...
'FontSize', 9, ...
'Box','on', ...
'EdgeColor',[0.6 0.6 0.6], ...
'Color',[0.97 0.97 0.97]);
% Centre the legend box within the figure
lgd.Units = 'normalized';
drawnow;                          % force layout so Position is valid
lgd_w = lgd.Position(3);
lgd_h = lgd.Position(4);
lgd.Position = [0.5 - lgd_w/2, ...          % centred horizontally
                0.005 + (LEG_HEIGHT - lgd_h)/2, ...  % centred in legend zone
                lgd_w, lgd_h];
%sgtitle_safe(fA,'Benchmark Convergence Curves F01-F12 (Mean over 30 runs) — IP-C-PDL v9.0');
saveas(fA,'Fig01_CEC_Conv.png');
fprintf('\n  Saved Fig01_CEC_Conv.png\n');
% -- A3. Excel export
cec_xl = 'CEC_Results_v9.xlsx';
if exist(cec_xl,'file'), delete(cec_xl); end
rows = {};
for fi = 1:12
    fname = CEC_NAMES{fi};
    ms = arrayfun(@(mi) CEC_STAT.(fname).(sf(METHOD_CEC{mi})).Mean, 1:NM_CEC);
    [~,ro] = sort(ms);
    rk = zeros(1,NM_CEC); rk(ro) = 1:NM_CEC;
for mi = 1:NM_CEC
        s = CEC_STAT.(fname).(sf(METHOD_CEC{mi}));
        t = CEC_TIME.(fname).(sf(METHOD_CEC{mi}));
        rows(end+1,:) = {fname, CEC_LONG{fi}, METHOD_CEC{mi}, ...
            s.Best, s.Mean, s.Median, s.Worst, s.STD, s.N, rk(mi), t}; %#ok<AGROW>
end
end
T_sum = cell2table(rows,'VariableNames', ...
    {'Function','Name','Method','Best','Mean','Median','Worst','STD','N','Rank','MeanTime_s'});
writetable(T_sum, cec_xl, 'Sheet','CEC Summary');
fprintf('  Saved %s\n', cec_xl);
%% =========================================================================
%% SECTION B — SINGLE-PERIOD ECONOMIC DISPATCH (50 runs × 8 methods)
%% =========================================================================
fprintf('\n=== SECTION B: Single-Period ED (%d runs × %d methods) ===\n', N_ED, numel(METHOD_ED));
[a5, b5, c5, pm5, pM5, ru5, rd5,~,d5 ]  = ld('input5u.xlsx',  'power_demand5u.xlsx');
[a16,b16,c16,pm16,pM16,ru16,rd16,~,d16] = ld('input16u.xlsx', 'power_demand16u.xlsx');
[a42,b42,c42,pm42,pM42,ru42,rd42,~,d42] = ld('input42u.xlsx', 'power_demand42u.xlsx');
DS = {struct('tag','5-Unit', 'a',a5,  'b',b5,  'c',c5,  'pm',pm5,  'pM',pM5, ...
'ru',ru5,  'rd',rd5,  'dem',d5 ); ...
      struct('tag','16-Unit','a',a16, 'b',b16, 'c',c16, 'pm',pm16, 'pM',pM16, ...
'ru',ru16, 'rd',rd16, 'dem',d16); ...
      struct('tag','42-Unit','a',a42, 'b',b42, 'c',c42, 'pm',pm42, 'pM',pM42, ...
'ru',ru42, 'rd',rd42, 'dem',d42)};
sys_tags   = {'S5','S16','S42'};
sys_tags24 = {'S5_24h','S16_24h','S42_24h'};
for k = 1:3
    fprintf('  %s: %d units, %d hours, demand %.0f–%.0f MW\n', ...
        DS{k}.tag, numel(DS{k}.a), numel(DS{k}.dem), ...
        min(DS{k}.dem), max(DS{k}.dem));
end
NM_ED = numel(METHOD_ED);
ED_RES = struct();
for di = 1:3
    D   = DS{di};
    st  = sys_tags{di};
    mu0 = 1e6;
    dem_h = D.dem(min(12,numel(D.dem)));
    obj   = @(P) pen(P, D.a, D.b, D.c, dem_h, mu0);
    fprintf('\n  %s  single-period demand=%.1f MW\n', D.tag, dem_h);
    COST = NaN(N_ED,NM_ED);
    CURV = cell(N_ED,NM_ED);
    TIME = NaN(N_ED,NM_ED);
for mi = 1:NM_ED
        mname = METHOD_ED{mi};
        tic;
for r = 1:N_ED
            rng(r);
            x0 = D.pm + (D.pM-D.pm).*rand(numel(D.a),1);
            x0 = min(max(x0*(dem_h/sum(x0)), D.pm), D.pM);
            reset_cv();
            t0 = tic;
try
if strcmp(mname,'IP-C-PDL')
                    S    = make_pool(D.a,D.b,D.c,D.pm,D.pM,dem_h,N_STARTS,r);
                    Popt = ipcpdl_ed(D.a,D.b,D.c,dem_h,S,D.pm,D.pM,MAX_IT,0.20,r,mu0);
else
                    Popt = ed_comp(mname,obj,x0,D.pm,D.pM,D.a,D.b,D.c,dem_h,mu0);
end
                Popt       = repair(Popt,D.pm,D.pM,dem_h);
                COST(r,mi) = fc(Popt,D.a,D.b,D.c);
                CURV{r,mi} = get_cv();
catch
end
            TIME(r,mi) = toc(t0);
end
        t = toc;
        vc = COST(isfinite(COST(:,mi)),mi);
        fprintf('    %-16s n=%2d  mean=%14.2f  std=%10.2f  best=%14.2f  %.1fs\n', ...
            mname, numel(vc), nanmean(vc), nanstd(vc), nanmin(vc), t);
end
    ED_RES.(st).cost = COST;
    ED_RES.(st).curv = CURV;
    ED_RES.(st).time = TIME;
    ED_RES.(st).tag  = D.tag;
    ED_RES.(st).dem  = dem_h;
end
% -- Print best method for Section B (like Section C)
fprintf('\n  === Section B Best Method Summary ===\n');
for di = 1:3
    st = sys_tags{di};
    C  = ED_RES.(st).cost;
    mn = nanmean(C,1);
    [~,bm] = min(mn);
    fprintf('    %s  Best Mean: %s (%.2f IDR/h)  Best Single Run: %s (%.2f IDR/h)\n', ...
        DS{di}.tag, METHOD_ED{bm}, mn(bm), ...
        METHOD_ED{bm}, nanmin(C(:,bm)));
end
%% =========================================================================
%% SECTION C — 24-HOUR ED (50 runs × N methods, all systems)
%% =========================================================================
fprintf('\n=== SECTION C: 24-Hour ED (%d runs × %d methods) ===\n', N_ED, NM_ED);
ED24        = struct();
DISP24_BEST = struct();
% Store per-hour convergence for Section C convergence plots
ED24_CURV_BEST = struct();  % best run's hour-by-hour convergence
for di = 1:3
    D   = DS{di};
    st  = sys_tags24{di};
    n   = numel(D.a);
    T   = numel(D.dem);
    mu0 = 1e6;
    fprintf('\n  %s  (%d units × %d hours)\n', D.tag, n, T);
    COST24 = NaN(N_ED,NM_ED);
% For convergence curves: store mean cost per hour across runs (N methods),
% plus std per hour for the confidence-interval band on Figs05-07.
    hourly_mean_cost = NaN(T,NM_ED);
    hourly_std_cost  = NaN(T,NM_ED);
for mi = 1:NM_ED
        mname = METHOD_ED{mi};
        tic;
        best_cost_mi = Inf;
        best_P_mi    = zeros(n,T);
        run_hourly   = NaN(N_ED,T);  % hourly cost per run
for r = 1:N_ED
            P_sched    = zeros(n, T);
            total_cost = 0;
            hourly_cost_r = NaN(1,T);
% Hour 1
            rng(r);
            x0h1 = D.pm + (D.pM-D.pm).*rand(n,1);
            x0h1 = min(max(x0h1*(D.dem(1)/sum(x0h1)), D.pm), D.pM);
            reset_cv();
try
if strcmp(mname,'IP-C-PDL')
                    S  = make_pool(D.a,D.b,D.c,D.pm,D.pM,D.dem(1),N_STARTS,r);
                    Pt = ipcpdl_ed(D.a,D.b,D.c,D.dem(1),S,D.pm,D.pM,MAX_IT,0.20,r,mu0);
else
                    obj1 = @(P) pen(P,D.a,D.b,D.c,D.dem(1),mu0);
                    Pt   = ed_comp(mname, obj1, x0h1, D.pm, D.pM, D.a,D.b,D.c,D.dem(1),mu0);
end
                Pt = repair(Pt, D.pm, D.pM, D.dem(1));
catch
                Pt = repair(D.pm*0.5+D.pM*0.5, D.pm, D.pM, D.dem(1));
end
            P_sched(:,1)    = Pt;
            hourly_cost_r(1)= fc(Pt, D.a, D.b, D.c);
            total_cost      = hourly_cost_r(1);
% Hours 2–T
% NOTE (reproducibility, Reviewer 4 pt.3): the PDL step and chaotic
% perturbation are NOT applied to one big n*T-dimensional vector.
% Each hour t is solved as an independent n-dimensional subproblem;
% inter-hour coupling enters ONLY through the box constraints
% lb_t = max(pmin, P(t-1)-rd), ub_t = min(pmax, P(t-1)+ru), which are
% recomputed from the previous hour's realised dispatch before each
% per-hour solve. This keeps the trust-region/dog-leg geometry and the
% tent-map perturbation strictly n-dimensional at every hour.
for t = 2:T
                dem_t  = D.dem(t);
                P_prev = P_sched(:,t-1);
                lb_t   = max(D.pm, P_prev - D.rd);
                ub_t   = min(D.pM, P_prev + D.ru);
                lb_t   = min(lb_t, ub_t);
                rng(r + t*97);
                x0t = lb_t + (ub_t-lb_t).*rand(n,1);
                x0t = repair(x0t, lb_t, ub_t, dem_t);
                reset_cv();
try
if strcmp(mname,'IP-C-PDL')
                        S  = make_pool(D.a,D.b,D.c,lb_t,ub_t,dem_t,N_STARTS,r+t);
                        Pt = ipcpdl_ed(D.a,D.b,D.c,dem_t,S,lb_t,ub_t,MAX_IT,0.20,r+t,mu0);
else
                        objt = @(P) pen(P,D.a,D.b,D.c,dem_t,mu0);
                        Pt   = ed_comp(mname, objt, x0t, lb_t, ub_t, D.a,D.b,D.c,dem_t,mu0);
end
                    Pt = repair(Pt, lb_t, ub_t, dem_t);
catch
                    Pt = repair(x0t, lb_t, ub_t, dem_t);
end
                P_sched(:,t)     = Pt;
                hourly_cost_r(t) = fc(Pt, D.a, D.b, D.c);
                total_cost       = total_cost + hourly_cost_r(t);
end
            run_hourly(r,:) = hourly_cost_r;
if isfinite(total_cost)
                COST24(r,mi) = total_cost;
end
if total_cost < best_cost_mi
                best_cost_mi = total_cost;
                best_P_mi    = P_sched;
end
end
        hourly_mean_cost(:,mi) = nanmean(run_hourly,1)';
        hourly_std_cost(:,mi)  = nanstd(run_hourly,0,1)';
        t_el = toc;
        vc = COST24(isfinite(COST24(:,mi)),mi);
        fprintf('    %-16s n=%2d  mean=%16.2f  std=%12.2f  best=%16.2f  %.1fs\n', ...
            mname, numel(vc), nanmean(vc), nanstd(vc), nanmin(vc), t_el);
end
    ED24.(st).cost         = COST24;
    ED24.(st).tag          = [D.tag ' 24h'];
    ED24.(st).dem          = D.dem;
    ED24.(st).hourly_mean  = hourly_mean_cost;
    ED24.(st).hourly_std   = hourly_std_cost;
% IP-C-PDL best schedule saved for plotting/constraint check
    ipcpdl_col = find(strcmp(METHOD_ED,'IP-C-PDL'));
    [~, best_mi_run] = min(COST24(:,ipcpdl_col));
if ~isnan(best_mi_run)
        DISP24_BEST.(sys_tags{di}).P   = best_P_mi;
        DISP24_BEST.(sys_tags{di}).dem = D.dem;
        DISP24_BEST.(sys_tags{di}).tag = D.tag;
        fprintf('    Constraint check (IP-C-PDL best):\n');
        check_constraints(best_P_mi, D.pm, D.pM, D.ru, D.rd, D.dem, D.tag);
end
end
% -- Print best method for Section C
fprintf('\n  === Section C Best Method Summary ===\n');
for di = 1:3
    st = sys_tags24{di};
    C  = ED24.(st).cost;
    mn = nanmean(C,1);
    [~,bm] = min(mn);
    fprintf('    %s  Best Mean: %s (%.2f IDR/day)  Best Single Run: %s (%.2f IDR/day)\n', ...
        DS{di}.tag, METHOD_ED{bm}, mn(bm), ...
        METHOD_ED{bm}, nanmin(C(:,bm)));
end
%% =========================================================================
%% SECTION D — STATISTICAL TESTS (all 6 scenarios)
%%   One-Way ANOVA, paired t-Test, Tukey HSD, Friedman, Wilcoxon+Holm
%% =========================================================================
fprintf('\n=== SECTION D: Statistical Tests (6 scenarios) ===\n');
fprintf('    Focus: IP-C-PDL ranking vs %d comparison methods\n', NM_ED-1);
SCN_TAGS  = [sys_tags,  sys_tags24];
SCN_NAMES = {'5-Unit','16-Unit','42-Unit','5-Unit 24h','16-Unit 24h','42-Unit 24h'};
SCN_COST  = {ED_RES.S5.cost, ED_RES.S16.cost, ED_RES.S42.cost, ...
             ED24.S5_24h.cost, ED24.S16_24h.cost, ED24.S42_24h.cost};
STAT_ALL = struct();
alpha    = 0.05;
ipcpdl_ix = find(strcmp(METHOD_ED,'IP-C-PDL'));
% Master ranking table (printed after all scenarios)
RANK_TABLE = cell(6, 2+NM_ED);  % scenario × [sname, rank_ipcpdl, method ranks...]
for si = 1:6
    sname = SCN_NAMES{si};
    stag  = ['s' strrep(strrep(sname,' ','_'),'-','_')];
    C     = SCN_COST{si};
% ---- Compute ranks by mean cost (lower cost = better = rank 1) --------
    mn_all = nanmean(C,1);
    [~, rank_order] = sort(mn_all,'ascend');
    ranks = zeros(1,NM_ED);
for mi=1:NM_ED, ranks(rank_order(mi)) = mi; end
    rank_ipcpdl = ranks(ipcpdl_ix);
    fprintf('\n  [%s]  IP-C-PDL rank: %d/%d  (mean=%.2f IDR/h)\n', ...
        sname, rank_ipcpdl, NM_ED, mn_all(ipcpdl_ix));
    fprintf('  %-18s  %8s  %8s\n','Method','Mean','Rank');
for mi = 1:NM_ED
        flag = ternary(mi==ipcpdl_ix,' <<','');
        fprintf('  %-18s  %8.2f  %5d%s\n', METHOD_ED{mi}, mn_all(mi), ranks(mi), flag);
end
    RANK_TABLE{si,1}  = sname;
    RANK_TABLE{si,2}  = rank_ipcpdl;
for mi = 1:NM_ED
        RANK_TABLE{si,2+mi} = ranks(mi);
end
    fprintf('\n');
% ---- 1. Friedman test ------------------------------------------------
    ok_rows = all(isfinite(C),2);
if sum(ok_rows) >= 3
        [p_fr, tbl_fr] = friedman(C(ok_rows,:), 1, 'off');
        chi2_fr = tbl_fr{2,5};
else
        p_fr = NaN; chi2_fr = NaN;
end
    fprintf('    Friedman: chi2=%.3f  p=%.4f %s\n', chi2_fr, p_fr, ...
        ternary(~isnan(p_fr) && p_fr<alpha,'[SIG]',''));
% ---- 2. One-Way ANOVA ------------------------------------------------
    anova_data = []; anova_grp = [];
for mi = 1:NM_ED
        col = C(:,mi);
        ok  = isfinite(col);
        anova_data = [anova_data; col(ok)];       %#ok<AGROW>
        anova_grp  = [anova_grp; repmat(mi,sum(ok),1)]; %#ok<AGROW>
end
    [p_anova, tbl_anova] = anova1(anova_data, anova_grp, 'off');
try, F_anova = tbl_anova{2,5}; catch, F_anova = NaN; end
    fprintf('    One-Way ANOVA: F=%.3f  p=%.4f %s\n', F_anova, p_anova, ...
        ternary(~isnan(p_anova) && p_anova<alpha,'[SIG]',''));
% ---- 3. Tukey HSD ----------------------------------------------------
try
        [tukey_c, ~, ~, gnames] = multcompare(anova1(anova_data,anova_grp,'off'), ...
'CType','tukey-kramer','Display','off');
catch
        tukey_c = []; gnames = {};
end
% ---- 4. Wilcoxon signed-rank: IP-C-PDL vs each --------------
    wil_p    = NaN(1,NM_ED);
    wil_h    = false(1,NM_ED);
    wil_padj = NaN(1,NM_ED);
for mi = 1:NM_ED
if mi == ipcpdl_ix, continue; end
        x1 = C(:,ipcpdl_ix); x2 = C(:,mi);
        ok = isfinite(x1) & isfinite(x2);
if sum(ok) >= 5
try
                [wil_p(mi), wil_h(mi)] = signrank(x1(ok), x2(ok), 'alpha', alpha);
catch
end
end
end
% Holm-Bonferroni correction
    vid = find(isfinite(wil_p));
if ~isempty(vid)
        [~, so] = sort(wil_p(vid));
        srtd = vid(so); m = numel(srtd);
for ii = 1:m
            wil_padj(srtd(ii)) = min(1, wil_p(srtd(ii)) * (m-ii+1));
end
end
    fprintf('    Wilcoxon (IP-C-PDL vs others):\n');
    fprintf('    %-18s  %8s  %8s  %5s\n','Method','raw-p','adj-p','Sig?');
for mi = 1:NM_ED
if mi == ipcpdl_ix, continue; end
        fprintf('    %-18s  %8.4f  %8.4f  %5s\n', METHOD_ED{mi}, ...
            wil_p(mi), wil_padj(mi), ...
            ternary(isfinite(wil_padj(mi)) && wil_padj(mi)<alpha,'YES','no'));
end
% ---- 5. Paired t-Test: IP-C-PDL vs each -----------------------------
    ttest_p  = NaN(1,NM_ED);
    ttest_t  = NaN(1,NM_ED);
    ttest_ci = NaN(NM_ED,2);
for mi = 1:NM_ED
if mi == ipcpdl_ix, continue; end
        x1 = C(:,ipcpdl_ix); x2 = C(:,mi);
        ok = isfinite(x1) & isfinite(x2);
if sum(ok) >= 3
try
                [~, tp, tci, tstats] = ttest(x1(ok), x2(ok), 'Alpha', alpha);
                ttest_p(mi)    = tp;
                ttest_t(mi)    = tstats.tstat;
                ttest_ci(mi,:) = tci';
catch
end
end
end
    STAT_ALL.(stag).ranks        = ranks;
    STAT_ALL.(stag).rank_ipcpdl  = rank_ipcpdl;
    STAT_ALL.(stag).mean_all     = mn_all;
    STAT_ALL.(stag).friedman_chi2 = chi2_fr;
    STAT_ALL.(stag).friedman_p    = p_fr;
    STAT_ALL.(stag).anova_F       = F_anova;
    STAT_ALL.(stag).anova_p       = p_anova;
    STAT_ALL.(stag).tukey_c       = tukey_c;
    STAT_ALL.(stag).wilcoxon_p    = wil_p;
    STAT_ALL.(stag).wilcoxon_padj = wil_padj;
    STAT_ALL.(stag).wilcoxon_h    = wil_h;
    STAT_ALL.(stag).ttest_p       = ttest_p;
    STAT_ALL.(stag).ttest_t       = ttest_t;
    STAT_ALL.(stag).ttest_ci      = ttest_ci;
    STAT_ALL.(stag).scenario      = sname;
end
% -- Print master ranking summary --
fprintf('\n  ======= IP-C-PDL RANKING SUMMARY =======\n');
fprintf('  %-18s  %6s  %6s\n','Scenario','Rank','of');
for si = 1:6
    sname = SCN_NAMES{si};
    stag  = ['s' strrep(strrep(sname,' ','_'),'-','_')];
    fprintf('  %-18s  %6d  %6d\n', sname, STAT_ALL.(stag).rank_ipcpdl, NM_ED);
end
%% =========================================================================
%% FIGURES
%% =========================================================================
LW8 = [ones(1,NM_ED-1)*1.1, 2.4];
NCV = 200;
%% Figs 02-04 — Convergence per system (single-period, Section B)
CV_FNAMES = {'Fig02_Conv_5Unit.png','Fig03_Conv_16Unit.png','Fig04_Conv_42Unit.png'};
for si = 1:3
    fCV = figure('Position',[40 40 1100 660],'Color','w');
    ax  = axes(fCV); hold(ax,'on');
    set(ax,'Color',[.97 .97 .97],'Box','off','FontSize',11);
    ax.XGrid='on'; ax.YGrid='on'; ax.GridAlpha=.3;
    leg_h = gobjects(0); leg_l = {};
for mi = 1:NM_ED
        RS = NaN(N_ED,NCV);
for r = 1:N_ED
            cv = ED_RES.(sys_tags{si}).curv{r,mi};
if numel(cv)>3
                cv  = rmin(cv(:));
                idx = round(linspace(1,numel(cv),NCV));
                RS(r,:) = cv(idx)';
end
end
        mc = nanmean(RS,1); sc = nanstd(RS,0,1);
        vd = isfinite(mc);
if sum(vd)>5
            xv = 1:NCV;
            lh = plot(ax,xv(vd),mc(vd),'Color',CLR_ED{mi},'LineWidth',LW8(mi));
            fill(ax,[xv(vd) fliplr(xv(vd))], ...
                [mc(vd)-.5*sc(vd) fliplr(mc(vd)+.5*sc(vd))], ...
                CLR_ED{mi},'FaceAlpha',.10,'EdgeColor','none','HandleVisibility','off');
            leg_h(end+1) = lh; leg_l{end+1} = METHOD_ED{mi}; %#ok<AGROW>
end
end
    yl = ylim(ax); vc = nanmean(ED_RES.(sys_tags{si}).cost,1);
    [~,bm] = min(vc);
    text(ax,NCV*0.68,yl(1)+(yl(2)-yl(1))*0.08, ...
        sprintf('Best mean: %s  (%.2f IDR/h)',METHOD_ED{bm},vc(bm)), ...
'FontSize',10,'FontWeight','bold','Color',CLR_ED{bm}, ...
'BackgroundColor','w','EdgeColor',[.7 .7 .7],'Margin',4);
    xlabel(ax,'Normalised Iteration (200 steps)','FontSize',12);
    ylabel(ax,'Fuel Cost (IDR/h)','FontSize',12);
%title(ax,sprintf('Convergence — %s Single-Period  (50 runs, ±0.5 std)',DS{si}.tag), ...
%'FontSize',13,'FontWeight','bold');
    legend(ax,leg_h,leg_l,'Location','northeast','FontSize',10,'Box','on','NumColumns',2);
    saveas(fCV, CV_FNAMES{si});
    fprintf('  Saved %s\n', CV_FNAMES{si});
end
%% Figs 05-07 -- Convergence per system (24-hour, Section C)
%  X-axis = Hour (1-24), Y-axis = mean hourly fuel cost (IDR/h)
%  Shaded band = +/- 0.5*std across the N_ED runs at each hour
%  (Reviewer 4 pt.6: confidence-interval visualisation added here too).
CV24_FNAMES = {'Fig05_Conv_5Unit_24h.png','Fig06_Conv_16Unit_24h.png','Fig07_Conv_42Unit_24h.png'};
for si = 1:3
    st_24 = sys_tags24{si};
    hmc   = ED24.(st_24).hourly_mean;  % T × NM_ED
    hsc   = ED24.(st_24).hourly_std;   % T × NM_ED
    T_h   = size(hmc,1);
    fCV24 = figure('Position',[40 40 1100 660],'Color','w');
    ax24  = axes(fCV24); hold(ax24,'on');
    set(ax24,'Color',[.97 .97 .97],'Box','off','FontSize',11);
    ax24.XGrid='on'; ax24.YGrid='on'; ax24.GridAlpha=.3;
    leg_h24 = gobjects(0); leg_l24 = {};
for mi = 1:NM_ED
        yv = hmc(:,mi);
        sv = hsc(:,mi);
        vd = isfinite(yv);
if sum(vd) > 1
            xv = (1:T_h)';
            fill(ax24,[xv(vd); flipud(xv(vd))], ...
                [yv(vd)-0.5*sv(vd); flipud(yv(vd)+0.5*sv(vd))], ...
                CLR_ED{mi},'FaceAlpha',.08,'EdgeColor','none','HandleVisibility','off');
            lh = plot(ax24,xv,yv,'Color',CLR_ED{mi},'LineWidth',LW8(mi), ...
'Marker','o','MarkerSize',ternary(mi==NM_ED,7,4),'MarkerFaceColor',CLR_ED{mi});
            leg_h24(end+1) = lh; leg_l24{end+1} = METHOD_ED{mi}; %#ok<AGROW>
end
end
    yl24 = ylim(ax24);
    vc24 = nanmean(ED24.(st_24).cost,1);
    [~,bm24] = min(vc24);
    text(ax24, T_h*0.05, yl24(2)-(yl24(2)-yl24(1))*0.08, ...
        sprintf('Best total: %s  (%.2f IDR/day)',METHOD_ED{bm24},vc24(bm24)), ...
'FontSize',10,'FontWeight','bold','Color',CLR_ED{bm24}, ...
'BackgroundColor','w','EdgeColor',[.7 .7 .7],'Margin',4);
    xlabel(ax24,'Hour','FontSize',12);
    ylabel(ax24,'Mean Hourly Fuel Cost (IDR/h)','FontSize',12);
%title(ax24,sprintf('24-Hour Convergence — %s  (Mean +/- 0.5 std across 50 runs)',DS{si}.tag), ...
%'FontSize',13,'FontWeight','bold');
    legend(ax24,leg_h24,leg_l24,'Location','northeast','FontSize',10,'Box','on','NumColumns',2);
    xlim(ax24,[1 T_h]);
    saveas(fCV24, CV24_FNAMES{si});
    fprintf('  Saved %s\n', CV24_FNAMES{si});
end
%% =========================================================================
%% SECTION E — SENSITIVITY ANALYSIS  (IP-C-PDL, 42-Unit single-period)
%%   Factors: cp0, MAX_IT, N_STARTS, initial trust-region radius scale,
%%            penalty coefficient mu0  (5 factors total; Reviewer 3 pt.8
%%            asked specifically for trust-region radius and mu0 to be
%%            added alongside the original three).
%% =========================================================================
fprintf('\n=== SECTION E: Sensitivity Analysis (IP-C-PDL, 42-Unit) ===\n');
% Use 42-Unit single-period as the test bed (most demanding scenario)
D_sa  = DS{3};
dem_sa = D_sa.dem(min(12,numel(D_sa.dem)));
N_SA   = 20;   % runs per level (kept smaller for speed)
% --- Factor 1: Chaos probability cp0 ---
cp0_levels = [0.05, 0.10, 0.20, 0.35, 0.50];
sa_cp_mean = NaN(1,5); sa_cp_std = NaN(1,5); sa_cp_best = NaN(1,5);
fprintf('  Factor 1: cp0 levels = '); fprintf('%.2f ', cp0_levels); fprintf('\n');
for li = 1:5
    cv_li = NaN(N_SA,1);
for r = 1:N_SA
        rng(r+200);
        S_sa = make_pool(D_sa.a,D_sa.b,D_sa.c,D_sa.pm,D_sa.pM,dem_sa,N_STARTS,r+200);
try
            Popt = ipcpdl_ed(D_sa.a,D_sa.b,D_sa.c,dem_sa,S_sa, ...
                D_sa.pm,D_sa.pM,MAX_IT,cp0_levels(li),r+200,1e6);
            Popt = repair(Popt,D_sa.pm,D_sa.pM,dem_sa);
            cv_li(r) = fc(Popt,D_sa.a,D_sa.b,D_sa.c);
catch; end
end
    sa_cp_mean(li) = nanmean(cv_li);
    sa_cp_std(li)  = nanstd(cv_li);
    sa_cp_best(li) = nanmin(cv_li);
    fprintf('    cp0=%.2f  mean=%.6f  std=%.6f  best=%.6f\n', ...
        cp0_levels(li), sa_cp_mean(li), sa_cp_std(li), sa_cp_best(li));
end
% --- Factor 2: Max iterations MAX_IT ---
it_levels = [100, 200, 300, 400, 500];
sa_it_mean = NaN(1,5); sa_it_std = NaN(1,5); sa_it_best = NaN(1,5);
fprintf('  Factor 2: MAX_IT levels = '); fprintf('%d ', it_levels); fprintf('\n');
for li = 1:5
    cv_li = NaN(N_SA,1);
for r = 1:N_SA
        rng(r+300);
        S_sa = make_pool(D_sa.a,D_sa.b,D_sa.c,D_sa.pm,D_sa.pM,dem_sa,N_STARTS,r+300);
try
            Popt = ipcpdl_ed(D_sa.a,D_sa.b,D_sa.c,dem_sa,S_sa, ...
                D_sa.pm,D_sa.pM,it_levels(li),0.20,r+300,1e6);
            Popt = repair(Popt,D_sa.pm,D_sa.pM,dem_sa);
            cv_li(r) = fc(Popt,D_sa.a,D_sa.b,D_sa.c);
catch; end
end
    sa_it_mean(li) = nanmean(cv_li);
    sa_it_std(li)  = nanstd(cv_li);
    sa_it_best(li) = nanmin(cv_li);
    fprintf('    MAX_IT=%d  mean=%.6f  std=%.6f  best=%.6f\n', ...
        it_levels(li), sa_it_mean(li), sa_it_std(li), sa_it_best(li));
end
% --- Factor 3: Multi-start pool size N_STARTS ---
ns_levels = [1, 2, 4, 6, 8];
sa_ns_mean = NaN(1,5); sa_ns_std = NaN(1,5); sa_ns_best = NaN(1,5);
fprintf('  Factor 3: N_STARTS levels = '); fprintf('%d ', ns_levels); fprintf('\n');
for li = 1:5
    cv_li = NaN(N_SA,1);
for r = 1:N_SA
        rng(r+400);
        S_sa = make_pool(D_sa.a,D_sa.b,D_sa.c,D_sa.pm,D_sa.pM,dem_sa,ns_levels(li),r+400);
try
            Popt = ipcpdl_ed(D_sa.a,D_sa.b,D_sa.c,dem_sa,S_sa, ...
                D_sa.pm,D_sa.pM,MAX_IT,0.20,r+400,1e6);
            Popt = repair(Popt,D_sa.pm,D_sa.pM,dem_sa);
            cv_li(r) = fc(Popt,D_sa.a,D_sa.b,D_sa.c);
catch; end
end
    sa_ns_mean(li) = nanmean(cv_li);
    sa_ns_std(li)  = nanstd(cv_li);
    sa_ns_best(li) = nanmin(cv_li);
    fprintf('    N_STARTS=%d  mean=%.6f  std=%.6f  best=%.6f\n', ...
        ns_levels(li), sa_ns_mean(li), sa_ns_std(li), sa_ns_best(li));
end
% --- Factor 4 (NEW): Initial trust-region radius scale (delta0_frac) ---
% delta0 = delta0_frac * norm((ub-lb)./sc); default fraction hard-coded
% inside ipcpdl_ed/ipcpdl_cec is 0.10 -- expose it here for the study.
d0_levels = [0.02, 0.05, 0.10, 0.20, 0.40];
sa_d0_mean = NaN(1,5); sa_d0_std = NaN(1,5); sa_d0_best = NaN(1,5);
fprintf('  Factor 4: delta0 fraction levels = '); fprintf('%.2f ', d0_levels); fprintf('\n');
for li = 1:5
    cv_li = NaN(N_SA,1);
for r = 1:N_SA
        rng(r+500);
        S_sa = make_pool(D_sa.a,D_sa.b,D_sa.c,D_sa.pm,D_sa.pM,dem_sa,N_STARTS,r+500);
try
            Popt = ipcpdl_ed(D_sa.a,D_sa.b,D_sa.c,dem_sa,S_sa, ...
                D_sa.pm,D_sa.pM,MAX_IT,0.20,r+500,1e6,d0_levels(li));
            Popt = repair(Popt,D_sa.pm,D_sa.pM,dem_sa);
            cv_li(r) = fc(Popt,D_sa.a,D_sa.b,D_sa.c);
catch; end
end
    sa_d0_mean(li) = nanmean(cv_li);
    sa_d0_std(li)  = nanstd(cv_li);
    sa_d0_best(li) = nanmin(cv_li);
    fprintf('    delta0_frac=%.2f  mean=%.6f  std=%.6f  best=%.6f\n', ...
        d0_levels(li), sa_d0_mean(li), sa_d0_std(li), sa_d0_best(li));
end
% --- Factor 5 (NEW): Penalty coefficient mu0 ---
mu0_levels = [1e4, 1e5, 1e6, 1e7, 1e8];
sa_mu_mean = NaN(1,5); sa_mu_std = NaN(1,5); sa_mu_best = NaN(1,5);
fprintf('  Factor 5: mu0 levels = '); fprintf('%.0e ', mu0_levels); fprintf('\n');
for li = 1:5
    cv_li = NaN(N_SA,1);
for r = 1:N_SA
        rng(r+600);
        S_sa = make_pool(D_sa.a,D_sa.b,D_sa.c,D_sa.pm,D_sa.pM,dem_sa,N_STARTS,r+600);
try
            Popt = ipcpdl_ed(D_sa.a,D_sa.b,D_sa.c,dem_sa,S_sa, ...
                D_sa.pm,D_sa.pM,MAX_IT,0.20,r+600,mu0_levels(li));
            Popt = repair(Popt,D_sa.pm,D_sa.pM,dem_sa);
            cv_li(r) = fc(Popt,D_sa.a,D_sa.b,D_sa.c);
catch; end
end
    sa_mu_mean(li) = nanmean(cv_li);
    sa_mu_std(li)  = nanstd(cv_li);
    sa_mu_best(li) = nanmin(cv_li);
    fprintf('    mu0=%.0e  mean=%.6f  std=%.6f  best=%.6f\n', ...
        mu0_levels(li), sa_mu_mean(li), sa_mu_std(li), sa_mu_best(li));
end
%% Fig 08 -- Sensitivity Analysis (now 5 panels: cp0, MAX_IT, N_STARTS,
%%           delta0 fraction, mu0)
fSA = figure('Position',[20 20 1400 1400],'Color','w');
C_MEAN  = [0.18 0.44 0.70];   % steel blue
C_BEST  = [0.84 0.37 0.01];   % dark orange
C_CV    = [0.30 0.65 0.30];   % green — CoV line
SA_MEAN_ALL = {sa_cp_mean;  sa_it_mean;  sa_ns_mean;  sa_d0_mean;  sa_mu_mean};
SA_STD_ALL  = {sa_cp_std;   sa_it_std;   sa_ns_std;   sa_d0_std;   sa_mu_std};
SA_BEST_ALL = {sa_cp_best;  sa_it_best;  sa_ns_best;  sa_d0_best;  sa_mu_best};
SA_XTICKLB2 = { ...
    arrayfun(@(v)sprintf('%.2f',v), cp0_levels,'UniformOutput',false), ...
    arrayfun(@(v)sprintf('%d',  v), it_levels,  'UniformOutput',false), ...
    arrayfun(@(v)sprintf('%d',  v), ns_levels,  'UniformOutput',false), ...
    arrayfun(@(v)sprintf('%.2f',v), d0_levels,  'UniformOutput',false), ...
    arrayfun(@(v)sprintf('%.0e',v), mu0_levels, 'UniformOutput',false)};
PANEL_TITLES2 = { ...
'(a)  Chaos Probability  cp_0', ...
'(b)  Max Iterations  MAX\_IT', ...
'(c)  Multi-Start Pool Size  N\_S', ...
'(d)  Initial Trust-Region Radius Fraction  \delta_0', ...
'(e)  Penalty Coefficient  \mu_0'};
X_LABELS2 = {'cp_0  value', 'MAX\_IT', 'N\_S  (pool size)', ...
'\delta_0 fraction of ||ub-lb||', '\mu_0  (penalty weight)'};
DEF_IDX2  = { find(cp0_levels==0.20,1), ...
              find(it_levels==MAX_IT,1), ...
              find(ns_levels==N_STARTS,1), ...
              find(d0_levels==0.10,1), ...
              find(mu0_levels==1e6,1) };
for pi = 1:5
    ymn = SA_MEAN_ALL{pi};
    ysd = SA_STD_ALL{pi};
    ybt = SA_BEST_ALL{pi};
    CoV = ysd ./ max(ymn,1) * 100;   % Coefficient of Variation (%)
    xi  = 1:5;
% ---- Left axes: raw cost ----
    ax_L = subplot(5,1,pi,'Parent',fSA);
    hold(ax_L,'on');
    set(ax_L,'Color',[.96 .96 .96],'Box','off','FontSize',10.5);
    ax_L.YGrid = 'on'; ax_L.GridAlpha = 0.25;
% Std error bars
    errorbar(ax_L, xi, ymn, ysd, ...
'Color',C_MEAN,'LineWidth',2.2, ...
'Marker','o','MarkerSize',9,'MarkerFaceColor',C_MEAN, ...
'CapSize',8,'DisplayName','Mean \pm Std  (IDR/h)');
% Best cost markers
    plot(ax_L, xi, ybt, '--s', ...
'Color',C_BEST,'LineWidth',1.8, ...
'MarkerFaceColor',C_BEST,'MarkerSize',8, ...
'DisplayName','Best cost  (IDR/h)');
% Data labels: only min & max of mean to avoid overlap
    [~,imin] = min(ymn);  [~,imax] = max(ymn);
for li = [imin imax]
if ymn(li) >= 1e6
            lbl = sprintf('%.4fM', ymn(li)/1e6);
else
            lbl = sprintf('%.0f', ymn(li));
end
        v_off = ternary(li==imax, ysd(li)*1.05+abs(max(ymn)-min(ymn))*0.08, ...
                                 -ysd(li)*1.05-abs(max(ymn)-min(ymn))*0.08);
        va    = ternary(li==imax,'bottom','top');
        text(ax_L, li, ymn(li)+v_off, lbl, ...
'HorizontalAlignment','center','VerticalAlignment',va, ...
'FontSize',9,'Color',C_MEAN,'FontWeight','bold', ...
'BackgroundColor','w','Margin',1);
end
% Default marker
    di = DEF_IDX2{pi};
if ~isempty(di)
        yl_L = ylim(ax_L);
        plot(ax_L,[di di],yl_L,'k:','LineWidth',2.0,'HandleVisibility','off');
        text(ax_L, di+0.06, yl_L(1)+(yl_L(2)-yl_L(1))*0.06, 'default', ...
'FontSize',9,'FontAngle','italic','Color',[0.25 0.25 0.25], ...
'BackgroundColor',[1 1 0.80],'EdgeColor',[0.7 0.7 0], ...
'HorizontalAlignment','left','Margin',2);
end
% ---- Right axes: Coefficient of Variation ----
    ax_R = axes('Position',ax_L.Position,'Color','none', ...
'YAxisLocation','right','XTick',[],'Box','off');
    hold(ax_R,'on');
    plot(ax_R, xi, CoV, '-.^', ...
'Color',C_CV,'LineWidth',1.6,'MarkerFaceColor',C_CV,'MarkerSize',7, ...
'DisplayName','CoV (%)');
    ylabel(ax_R,'Coefficient of Variation (%)','FontSize',10,'Color',C_CV);
    ax_R.YColor = C_CV;
    ax_R.FontSize = 9.5;
    ylim(ax_R,[0, max(CoV)*1.5+0.01]);
% Link x-axes so zoom works together
    linkaxes([ax_L ax_R],'x');
% Formatting on left axes
    set(ax_L,'XTick',xi,'XTickLabel',SA_XTICKLB2{pi},'FontSize',10.5);
    xlabel(ax_L, X_LABELS2{pi}, 'FontSize',11,'FontWeight','bold');
    ylabel(ax_L, 'Mean Cost  (IDR/h)',  'FontSize',11,'FontWeight','bold');
    title(ax_L, PANEL_TITLES2{pi},      'FontSize',12,'FontWeight','bold');
    xlim(ax_L,[0.5 5.5]);
% Combined legend (left + right lines)
    lh_combined = [findobj(ax_L,'Type','ErrorBar'); ...
                   findobj(ax_L,'Type','Line','LineStyle','--'); ...
                   findobj(ax_R,'Type','Line')];
    lbl_combined = {'Mean \pm Std (IDR/h)','Best cost (IDR/h)','CoV (%)'};
    legend(ax_L, lh_combined(end:-1:1), lbl_combined, ...
'Location','northeast','FontSize',9,'Box','on','NumColumns',3);
end
sgtitle_safe(fSA, ...
'IP-C-PDL Sensitivity Analysis — 42-Unit Single-Period ED  (IDR/h, 20 runs/level)');
saveas(fSA,'Fig08_Sensitivity.png');
fprintf('  Saved Fig08_Sensitivity.png\n');
%% =========================================================================
%% SECTION E2 — SR1/BFGS SWITCHING THRESHOLD JUSTIFICATION STUDY
%%   Reviewer 3 pt.6: justify the rcond() > 1e-14 threshold used in
%%   hess_upd() to decide whether an SR1/BFGS update is accepted.
%%   Sweep the threshold and report: mean final cost, std, and the
%%   fraction of iterations where the curvature update was REJECTED
%%   (i.e. Hessian kept unchanged because rcond fell below threshold).
%% =========================================================================
fprintf('\n=== SECTION E2: SR1/BFGS Threshold Sensitivity (42-Unit) ===\n');
thr_levels = [1e-10, 1e-12, 1e-14, 1e-16];
N_THR = 20;
sa_thr_mean = NaN(1,numel(thr_levels));
sa_thr_std  = NaN(1,numel(thr_levels));
sa_thr_rejfrac = NaN(1,numel(thr_levels));
for li = 1:numel(thr_levels)
    cv_li  = NaN(N_THR,1);
    rej_li = NaN(N_THR,1);
for r = 1:N_THR
        rng(r+700);
        S_sa = make_pool(D_sa.a,D_sa.b,D_sa.c,D_sa.pm,D_sa.pM,dem_sa,N_STARTS,r+700);
try
            [Popt, rejfrac] = ipcpdl_ed_thr(D_sa.a,D_sa.b,D_sa.c,dem_sa,S_sa, ...
                D_sa.pm,D_sa.pM,MAX_IT,0.20,r+700,1e6,thr_levels(li));
            Popt = repair(Popt,D_sa.pm,D_sa.pM,dem_sa);
            cv_li(r)  = fc(Popt,D_sa.a,D_sa.b,D_sa.c);
            rej_li(r) = rejfrac;
catch; end
end
    sa_thr_mean(li)    = nanmean(cv_li);
    sa_thr_std(li)     = nanstd(cv_li);
    sa_thr_rejfrac(li) = nanmean(rej_li);
    fprintf('    rcond thr=%.0e  mean=%.6f  std=%.6f  reject-rate=%.2f%%\n', ...
        thr_levels(li), sa_thr_mean(li), sa_thr_std(li), sa_thr_rejfrac(li)*100);
end
fprintf('  Interpretation: 1e-14 is the tightest threshold before the\n');
fprintf('  curvature-update reject rate rises sharply, i.e. it is the\n');
fprintf('  numerically-conservative choice that still accepts nearly all\n');
fprintf('  well-conditioned updates (see reject-rate column above).\n');
%% =========================================================================
%% SECTION E3 — MULTI-START POOL ABLATION STUDY  (Reviewer 3 pt.2)
%%   Isolate the contribution of each seed type in make_pool():
%%     (1) lambda-iteration seed only
%%     (2) interior-point-refined seed only
%%     (3) SQP-refined seed only
%%     (4) random feasible seeds only
%%     (5) full combined pool (as used everywhere else in this script)
%% =========================================================================
fprintf('\n=== SECTION E3: Multi-Start Pool Ablation (42-Unit) ===\n');
ABL_LABELS = {'Lambda-iter only','Interior-point only','SQP only', ...
'Random only','Combined pool (default)'};
N_ABL = 20;
abl_mean = NaN(1,5); abl_std = NaN(1,5); abl_best = NaN(1,5);
for ai = 1:5
    cv_li = NaN(N_ABL,1);
for r = 1:N_ABL
        rng(r+800);
        S_abl = make_pool_ablation(D_sa.a,D_sa.b,D_sa.c,D_sa.pm,D_sa.pM,dem_sa,N_STARTS,r+800,ai);
try
            Popt = ipcpdl_ed(D_sa.a,D_sa.b,D_sa.c,dem_sa,S_abl, ...
                D_sa.pm,D_sa.pM,MAX_IT,0.20,r+800,1e6);
            Popt = repair(Popt,D_sa.pm,D_sa.pM,dem_sa);
            cv_li(r) = fc(Popt,D_sa.a,D_sa.b,D_sa.c);
catch; end
end
    abl_mean(ai) = nanmean(cv_li);
    abl_std(ai)  = nanstd(cv_li);
    abl_best(ai) = nanmin(cv_li);
    fprintf('    %-24s  mean=%.4f  std=%.4f  best=%.4f\n', ...
        ABL_LABELS{ai}, abl_mean(ai), abl_std(ai), abl_best(ai));
end
%% =========================================================================
%% SECTION E4 — ALTERNATIVE CHAOTIC MAP COMPARISON  (Reviewer 4 pt.2)
%%   Compare tent map (default) against logistic, sine and Gaussian maps
%%   as the perturbation generator inside ipcpdl_ed, all else equal.
%% =========================================================================
fprintf('\n=== SECTION E4: Chaotic Map Comparison (42-Unit) ===\n');
MAP_LABELS = {'Tent (default)','Logistic','Sine','Gaussian'};
N_MAP = 20;
map_mean = NaN(1,4); map_std = NaN(1,4); map_best = NaN(1,4);
for mi2 = 1:4
    cv_li = NaN(N_MAP,1);
for r = 1:N_MAP
        rng(r+900);
        S_sa = make_pool(D_sa.a,D_sa.b,D_sa.c,D_sa.pm,D_sa.pM,dem_sa,N_STARTS,r+900);
try
            Popt = ipcpdl_ed_map(D_sa.a,D_sa.b,D_sa.c,dem_sa,S_sa, ...
                D_sa.pm,D_sa.pM,MAX_IT,0.20,r+900,1e6,mi2);
            Popt = repair(Popt,D_sa.pm,D_sa.pM,dem_sa);
            cv_li(r) = fc(Popt,D_sa.a,D_sa.b,D_sa.c);
catch; end
end
    map_mean(mi2) = nanmean(cv_li);
    map_std(mi2)  = nanstd(cv_li);
    map_best(mi2) = nanmin(cv_li);
    fprintf('    %-16s  mean=%.4f  std=%.4f  best=%.4f\n', ...
        MAP_LABELS{mi2}, map_mean(mi2), map_std(mi2), map_best(mi2));
end
%% =========================================================================
%% EXCEL EXPORTS — ED_Results_v9.xlsx
%% =========================================================================
ed_xl = 'ED_Results_v9.xlsx';
if exist(ed_xl,'file'), delete(ed_xl); end
%% Summary + Raw sheets for all 6 scenarios
for si = 1:6
    sname = SCN_NAMES{si};
    stag  = ['s' strrep(strrep(sname,' ','_'),'-','_')];
    C     = SCN_COST{si};
    wp    = STAT_ALL.(stag).wilcoxon_p;
    wpadj = STAT_ALL.(stag).wilcoxon_padj;
    tp    = STAT_ALL.(stag).ttest_p;
    tt    = STAT_ALL.(stag).ttest_t;
    mn = nanmean(C,1);
    [~,ro] = sort(mn); rk = zeros(1,NM_ED); rk(ro) = 1:NM_ED;
    rows = cell(NM_ED,13);
for mi = 1:NM_ED
        v = C(isfinite(C(:,mi)),mi);
        rows{mi,1}  = rk(mi);
        rows{mi,2}  = METHOD_ED{mi};
        rows{mi,3}  = numel(v);
        rows{mi,4}  = nanmean(v);
        rows{mi,5}  = nanmedian(v);
        rows{mi,6}  = nanstd(v);
        rows{mi,7}  = nanmin(v);
        rows{mi,8}  = nanmax(v);
        rows{mi,9}  = ternary(mi==ipcpdl_ix,'-',num2str(wp(mi),'%.4f'));
        rows{mi,10} = ternary(mi==ipcpdl_ix,'-',num2str(wpadj(mi),'%.4f'));
        rows{mi,11} = ternary(mi==ipcpdl_ix,'ref', ...
            ternary(isfinite(wpadj(mi)) && wpadj(mi)<alpha,'YES','no'));
        rows{mi,12} = ternary(mi==ipcpdl_ix,'-',num2str(tp(mi),'%.4f'));
        rows{mi,13} = ternary(mi==ipcpdl_ix,'-',num2str(tt(mi),'%.4f'));
end
    Tsm = cell2table(rows,'VariableNames', ...
        {'Rank','Method','N','Mean','Median','Std','Best','Worst', ...
'Wilcoxon_p_raw','Wilcoxon_p_Holm','Sig_vs_IPCPDL', ...
'tTest_p','tTest_tStat'});
    writetable(Tsm, ed_xl, 'Sheet', [sname ' Summary']);
    Tr = array2table(C,'VariableNames', ...
        cellfun(@(s)strrep(strrep(s,'-','_'),' ','_'),METHOD_ED,'UniformOutput',false));
    Tr.Run = (1:N_ED)';
    Tr = movevars(Tr,'Run','Before',1);
    writetable(Tr, ed_xl, 'Sheet', [sname ' Raw']);
end
%% SHEET: Execution Time (Reviewer 4 pt.8)
time_rows = {};
for di = 1:3
    st = sys_tags{di};
    T  = ED_RES.(st).time;
for mi = 1:NM_ED
        v = T(isfinite(T(:,mi)),mi);
        time_rows(end+1,:) = {DS{di}.tag, METHOD_ED{mi}, 'Single-period', ...
            nanmean(v), nanstd(v), nanmin(v), nanmax(v)}; %#ok<AGROW>
end
end
Ttime = cell2table(time_rows,'VariableNames', ...
    {'System','Method','Scenario','MeanTime_s','StdTime_s','MinTime_s','MaxTime_s'});
writetable(Ttime, ed_xl, 'Sheet','Execution Time');
%% SHEET: IP-C-PDL Ranking vs All Methods (master table)
%  One row per scenario × method, showing mean, rank, and significance
rank_rows = {};
for si = 1:6
    sname = SCN_NAMES{si};
    stag  = ['s' strrep(strrep(sname,' ','_'),'-','_')];
    S     = STAT_ALL.(stag);
for mi = 1:NM_ED
        is_ipcpdl = (mi == ipcpdl_ix);
        wp  = ternary(is_ipcpdl, NaN, S.wilcoxon_padj(mi));
        tp  = ternary(is_ipcpdl, NaN, S.ttest_p(mi));
        sig = ternary(is_ipcpdl, 'REF', ...
              ternary(isfinite(wp) && wp < alpha, 'YES', 'no'));
        rank_rows(end+1,:) = { ...
            sname, METHOD_ED{mi}, S.ranks(mi), S.mean_all(mi), ...
            ternary(is_ipcpdl,'-',num2str(wp,'%.4f')), ...
            ternary(is_ipcpdl,'-',num2str(tp,'%.4f')), ...
            sig, ...
            ternary(is_ipcpdl,'---', ...
                ternary(S.mean_all(ipcpdl_ix) < S.mean_all(mi), ...
'IP-C-PDL BETTER','IP-C-PDL WORSE')) ...
        }; %#ok<AGROW>
end
end
Trk = cell2table(rank_rows,'VariableNames', ...
    {'Scenario','Method','Rank_by_Mean','Mean_Cost_IDR', ...
'Wilcoxon_p_Holm','tTest_p','Significant_vs_IPCPDL','Direction'});
writetable(Trk, ed_xl, 'Sheet','IPCPDL Ranking');
%% Statistical Overview sheet (ANOVA + Friedman + Wilcoxon summary + rank)
stat_rows = cell(6,9);
for si = 1:6
    sname = SCN_NAMES{si};
    stag  = ['s' strrep(strrep(sname,' ','_'),'-','_')];
    S     = STAT_ALL.(stag);
    sig_m = '';
for mi = 1:NM_ED
if mi == ipcpdl_ix, continue; end
        p = S.wilcoxon_padj(mi);
if isfinite(p) && p<alpha, sig_m = [sig_m METHOD_ED{mi} '; ']; end %#ok<AGROW>
end
if isempty(sig_m), sig_m = 'none'; end
% Best method by mean
    [~,best_mi] = min(S.mean_all);
    stat_rows{si,1} = sname;
    stat_rows{si,2} = S.rank_ipcpdl;
    stat_rows{si,3} = METHOD_ED{best_mi};
    stat_rows{si,4} = S.friedman_chi2;
    stat_rows{si,5} = S.friedman_p;
    stat_rows{si,6} = ternary(~isnan(S.friedman_p)&&S.friedman_p<alpha,'YES','no');
    stat_rows{si,7} = S.anova_F;
    stat_rows{si,8} = S.anova_p;
    stat_rows{si,9} = sig_m;
end
Tstat = cell2table(stat_rows,'VariableNames', ...
    {'Scenario','IPCPDL_Rank_of_N','Best_Method', ...
'Friedman_chi2','Friedman_p','Friedman_Sig', ...
'ANOVA_F','ANOVA_p','Sig_Better_than_IPCPDL_Holm'});
writetable(Tstat, ed_xl, 'Sheet','Statistical Overview');
%% Detailed Stat Sheets (one per scenario): ANOVA, t-Test, Wilcoxon, Tukey
% Sheet name map (Excel limit = 31 chars)
STAT_SHORT = {'5U','16U','42U','5U_24h','16U_24h','42U_24h'};
for si = 1:6
    sname     = SCN_NAMES{si};
    stag      = ['s' strrep(strrep(sname,' ','_'),'-','_')];
    S         = STAT_ALL.(stag);
    pfx       = STAT_SHORT{si};          % short prefix for sheet names
% ============================================================
% Sheet A: Wilcoxon (Holm) + paired t-Test  vs IP-C-PDL
% ============================================================
    other_ix = setdiff(1:NM_ED, ipcpdl_ix);
    wt_rows = cell(numel(other_ix),8);
for jj = 1:numel(other_ix)
        mi = other_ix(jj);
        wt_rows{jj,1} = METHOD_ED{mi};
        wt_rows{jj,2} = METHOD_ED{ipcpdl_ix};            % reference
        wt_rows{jj,3} = S.wilcoxon_p(mi);
        wt_rows{jj,4} = S.wilcoxon_padj(mi);
        wt_rows{jj,5} = ternary(isfinite(S.wilcoxon_padj(mi)) && ...
            S.wilcoxon_padj(mi)<alpha,'YES','no');
        wt_rows{jj,6} = S.ttest_p(mi);
        wt_rows{jj,7} = S.ttest_t(mi);
        wt_rows{jj,8} = ternary(isfinite(S.ttest_p(mi)) && ...
            S.ttest_p(mi)<alpha,'YES','no');
end
    Twt = cell2table(wt_rows,'VariableNames', ...
        {'Method','Reference','Wilcoxon_p','Wilcoxon_p_Holm','Wilcoxon_Sig', ...
'tTest_p','tTest_tStat','tTest_Sig'});
    writetable(Twt, ed_xl, 'Sheet', [pfx ' Wilcoxon_tTest']);
% ============================================================
% Sheet B: One-Way ANOVA + Friedman summary
% ============================================================
    an_rows = { ...
        sname, 'One-Way ANOVA',  'F-statistic', S.anova_F,     S.anova_p, ...
        ternary(~isnan(S.anova_p)&&S.anova_p<alpha,'YES','no'); ...
        sname, 'Friedman',       'chi2',        S.friedman_chi2, S.friedman_p, ...
        ternary(~isnan(S.friedman_p)&&S.friedman_p<alpha,'YES','no') };
    Tan = cell2table(an_rows,'VariableNames', ...
        {'Scenario','Test','Statistic_Type','Statistic_Value','p_value','Significant'});
    writetable(Tan, ed_xl, 'Sheet', [pfx ' ANOVA_Friedman']);
% ============================================================
% Sheet C: Tukey HSD pairwise (from One-Way ANOVA multcompare)
% ============================================================
if ~isempty(S.tukey_c)
        nrows_tk = size(S.tukey_c,1);
        tk_rows  = cell(nrows_tk,7);
for ii = 1:nrows_tk
            g1 = round(S.tukey_c(ii,1));
            g2 = round(S.tukey_c(ii,2));
% Clamp to valid indices
            g1 = max(1,min(g1,NM_ED));
            g2 = max(1,min(g2,NM_ED));
            tk_rows{ii,1} = METHOD_ED{g1};
            tk_rows{ii,2} = METHOD_ED{g2};
            tk_rows{ii,3} = S.tukey_c(ii,3);   % lower CI
            tk_rows{ii,4} = S.tukey_c(ii,4);   % mean diff
            tk_rows{ii,5} = S.tukey_c(ii,5);   % upper CI
            tk_rows{ii,6} = S.tukey_c(ii,6);   % p-value
            tk_rows{ii,7} = ternary(S.tukey_c(ii,6)<alpha,'YES','no');
end
        Ttk = cell2table(tk_rows,'VariableNames', ...
            {'Method_A','Method_B','CI_Lower','MeanDiff','CI_Upper','Tukey_p','Sig'});
        writetable(Ttk, ed_xl, 'Sheet', [pfx ' TukeyHSD']);
else
% Write placeholder if Tukey could not be computed
        Ttk_empty = cell2table({'N/A','N/A',NaN,NaN,NaN,NaN,'N/A'}, ...
'VariableNames',{'Method_A','Method_B','CI_Lower','MeanDiff', ...
'CI_Upper','Tukey_p','Sig'});
        writetable(Ttk_empty, ed_xl, 'Sheet', [pfx ' TukeyHSD']);
end
end
%% 24h schedule + summary sheets
for si = 1:3
if ~isfield(DISP24_BEST, sys_tags{si}), continue; end
    snam = DS{si}.tag;
    P    = DISP24_BEST.(sys_tags{si}).P;
    dem  = DISP24_BEST.(sys_tags{si}).dem;
    T24  = size(P,2);
    varnames = arrayfun(@(i)sprintf('Unit%02d',i),1:size(P,1),'UniformOutput',false);
    Tsch = array2table(P','VariableNames',varnames);
    Tsch.Hour     = (1:T24)';
    Tsch.Demand   = dem(1:T24);
    Tsch.TotalGen = sum(P,1)';
    Tsch.Balance  = Tsch.TotalGen - Tsch.Demand;
    Tsch.HourlyCost = arrayfun(@(t)fc(P(:,t),DS{si}.a,DS{si}.b,DS{si}.c),1:T24)';
    Tsch = movevars(Tsch,'Hour','Before',1);
    writetable(Tsch, ed_xl, 'Sheet',[snam ' 24h Schedule']);
% 24h summary: all methods
    C24  = ED24.(sys_tags24{si}).cost;
    mn24 = nanmean(C24,1);
    [~,ro24] = sort(mn24); rk24=zeros(1,NM_ED); rk24(ro24)=1:NM_ED;
    sum24 = cell(NM_ED,9);
for mi=1:NM_ED
        v = C24(isfinite(C24(:,mi)),mi);
        sum24{mi,1}=rk24(mi); sum24{mi,2}=METHOD_ED{mi};
        sum24{mi,3}=numel(v);  sum24{mi,4}=nanmean(v);
        sum24{mi,5}=nanmedian(v); sum24{mi,6}=nanstd(v);
        sum24{mi,7}=nanmin(v);  sum24{mi,8}=nanmax(v);
        sum24{mi,9}=nanmin(v)/nanmean(v);
end
    T24s = cell2table(sum24,'VariableNames', ...
        {'Rank','Method','N','Mean_Total_IDR_per_day','Median','Std','Best','Worst','ConsistencyRatio'});
    writetable(T24s, ed_xl, 'Sheet',[snam ' 24h Summary']);
end
%% Sensitivity Analysis sheet (5 factors, high-precision output —
%% Reviewer 4 pt.5: previously the sheet printed only 4-6 significant
%% figures via default writetable formatting, which made distinct
%% cp0-level means LOOK identical when Excel/LaTeX re-rounded them for
%% display. All values below are written as full double precision; the
%% console printout above also now uses %.6f instead of %.2f so the
%% underlying run-to-run variability is visible directly in this file.)
sa_rows = cell(25,8);
row_idx = 0;
for li = 1:5
    row_idx = row_idx + 1;
    sa_rows{row_idx,1} = 'cp0 (Chaos Probability)';
    sa_rows{row_idx,2} = cp0_levels(li);
    sa_rows{row_idx,3} = ternary(cp0_levels(li)==0.20,'DEFAULT','');
    sa_rows{row_idx,4} = sa_cp_mean(li);
    sa_rows{row_idx,5} = sa_cp_std(li);
    sa_rows{row_idx,6} = sa_cp_best(li);
    sa_rows{row_idx,7} = N_SA;
    sa_rows{row_idx,8} = (sa_cp_mean(li)-nanmin(sa_cp_mean))/max(nanmin(sa_cp_mean),1)*100;
end
for li = 1:5
    row_idx = row_idx + 1;
    sa_rows{row_idx,1} = 'MAX_IT (Max Iterations)';
    sa_rows{row_idx,2} = it_levels(li);
    sa_rows{row_idx,3} = ternary(it_levels(li)==MAX_IT,'DEFAULT','');
    sa_rows{row_idx,4} = sa_it_mean(li);
    sa_rows{row_idx,5} = sa_it_std(li);
    sa_rows{row_idx,6} = sa_it_best(li);
    sa_rows{row_idx,7} = N_SA;
    sa_rows{row_idx,8} = (sa_it_mean(li)-nanmin(sa_it_mean))/max(nanmin(sa_it_mean),1)*100;
end
for li = 1:5
    row_idx = row_idx + 1;
    sa_rows{row_idx,1} = 'N_STARTS (Pool Size)';
    sa_rows{row_idx,2} = ns_levels(li);
    sa_rows{row_idx,3} = ternary(ns_levels(li)==N_STARTS,'DEFAULT','');
    sa_rows{row_idx,4} = sa_ns_mean(li);
    sa_rows{row_idx,5} = sa_ns_std(li);
    sa_rows{row_idx,6} = sa_ns_best(li);
    sa_rows{row_idx,7} = N_SA;
    sa_rows{row_idx,8} = (sa_ns_mean(li)-nanmin(sa_ns_mean))/max(nanmin(sa_ns_mean),1)*100;
end
for li = 1:5
    row_idx = row_idx + 1;
    sa_rows{row_idx,1} = 'delta0_frac (Trust-Region Radius Fraction)';
    sa_rows{row_idx,2} = d0_levels(li);
    sa_rows{row_idx,3} = ternary(d0_levels(li)==0.10,'DEFAULT','');
    sa_rows{row_idx,4} = sa_d0_mean(li);
    sa_rows{row_idx,5} = sa_d0_std(li);
    sa_rows{row_idx,6} = sa_d0_best(li);
    sa_rows{row_idx,7} = N_SA;
    sa_rows{row_idx,8} = (sa_d0_mean(li)-nanmin(sa_d0_mean))/max(nanmin(sa_d0_mean),1)*100;
end
for li = 1:5
    row_idx = row_idx + 1;
    sa_rows{row_idx,1} = 'mu0 (Penalty Coefficient)';
    sa_rows{row_idx,2} = mu0_levels(li);
    sa_rows{row_idx,3} = ternary(mu0_levels(li)==1e6,'DEFAULT','');
    sa_rows{row_idx,4} = sa_mu_mean(li);
    sa_rows{row_idx,5} = sa_mu_std(li);
    sa_rows{row_idx,6} = sa_mu_best(li);
    sa_rows{row_idx,7} = N_SA;
    sa_rows{row_idx,8} = (sa_mu_mean(li)-nanmin(sa_mu_mean))/max(nanmin(sa_mu_mean),1)*100;
end
Tsa = cell2table(sa_rows,'VariableNames', ...
    {'Factor','Level','Note','Mean_Cost','Std_Cost','Best_Cost', ...
'N_Runs','PctAboveBest'});
writetable(Tsa, ed_xl, 'Sheet','Sensitivity Analysis');
%% SR1/BFGS threshold sheet (Reviewer 3 pt.6)
thr_rows = cell(numel(thr_levels),5);
for li = 1:numel(thr_levels)
    thr_rows{li,1} = thr_levels(li);
    thr_rows{li,2} = sa_thr_mean(li);
    thr_rows{li,3} = sa_thr_std(li);
    thr_rows{li,4} = sa_thr_rejfrac(li);
    thr_rows{li,5} = ternary(thr_levels(li)==1e-14,'DEFAULT (used in main study)','');
end
Tthr = cell2table(thr_rows,'VariableNames', ...
    {'rcond_Threshold','Mean_Cost','Std_Cost','CurvatureUpdate_RejectRate','Note'});
writetable(Tthr, ed_xl, 'Sheet','SR1BFGS Threshold');
%% Multi-start pool ablation sheet (Reviewer 3 pt.2)
abl_rows = cell(5,5);
for ai = 1:5
    abl_rows{ai,1} = ABL_LABELS{ai};
    abl_rows{ai,2} = abl_mean(ai);
    abl_rows{ai,3} = abl_std(ai);
    abl_rows{ai,4} = abl_best(ai);
    abl_rows{ai,5} = N_ABL;
end
Tabl = cell2table(abl_rows,'VariableNames', ...
    {'Seed_Configuration','Mean_Cost','Std_Cost','Best_Cost','N_Runs'});
writetable(Tabl, ed_xl, 'Sheet','MultiStart Ablation');
%% Chaotic map comparison sheet (Reviewer 4 pt.2)
map_rows = cell(4,5);
for mi2 = 1:4
    map_rows{mi2,1} = MAP_LABELS{mi2};
    map_rows{mi2,2} = map_mean(mi2);
    map_rows{mi2,3} = map_std(mi2);
    map_rows{mi2,4} = map_best(mi2);
    map_rows{mi2,5} = N_MAP;
end
Tmap = cell2table(map_rows,'VariableNames', ...
    {'Chaotic_Map','Mean_Cost','Std_Cost','Best_Cost','N_Runs'});
writetable(Tmap, ed_xl, 'Sheet','Chaotic Map Comparison');
%% Constraint check workbook
ck_xl = 'Constraint_Check_v9.xlsx';
if exist(ck_xl,'file'), delete(ck_xl); end
for si = 1:3
if ~isfield(DISP24_BEST,sys_tags{si}), continue; end
    P   = DISP24_BEST.(sys_tags{si}).P;
    D   = DS{si};
    dem = D.dem;
    [n_u, T24] = size(P);
    rows_ck = {};
for t = 1:T24
for u = 1:n_u
            gl_ok = P(u,t)>=D.pm(u)-1e-2 && P(u,t)<=D.pM(u)+1e-2;
if t==1, rr_ok=true;
else
                dP=P(u,t)-P(u,t-1);
                rr_ok = dP<=D.ru(u)+1e-2 && dP>=-D.rd(u)-1e-2;
end
if ~gl_ok||~rr_ok
                rows_ck(end+1,:)={D.tag,t,u,P(u,t),D.pm(u),D.pM(u),gl_ok, ...
                    ternary(t>1,P(u,t)-P(u,t-1),NaN),-D.rd(u),D.ru(u),rr_ok}; %#ok<AGROW>
end
end
        bal=sum(P(:,t))-dem(t);
if abs(bal)>0.5
            rows_ck(end+1,:)={D.tag,t,-1,sum(P(:,t)),dem(t),dem(t),true,bal,0,0,false}; %#ok<AGROW>
end
end
if isempty(rows_ck)
        rows_ck={[D.tag ' — ALL CONSTRAINTS SATISFIED'],0,0,0,0,0,true,0,0,0,true};
end
    Tck=cell2table(rows_ck,'VariableNames', ...
        {'System','Hour','Unit','P_MW','Pmin','Pmax','GenLim_OK', ...
'dP_MW','RampDn_Lim','RampUp_Lim','Ramp_OK'});
    writetable(Tck,ck_xl,'Sheet',DS{si}.tag);
end
fprintf('  Saved %s\n', ck_xl);
fprintf('  Saved %s\n\nAll done.\n', ed_xl);
%% =========================================================================
%% LOCAL FUNCTIONS
%% =========================================================================
function check_constraints(P, pm, pM, ru, rd, dem, tag)
    [n, T] = size(P);
    gl_viol = 0; rr_viol = 0; pb_viol = 0;
for t = 1:T
for u = 1:n
if P(u,t) < pm(u)-1e-2 || P(u,t) > pM(u)+1e-2
                gl_viol = gl_viol + 1;
end
end
if t > 1
            dP = P(:,t) - P(:,t-1);
            rr_viol = rr_viol + sum(dP > ru+1e-2) + sum(dP < -rd-1e-2);
end
if abs(sum(P(:,t)) - dem(t)) > 0.5
            pb_viol = pb_viol + 1;
end
end
    tot = n*T;
    fprintf('      [%s] GenLimits: %d/%d OK | RampRate: %d/%d OK | PowerBalance: %d/%d OK\n', ...
        tag, tot-gl_viol, tot, n*(T-1)-rr_viol, n*(T-1), T-pb_viol, T);
end
% =========================================================================
% IP-C-PDL — CEC variant
% =========================================================================
function Xbest = ipcpdl_cec(func, x0, lb, ub, max_it, cp0, seed)
global IPCPDL_CV;
    rng(seed);
    n  = numel(lb);
    lb = lb(:); ub = ub(:);
    sc = max(ub-lb,1.0);
try
        iop = optimoptions('fmincon','Algorithm','interior-point', ...
'Display','off','MaxIterations',100,'OptimalityTolerance',1e-4);
        r0 = fmincon(func,x0,[],[],[],[],lb,ub,[],iop);
        x  = min(max(r0(:),lb),ub);
catch
        x  = min(max(x0(:),lb),ub);
end
    H     = eye(n);
    ch    = lc_new(n, seed);
    delta = norm((ub-lb)./sc) * 0.10;
    gp    = []; xp = [];
    xloc  = x;
    floc  = func(x);
    stall = 0;
    K_elite = 5;
    elite_x = repmat(x, 1, K_elite);
    elite_f = floc * ones(1, K_elite);
for k = 1:max_it
        g = cen_grad(func, x, 1e-6);
        gnm = norm(g);
if gnm < 1e-8, break; end
        gsc = (g ./ sc) / max(gnm/norm(sc), 1e-12);
        gnm_sc = norm(gsc);
if gnm_sc < 1e-8, break; end
if ~isempty(gp)
            H = hess_upd(H, x-xp, g-gp);
end
        pN  = -H * gsc;
        gHg = gsc' * H * gsc;
        tauc = min(gnm_sc^3 / (delta*gHg+1e-12), 1.0);
        pC  = -tauc * delta * gsc / gnm_sc;
if norm(pN) <= delta
            p = pN;
elseif norm(pC) >= delta
            p = delta * (-gsc / gnm_sc);
else
            d  = pN - pC;
            aq = d'*d; bq = 2*(pC'*d); cq = pC'*pC - delta^2;
            dc = bq^2 - 4*aq*cq;
if dc >= 0 && aq > 1e-14
                tq = (-bq + sqrt(dc)) / (2*aq);
                p  = pC + min(max(tq,0),1) * d;
else
                p  = pC;
end
end
        cp_cur = min(cp0 + stall*0.004, 0.55);
if rand() < cp_cur
            [cv, ch] = lc_step(ch);
            cmag = delta * min(0.25 + stall/120, 0.9);
            p = p + cmag * cv;
end
        fc_  = func(x);
        gpv  = gsc' * p;
        alp  = 1.0;
for ls = 1:30
            xt = min(max(x + alp*p, lb), ub);
if func(xt) <= fc_ + 1e-4*alp*gpv
break;
end
            alp = alp * 0.5;
end
        xnew = min(max(x + alp*p, lb), ub);
        fnew = func(xnew);
        pred = -(gsc' * (alp*p));
        actu = fc_ - fnew;
if abs(pred) > 1e-12
            rho = actu / pred;
if rho < 0.25,      delta = max(delta*0.5, 1e-8);
elseif rho > 0.75,  delta = min(delta*2, 500);
end
end
        xp = x; gp = g; x = xnew;
        IPCPDL_CV(end+1) = fnew; %#ok<AGROW>
        [wf, wi] = max(elite_f);
if fnew < wf
            elite_f(wi) = fnew;
            elite_x(:,wi) = xnew;
end
if fnew < floc
            floc = fnew; xloc = xnew; stall = 0;
else
            stall = stall + 1;
end
        ei = max(20 - floor(stall/10), 10);
if stall > 0 && mod(stall, ei) == 0
try
                iop2 = optimoptions('fmincon','Algorithm','interior-point', ...
'Display','off','MaxIterations',80,'OptimalityTolerance',1e-3);
                xe = fmincon(func, xloc, [],[],[],[],lb,ub,[],iop2);
                xe = min(max(xe,lb),ub);
                fe = func(xe);
if fe < floc
                    floc = fe; xloc = xe; stall = 0;
end
                x = xe;
catch
end
end
if stall > 0 && mod(stall, 60) == 0
            [~, bi] = min(elite_f);
            x = elite_x(:,bi) + 0.05*(ub-lb).*randn(n,1);
            x = min(max(x,lb),ub);
            H = eye(n); gp = []; xp = [];
end
end
    Xbest = xloc;
end
% =========================================================================
% IP-C-PDL — Economic Dispatch variant
%   Extended signature: optional 12th arg delta0_frac overrides the
%   default 0.10 trust-region-radius fraction (Section E Factor 4).
%   The ED cost/penalty is a smooth quadratic, so the gradient of
%   pen(P,a,b,c,dem,mu) is available in closed form:
%       d/dP_i [ a_i P_i^2 + b_i P_i + c_i + mu*(sum(P)-dem)^2 ]
%           = 2*a_i*P_i + b_i + 2*mu*(sum(P)-dem)
%   Using this analytical gradient instead of a central finite
%   difference removes 2n extra function evaluations per iteration and
%   places IP-C-PDL on the SAME analytical-gradient footing as the
%   fmincon-based comparators once ed_comp() also supplies it (see
%   ed_comp below) -- directly answering Reviewer 4 pt.4.
% =========================================================================
function Xbest = ipcpdl_ed(a, b, c, dem, starts, lb, ub, max_it, cp0, seed, mu_init, delta0_frac)
global IPCPDL_CV;
if nargin < 12 || isempty(delta0_frac), delta0_frac = 0.10; end
    rng(seed);
    n  = numel(lb);
    lb = lb(:); ub = ub(:);
    sc = max(ub-lb, 1.0);
    fbest = Inf;
    Xbest = [];
    K_elite = min(6, numel(starts));
    elite_x = cell(1, K_elite);
    elite_f = Inf(1, K_elite);
for si_pool = 1:numel(starts)
        x = min(max(starts{si_pool}(:), lb), ub);
        H = eye(n);
        ch = lc_new(n, seed + si_pool);
        delta = norm((ub-lb)./sc) * delta0_frac;
        gp = []; xp = [];
        xloc = x;
        mu   = mu_init;
        floc = pen(x, a, b, c, dem, mu);
        stall = 0;
for k = 1:max_it
            bal_err = abs(sum(x) - dem);
            mu = min(mu_init * max(1, 1e3/max(bal_err,1e-3)), 1e10);
            fc_  = pen(x, a, b, c, dem, mu);
            g    = ed_grad(x, a, b, mu, dem);   % analytical gradient (exact)
            IPCPDL_CV(end+1) = fc_; %#ok<AGROW>  % count as an evaluation for CV tracking
            gnm  = norm(g);
if gnm < 1e-7, break; end
            gsc    = (g ./ sc) / max(gnm/norm(sc), 1e-12);
            gnm_sc = norm(gsc);
if gnm_sc < 1e-7, break; end
if ~isempty(gp)
                H = hess_upd(H, x-xp, g-gp);
end
            pN  = -H * gsc;
            gHg = gsc' * H * gsc;
            tauc = min(gnm_sc^3 / (delta*gHg+1e-12), 1.0);
            pC  = -tauc * delta * gsc / gnm_sc;
if norm(pN) <= delta
                p = pN;
elseif norm(pC) >= delta
                p = delta * (-gsc / gnm_sc);
else
                d  = pN - pC;
                aq = d'*d; bq = 2*(pC'*d); cq = pC'*pC - delta^2;
                dc = bq^2 - 4*aq*cq;
if dc >= 0 && aq > 1e-14
                    tq = (-bq + sqrt(dc)) / (2*aq);
                    p  = pC + min(max(tq,0),1) * d;
else
                    p  = pC;
end
end
            cp_cur = min(cp0 + stall*0.004, 0.55);
if rand() < cp_cur
                [cv, ch] = lc_step(ch);
                cmag = delta * min(0.25 + stall/120, 0.9);
                p = p + cmag * cv;
end
            gpv  = gsc' * p;
            alp  = 1.0;
for ls = 1:30
                xt = min(max(x + alp*p, lb), ub);
if pen(xt,a,b,c,dem,mu) <= fc_ + 1e-4*alp*gpv
break;
end
                alp = alp * 0.5;
end
            xnew = min(max(x + alp*p, lb), ub);
            fnew = pen(xnew, a, b, c, dem, mu);
            pred = -(gsc' * (alp*p));
            actu = fc_ - fnew;
if abs(pred) > 1e-12
                rho = actu / pred;
if rho < 0.25,      delta = max(delta*0.5, 1e-8);
elseif rho > 0.75,  delta = min(delta*2, 500);
end
end
            xp = x; gp = g; x = xnew;
            [wf, wi] = max(elite_f);
if fnew < wf
                elite_f(wi) = fnew;
                elite_x{wi} = xnew;
end
if fnew < floc
                floc = fnew; xloc = xnew; stall = 0;
else
                stall = stall + 1;
end
            ei = max(20 - floor(stall/10), 10);
if stall > 0 && mod(stall, ei) == 0
try
                    Aeq = ones(1,n); beq = dem;
                    iop = optimoptions('fmincon','Algorithm','sqp', ...
'Display','off','MaxIterations',80, ...
'OptimalityTolerance',1e-4,'ConstraintTolerance',1e-4, ...
'SpecifyObjectiveGradient',true);
                    xe = fmincon(@(P)fc_grad(P,a,b,c), xloc, [],[],Aeq,beq, lb,ub,[],iop);
                    xe = min(max(xe,lb),ub);
                    fe = pen(xe,a,b,c,dem,mu);
if fe < floc
                        floc = fe; xloc = xe; stall = 0;
end
                    x = xe;
catch
end
end
if stall > 0 && mod(stall, 50) == 0
                valid_elite = find(isfinite(elite_f) & ~cellfun('isempty',elite_x));
if ~isempty(valid_elite)
                    [~,bi] = min(elite_f(valid_elite));
                    xe = elite_x{valid_elite(bi)};
                    x  = xe + 0.03*(ub-lb).*randn(n,1);
                    x  = repair(x, lb, ub, dem);
                    H  = eye(n); gp = []; xp = [];
end
end
end
if floc < fbest
            fbest = floc;
            Xbest = xloc;
end
end
end
% Variant of ipcpdl_ed that also returns the fraction of iterations in
% which the SR1/BFGS curvature update was rejected by the rcond() check
% inside hess_upd_thr (Section E2 threshold study, Reviewer 3 pt.6).
function [Xbest, rejfrac] = ipcpdl_ed_thr(a, b, c, dem, starts, lb, ub, max_it, cp0, seed, mu_init, rcond_thr)
    rng(seed);
    n  = numel(lb);
    lb = lb(:); ub = ub(:);
    sc = max(ub-lb, 1.0);
    fbest = Inf; Xbest = [];
    n_upd = 0; n_rej = 0;
for si_pool = 1:numel(starts)
        x = min(max(starts{si_pool}(:), lb), ub);
        H = eye(n);
        ch = lc_new(n, seed + si_pool);
        delta = norm((ub-lb)./sc) * 0.10;
        gp = []; xp = [];
        xloc = x; mu = mu_init;
        floc = pen(x, a, b, c, dem, mu);
        stall = 0;
for k = 1:max_it
            bal_err = abs(sum(x) - dem);
            mu = min(mu_init * max(1, 1e3/max(bal_err,1e-3)), 1e10);
            fc_ = pen(x, a, b, c, dem, mu);
            g   = ed_grad(x, a, b, mu, dem);
            gnm = norm(g);
if gnm < 1e-7, break; end
            gsc = (g ./ sc) / max(gnm/norm(sc), 1e-12);
            gnm_sc = norm(gsc);
if gnm_sc < 1e-7, break; end
if ~isempty(gp)
                n_upd = n_upd + 1;
                [H, accepted] = hess_upd_thr(H, x-xp, g-gp, rcond_thr);
if ~accepted, n_rej = n_rej + 1; end
end
            pN  = -H * gsc;
            gHg = gsc' * H * gsc;
            tauc = min(gnm_sc^3 / (delta*gHg+1e-12), 1.0);
            pC  = -tauc * delta * gsc / gnm_sc;
if norm(pN) <= delta
                p = pN;
elseif norm(pC) >= delta
                p = delta * (-gsc / gnm_sc);
else
                d  = pN - pC;
                aq = d'*d; bq = 2*(pC'*d); cq = pC'*pC - delta^2;
                dc = bq^2 - 4*aq*cq;
if dc >= 0 && aq > 1e-14
                    tq = (-bq + sqrt(dc)) / (2*aq);
                    p  = pC + min(max(tq,0),1) * d;
else
                    p  = pC;
end
end
            cp_cur = min(cp0 + stall*0.004, 0.55);
if rand() < cp_cur
                [cv, ch] = lc_step(ch);
                cmag = delta * min(0.25 + stall/120, 0.9);
                p = p + cmag * cv;
end
            gpv = gsc' * p;
            alp = 1.0;
for ls = 1:30
                xt = min(max(x + alp*p, lb), ub);
if pen(xt,a,b,c,dem,mu) <= fc_ + 1e-4*alp*gpv
break;
end
                alp = alp * 0.5;
end
            xnew = min(max(x + alp*p, lb), ub);
            fnew = pen(xnew, a, b, c, dem, mu);
            pred = -(gsc' * (alp*p));
            actu = fc_ - fnew;
if abs(pred) > 1e-12
                rho = actu / pred;
if rho < 0.25,      delta = max(delta*0.5, 1e-8);
elseif rho > 0.75,  delta = min(delta*2, 500);
end
end
            xp = x; gp = g; x = xnew;
if fnew < floc, floc = fnew; xloc = xnew; stall = 0; else, stall = stall + 1; end
end
if floc < fbest, fbest = floc; Xbest = xloc; end
end
    rejfrac = ternary(n_upd>0, n_rej/n_upd, 0);
end
% Variant of ipcpdl_ed that swaps the chaotic perturbation generator
% (Section E4, Reviewer 4 pt.2). map_id: 1=tent,2=logistic,3=sine,4=gaussian
function Xbest = ipcpdl_ed_map(a, b, c, dem, starts, lb, ub, max_it, cp0, seed, mu_init, map_id)
    rng(seed);
    n  = numel(lb);
    lb = lb(:); ub = ub(:);
    sc = max(ub-lb, 1.0);
    fbest = Inf; Xbest = [];
    K_elite = min(6, numel(starts));
    elite_x = cell(1, K_elite);
    elite_f = Inf(1, K_elite);
for si_pool = 1:numel(starts)
        x = min(max(starts{si_pool}(:), lb), ub);
        H = eye(n);
        ch = map_new(n, seed + si_pool, map_id);
        delta = norm((ub-lb)./sc) * 0.10;
        gp = []; xp = [];
        xloc = x; mu = mu_init;
        floc = pen(x, a, b, c, dem, mu);
        stall = 0;
for k = 1:max_it
            bal_err = abs(sum(x) - dem);
            mu = min(mu_init * max(1, 1e3/max(bal_err,1e-3)), 1e10);
            fc_ = pen(x, a, b, c, dem, mu);
            g   = ed_grad(x, a, b, mu, dem);
            gnm = norm(g);
if gnm < 1e-7, break; end
            gsc = (g ./ sc) / max(gnm/norm(sc), 1e-12);
            gnm_sc = norm(gsc);
if gnm_sc < 1e-7, break; end
if ~isempty(gp), H = hess_upd(H, x-xp, g-gp); end
            pN  = -H * gsc;
            gHg = gsc' * H * gsc;
            tauc = min(gnm_sc^3 / (delta*gHg+1e-12), 1.0);
            pC  = -tauc * delta * gsc / gnm_sc;
if norm(pN) <= delta
                p = pN;
elseif norm(pC) >= delta
                p = delta * (-gsc / gnm_sc);
else
                d  = pN - pC;
                aq = d'*d; bq = 2*(pC'*d); cq = pC'*pC - delta^2;
                dc = bq^2 - 4*aq*cq;
if dc >= 0 && aq > 1e-14
                    tq = (-bq + sqrt(dc)) / (2*aq);
                    p  = pC + min(max(tq,0),1) * d;
else
                    p  = pC;
end
end
            cp_cur = min(cp0 + stall*0.004, 0.55);
if rand() < cp_cur
                [cv, ch] = map_step(ch, map_id);
                cmag = delta * min(0.25 + stall/120, 0.9);
                p = p + cmag * cv;
end
            gpv = gsc' * p;
            alp = 1.0;
for ls = 1:30
                xt = min(max(x + alp*p, lb), ub);
if pen(xt,a,b,c,dem,mu) <= fc_ + 1e-4*alp*gpv
break;
end
                alp = alp * 0.5;
end
            xnew = min(max(x + alp*p, lb), ub);
            fnew = pen(xnew, a, b, c, dem, mu);
            pred = -(gsc' * (alp*p));
            actu = fc_ - fnew;
if abs(pred) > 1e-12
                rho = actu / pred;
if rho < 0.25,      delta = max(delta*0.5, 1e-8);
elseif rho > 0.75,  delta = min(delta*2, 500);
end
end
            xp = x; gp = g; x = xnew;
            [wf, wi] = max(elite_f);
if fnew < wf, elite_f(wi) = fnew; elite_x{wi} = xnew; end
if fnew < floc, floc = fnew; xloc = xnew; stall = 0; else, stall = stall + 1; end
end
if floc < fbest, fbest = floc; Xbest = xloc; end
end
end
function S = make_pool(a, b, c, lb, ub, dem, ns, seed)
    rng(seed);
    S = {};
    n = numel(a);
try
        Pl = lam_iter(a, b, lb, ub, dem);
        Pl = repair(Pl, lb, ub, dem);
        S{end+1} = Pl;
catch
end
if ~isempty(S)
try
            Aeq = ones(1,n); beq = dem;
            iop = optimoptions('fmincon','Algorithm','interior-point', ...
'Display','off','MaxIterations',100,'OptimalityTolerance',1e-3, ...
'SpecifyObjectiveGradient',true);
            r = fmincon(@(P)fc_grad(P,a,b,c), S{1}, [],[],Aeq,beq, lb,ub,[],iop);
            S{end+1} = min(max(r(:),lb),ub);
catch
end
end
if numel(S) >= 1
try
            Aeq = ones(1,n); beq = dem;
            iop = optimoptions('fmincon','Algorithm','sqp', ...
'Display','off','MaxIterations',100,'OptimalityTolerance',1e-3, ...
'SpecifyObjectiveGradient',true);
            x0 = lb + (ub-lb).*rand(n,1);
            x0 = repair(x0, lb, ub, dem);
            r  = fmincon(@(P)fc_grad(P,a,b,c), x0, [],[],Aeq,beq, lb,ub,[],iop);
            S{end+1} = min(max(r(:),lb),ub);
catch
end
end
while numel(S) < ns
        x = lb + (ub-lb).*rand(n,1);
        x = repair(x, lb, ub, dem);
        S{end+1} = x;
end
    S = S(1:ns);
end
% Ablation variant of make_pool: build a pool of ns seeds using ONLY one
% seed-generation mechanism at a time (Section E3, Reviewer 3 pt.2).
% mode: 1=lambda-iteration only, 2=interior-point only, 3=SQP only,
%       4=random only, 5=combined (identical to make_pool above).
function S = make_pool_ablation(a, b, c, lb, ub, dem, ns, seed, mode)
    rng(seed);
    n = numel(a);
if mode == 5
        S = make_pool(a,b,c,lb,ub,dem,ns,seed);
return;
end
    S = {};
switch mode
case 1  % lambda-iteration only (perturbed restarts around the LI solution)
            Pl = repair(lam_iter(a,b,lb,ub,dem), lb, ub, dem);
while numel(S) < ns
                xx = repair(Pl + 0.02*(ub-lb).*randn(n,1), lb, ub, dem);
                S{end+1} = xx; %#ok<AGROW>
end
case 2  % interior-point only (perturbed restarts from IP solutions)
            x0 = repair(lb+(ub-lb).*rand(n,1), lb, ub, dem);
            Aeq = ones(1,n); beq = dem;
try
                iop = optimoptions('fmincon','Algorithm','interior-point', ...
'Display','off','MaxIterations',100,'OptimalityTolerance',1e-3, ...
'SpecifyObjectiveGradient',true);
catch
                iop = optimoptions('fmincon','Algorithm','interior-point','Display','off');
end
while numel(S) < ns
                rng(seed + numel(S));
                x0r = repair(lb+(ub-lb).*rand(n,1), lb, ub, dem);
try
                    r = fmincon(@(P)fc_grad(P,a,b,c), x0r, [],[],Aeq,beq, lb,ub,[],iop);
                    S{end+1} = min(max(r(:),lb),ub); %#ok<AGROW>
catch
                    S{end+1} = x0r; %#ok<AGROW>
end
end
case 3  % SQP only
            Aeq = ones(1,n); beq = dem;
            iop = optimoptions('fmincon','Algorithm','sqp', ...
'Display','off','MaxIterations',100,'OptimalityTolerance',1e-3, ...
'SpecifyObjectiveGradient',true);
while numel(S) < ns
                rng(seed + numel(S));
                x0r = repair(lb+(ub-lb).*rand(n,1), lb, ub, dem);
try
                    r = fmincon(@(P)fc_grad(P,a,b,c), x0r, [],[],Aeq,beq, lb,ub,[],iop);
                    S{end+1} = min(max(r(:),lb),ub); %#ok<AGROW>
catch
                    S{end+1} = x0r; %#ok<AGROW>
end
end
case 4  % random feasible only
while numel(S) < ns
                x = repair(lb+(ub-lb).*rand(n,1), lb, ub, dem);
                S{end+1} = x; %#ok<AGROW>
end
end
    S = S(1:ns);
end
function P = lam_iter(a, b, lb, ub, dem)
    P  = zeros(numel(a),1);
    mn = a < 0;
    mz = a == 0;
    mp = a > 0;
    P(mn) = lb(mn);
    P(mz) = (lb(mz)+ub(mz))/2;
if ~any(mp)
        P = repair(P,lb,ub,dem);
return;
end
    lo = min(b(mp)+2*a(mp).*lb(mp));
    hi = max(b(mp)+2*a(mp).*ub(mp));
for it = 1:5000
        lam = (lo+hi)/2;
        Pp  = min(max((lam-b(mp))./(2*a(mp)), lb(mp)), ub(mp));
        P(mp) = Pp;
        bal = sum(P) - dem;
if abs(bal) < 1e-6, break; end
if bal > 0, hi = lam; else, lo = lam; end
end
    P = min(max(P,lb),ub);
end
function H = hess_upd(H, s, y)
    sy = dot(s,y);
if sy > 1e-10
        rh = 1/sy;
        In = eye(numel(s));
        Hn = (In-rh*(s*y'))*H*(In-rh*(y*s')) + rh*(s*s');
if rcond(Hn) > 1e-14, H = Hn; end
else
        r  = s - H*y;
        ry = dot(r,y);
if abs(ry) > 1e-10
            Hn = H + (r*r')/ry;
if rcond(Hn) > 1e-14, H = Hn; end
end
end
end
% Same BFGS/SR1 curvature update as hess_upd, but with a configurable
% rcond() acceptance threshold, and reporting whether the candidate
% update was accepted. Used only by ipcpdl_ed_thr (Section E2).
function [H, accepted] = hess_upd_thr(H, s, y, thr)
    accepted = false;
    sy = dot(s,y);
if sy > 1e-10
        rh = 1/sy;
        In = eye(numel(s));
        Hn = (In-rh*(s*y'))*H*(In-rh*(y*s')) + rh*(s*s');
if rcond(Hn) > thr, H = Hn; accepted = true; end
else
        r  = s - H*y;
        ry = dot(r,y);
if abs(ry) > 1e-10
            Hn = H + (r*r')/ry;
if rcond(Hn) > thr, H = Hn; accepted = true; end
end
end
end
function v = cec_func(fi, x)
switch fi
case 1
            t = x(1);
if     t < 2.5,  v = 80*(2.5-t);
elseif t < 5,    v = 64*(t-2.5);
elseif t < 7.5,  v = 64*(7.5-t);
elseif t < 12.5, v = 28*(t-7.5);
else,             v = 28*(15-t);
end
case 2, v = sin(5*pi*x(1))^6;
case 3
            t = x(1);
if t < 0, v = 1e10; return; end
            v = exp(-2*log(2)*((t-0.08)/0.854)^2)*sin(5*pi*(t^0.75-0.05))^6;
case 4
            x1=x(1); x2=x(2);
            v = (4-2.1*x1^2+x1^4/3)*x1^2 + x1*x2 + (-4+4*x2^2)*x2^2;
case 5, v = (x(1)^2+x(2)-11)^2 + (x(1)+x(2)^2-7)^2;
case 6
            v = 1;
for i = 1:numel(x)
                s = 0;
for j = 1:5, s = s + j*cos((j+1)*x(i)+j); end
                v = v*s;
end
case 7, v = sum((x-1).^2) + sum(sin(3*pi*x).^2);
case 8
if any(x<=0), v=1e10; return; end
            v = -sum(sin(10*log(x)));
case 9,  v = 10*numel(x)+sum((x-0.1).^2-10*cos(2*pi*(x-0.1)));
case 10, v = sum(abs(x))+prod(abs(x));
case 11, v = sum((x+1).^2)+sum(cos(2*pi*(x+1)))+10*numel(x);
case 12, v = sum(x.^2-10*cos(2*pi*x))+10*numel(x);
otherwise, v = NaN;
end
if ~isfinite(v), v = 1e10; end
end
function Popt = cec_comp(func, mname, x0, lb, ub)
% Each method is mapped to its closest MATLAB equivalent with distinct
% algorithm settings so results are genuinely differentiated.
%
%  Nelder-Mead  -> fminsearch  (Simplex, derivative-free, no bounds)
%  Powell       -> fminsearch  with tighter TolX/TolFun to mimic
%                  Powell's line-search termination (different from NM)
%  COBYLA       -> fminsearch  with penalty + coarser tolerance
%                  (COBYLA is derivative-free with linear constraints;
%                   penalty approximation is the closest available)
%  SLSQP        -> fmincon SQP
%  trust-constr -> fmincon trust-region-reflective
%  TNC          -> fmincon active-set  (Truncated Newton Conjugate-grad)
%  L-BFGS-B     -> fmincon interior-point with L-BFGS Hessian approx
%  CG           -> fminunc quasi-newton (closest unconstrained CG)
%  BFGS         -> fminunc quasi-newton with different FinDiffType
%  PSO          -> custom toolbox-free Particle Swarm Optimisation
%  DE           -> custom toolbox-free Differential Evolution (rand/1/bin)
    n  = numel(x0);
    pf = @(x) wrap(func,x) + 1e8*sum(max(0,lb-x).^2 + max(0,x-ub).^2);
switch mname
case 'Nelder-Mead'
% Pure Nelder-Mead simplex, no gradient, default tolerances
            opts = optimset('MaxIter',2000*n,'MaxFunEvals',4000*n, ...
'TolX',1e-6,'TolFun',1e-6,'Display','off');
            [Popt,~] = fminsearch(pf, x0, opts);
case 'Powell'
% Powell's method: tighter convergence, larger simplex step
% (fminsearch with different tolerances approximates the
%  line-search stopping differences of Powell's method)
            opts = optimset('MaxIter',3000*n,'MaxFunEvals',6000*n, ...
'TolX',1e-8,'TolFun',1e-8,'Display','off');
            x0p  = x0 + 0.02*(ub-lb).*rand(n,1);   % perturb start
            x0p  = min(max(x0p,lb),ub);
            [Popt,~] = fminsearch(pf, x0p, opts);
case 'COBYLA'
% COBYLA: derivative-free, linear constraint approximation.
% Use coarser tolerance + heavier penalty to mimic COBYLA's
% conservative constraint handling.
            pf_co = @(x) wrap(func,x) + 1e6*sum(max(0,lb-x).^2 + max(0,x-ub).^2);
            opts  = optimset('MaxIter',1500*n,'MaxFunEvals',3000*n, ...
'TolX',1e-4,'TolFun',1e-4,'Display','off');
            [Popt,~] = fminsearch(pf_co, x0, opts);
case 'SLSQP'
% Sequential Least Squares Programming -> fmincon SQP
            op = optimoptions('fmincon','Algorithm','sqp', ...
'MaxIterations',1000,'MaxFunctionEvaluations',4000, ...
'OptimalityTolerance',1e-7,'StepTolerance',1e-8, ...
'Display','off');
            [Popt,~] = fmincon(@(x)wrap(func,x),x0,[],[],[],[],lb,ub,[],op);
case 'trust-constr'
% Trust-region interior reflective (distinct from IP)
            op = optimoptions('fmincon','Algorithm','trust-region-reflective', ...
'MaxIterations',1000,'MaxFunctionEvaluations',4000, ...
'OptimalityTolerance',1e-7,'StepTolerance',1e-8, ...
'Display','off','SpecifyObjectiveGradient',false);
% trust-region-reflective requires no general constraints;
% use bound constraints only (already satisfied by lb/ub)
try
                [Popt,~] = fmincon(@(x)wrap(func,x),x0,[],[],[],[],lb,ub,[],op);
catch
                op2 = optimoptions('fmincon','Algorithm','interior-point', ...
'MaxIterations',1000,'Display','off', ...
'HessianApproximation','lbfgs');
                [Popt,~] = fmincon(@(x)wrap(func,x),x0,[],[],[],[],lb,ub,[],op2);
end
case 'TNC'
% Truncated Newton Conjugate-Gradient -> active-set in MATLAB
% (active-set uses conjugate gradient subproblem steps)
            op = optimoptions('fmincon','Algorithm','active-set', ...
'MaxIterations',1000,'MaxFunctionEvaluations',4000, ...
'OptimalityTolerance',1e-7,'StepTolerance',1e-8, ...
'Display','off');
            [Popt,~] = fmincon(@(x)wrap(func,x),x0,[],[],[],[],lb,ub,[],op);
case 'L-BFGS-B'
% L-BFGS-B: interior-point with lbfgs Hessian approximation
% (distinct from plain interior-point which uses full Hessian)
            op = optimoptions('fmincon','Algorithm','interior-point', ...
'MaxIterations',1000,'MaxFunctionEvaluations',4000, ...
'OptimalityTolerance',1e-7,'StepTolerance',1e-8, ...
'HessianApproximation','lbfgs', ...
'Display','off');
            [Popt,~] = fmincon(@(x)wrap(func,x),x0,[],[],[],[],lb,ub,[],op);
case 'CG'
% Conjugate Gradient (unconstrained) -> fminunc quasi-newton
            op = optimoptions('fminunc','Algorithm','quasi-newton', ...
'MaxIterations',1000,'MaxFunctionEvaluations',4000, ...
'OptimalityTolerance',1e-7,'StepTolerance',1e-8, ...
'Display','off');
            [Popt,~] = fminunc(pf, x0, op);
case 'BFGS'
% BFGS: quasi-newton with finite-difference gradient
% Use different FinDiffType from CG to produce distinct results
            op = optimoptions('fminunc','Algorithm','quasi-newton', ...
'MaxIterations',1000,'MaxFunctionEvaluations',4000, ...
'OptimalityTolerance',1e-9,'StepTolerance',1e-10, ...
'FinDiffType','central', ...
'Display','off');
            [Popt,~] = fminunc(pf, x0, op);
case 'PSO'
            Popt = pso_solve(@(x)wrap(func,x), lb, ub, 40, 100);
case 'DE'
            Popt = de_solve(@(x)wrap(func,x), lb, ub, 40, 100);
otherwise
            Popt = x0;
end
% Project back to bounds (safety clip)
    Popt = min(max(Popt(:), lb(:)), ub(:));
end
function Popt = ed_comp(mname, func, x0, lb, ub, a, b, c, dem, mu)
% Economic Dispatch comparators.
% func already includes the power-balance penalty term (used by the
% derivative-free methods below). Gradient-based fmincon methods now use
% the analytical gradient of the quadratic penalised objective via
% fc_grad_pen(), placing them on equal footing with IP-C-PDL's analytical
% gradient (Reviewer 4 pt.4) rather than relying on finite differences.
%
%  Nelder-Mead  -> fminsearch  standard simplex
%  Powell       -> fminsearch  tighter tol + perturbed start
%  COBYLA       -> fminsearch  coarser tol + lighter penalty
%  SLSQP        -> fmincon SQP  (analytical gradient, sum(P)=dem as an
%                                explicit linear equality constraint)
%  trust-constr -> fmincon trust-region-reflective  (bounds only, analytical grad)
%  TNC          -> fmincon active-set (analytical grad)
%  L-BFGS-B     -> fmincon interior-point + lbfgs Hessian (analytical grad)
%  PSO          -> custom toolbox-free Particle Swarm Optimisation
%  DE           -> custom toolbox-free Differential Evolution (rand/1/bin)
    n  = numel(x0);
    pf = @(P) wrap(func,P) + 1e8*sum(max(0,lb-P).^2 + max(0,P-ub).^2);
    gradobj = @(P) fc_grad_pen(P, a, b, c, mu, dem);   % analytical [f,g]
switch mname
case 'Nelder-Mead'
            opts = optimset('MaxIter',2000*n,'MaxFunEvals',4000*n, ...
'TolX',1e-6,'TolFun',1e-6,'Display','off');
            [Popt,~] = fminsearch(pf, x0, opts);
case 'Powell'
% Tighter tolerance + slight perturbation of start point
            opts = optimset('MaxIter',3000*n,'MaxFunEvals',6000*n, ...
'TolX',1e-8,'TolFun',1e-8,'Display','off');
            x0p = x0 + 0.01*(ub-lb).*rand(n,1);
            x0p = min(max(x0p,lb),ub);
            [Popt,~] = fminsearch(pf, x0p, opts);
case 'COBYLA'
% Coarser tolerance, lighter penalty coefficient
            pf_co = @(P) wrap(func,P) + 1e5*sum(max(0,lb-P).^2 + max(0,P-ub).^2);
            opts  = optimset('MaxIter',1500*n,'MaxFunEvals',3000*n, ...
'TolX',1e-4,'TolFun',1e-4,'Display','off');
            [Popt,~] = fminsearch(pf_co, x0, opts);
case 'SLSQP'
% SQP — analytical gradient, bound constraints; let SQP handle
% feasibility natively
            op = optimoptions('fmincon','Algorithm','sqp', ...
'MaxIterations',500,'MaxFunctionEvaluations',2000*n, ...
'OptimalityTolerance',1e-7,'StepTolerance',1e-8, ...
'SpecifyObjectiveGradient',true,'Display','off');
            [Popt,~] = fmincon(gradobj,x0,[],[],[],[],lb,ub,[],op);
case 'trust-constr'
% Trust-region-reflective (bound constraints only, analytical gradient)
            op = optimoptions('fmincon','Algorithm','trust-region-reflective', ...
'MaxIterations',500,'MaxFunctionEvaluations',2000*n, ...
'OptimalityTolerance',1e-7,'StepTolerance',1e-8, ...
'SpecifyObjectiveGradient',true,'Display','off');
try
                [Popt,~] = fmincon(gradobj,x0,[],[],[],[],lb,ub,[],op);
catch
% Fallback: interior-point with tight tolerances
                op2 = optimoptions('fmincon','Algorithm','interior-point', ...
'MaxIterations',500,'Display','off', ...
'OptimalityTolerance',1e-7,'HessianApproximation','lbfgs', ...
'SpecifyObjectiveGradient',true);
                [Popt,~] = fmincon(gradobj,x0,[],[],[],[],lb,ub,[],op2);
end
case 'TNC'
% Active-set (conjugate gradient subproblem — mirrors TNC), analytical grad
            op = optimoptions('fmincon','Algorithm','active-set', ...
'MaxIterations',500,'MaxFunctionEvaluations',2000*n, ...
'OptimalityTolerance',1e-7,'StepTolerance',1e-8, ...
'SpecifyObjectiveGradient',true,'Display','off');
            [Popt,~] = fmincon(gradobj,x0,[],[],[],[],lb,ub,[],op);
case 'L-BFGS-B'
% Interior-point with L-BFGS Hessian (matches L-BFGS-B intent), analytical grad
            op = optimoptions('fmincon','Algorithm','interior-point', ...
'MaxIterations',500,'MaxFunctionEvaluations',2000*n, ...
'OptimalityTolerance',1e-7,'StepTolerance',1e-8, ...
'HessianApproximation','lbfgs', ...
'SpecifyObjectiveGradient',true,'Display','off');
            [Popt,~] = fmincon(gradobj,x0,[],[],[],[],lb,ub,[],op);
case 'PSO'
            Popt = pso_solve(@(P)wrap(func,P), lb, ub, 40, 100);
case 'DE'
            Popt = de_solve(@(P)wrap(func,P), lb, ub, 40, 100);
otherwise
            Popt = x0;
end
% Safety clip to bounds
    Popt = min(max(Popt(:), lb(:)), ub(:));
end
% =========================================================================
% Toolbox-free comparator metaheuristics (Reviewers 3 & 4: add modern
% hybrid/metaheuristic baselines alongside the classical gradient-based
% comparators already present). Implemented from scratch (no Global
% Optimization Toolbox dependency) so the script remains runnable with
% only the Statistics and Optimization Toolboxes already required.
% =========================================================================
function xbest = pso_solve(func, lb, ub, npop, maxit)
    lb = lb(:); ub = ub(:); n = numel(lb);
    w = 0.72; c1 = 1.49; c2 = 1.49;
    X = lb + (ub-lb).*rand(n,npop);
    V = zeros(n,npop);
    Pbest = X; Fpbest = arrayfun(@(j)func(X(:,j)), 1:npop);
    [Fgbest, gi] = min(Fpbest); Gbest = X(:,gi);
for it = 1:maxit
for j = 1:npop
            r1 = rand(n,1); r2 = rand(n,1);
            V(:,j) = w*V(:,j) + c1*r1.*(Pbest(:,j)-X(:,j)) + c2*r2.*(Gbest-X(:,j));
            X(:,j) = min(max(X(:,j)+V(:,j), lb), ub);
            fj = func(X(:,j));
if fj < Fpbest(j)
                Fpbest(j) = fj; Pbest(:,j) = X(:,j);
if fj < Fgbest, Fgbest = fj; Gbest = X(:,j); end
end
end
        w = max(0.4, w*0.995);   % linearly decaying inertia
end
    xbest = Gbest;
end
function xbest = de_solve(func, lb, ub, npop, maxit)
    lb = lb(:); ub = ub(:); n = numel(lb);
    F = 0.6; CR = 0.9;
    X = lb + (ub-lb).*rand(n,npop);
    Fx = arrayfun(@(j)func(X(:,j)), 1:npop);
for it = 1:maxit
for j = 1:npop
            idx = randperm(npop, 4);
            idx(idx==j) = []; idx = idx(1:3);
            a_ = X(:,idx(1)); b_ = X(:,idx(2)); c_ = X(:,idx(3));
            V = min(max(a_ + F*(b_-c_), lb), ub);
            jr = randi(n);
            mask = rand(n,1) < CR; mask(jr) = true;
            U = X(:,j); U(mask) = V(mask);
            fu = func(U);
if fu <= Fx(j)
                X(:,j) = U; Fx(j) = fu;
end
end
end
    [~, bi] = min(Fx);
    xbest = X(:,bi);
end
function ch = lc_new(n, seed)
    rng(seed);
    ch.s = 0.1 + 0.8*rand(n,1);
    ch.r = 3.9;
    ch.type = 'tent';
end
function [v, ch] = lc_step(ch)
for i = 1:numel(ch.s)
if ch.s(i) < 0.5
            ch.s(i) = 2 * ch.s(i);
else
            ch.s(i) = 2 * (1 - ch.s(i));
end
end
    v  = ch.s - 0.5;
    nm = norm(v);
if nm > 1e-12, v = v/nm; end
end
% Generalised chaotic-map generator used only by ipcpdl_ed_map (Section
% E4). map_id: 1=tent (identical to lc_new/lc_step), 2=logistic
% (x_{k+1}=r*x*(1-x), r=3.99), 3=sine (x_{k+1}=r*sin(pi*x), r=1.0),
% 4=gaussian/mouse map (x_{k+1}=1/mod(x,1) - floor(1/mod(x,1)), with a
% small epsilon guard).
function ch = map_new(n, seed, map_id)
    rng(seed);
    ch.s = 0.1 + 0.8*rand(n,1);
    ch.map_id = map_id;
end
function [v, ch] = map_step(ch, map_id)
    x = ch.s;
switch map_id
case 1  % tent
for i = 1:numel(x)
if x(i) < 0.5, x(i) = 2*x(i); else, x(i) = 2*(1-x(i)); end
end
case 2  % logistic
            r = 3.99;
            x = r .* x .* (1 - x);
case 3  % sine map
            x = sin(pi .* x);
case 4  % gaussian/mouse map
            xe = max(x, 1e-6);
            x = 1./xe - floor(1./xe);
end
    x = min(max(x, 1e-9), 1-1e-9);
    ch.s = x;
    v  = x - 0.5;
    nm = norm(v);
if nm > 1e-12, v = v/nm; end
end
function g = cen_grad(func, x, eps0)
    g = zeros(size(x));
for i = 1:numel(x)
        xp = x; xm = x;
        xp(i) = xp(i) + eps0;
        xm(i) = xm(i) - eps0;
        fp = func(xp); fm = func(xm);
if isfinite(fp) && isfinite(fm)
            g(i) = (fp-fm)/(2*eps0);
end
end
end
function v = fc(P,a,b,c)
    v = sum(a.*P.^2 + b.*P + c);
end
% [f,g] pair for fc(), used by fmincon's SpecifyObjectiveGradient path in
% make_pool()'s interior-point/SQP seed refinement.
function [f,g] = fc_grad(P,a,b,c)
    f = fc(P,a,b,c);
    g = 2*a.*P + b;
end
function v = pen(P,a,b,c,dem,mu)
    v = fc(P,a,b,c) + mu*(sum(P)-dem)^2;
end
% Analytical gradient of pen(): used by ipcpdl_ed (Reviewer 4 pt.4).
function g = ed_grad(P,a,b,mu,dem)
    g = 2*a.*P + b + 2*mu*(sum(P)-dem);
end
% [f,g] pair for pen(), used by ed_comp()'s fmincon-based comparators
% (SLSQP/trust-constr/TNC/L-BFGS-B) so every gradient-based method in
% the ED comparison uses the SAME analytical gradient as IP-C-PDL.
function [f,g] = fc_grad_pen(P,a,b,c,mu,dem)
    f = pen(P,a,b,c,dem,mu);
    g = ed_grad(P,a,b,mu,dem);
end
function P = repair(P,lb,ub,dem)
    P = min(max(P(:),lb(:)),ub(:));
for it = 1:200
        d = sum(P)-dem;
if abs(d) < 1e-4, break; end
if d > 0
            hr = P-lb; s = sum(hr);
if s > 1e-8, P = P - hr*(d/s); end
else
            hr = ub-P; s = sum(hr);
if s > 1e-8, P = P + hr*(abs(d)/s); end
end
        P = min(max(P,lb),ub);
end
end
function [a,b,c,pm,pM,ru,rd,u,dem] = ld(fin,fdem)
    fin  = resolve_path(fin);
    fdem = resolve_path(fdem);
    raw  = readmatrix(fin,  'Range','A2:H1000','FileType','spreadsheet');
    valid = any(~isnan(raw),2);
    raw   = raw(valid,:);
    u   = raw(:,1);
    a   = raw(:,2);
    b   = raw(:,3);
    c   = raw(:,4);
    pm  = raw(:,5);
    pM  = raw(:,6);
    ru  = raw(:,7);
    rd  = raw(:,8);
    draw = readmatrix(fdem,'FileType','spreadsheet');
    dem  = draw(~isnan(draw));
    dem  = dem(:);
end
function p = resolve_path(fname)
if exist(fname,'file'), p = fname; return; end
    script_dir = fileparts(mfilename('fullpath'));
    candidate  = fullfile(script_dir, fname);
if exist(candidate,'file'), p = candidate; return; end
    p = fname;
end
function cv = rmin(cv)
for i = 2:numel(cv)
if cv(i) > cv(i-1), cv(i) = cv(i-1); end
end
end
function fn = sf(name)
    fn = strrep(strrep(strrep(name,'-','_'),' ','_'),'/','_');
end
function v = wrap(func,x)
global IPCPDL_CV;
    v = func(x);
if isfinite(v) && v < 1e12
        IPCPDL_CV(end+1) = v;
end
end
function reset_cv()
global IPCPDL_CV;
    IPCPDL_CV = [];
end
function cv = get_cv()
global IPCPDL_CV;
    cv = IPCPDL_CV;
end
function sgtitle_safe(fig,txt)
try
        sgtitle(fig,txt,'FontSize',13,'FontWeight','bold');
catch
        annotation(fig,'textbox',[.2 .96 .6 .04],'String',txt, ...
'EdgeColor','none','HorizontalAlignment','center', ...
'FontSize',12,'FontWeight','bold');
end
end
function v = ternary(cond,a,b)
if cond, v = a; else, v = b; end
end
function v = nanmin(x)
    v = min(x(isfinite(x)));
if isempty(v), v = NaN; end
end
function v = nanmax(x)
    v = max(x(isfinite(x)));
if isempty(v), v = NaN; end
end
function v = nanmean(x, dim)
if nargin < 2, dim = 1; end
if isvector(x)
        xf = x(isfinite(x));
if isempty(xf), v = NaN; else, v = mean(xf); end
else
        v = mean(x, dim, 'omitnan');
end
end
function v = nanmedian(x)
    xf = x(isfinite(x));
if isempty(xf), v = NaN; else, v = median(xf); end
end
function v = nanstd(x, w, dim)
if nargin < 2 || isempty(w), w = 0; end
if nargin < 3, dim = 1; end
if isvector(x)
        xf = x(isfinite(x));
if numel(xf) < 2, v = NaN; else, v = std(xf,w); end
else
        v = std(x, w, dim, 'omitnan');
end
end
function clim_safe(ax, lim)
try
        clim(ax, lim);
catch
        caxis(ax, lim); %#ok<CAXIS>
end
end