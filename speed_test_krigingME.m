function [tCPU_orig, tCPU_stg] = speed_test_krigingME(softDataParam, hardDataParam, krigingParam, nk_vec, speed_plot)

    if nargin < 1
        softDataParam.nMS = 10000;
        softDataParam.nME = 100;
        softDataParam.nan_target_ratio = 0.1;
        softDataParam.nsmax = 10;
    end

    if nargin < 2
        hardDataParam.nMS = 200;
        hardDataParam.nME = 1000;
        hardDataParam.nan_target_ratio = 0.8;
        hardDataParam.nhmax = 10;
    end

    if nargin < 3
        krigingParam.covmodel = {'exponentialC'};
        krigingParam.covparam = {[0.4223 0.1200 4]};
        krigingParam.sill = 0.6812;
        krigingParam.order = NaN;
        krigingParam.dmax = [100 1000 0.2];
    end

    if nargin < 4
        nk_vec = [1, 10, 100, 1000, 5000, 10000];
    end

    if nargin < 5
        speed_plot = 1;
    end

    % % control params
    % nMS_h = 200; % 100
    % nME_h = 1000; %1000
    % nan_target_ratio_h = 0.8;
    nMS_h = hardDataParam.nMS;
    nME_h = hardDataParam.nME;
    nan_target_ratio_h = hardDataParam.nan_target_ratio;
    nhmax = hardDataParam.nhmax;

    % nMS_s = 10000;
    % nME_s = 100;
    % nan_target_ratio_s = 0.1;
    nMS_s = softDataParam.nMS;
    nME_s = softDataParam.nME;
    nan_target_ratio_s = softDataParam.nan_target_ratio;
    nsmax = softDataParam.nsmax;

    % nk_vec = [1, 10, 100, 1000, 5000, 10000];

    % harddata params
    sMS_h = 100*(rand(nMS_h,2)-0.5);
    tME_h = 1:nME_h;
    tME_h = tME_h + 20;

    Zh = 100*rand(nMS_h, nME_h);
    Zh(Zh < nan_target_ratio_h*100) = NaN;
    Zh(1,:) = 0;
    Zh(:,1) = 0;
    Zisnotnan_h = ~isnan(Zh);
    nanratio_h = sum(~Zisnotnan_h(:))/length(Zisnotnan_h(:));
    index_stg_to_stv_h = cumsum(Zisnotnan_h(:));

    [ph_stg, zh_stg] = valstg2stv(Zh, sMS_h, tME_h);
    ph = ph_stg(~isnan(zh_stg), :);
    zh = zh_stg(~isnan(zh_stg));

    harddata.sMS = sMS_h;
    harddata.tME = tME_h;
    harddata.Z = Zh;
    harddata.Zisnotnan = Zisnotnan_h;
    harddata.nanratio = nanratio_h;
    harddata.index_stg_to_stv = index_stg_to_stv_h;
    harddata.p = ph;
    harddata.z = zh;
    % [harddata.p, harddata.z] = valstg2stv(Zh, sMS_h, tME_h);

    % softdata params
    sMS_s = 100*(rand(nMS_s,2)-0.5);

    tME_s = 1:nME_s;
    tME_s = tME_s + 20;

    Xms = 100*rand(nMS_s, nME_s);
    Xvs = 100*rand(nMS_s, nME_s);
    Xms(Xms < nan_target_ratio_s*100) = NaN;
    Xms(1,:) = 0;
    Xms(:,1) = 0;
    Xisnotnan_s = ~isnan(Xms);
    nanratio_s = sum(~Xisnotnan_s(:))/length(Xisnotnan_s(:));
    index_stg_to_stv_s = cumsum(Xisnotnan_s(:));
    [ps_stg, xs_stg] = valstg2stv(Xms, sMS_s, tME_s);
    ps = ps_stg(~isnan(xs_stg), :);
    xs = xs_stg(~isnan(xs_stg));

    % at the same index as Xms, Xvs is also NaN
    Xvs(~Xisnotnan_s) = NaN;
    [~, vs_stg] = valstg2stv(Xvs, sMS_s, tME_s);
    vs = vs_stg(~isnan(vs_stg));

    softdata.sMS = sMS_s;
    softdata.tME = tME_s;
    softdata.Xms = Xms;
    softdata.Xvs = Xvs;
    softdata.Zisnotnan = Xisnotnan_s;
    softdata.nanratio = nanratio_s;
    softdata.index_stg_to_stv = index_stg_to_stv_s;
    % [softdata.p, softdata.z] = valstg2stv(Xms, sMS_s, tME_s);
    softdata.p = ps;
    softdata.z = xs;
    softdata.vs = vs;

    % kriging params
    covmodel = krigingParam.covmodel;
    covparam = krigingParam.covparam;
    % sill = krigingParam.sill;
    order = krigingParam.order;
    dmax = krigingParam.dmax;

    sk = 100*(rand(max(nk_vec),2)-0.5);
    tk = rand(max(nk_vec),1);
    tk = tk + 20;
    pk_all = [sk tk];

    tCPU_stg = zeros(length(nk_vec), 1);
    tCPU_orig = zeros(length(nk_vec), 1);

    for i = 1:length(nk_vec)
        nk = nk_vec(i);
        pk = pk_all(1:nk, :);
        tic;
        [zk_stg, vk_stg] = krigingME_stg(pk, harddata, softdata, covmodel, covparam, nhmax, nsmax, dmax, order);
        tCPU_stg(i) = toc;
        
        tic;
        [zk, vk] = krigingME(pk, ph, ps, zh, xs, vs, covmodel, covparam, nhmax, nsmax, dmax, order);
        tCPU_orig(i) = toc;

        if ~isequal(zk, zk_stg) 
            warning('zk are not equal');
        end
        if ~isequal(vk, vk_stg)
            warning('vk are not equal');
        end
    end

    % calculate the speedup
    tCPU_speedup = tCPU_orig./tCPU_stg;

    if speed_plot
        % plot the results
        figure;
        hold on;
        h_orig = plot(nk_vec, tCPU_orig, 'g-o');
        h_stg = plot(nk_vec, tCPU_stg, 'b-o');
        xlabel('nk');
        ylabel('CPU time (s)');
        yyaxis right;
        h_speedup = plot(nk_vec, tCPU_speedup, 'r--o');
        ylabel('krigingME\_stg speedup');
        legend([h_orig, h_stg, h_speedup], {'krigingME', 'krigingME\_stg', 'speedup'});
        title('krigingME vs krigingME\_stg speed test');
        hold off;
    end

    
end