% this code is used to run the krigingME_stg function and run the profiler 
clear;
rng(1);

% control params
nMS_h = 200; % 100
nME_h = 1000; %1000
nan_target_ratio_h = 0.8;

nMS_s = 10000;
nME_s = 100;
nan_target_ratio_s = 0.1;

nk_vec = [10000];

% harddata params
sMS_h = 100*(rand(nMS_h,2)-0.5);
tME_h = 1:nME_h;
tME_h = tME_h + 20;

nhmax = 10;
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

% softdata params
sMS_s = 100*(rand(nMS_s,2)-0.5);

tME_s = 1:nME_s;
tME_s = tME_s + 20;

nsmax = 10;

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
covmodel = {'exponentialC'};
covparam = {[0.4223 0.1200 4]};
sill = 0.6812;
order = NaN;
dmax = [100 1000 0.2];

sk = 100*(rand(max(nk_vec),2)-0.5);
tk = rand(max(nk_vec),1);
tk = tk + 20;
pk_all = [sk tk];

for i = 1:length(nk_vec)
    nk = nk_vec(i);
    pk = pk_all(1:nk, :);
    [zk_stg, vk_stg] = krigingME_stg(pk, harddata, softdata, covmodel, covparam, nhmax, nsmax, dmax, order);
end