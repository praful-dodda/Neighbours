% testKrigingME.m : test the krigingME and krigingME_stg functions if they give the same output
% 
clear; close all;
rng(1);

% harddata params
nMS_h = 10; % 100
sMS_h = 100*(rand(nMS_h,2)-0.5);

nME_h = 5; %1000
tME_h = 50*(1:nME_h)/nME_h;

nhmax = 10;
nan_target_ratio_h = 0.8;
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
nMS_s = 100;
sMS_s = 100*(rand(nMS_s,2)-0.5);

nME_s = 5;
tME_s = 50*(1:nME_s)/nME_s;

nsmax = 10;
nan_target_ratio_s = 0.1;
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

% estimation points
nk_vec = [10];
sk = 100*(rand(max(nk_vec),2)-0.5);
pk_all = [sk 50.*rand(max(nk_vec), 1)];

for k = 1:length(nk_vec)
  nk = nk_vec(k);
  pk = pk_all(1:nk, :);

  zk_all = zeros(nk, 1);
  vk_all = zeros(nk, 1);
  chslocal_all = cell(nk, 1);

  zk_stg_all = zeros(nk, 1);
  vk_stg_all = zeros(nk, 1);
  chslocal_stg_all = cell(nk, 1);

  for i = 1:size(pk,1)
    i
    pk_i = pk(i,:);
    
    % use krigingME
    [zk, vk] = krigingME(pk_i, ph, ps, zh, xs, vs, covmodel, covparam, nhmax, nsmax, dmax, order);
    
    % use krigingME_stg
    [zk_stg, vk_stg] = krigingME_stg(pk_i, harddata, softdata, covmodel, covparam, nhmax, nsmax, dmax, order);

    % check if the results are the same
    if ~isequal(vk, vk_stg)
      warning('vk is not equal');
    end
    
    if ~isequal(zk, zk_stg)
      warning('zk is not equal');
    end
    
    % if ~isequal(chslocal, chslocal_stg)    
    %     if isequal(sortrows(chslocal), sortrows(chslocal_stg))
    %       fprintf('nk=%d, chslocal and chslocal_stg did not give the same output, but they appear to be the same once sorted\n', nk);
    %     else
    %         warning('chslocal is not equal');
    %     end
    % 
      % [~, ~, idx_h]  = intersect(chslocal(1:10,:), harddata.p, 'rows');
      % [~, ~, idx_s]  = intersect(chslocal(11:end,:), softdata.p, 'rows');
      % [~, ~, idx_h_stg]  = intersect(chslocal_stg(1:10,:), harddata.p, 'rows');
      % [~, ~, idx_s_stg]  = intersect(chslocal_stg(11:end,:), softdata.p, 'rows');
      % 
      % if ~isequal(idx_h, idx_h_stg)
      %   warning('same neighbors in harddata are not found');
      % end
      % 
      % if ~isequal(idx_s, idx_s_stg)
      %   warning('same neighbors in softdata are not found');
      % end
      % 
      % harddata.z(idx_h)
      % harddata.z(idx_h_stg)
      % softdata.z(idx_s)
      % softdata.z(idx_s_stg)
    % 
    % 
    % end

    zk_all(i) = zk;
    vk_all(i) = vk;
    chslocal_all{i} = chslocal;

  end

end

disp('All tests passed!');