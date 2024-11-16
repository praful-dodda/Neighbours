% this script is used to test the krigingME_stg and krigingME functions
clear; close all;
%% --------
% krigInputs.mat has the following variables: 
% BMEparam = 
%   struct with fields:
%       nhmax: 100
%       nsmax: 3
%       order: NaN
%        dmax: [10 100 0.0296]
%      maxpts: 1000000
%        rEps: 0.0500
%        nMom: 2
%     options: [0 1.0000e-04 1000000 0.0500 0 25 1.0000e-03 2 0 0 0 0 0 100 0 0 0 0 0 0.6800 0 0 0 0 0 0 0 0 0]
% 
% KG =  struct with fields:

%     covmodel: {'exponentialC/exponentialC'  'exponentialC/exponentialC'  'exponentialC/exponentialC'}
%     covparam: {[0.4223 0.1200 4]  [0.0341 0.1200 450]  [0.2248 15 450]}
%         sill: [0.6812]
%        order: [NaN]
% 
% KS = struct with fields:
% harddata: [1x1 struct]
% softdata: [1x1 struct]
% softpdftype: 2
%       nl: []
%     limi: []
% probdens: []
% KS.harddata = struct with fields:
%   sMSh: [182x2 double]
%   tMEh: [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 ... ] (1x4144 double)
%     xh: [182x4144 double]
% 
% KS.softdata = struct with fields:
%   sMSs: [34840x2 double]
%   tMEs: [2400 2401 2402 2403 2404 2405 2406 2407 2408 2409 2410 2411 2412 2413 2414 2415 2416 2417 2418 2419 2420 2421 ... ] (1x41 double)
%    xms: [34840x41 double]
%    xvs: [34840x41 double]
%% --------

% load krigInputs.mat
load krigInputs.mat 

% estimation point(s)
load point_data.mat
nk = 1;
ck = point_data(1:min(nk, size(point_data, 2)), :);
pk = [-90.1, 27.1, 2404];

% harddata and softdata
harddata.sMS = KS.harddata.sMSh;
harddata.tME = KS.harddata.tMEh;
harddata.z = KS.harddata.Xh;
harddata.Zisnotnan = ~isnan(harddata.z);
harddata.nanratio = sum(~harddata.Zisnotnan(:))/length(harddata.Zisnotnan(:));
harddata.index_stg_to_stv = cumsum(harddata.Zisnotnan(:));
[harddata.p, ~] = valstg2stv(harddata.z, harddata.sMS, harddata.tME);

softdata.sMS = KS.softdata.sMSs;
softdata.tME = KS.softdata.tMEs;
softdata.z = KS.softdata.Xms;
softdata.Xvs = KS.softdata.Xvs;

softdata.Zisnotnan = ~isnan(softdata.z);
softdata.nanratio = sum(~softdata.Zisnotnan(:))/length(softdata.Zisnotnan(:));
softdata.index_stg_to_stv = cumsum(softdata.Zisnotnan(:));
[softdata.p, ~] = valstg2stv(softdata.z, softdata.sMS, softdata.tME);

% BME parameters
covmodel = KG.covmodel;
covparam = KG.covparam;
nhmax = BMEparam.nhmax;
nsmax = 10;
dmax = BMEparam.dmax;
order = BMEparam.order;
options = BMEparam.options;

% use krigingME_stg 
[zk_stg, vk_stg] = krigingME_stg(pk, harddata, softdata, covmodel, covparam, nhmax, nsmax, dmax, order, options);

% Display the results
disp('Estimated values (zk_stg) using krigingME_stg:');
disp(zk_stg);
disp('Kriging variances (vk_stg) using krigingME_stg:');
disp(vk_stg);

% harddata and softdata
[ch, zh] = valstg2stv(KS.harddata.Xh, KS.harddata.sMSh, KS.harddata.tMEh);
[cs, zs] = valstg2stv(KS.softdata.Xms, KS.softdata.sMSs, KS.softdata.tMEs);
[~, vs] = valstg2stv(KS.softdata.Xms, KS.softdata.sMSs, KS.softdata.tMEs);

% use krigingME
[zk, vk] = krigingME(pk, ch, cs, zh, zs, vs, covmodel, covparam, nhmax, nsmax, dmax, order, options);

% Display the results
disp('Estimated values (zk) using krigingME:');
disp(zk);
disp('Kriging variances (vk) using krigingME:');
disp(vk);


