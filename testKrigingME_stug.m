% this code is used to test the KrigingME_stg function to test if the function is working as expected
% run('../../BMELIB2.0c_MATLAB2009a/startup.m');
% addpath('../../analysis/');

load('krigInputs.mat', 'KG', 'KS', 'BMEparam');

sk = KS.harddata.sMSh(:,1:2);
tMEest = [2408];

BMEmeans = cell(length(tMEest),1);
BMEvars = cell(length(tMEest),1);

for i=1:length(tMEest)
    tk = tMEest(i);
    pk = [sk ones(size(sk,1),1)*tk];

    % use the new function
    [BMEm,BMEv] = krigingME_stg(pk,KS.harddata,KS.softdata,KG.covmodel,KG.covparam,BMEparam.nhmax,BMEparam.nsmax,BMEparam.dmax,BMEparam.order,BMEparam.options);
    BMEmeans{i} = BMEm;
    BMEvars{i} = BMEv;
end

% using old function
% harddata
ph_s = KS.harddata.sMSh(:,1:2);
ph_t = KS.harddata.tMEh;
Xh = KS.harddata.xh;
[ph, xh] = valstg2stv(Xh, ph_s, ph_t);

% softdata
ps_s = KS.softdata.sMSs(:,1:2);
ps_t = KS.softdata.tMEs;
% Xs = KS.softdata.xms;
[ps, xms] = valstg2stv(KS.softdata.xms, ps_s, ps_t);
[~, xvs] = valstg2stv(KS.softdata.xvs, ps_s, ps_t);

BMEoldmeans = cell(length(tMEest),1);
BMEoldvars = cell(length(tMEest),1);

for i=1:length(tMEest)
    tk = tMEest(i);
    pk = [sk ones(size(sk,1),1)*tk];

    % use the old function
    [BMEm,BMEv] = krigingME(pk,ph,ps,xh,xms,xvs,KG.covmodel,KG.covparam,BMEparam.nhmax,BMEparam.nsmax,BMEparam.dmax,KG.order,BMEparam.options);
    BMEoldmeans{i} = BMEm;
    BMEoldvars{i} = BMEv;
end

% check if the results are the same
assert(isequal(BMEmeans, BMEoldmeans), 'means are not equal');