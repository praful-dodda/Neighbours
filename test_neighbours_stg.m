% this script is used to test the function neighbours_stg.m

clear; close all;

% load krigInputs.mat
load krigInputs.mat 

% estimation point(s)
% load point_data.mat

p0 = [-90.1, 27.1, 2404]; % July, 31st 2012

% harddata and softdata for `stg`
harddata.sMS = KS.harddata.sMSh;
harddata.tME = KS.harddata.tMEh;
% harddata.z = KS.harddata.xh;
harddata.Zisnotnan = ~isnan(harddata.z);
harddata.nanratio = sum(~harddata.Zisnotnan(:))/length(harddata.Zisnotnan(:));
harddata.index_stg_to_stv = cumsum(harddata.Zisnotnan(:));
[harddata.p, harddata.z] = valstg2stv(KS.harddata.xh, harddata.sMS, harddata.tME);

softdata.sMS = KS.softdata.sMSs;
softdata.tME = KS.softdata.tMEs;
% softdata.z = KS.softdata.xms;
softdata.xvs = KS.softdata.xvs;

softdata.Zisnotnan = ~isnan(softdata.z);
softdata.nanratio = sum(~softdata.Zisnotnan(:))/length(softdata.Zisnotnan(:));
softdata.index_stg_to_stv = cumsum(softdata.Zisnotnan(:));
[softdata.p, softdata.z] = valstg2stv(KS.softdata.xms, softdata.sMS, softdata.tME);

% harddata and softdata for `stv`
[ch, zh] = valstg2stv(KS.harddata.xh, KS.harddata.sMSh, KS.harddata.tMEh);
[cs, zs] = valstg2stv(KS.softdata.xms, KS.softdata.sMSs, KS.softdata.tMEs);
[~, vs] = valstg2stv(KS.softdata.xms, KS.softdata.sMSs, KS.softdata.tMEs);

nhmax = 10; nsmax = 5; dmax = BMEparam.dmax;

% neighbours.m for reference
% soft data neighbours
[psub_orig,zsub_orig,dsub_orig,nsub_orig,index_orig]=neighbours(p0,cs,zs,nsmax,dmax);
neigh.softdata.psub_orig = psub_orig;
neigh.softdata.zsub_orig = zsub_orig;
neigh.softdata.dsub_orig = dsub_orig;
neigh.softdata.nsub_orig = nsub_orig;
neigh.softdata.index_orig = index_orig;

% hard data neighbours
[psub_orig,zsub_orig,dsub_orig,nsub_orig,index_orig]=neighbours(p0,ch,zh,nhmax,dmax);
neigh.harddata.psub_orig = psub_orig;
neigh.harddata.zsub_orig = zsub_orig;
neigh.harddata.dsub_orig = dsub_orig;
neigh.harddata.nsub_orig = nsub_orig;
neigh.harddata.index_orig = index_orig;

% neighbours_stg.m
% soft data neighbours
[psub_stg,zsub_stg,dsub_stg,nsub_stg,index_stg]=neighbours_stg(p0,softdata,nsmax,dmax);
neigh.softdata.psub_stg = psub_stg;
neigh.softdata.zsub_stg = zsub_stg;
neigh.softdata.dsub_stg = dsub_stg;
neigh.softdata.nsub_stg = nsub_stg;
neigh.softdata.index_stg = index_stg;

% hard data neighbours
[psub_stg,zsub_stg,dsub_stg,nsub_stg,index_stg]=neighbours_stg(p0,harddata,nhmax,dmax);
neigh.harddata.psub_stg = psub_stg;
neigh.harddata.zsub_stg = zsub_stg;
neigh.harddata.dsub_stg = dsub_stg;
neigh.harddata.nsub_stg = nsub_stg;
neigh.harddata.index_stg = index_stg;

% check if the results are the same
% softdata
%          psub_orig
if ~isequal(neigh.softdata.psub_orig, neigh.softdata.psub_stg)
    warning('neighbours and neighbours_stg did not give the same output for softdata.psub');
end
%          zsub_orig
if ~isequal(neigh.softdata.zsub_orig(~isnan(neigh.softdata.zsub_orig)), neigh.softdata.zsub_stg(~isnan(neigh.softdata.zsub_stg)))
    warning('neighbours and neighbours_stg did not give the same output for softdata.zsub');
end
%          dsub_orig
if ~isequal(neigh.softdata.dsub_orig, neigh.softdata.dsub_stg)
    warning('neighbours and neighbours_stg did not give the same output for softdata.dsub');
end
%          nsub_orig
if ~isequal(neigh.softdata.nsub_orig, neigh.softdata.nsub_stg)
    warning('neighbours and neighbours_stg did not give the same output for softdata.nsub');
end
%          index_orig
if ~isequal(neigh.softdata.index_orig, neigh.softdata.index_stg)
    warning('neighbours and neighbours_stg did not give the same output for softdata.index');
end

% harddata
if ~isequal(neigh.harddata.psub_orig, neigh.harddata.psub_stg)
    warning('neighbours and neighbours_stg did not give the same output for harddata.psub');
end
if ~isequal(neigh.harddata.zsub_orig(~isnan(neigh.harddata.zsub_orig)), neigh.harddata.zsub_stg(~isnan(neigh.harddata.zsub_stg)))
    warning('neighbours and neighbours_stg did not give the same output for harddata.zsub');
end
if ~isequal(neigh.harddata.dsub_orig, neigh.harddata.dsub_stg)
    warning('neighbours and neighbours_stg did not give the same output for harddata.dsub');
end
if ~isequal(neigh.harddata.nsub_orig, neigh.harddata.nsub_stg)
    warning('neighbours and neighbours_stg did not give the same output for harddata.nsub');
end
if ~isequal(neigh.harddata.index_orig, neigh.harddata.index_stg)
    warning('neighbours and neighbours_stg did not give the same output for harddata.index');
end
