% testNeighbours_nonan : speed the neighbours functions for speed
%

nMS=10000;  % try nMS=100 and nMS=10000
nMEvec=[1 10 100 250 500 1000 2000];
nmax=100; % try nMS=2, nMS=10 and nMS=200
dmax=[1000 10000 0.1];

sMS=100*(rand(nMS,2)-0.5);
tCPU_orig=nan(1,length(nMEvec));
tCPU_stv=nan(1,length(nMEvec));
tCPU_stg=nan(1,length(nMEvec));
tCPU_stg_nonan = nan(1, length(nMEvec));

for jj=1:length(nMEvec)

  % Generate the data
  nME=nMEvec(jj);
  tME=1:nME;
  p0=[0 0 ceil(tME(end)/2)];
  Z=100*rand(nMS,nME);
  if sum(~isnan(Z(:)))==0
    Z(1,1)=0;
  end
  [p,z]=valstg2stv(Z,sMS,tME);
  
  % use neighbours.m for reference
  tic;
  [psub_orig,zsub_orig,dsub_orig,nsub_orig,index_orig]=neighbours(p0,p,z,nmax,dmax);
  tCPU_orig(jj)=toc;

  % test neighbours_stv.m
  tic;
  [psub_stv,zsub_stv,dsub_stv,nsub_stv,index_stv]=neighbours_stv(p0,p,z,nmax,dmax);
  tCPU_stv(jj)=toc;
  if ~isequal(index_orig,index_stv) || ~isequal(psub_orig,psub_stv) ||  ~isequal(zsub_orig,zsub_stv) ||  ~isequal(dsub_orig,dsub_stv)
    error('neighbours and neighbours_stv did not give the same output');
  end

  % test neighbours_stg_nonan.m
  data.Z=Z;  
  data.sMS=sMS;
  data.tME=tME;
  data.z=z;
  data.p=p;
  
  tic;
  [psub_stg,zsub_stg,dsub_stg,nsub_stg,index_stg]=neighbours_stg(p0,data,nmax,dmax);
  tCPU_stg(jj)=toc;
  if ~isequal(index_stv,index_stg) || ~isequal(psub_stv,psub_stg) || ~isequal(zsub_stv,zsub_stg) ||  ~isequal(dsub_stv,dsub_stg)
    error('neighbours_stv and neighbours_stg did not give the same output');
  end

  tic;
  [psub_stg_nonan,zsub_stg_nonan,dsub_stg_nonan,nsub_stg_nonan,index_stg_nonan]=neighbours_stg_nonan(p0,data,nmax,dmax);
  tCPU_stg_nonan(jj)=toc;
  if ~isequal(index_stv,index_stg_nonan) || ~isequal(psub_stv,psub_stg_nonan) || ~isequal(zsub_stv,zsub_stg_nonan) ||  ~isequal(dsub_stv,dsub_stg_nonan)
    error('neighbours_stv and neighbours_stg_nonan did not give the same output');
  end
end

figure
hold on
h_orig=plot(nMEvec,tCPU_orig,'o-g','LineWidth',1);
h_stv=plot(nMEvec,tCPU_stv,'o-b','LineWidth',1);
h_stg=plot(nMEvec,tCPU_stg,'o-r','LineWidth',1);
h_stg_nonan=plot(nMEvec,tCPU_stg_nonan,'o-k','LineWidth',1);
xlabel('nME');
ylabel('CPU time (sec)');
yyaxis right;
h_speedup=plot(nMEvec,tCPU_stv./tCPU_stg, 'LineWidth',1);
h_speedup_nonan=plot(nMEvec,tCPU_stv./tCPU_stg_nonan,'LineWidth',1);
ylabel('stv to stg speedup');
legend([h_orig h_stv h_stg h_stg_nonan h_speedup h_speedup_nonan],{'CPU time orig','CPU time stv', 'CPU time stg', 'CPU time stg nonan','speedup', 'speedup nonan'});
title(sprintf('CPU time and speedup, nMS=%d, nmax=%d, nonan',nMS,nmax));