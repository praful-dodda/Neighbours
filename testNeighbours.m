% testNeighbours : speed the neighbours functions for speed
%

nMS=10000;  % try nMS=100 and nMS=10000
nMEvec=[5 10 100 250 500 1000 2000];
nmax=200; % try nMS=2, nMS=10 and nMS=200
nanratio_target=0.5;
dmax=[1000 10000 0.1];

rng(1);
sMS=100*(rand(nMS,2)-0.5);
tCPU_orig=nan(1,length(nMEvec));
tCPU_stv=nan(1,length(nMEvec));
tCPU_stg=nan(1,length(nMEvec));

for jj=1:length(nMEvec)

  % generate the data
  nME=nMEvec(jj);
  tME=1:nME;
  p0=[0 0 ceil(tME(end)/2)];
  Z=100*rand(nMS,nME);
  Z(Z<nanratio_target*100)=nan;
  Z(1,:)=0;
  Z(:,1)=0;
  Zisnotnan=~isnan(Z);
  nanratio=sum(~Zisnotnan(:))/length(Zisnotnan(:));
  index_stg_to_stv=cumsum(Zisnotnan(:));
  [p_stg,z_stg]=valstg2stv(Z,sMS,tME);
  p=p_stg(~isnan(z_stg),:);
  z=z_stg(~isnan(z_stg));

  % use neighbours.m for reference
  tic;
  [psub_orig,zsub_orig,dsub_orig,nsub_orig,index_orig]=neighbours(p0,p,z,nmax,dmax);
  tCPU_orig(jj)=toc;

  % test neighbours_stv.m
  tic;
  [psub_stv,zsub_stv,dsub_stv,nsub_stv,index_stv]=neighbours_stv(p0,p,z,nmax,dmax);
  tCPU_stv(jj)=toc;
  if ~isequal(index_orig,index_stv) || ~isequal(psub_orig,psub_stv) ||  ~isequal(zsub_orig(~isnan(zsub_orig)),zsub_stv(~isnan(zsub_stv))) ||  ~isequal(dsub_orig,dsub_stv)
    error('neighbours and neighbours_stv did not give the same output');
  end

  % test neighbours_stg_nonan.m
  data.Zisnotnan=Zisnotnan; 
  data.nanratio=nanratio;
  data.sMS=sMS;
  data.tME=tME;
  data.z=z;
  data.p=p;
  data.index_stg_to_stv=index_stg_to_stv;

  tic;
  [psub_stg,zsub_stg,dsub_stg,nsub_stg,index_stg]=neighbours_stg(p0,data,nmax,dmax);
  tCPU_stg(jj)=toc;
  if ~isequal(index_stv,index_stg) || ~isequal(psub_stv,psub_stg) ||  ~isequal(zsub_stv,zsub_stg) ||  ~isequal(dsub_stv,dsub_stg)
    if isequal(sort(index_stv),sort(index_stg)) || isequal(sort(zsub_stv),sort(zsub_stg))
      fprintf('nME=%d,neighbours_stv and neighbours_stg did not give the same output, but they appear to be the same once sorted\n',nME);
    else
      fprintf('neighbours_stv and neighbours_stg did not give the same output, but they may be close. See below\n');
      dstOutput=[dsub_stv(:,1)+dmax(3)*dsub_stv(:,2) dsub_stg(:,1)+dmax(3)*dsub_stg(:,2) dsub_stv dsub_stg];
      disp('    dst_stv   dst_stg   ds_stv    dt_stv    ds_stg    dt_stg')
      disp(dstOutput(1:min(end,10),:));
    end
  end
end

figure
hold on
h_orig=plot(nMEvec,tCPU_orig,'o-g','LineWidth',1);
h_stv=plot(nMEvec,tCPU_stv,'o-b','LineWidth',1);
h_stg=plot(nMEvec,tCPU_stg,'o-r','LineWidth',1);
xlabel('nME');
ylabel('CPU time (sec)');
yyaxis right;
h_speedup=plot(nMEvec,tCPU_stv./tCPU_stg,'LineWidth',1);
ylabel('stv to stg speedup');
legend([h_orig h_stv h_stg h_speedup],{'CPU time orig','CPU time stv','CPU time stg','stv to stg speedup'});
title(sprintf('CPU time and speedup, nMS=%d, nmax=%d, nanratio=%.3f',nMS,nmax,nanratio));



