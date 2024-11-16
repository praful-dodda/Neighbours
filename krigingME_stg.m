function [zk, vk] = krigingME_stg(pk,harddata,softdata,covmodel,covparam,nhmax,nsmax,dmax,order,options)

    % krigingME                 - prediction using kriging with measurement errors using space-time grid data
    % 
    % Account for uncertainties by using a variation of kriging.
    % The uncertainties at the data locations are modeled as random
    % noises, and the resulting data can be considered as probabilitic
    % soft data. These probabilistic soft data are assumed to be
    % completely characterized by their mean and their variance. As
    % these are only second-order moments related information, a linear
    % kriging algorithm can still process adequately any mixture of these
    % soft data with hard data. 
    % 
    % SYNTAX :
    %
    % [zk,vk]=krigingME_stg(pk,harddata,softdata,covmodel,covparam,nhmax,nsmax,dmax,order,options);
    %
    % INPUT :
    % 
    
    if nargin < 10
        options(1) = 0;
    end

    % noindex = 1;

    npk = size(pk,1);
    % nsMSh = size(harddata.sMSh,1);
    % nsMSs = size(softdata.sMSs,1);

    if options(1) == 1
        num2strnk = num2str(npk);
    end

    zk = NaN(npk, 1);
    vk = NaN(npk,1);

    for i = 1:npk
        pk0 = pk(i,:);

        [phlocal, zhlocal, ~, sumnhlocal, ~] = neighbours_stg(pk0, harddata, nhmax, dmax);
        [pslocal, zslocal, ~, sumnslocal, index] = neighbours_stg(pk0, softdata, nsmax, dmax);

        vslocal = softdata.Xvs(index);

        Khh = coord2K(phlocal, phlocal, covmodel, covparam);
        Kss = coord2K(pslocal, pslocal, covmodel, covparam);
        Khs = coord2K(phlocal, pslocal, covmodel, covparam);

        Kss=Kss+diag(vslocal);                               % add the error variances on the diagonal
        kh=coord2K(phlocal,pk0,covmodel,covparam);           % built the right-hand side vector for hard data
        ks=coord2K(pslocal,pk0,covmodel,covparam);           % built the right-hand side vector for soft data

        K=[[Khh,Khs];[Khs',Kss]];                            % built the composite left-hand side matrix
        k=[kh;ks];                                           % built the composite right-hand side matrix

        chslocal=[phlocal;pslocal];
        [X,x]=krigconstr(chslocal,pk0,order);
        index=findpairs(pk0,pslocal);

        if isempty(index)
            k0=coord2K(pk0,pk0,covmodel,covparam);             % compute the variance at pk0
        end

        if (sumnhlocal+sumnslocal)>0
            nx=size(X,2);
            Kadd=[[K,X];[X',zeros(nx)]];
            kadd=[k;x];
            lam=Kadd\kadd;                                     % compute the kriging weights lam
            lam=lam(1:(sumnhlocal+sumnslocal));                % remove the Lagrangians from the solution
            lamt=lam';
            zk(i)=lamt*[zhlocal;zslocal];                      % compute the kriging estimates zk(i)
            vk(i)=k0-2*lamt*k+lamt*K*lam;                      % compute the kriging variance vk(i)
        end

        if options(1) == 1
            disp([num2str(i),'/',num2strnk]);
        end
    end
end