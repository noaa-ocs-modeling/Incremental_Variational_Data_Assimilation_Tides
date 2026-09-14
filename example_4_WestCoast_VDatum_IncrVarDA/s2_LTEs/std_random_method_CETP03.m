function [norm01]=std_random_method_CETP03(dt,half_dt_iter,D5,D12,lnodes,whichnode,nsta,sAenodr,meshele,meshnod,pathout)
total_rand_vector=500;
total_rand_vector=1000;
fprintf('std_random_method_CETP03\n')
cent01=-dt*D12+D5;
%    tic
% 6) norm01: diagnal element for standardization of the error covariance
    rep=zeros(lnodes,total_rand_vector);
    %rep(:,1:total_rand_vector)=normrnd(0,1,[lnodes,total_rand_vector]);% rep_ori=rep;
    rep(:,1:total_rand_vector)=randn(lnodes,total_rand_vector);% rep_ori=rep; Changed by Rachel 09/29/2021
    cent03=rep; clear rep meshele meshnod
  cent03=sqrt(sAenodr)*cent03; % this apply sqrt of the w^(-1)
    for kkk=1:half_dt_iter; 
%%%        cent03=cent01\(D5*cent03);
%         cent02=D5*cent03;
%         cent03=cent01\cent02;
        cent03=cent01\(D5*cent03);
    end;
      clear cent01;
      norm01=sqrt(mean(cent03(:,1:total_rand_vector).*cent03(:,1:total_rand_vector),2));


% 7) exact standard deviation at the observation point and update 
    rep=zeros(lnodes,nsta);
    for sta=1:nsta;rep(whichnode(sta),sta)=1;end;
    cent03=rep; clear rep
%%%    cent03=snorm01r*cent03;% normalization
        
    cent01=-dt*D12+D5; clear D12
%    tic
fprintf(1,'Iterating cent03...\n');
    % adjoint 1st
  for kkk=1:half_dt_iter;
%     cent02=(cent01')\cent03;
%     cent03=(D5')*cent02;
cent03=(D5')*((cent01')\cent03);
fprintf('%d of %d\n',kkk,half_dt_iter);

  end;

  cent03=sAenodr*cent03; clear sAenodr
  fprintf(1,'\nIterating cent03...\n');

  for kkk=1:half_dt_iter;
%     cent02=D5*cent03;
%     cent03=cent01\cent02;
    fprintf('%d of %d\n',kkk,half_dt_iter);
    cent03=cent01\(D5*cent03);

  end;
  clear D5 cent01
%%%    cent03=snorm01r*cent03;% normalization
    for sta=1:nsta;norm01(whichnode(sta))=sqrt(cent03(whichnode(sta),sta));end;

    snorm01r=sparse([1:lnodes],[1:lnodes],1./norm01,lnodes,lnodes);
   cent03=snorm01r*cent03*snorm01r(whichnode(1:nsta),whichnode(1:nsta));
 
 eval(['save ' pathout...
     '/std_random_snorm02r_non_const_diff_sqrtgh_60000.mat cent03 whichnode nsta norm01']);
  fprintf('End of std_random_method_CETP03\n')
