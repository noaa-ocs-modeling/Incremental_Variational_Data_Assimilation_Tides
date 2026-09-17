function [dxb,dxf]=DA_diff_incr_BC_TP(dd,dt,half_dt_iter,sAenodr,D5,D12,norm01,lnodes,opnod,kki,HH,HT,A11,nA12,Af1,uncerB,uncer,R);
%        	[dxb,dxf]=DA_diff_incr_BC_TP(dd,dt,half_dt_iter,A6,D12,norm01,lnodes,opnod,kki,HH,HT,A11,nA12,Af1,uncerB,uncer,R);
            % what is original boundary condition
    fprintf('DA_diff_incr_BC_TP.m\n')
            lopnod=length(opnod);
    fprintf(1,'sum(abs(dd)=%f; %d stations \n',sum(abs(dd)),length(dd));
    
    
    cent01=(A11')\HT; clear HT

%     kent02_up=nA12'*((A11')\HT);
%     kent02_dn=Af1'*((A11')\HT);
    kent02_up=nA12'*cent01;
    kent02_dn=Af1'*cent01; clear cent01
    
    AAf=[nA12,Af1];clear nA12 Af1
       fprintf('AAf done\n');
%     HWT=AAf'*((A11')\HT);
%     kent02=HWT;

  %      HB=sparse(lnodes,lopnod);for i=1:lopnod;HB(opnod(i),i)=uncerB(i);end;HB=HB';       
        I=opnod;
      J=1:lopnod;
    HB=sparse (I,J,uncerB,lnodes,lopnod);
    HB=HB';
      
        
%         Lambda=sparse(lnodes,lnodes);
%     for i=1:lnodes;Lambda(i,i)=1/norm01(i);end;
    I=1:lnodes;
   X=1./norm01;
    Lambda=sparse (I,I,X,lnodes,lnodes);
%------------------------------
%     SigmaB=sparse(lnodes,lnodes);
%     for i=1:lnodes;SigmaB(i,i)=1;end;
      X=X*0+1;
    SigmaB=sparse (I,I,X,lnodes,lnodes);
  
    
  %--------------------
%         SigmaI=sparse(lnodes,lnodes);
%     for i=1:lnodes;SigmaI(i,i)=uncer(i);end;
       SigmaI=sparse (I,I,uncer,lnodes,lnodes);
 
%        D5=A6;
        
    cent01=-dt*D12+D5;clear D12
%    tic
fprintf(1,'Iterating cent03...\n');

    cent03=Lambda'*(SigmaB'*HB'*kent02_up); 
    % adjoint 1st
  for kkk=1:half_dt_iter;
%     cent02=(cent01')\cent03;
%     cent03=(D5')*cent02;
        fprintf('%d of %d\n',kkk,half_dt_iter);
    cent03=(D5')*((cent01')\cent03);

  end;
  
  
  cent03=sAenodr*cent03;
  

  for kkk=1:half_dt_iter;
      fprintf('%d of %d\n',kkk,half_dt_iter);

%     cent02=D5*cent03;
%     cent03=cent01\cent02;
        cent03=cent01\(D5*cent03);

  end;
  
  cent04_up=HB*SigmaB*(Lambda*cent03); clear SigmaB HB
  cent03=Lambda'*(SigmaI'*kent02_dn); 
    % adjoint 1st
  for kkk=1:half_dt_iter;
%     cent02=(cent01')\cent03;
%     cent03=(D5')*cent02;
    
    fprintf('%d of %d\n',kkk,half_dt_iter);
    cent03=(D5')*((cent01')\cent03);

  end;
  
  
  cent03=sAenodr*cent03; clear sAenodr
  

  for kkk=1:half_dt_iter;
%     cent02=D5*cent03;
%     cent03=cent01\cent02;
    
    fprintf('%d of %d\n',kkk,half_dt_iter);
    cent03=cent01\(D5*cent03);

  end;
  
  cent04_dn=SigmaI*(Lambda*cent03);clear SigmaI Lambda cent03 D5 cent01

  nsta=length(HH(:,1));
  kent04=zeros(lnodes+lopnod,nsta);% kent04=zeros(lnodes+lopnod,nsta);
  kent04(1:lopnod,1:nsta)=cent04_up; clear cent04_up
  kent04(lopnod+1:lopnod+lnodes,1:nsta)=cent04_dn;
   
rep=A11\(AAf*kent04); clear A11 AAf cent04_dn cent02
cent05=HH*rep+R; clear rep R HH
kent06=cent05\dd; clear dd cent05
dz=kent04*kent06; clear kent04 kent06
dxb=dz(1:lopnod);
dxf=dz(lopnod+1:lopnod+length(kki)+lopnod);
    fprintf('End of DA_diff_incr_BC_TP.m\n')

% dx=rep*kent06; % this is the end;
%     nFeq(kki)=Feq(kki)+dxf(1:length(kki));
%     nFeq(opnod)=Feq(opnod)+dxf(length(kki)+1:length(kki)+lopnod);
