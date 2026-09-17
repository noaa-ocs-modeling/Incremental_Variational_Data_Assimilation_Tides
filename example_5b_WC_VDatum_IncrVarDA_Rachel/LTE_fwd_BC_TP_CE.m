function [PSI]=LTE_fwd_BC_TP_CE(xbc,Feq,CE,nbc,kki,A11,nA12,Af1,nA61);
% here we all agree LTE(0,0,0)=0, and it is linear. so we would use this 
% function as the LTE forwards model. Any input forcing term can be
% replaced with empty array [], if that forcing term not apply. Here we
% use the if loop to go through different option
%       if     xbc=[]; nA12=[];
%       elseif Feq=[]; Af1 =[];
%       elseif  CE=[]; nA61=[];
%       end
    PSI=zeros(length(kki)+length(nbc),1);
if ~(isempty(xbc) & isempty(Feq) & isempty(CE))
    BFib=[xbc;Feq;CE];
    AAf=[nA12 Af1 nA61];
        PSI(kki)=A11\(AAf*BFib);
        if ~isempty(xbc);PSI(nbc)=xbc;else;PSI(nbc)=zeros(length(nbc),1);end;%        PSI(nbc)=B2;
end

