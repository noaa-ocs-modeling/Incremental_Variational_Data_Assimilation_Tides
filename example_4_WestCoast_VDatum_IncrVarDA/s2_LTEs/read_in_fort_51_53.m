function [HC5153, hc53]=read_in_fort_51_53(infile, pathout);

fid = fopen(infile,'rt');
% line = fgetl(fid);

NTIF = fscanf(fid,'%i\n',[1]);
fprintf('fort.53 has %d harmonic constituents\n',NTIF)
hc53=[];
for k=1:NTIF;
    %line = fgetl(fid);
    hc53(k).in(1,1:3)=fscanf(fid,'%f',3);
    temp=fscanf(fid,'%s',1);
    hc53(k).name=temp([1 3]);
end;
NSTN = fscanf(fid,'%i\n',[1]);
fprintf('fort.53 has %d\n',NSTN)
HC5153=zeros(NSTN,NTIF*2);% HC53=zeros(NSTN,(NTIF-1)*6);

for i=1:NSTN;
%    i = fscanf(fid,'%i\n',[1]);
    line = fgetl(fid);
    A=fscanf(fid,'%f\n',[2 NTIF]);
    mc=0;  
    for k=[1:NTIF]; 
        mc=mc+1; % mc=0;  for k=[1 2 4 5]; mc=mc+1;
        HC5153(i,2*(mc-1)+[1:2])=A(:,k)';
    end;
end;

fclose(fid);
eval(['save ' pathout '/HC5153 HC5153 hc53'])

