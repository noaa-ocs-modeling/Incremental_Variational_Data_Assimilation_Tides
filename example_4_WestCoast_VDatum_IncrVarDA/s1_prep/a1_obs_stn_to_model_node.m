% ok here estabilish obs stn to model 1-to-1 connection
% /disks/NASUSER/lshi/file_for_Neil
% ec2012_v3d_chk.grd	 read_input.m	     sea-level-rise.txt
% Mesh.grd		 sea_level_rise.dat
% obs_stn_to_model_node.m  sea_level_rise.txt
%-----------------------------
% Changes added by Rachel 09/29/2021
% output all files to pathout: out/
%
function obs_stn_to_model_node();
% read xlsx file
%infile='/disks/NASUSER/lshi/NOAA_COOPS_tidal_datum/Copy_of_USA_WestCoast_Complete_List_xGEOID17B_copy.xlsx';
infile='C:\Users\Liujuan.Tang\Documents\VDatum\WestCoast\tide_stations\USA_WestCoast_Complete_List_xGEOID17B.xlsx';
infile='C:\Users\Liujuan.Tang\Documents\VDatum\datasets\obs\tide_stations\indata\USA_WestCoast_Complete_List_xGEOID17B.xlsx';
pathout='out/';
pathin='../../../indata/';
% infile='indata/coops149HC_WestCoast_20180309';
% eval(['load ' infile])
% [num, txt, raw]= xlsread(infile,'Harmonics');
%         %num 5510x6 for west coast
% 
% aaaa=0;
% %cent01=xyp(:[3 1 2]
% cent01=num(:,[1 3 4 5 7 8]); %remove 2 column of nan for text, keep only the numerical values
% % reform num 5510 x 7
% num=zeros(length(cent01(:,1)),7);num(:,1)=[1:length(cent01(:,1))];num(:,[2:7])=cent01;clear cent01;
% allhc=num;
%--------------
%read station table data
T=readtable(infile,'sheet',3);
sid=unique(T{:,1});% how many tide stations
nid=length(sid);
fprintf('%d tide stations\n',nid);
nhc=max(T{:,5}); % not every tide station has nhc rows
fprintf('%d max HC data\n',nhc);

nrow=nid*nhc;
allhc=zeros(nrow,7);
allhc(:,1)=[1:nrow]';
corrok=1;
 if corrok==1
        eval(['load ' pathin 'obs268_westcoast_20210928 obs'])
    end
for i=1:nid
    iid=sid(i);
    loc=find(T{:,1}==iid);
    nloc=length(loc);
    
    fprintf('%d  %d found %d HCs\n',i,iid,nloc);
    for j=1:nloc
        jrow=loc(j);
        temp=T{jrow,5};
        
        allhc((i-1)*37+temp,2:end)=T{jrow,[1 3 4 5 7 8]};
    end
    loc=find(allhc((i-1)*37+[1:37],2)==0);
    if ~isempty(loc)
        fprintf('Assign station id for missing %d HCs\n,',length(loc))
        allhc((i-1)*37+loc,2)=iid;
    end
    %-----------
    %make correction to station coordinates
    if corrok==1
        sloc=find(obs.id==iid);
        if ~isempty(sloc)
            sloc=sloc(1);
        allhc((i-1)*37+[1:37],3)=obs.x(sloc);
        allhc((i-1)*37+[1:37],4)=obs.y(sloc);   
        end       
    end
end

% load /disks/NASUSER/lshi/NOAA_COOPS_tidal_datum/Harmonics_all_correct.txt;
% load /disks/NASUSER/lshi/NOAA_COOPS_tidal_datum/Harmonics_all_16July2013_tidyup.txt;allhc=Harmonics_all_16July2013_tidyup;
lon=allhc(:,3);
lat=allhc(:,4);
stnnam=allhc(:,2);

%comment out by Rachel
%allhc=allhc(find(lon>-123.1874 & lon < -121.434 & lat > 37.3793 & lat < 38.3199 & stnnam ~= 9414958),:);
%allhc=allhc(find( stnnam ~= 9414958),:); %5473 x 7


allhc1L=allhc(1:37:length(allhc(:,1)),:); %148 x7 37 harmonic constituents
save allhc1L_new.dat allhc1L -ASCII;
stnid=allhc1L(:,2);

%unique(stnid);

% this for database Harmonics_all_16July2013_tidyup.txt; stnseq=[16 1 2 3 4 12 5 11 6 7 9 8 10 14 24 21 19 20 22 13 18 17 23];% stnseq=[16 15 1 2 3 4 12 5 11 6 7 9 8 10 14 24 21 19 20 22 13 18 17 23];
% this is for database Copy_of_USA_WestCoast_Complete_List_xGEOID17B_copy.xlsx
% cent01=zeros(37*32,7);iii=0;
% for i=[2:24 26 27 28 30 32:36];
%     iii=iii+1;
%     cent01(37*(iii-1)+1:37*iii,:)=allhc(37*(i-1)+1:37*i,:);
% end
cent01=zeros(37*nid,7);iii=0;  %5473 by 7
for i=1:nid
    iii=iii+1;
%     if i==nid
%         nid
%     end
    cent01(37*(iii-1)+1:37*iii,:)=allhc(37*(i-1)+[1:37],:);
end


allhc=cent01;
allhc1L=allhc(1:37:length(allhc(:,1)),:);
eval(['save ' pathout 'allhc1L_new.dat allhc1L -ASCII']);

ax=allhc1L(:,3);ay=allhc1L(:,4);
m=length(ax); % m=length(merged(:,1));
fprintf('%d tide stations \n',m);
%----------------------
% T.CONST_NAME(1:37)
%   37×1 cell array
%     {'M2'  }
%     {'S2'  }
%     {'N2'  }
%     {'K1'  }
%     {'M4'  }
%     {'O1'  }
%     {'M6'  }
%     {'MK3' }
%     {'S4'  }
%     {'MN4' }
%     {'NU2' }
%     {'S6'  }
%     {'MU2' }
%     {'2N2' }
%     {'OO1' }
%     {'LAM2'}
%     {'S1'  }
%     {'M1'  }
%     {'J1'  }
%     {'MM'  }
%     {'SSA' }
%     {'SA'  }
%     {'MSF' }
%     {'MF'  }
%     {'RHO' }
%     {'Q1'  }
%     {'T2'  }
%     {'R2'  }
%     {'2Q1' }
%     {'P1'  }
%     {'2SM2'}
%     {'M3'  }
%     {'L2'  }
%     {'2MK3'}
%     {'K2'  }
%     {'M8'  }
%     {'MS4' }
coops_HC_list=T.CONST_NAME(1:37);
TPC_list_adcirc={'K1';'O1';'P1';'Q1';'M2';'S2';'N2';'K2'};% the list for ADCIRC

NTIF=length(TPC_list_adcirc);
STNAMP=zeros(m,NTIF);STNPHA=zeros(m,NTIF);
for i=1:NTIF
    %TPC=i;
    NHC=strmatch(TPC_list_adcirc{i},coops_HC_list);
    fprintf('%d : %s found on COOPS list order No. %d\n',i, TPC_list_adcirc{i},NHC)
    STNAMP(:,i)=allhc(NHC:37:length(allhc(:,1)),6);
    STNPHA(:,i)=allhc(NHC:37:length(allhc(:,1)),7);
end

merged=allhc1L;

testcase='R113_LTE1R97_CF5_flow1p2';
pathrun=['C:/Users/Liujuan.Tang/Documents/VDatum/WestCoast/runoutput/' testcase '/adcirc_output/'];
fn_grd=[pathrun 'fort.14'];

[meshele, meshnod, bathy] = xmgrid2quoddy01(fn_grd);
% figure(2);trimesh(meshele(:,2:4),meshnod(:,2),meshnod(:,3));

% load /disks/NASUSER/lshi/coastline_copy/world_ascii/8072.dat;
figure(1);clf;hold on;
% plot(X8072(:,1),X8072(:,2),'color',[0.5 0.5 0.5]);hold on;

ind=zeros(m,7);indist=zeros(m,1);
for j=1:m;
    ind(j,1)=j;ind(j,2)=merged(j,6);%here is HC;%ind(j,2)=merged(j,5);% here 4 means MHHW;% here 5 means MHW;% here 6 means MLW;% here 7 means MLLW;
      ind(j,4)=merged(j,2);% ind(j,4)=merged(j,1);% COOPS station number
      ind(j,5)=merged(j,3);ind(j,6)=merged(j,4);
    cent01=(meshnod(:,2)-merged(j,3))*cos(38.0*pi()/180.0)*1000.0*6400.0*2.0*pi()/360.0;
    cent02=(meshnod(:,3)-merged(j,4))*1000.0*6400.0*2.0*pi()/360.0;
    cent03=sqrt(cent01.*cent01+cent02.*cent02);
    cent04=min(cent03);indist(j,1)=cent04;
    ind(j,7)=cent04;
    ind(j,3)=find(cent03==cent04,1,'first');
end;
eval(['save ' pathout 'indist.txt indist -ASCII']);
%
stnloc=zeros(m,3);
stnloc(:,1)=360.0+ax;stnloc(:,2)=ay;stnloc(:,3)=[1:m]';
eval(['save ' pathout 'stnloc.dat stnloc -ASCII']);

cent01=ind(:,7);
aaa=ind(cent01<5000,[1:6]);
%aaa(20,3)=3232;% 2013 database, number sequence is 15, now it is 20;aaa(15,3)=3232;% this is station 9415020; change from 3399 to 3232
%STNAMP=STNAMP(cent01<5000,:);STNPHA=STNPHA(cent01<5000,:);

%Add DART bpr data
adddartok=1;
if adddartok==1
    eval(['load ' pathin 'bpr_v130121_cm_8_20180928.mat'])
    % Find DARTs for the West Coast grid
    loc=find(xyp(:,1)<52.25 & xyp(:,1)>26 & xyp(:,2)>-132.87 & xyp(:,2)<-115.05);
    ndart=length(loc);
    fprintf('%d DART stations found\n',ndart);
    xyp=xyp(loc,:);
    for j=1:NTIF
        iname=TPC_list_adcirc{j};
        eval([iname '=' iname '(loc,:);']);
    end
    for i=1:ndart
        ixy=xyp(i,1:2);
        dist = pos2dist_1(ixy(1),ixy(2),meshnod(:,3),meshnod(:,2),2)*1000; %  (m)
        [dist_min,eloc]=sort(dist);
        fprintf(1,'min = %6.0f (m)  \n',dist_min(1));
        snode(i)=eloc(1);
        aaa(m+i,1:6)=zeros(1,6);
        %all_dmin(i)=dist_min(1);
    end
    
    aaa(m+[1:ndart],1)=[m+[1:ndart]]';
    aaa(m+[1:ndart],3)=snode';
    aaa(m+[1:ndart],4)=90000+[1:ndart]';
    aaa(m+[1:ndart],5)=meshnod(snode,2);
    aaa(m+[1:ndart],6)=meshnod(snode,3);

    for j=1:NTIF
        iname=TPC_list_adcirc{j};
        eval(['STNAMP(m+[1:ndart],j)=' iname '(:,1)/100;']);
        eval(['STNPHA(m+[1:ndart],j)=' iname '(:,2);']);
    end
end
%stnloc(m+[1:ndart],1)=xyp(:,2)+360;stnloc(m+[1:ndart],2)=xyp(:,1);stnloc(m+[1:ndart],3)=snode';
STNLOC=zeros(length(aaa(:,1)),5);
STNLOC(:,1)=aaa(:,1);
STNLOC(:,2)=aaa(:,4);
STNLOC(:,3)=aaa(:,5);
STNLOC(:,4)=aaa(:,6);
STNLOC(:,5)=aaa(:,3);
STNCLX=complex(STNAMP.*cos(STNPHA*pi()/180),    -STNAMP.*sin(STNPHA*pi()/180));

%eval(['save ' pathout 'stnloc.dat stnloc -ASCII']);

%----------------------------------
save([pathout 'STNAMPPHA.mat'],'STNAMP','STNPHA', 'STNCLX','STNLOC','TPC_list_adcirc');
        cent01=STNAMP(:,2).*cos(STNPHA(:,2)*pi()/180);
        cent02=STNAMP(:,2).*sin(STNPHA(:,2)*pi()/180);
        cent03=complex(cent01,cent02);
 
% if  adddartok==1
           % stnseq1=[34 33 20 1 2 3 4 13 5 11 6 7 9 8 10 19 31 27 25 26 29 17 24 22 30];
  stnseq1=[1:149];
        stnseq2=[150:155];
        stnseq3=[150];

  %   else
      %  end        
        plot(cent03(stnseq1),'r');hold on;
       plot(cent03(stnseq1),'r*');hold on;grid on;
       plot(cent03(stnseq2),'b');plot(cent03(stnseq2),'b^')
       %plot(cent03(stnseq2),'b');
       plot(cent03(stnseq3),'g^')
%        plot(cent03(1,1),'b*');hold on;
        grid on;axis([-1 1 -0 1]);axis equal;






function [meshele, meshnod, bathy] = xmgrid2quoddy01(filename_xmg)

% [meshele, meshnod,bathy] = xmgrid2quoddy(filename_xmg);
% Convert the xmgredit5 mesh file into quoddy arrays.
%
%filename_xmg='/disks/NASWORK/cheops/GRIDBUILDING/meshes/phil13.grid'
%filename_xmg='/disks/NASWORK/cheops/GRIDBUILDING/tom2.grid'
%filename_xmg='/disks/NASWORK/cheops/GRIDBUILDING/tomapril2003.grid'

fid = fopen(filename_xmg,'rt');

line = fgetl(fid);

[nums] = fscanf(fid,'%i %i',[2]);
numele=nums(1); 
numnode=nums(2);
A = fscanf(fid,'%f',[4,numnode])';

meshnod= A(:,1:3);
bathy=A(:,[1 4]);

A = fscanf(fid,'%f',[5,numele])';

meshele=A(:,[1 3 4 5]);

fclose(fid);
clear A
