% this program is to prepare input files for LTEs
% We have done three test cases
% (1) Global (plug in)
% (2) San Francisco high resolution (plug in)
% (3) West Coast VDatum model
%
% For each test case: 
% input: (1) bathy file fort.14 or all.mat
%        (2) BPR harmonic data 
% internal tide dissipate input file need to use oceanmesh2d to prepare
clear

%load bpr_v130121_cm_4.mat
%load bpr_v130121_cm_4_20180928.mat
%load bpr_v130121_cm_8_20180928.mat
%load bpr_v130121_cm_13_20230406
%-------------
%recheck in 09/09/2026
load ../../../datasets/obs/DeepOceanData/check_bpr/out/bpr_v130121_cm_13_20230406.mat

casenumber=-1;% no longer use
ctestcase='global';
%ctestcase='SF';
%tlist={'M2';'S2';'K1';'O1';};
%tlist={'K2';'N2';'P1';'Q1';};
testyear=2026; %change laptop
switch ctestcase
    case 'global'
        % pathin='C:\Users\Liujuan.Tang\Documents\VDatum\data_assimilation\OceanMesh2D\OceanMesh2D-master1\OceanMesh2D-master\out_grids\';
        % gridfile='G5km6s_mEWc_20180927_depth.14';
        
        %pathin='C:\Users\Liujuan.Tang\Documents\VDatum\data_assimilation\OceanMesh2D\OceanMesh2D-Projection\out_grids\';
        %gridfile='GF5km_s1t6s_2.14';
        switch testyear
            case 2023
                pathin='C:\Users\Liujuan.Tang\Documents\VDatum\tides_assimilation\LTEs\input_data_GFEN2p5km300m_151\';
                gridfile='GFN2p5km_s1t6_300m_US1AK2HI1GU1_20181113_dn1_AE1.14';

            %----
            case 2026
                %20260912 test
                
                testcase='GFEN2p5km300m_151';
                pathin=['input_data_' testcase '/'];
                gridfile='GFN2p5km_s1t6_300m_US1AK2HI1GU1_20181113_dn1_AE1.14';
                gridfile='GFN2p5km_s1t6_300m_US1AK2HI1GU1_20181113_dn1_AE1_all.mat';
        end
    case 'SF'
        % pathin='C:\Users\Liujuan.Tang\Documents\VDatum\data_assimilation\OceanMesh2D\OceanMesh2D_20190319\OceanMesh2D-Projection\out_grids\';
        gridfile='SFWC300s3_mergall_dbc_20190417.14';
        testcase='SFWC300s3';
        pathin=['input_data_' testcase '/'];
        
    case 'WestCoast'
        % pathin='C:\Users\liujuan.tang\Documents\VDatum\WestCoast\runoutput\R89_c15_F0025R5_TAUp1_ELSMnp2_CRd3p3_HC\adcirc_output\';
        % gridfile='fort.14';
end
%----------------
%
switch testyear
    case 2023
switch casenumber
    case 1
        testcase='WCR87';
        %testcase='SFWC300s3'
%         testcase='GA2p5km151';
%         testcase='GA1p8km151';
%         testcase='GN2p5km_p9_151';
%         testcase='G300m_2p5km_151';
%         testcase='GN5km_151';
%         testcase='GF5km2_151';
%        testcase='GFTB5km_151';
%         testcase='GW4km_151';
%         testcase='GFN2p5km_151';
         testcase='GFEN2p5km300m_151';
%         testcase='GFWE2p5km300m_151';

    case 2
%-----------------------------
loc=find(xyp(:,1)>55 & xyp(:,2)>-20 & xyp(:,2)<20);
length(loc);
xyp=xyp(loc,:);
for j=1:4
    iname=tlist{j};
    eval([iname '=' iname '(loc,:);']);
end
testcase='G5k5s';

    case 3
        loc=randi([1 151],4,1);
        xyp=xyp(loc,:);
for j=1:4
    iname=tlist{j};
    eval([iname '=' iname '(loc,:);']);
end
testcase='G5kR4';
    case 4
        loc=randi([1 151],1,1);
        xyp=xyp(loc,:);
for j=1:4
    iname=tlist{j};
    eval([iname '=' iname '(loc,:);']);
end
testcase='G5kR1'; 

end
end
%----------------------

% used wrong file in the first try
switch gridfile(end-1:end)
    case '14'
        all=f_read_bathy_adcirc_grd(pathin,gridfile,0);
    case 'at'
        eval(['load ' pathin gridfile ' all'])
end

nsta=size(xyp,1);  % number of stations
snodeok=0;
snode=[];
if snodeok==0
p=all.nodes;  
p(:,[1 end])=[];
all_dmin=[];
for i=1:nsta
    ixy=xyp(i,1:2);
    dist = pos2dist_1(ixy(1),ixy(2),p(:,2),p(:,1),2)*1000; %  (m)
    [dist_min,eloc]=sort(dist);
    fprintf(1,'min = %6.0f (m)  \n',dist_min(1));
    snode(i)=eloc(1);
    all_dmin(i)=dist_min(1);
end
fprintf(1,'mean= %d; min=%d ; max = %d m\n',...
    mean(all_dmin),min(all_dmin),max(all_dmin));
save snode2 snode
else
    load snode2
end
%pp
all_dmin=all_dmin';
inloc=find(all_dmin<1e5);
nin=length(inloc);
fprintf('%d BPR found in the grid\n',nin);
%trimesh(all.elements(:,3:end),all.nodes(:,2),all.nodes(:,3),all.nodes(:,4));
hold on
%plot3(xyp(inloc,2),xyp(inloc,1),xyp(inloc,1)*0+40000,'ro')
% for i=1:nin
%     text();
% end
xyp=xyp(inloc,:);
snode=snode(inloc);

pathout=['input_data_' testcase '/'];
if ~exist(pathout,'dir')
    eval(['mkdir ' pathout])
end
    
%---------------------
% prepare station data files
%pathout='output/';    

%tlist={'M2';'S2';'K1';'O1';};
%tlist={'K2';'N2';'P1';'Q1';'M4';'MS4';'mu2';'nu2'};
tlist={'M2';'S2';'O1';'K1';
    'N2';'K2';'P1';'Q1';
    'M4';'MS4';%'MSF'; %'M6'
    'mu2';'nu2';'L2';};

temp=all.nodes(:,1)*0;
for j=1:8
    jsnode=snode;
    iname=tlist{j};
    eval(['ti=' iname '(inloc,:);']);
    loc=find(isnan(ti(:,1)));
    if ~isempty(loc)
    fprintf('%d Nan found; Removed\n',length(loc));
    %loc
    ti(loc,:)=[];
    jsnode(loc)=[];
    end
    
    for k=1:2
        ap='ap';
        outfilename=[pathout testcase '_B_' upper(iname) '_h' ap(k) '_01.txt'];
        fid=fopen(outfilename,'wt');
        fprintf(fid,'%s bpr data\n%d\n',testcase,length(ti));
        for i=1:size(ti,1)
            fprintf(fid,'%3d %5.1f %7d 0 0 0\n',i,ti(i,k),jsnode(i));
        end
        fclose(fid);
    end
    %----------------
    iname=upper(iname);
    temp=all.nodes(:,1)*0;
    temp(jsnode)=ti(:,1);
    outfile2=[pathout testcase  '_HeIn_' iname '.dat'];
    fid=fopen(outfile2,'wt');
    fprintf(fid,'%8.3f\n',temp);
    fclose(fid);
    
    temp=all.nodes(:,1)*0;
    temp(jsnode)=ti(:,2);
    
    outfile2=[pathout testcase '_HrIn_' iname '.dat'];
    fid=fopen(outfile2,'wt');
    fprintf(fid,'%8.3f\n',temp);
    fclose(fid);

%end
%--------------
% prepare inpsiLLfinal.dat
outfile=[pathout testcase '_inpsiLLfinal_' iname '.dat' ];
fid=fopen(outfile,'wt');
for i=1:size(ti,1)
    fprintf(fid,'%3d %5.1f %7d   0 %12.4f %12.4f\n',i,ti(i,1),snode(i),xyp(i,2),xyp(i,1));
end

fclose(fid);
end
%-----------------

%copy file 
%eval(['copyfile ' pathin gridfile(1:end-3) '* ' pathout])
%eval(['copyfile ' pathin gridfile ' ' pathout])

%------------------------
% prepare tide gage file