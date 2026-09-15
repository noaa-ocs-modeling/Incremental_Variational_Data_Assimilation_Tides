function LTEs_head()

%-%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% V02: Developed by Lei Shi;                    
%      Total run time on O1 is ~225 sec
%                                               Last updated on 05/01/2018
% V03c: Changes include: 
%      1. Re-structrured the codes to reduce computational time near half
%      2. Added ADCIRC CPP approach
%       Total run time is ~105 sec
%                                               Last updated on 05/11/2018
% V03d: Added spherical triangle area computation  
%                                               Last updated on 05/14/2018
% V03e: Added spherical Linear SWE  approach
%       Note for 3e version, Ae_method=1 is better than 2
%                                               Last updated on 05/15/2018
% V03f: Added spherical2 Linear SWE  approach
%       Corrected one error in V03e  Ae_method=2 when computing Ae
%                                               Last updated on 05/16/2018
% V03g: Lei added K1, M2, S2;                   Last updated on 05/18/2018
% V03h: Added triplets approach to assemble matrix A (50 times faster). 
%       Total run time on O1 is ~54 sec 
%       or 34 sec w/o writing grid.14          Last updated on 05/31/2018
% V03i  Added tidal astronomic potential�forcing�terms  
%       still under testing                    Last updated on 06/04/2018
% V03j  Add global tide model test             
%       This vesrion has a phase jump between -180 and 180 in Pacific
%       Total run time on M2 global grid is 59 sec.
%                                              Last updated on 09/28/2018
% V03k  Solved the bounary jump between -180 & 180
%                                              Last updated on 10/11/2018
% V03l Added global grid without open boundary
%                                               Last updated 11/01/2018
% V03m Added inverse algorithm developed by Lei
%                                               Last updated 12/12/2018
% V03np Added Local Projection developed by Lei
%                                               Last updated 02/06/2019
% V03q Added Matlab distance calculation.
%                                               Last updated 02/06/2019
% V03r Added internal wave energy dissipation
%                                               Last updated 03/04/2019
% v03s Added new local projection and velocity Components developed 
% by Lei Shi
%                                               Last updated 03/2019
% v03t Speeded up Velocity computation. Cleaned code a bit.
%                                               Last updated 04/04/2019
%-%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% note by Lei, three features of the code
%   1. boundary condition optimization
%   2. inverse code for representor or basis function for interpolation
%   3. use a local coordinate in finite element method
%   4. diffusion operator
%   another thing need to pay attention to is the a few parameters for
%   example, is it necessary to use the diffusion operator? what is the
%   observation error covariane and backgroud error covariance matrix?
% the purpose here is to have more easy way for testing. 
% above is just something plan to do, not done yet
% now we are going to interpolate TSS using LTE by plugin. and will using
% inverse afterwards. 
% plugin: need a boundary condition
%   1) similar boundary as TCARI (Neumann boundary condition)
%   2) calculate the open boundary condition (Dirichlet boundary condition)
% this particular program will deal with neumann boundary condition
%-%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% note by Lei, recent change to the code 9/6/2019
%   1. to integrate the boundary optimization and model error adjustment
%   into single calculation.
%   2. use a random approach to approximate the background error
%-%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

testcase='Bohai';
testcase='SFB';
testcase='WC'; % that's us west coast
% testcase='GF5km2_151' % 5km global grid
% testcase='case_shen';% testcase='broken_square';
% pathin=['../input/'];
% fn_grd   = '/disks/NASUSER/lshi/file_for_Ed/LTEs/input_data/bohai_1km_2.grd';%fn_grd=[pathin 'shinfinal.14'];% fn_grd=[pathin 'fort.14.tmp'];% fn_grd=[pathin 'fort_inlet.14'];fn_grd=[pathin 'shinfinal.14'];
% fn_grd   = '/disks/NASUSER/lshi/file_for_Ed/LTEs/unified/test_bc_in/code_adcirc/results_02/fort.14';% fn_grd='/disks/NASUSER/lshi/file_for_Ed/LTEs/unified/test_diff/case_broken_square/input/fort.14.tmp';
% fn_grd   ='/disks/NASUSER/lshi/file_for_Ed/LTEs/unified/test_bc_in/code_2020_04_25/sfbofs/SFB_input/fort.14_land_bc_change_to_0_1';
fn_grd   ='/disks/NASUSER/lshi/file_for_Ed/LTEs/unified/test_bc_in/code_2020_04_25/wcadcirc/data/fort.14';
%-----------------------------------
fn_err   = ''; % '../input/data_pts_input.err';

CITE=0; % 0 = no internal tide energy conversion
        % 1 = Local method        % 2 = Nonlocal method (not available yet)
        
NLI_list={'LTE_Local';'NL_Import'};% option for linear or non-linear model as the forwards model, sometimes for linear model, code can be modified to iterate multiple times (outer loop)
NLI=NLI_list{2};                   % while for non-linear model as forwards model, the calculation of the incremental is only once before run the non-linear model again (shown in the flow chart) 

if strcmp(NLI,'NL_Import');
    NLI_MFO_list={'ADCIRC'}; % this should include a list of non-linear model, now we only read-in ADCIRC output
    NLI_MFO=NLI_MFO_list{1}; 
elseif strcmp(NLI,'LTE_Local');
    NLI_MFO_list={'LTE_MFO_BC';'LTE_MFO_Potential';'LTE_MFO_CE';...
        'LTE_MFO_BC_Potential';'LTE_MFO_BC_CE';'LTE_MFO_Potential_CE';'LTE_MFO_BC_Potential_CE'}; % this is list of LTE forwards model; forcing include BC, incremental adjustment to continuity equation, tidal potential etc.
    NLI_MFO=NLI_MFO_list{4}; NLI_MFO=NLI_MFO_list{7};
end;

ERR_list={'ERR_internal';'ERR_import'}; 
ERR=ERR_list{1};                        % ERR even it is a input in the code, but it never used
        
DAS_list={'DA_NO'; 'DA_BC';'DA_In';'DA_BC_In';...
    'DA_CE_dyn1';'DA_BC_CE_dyn1'};
DAS=DAS_list{4}; % DAS=DAS_list{4}; % only 2 option is used here, DAS=DAS_list{4} or DAS=DAS_list{2};
ICS_list={'Cartesian'; 'Spherical2';'Local2'}; % coordinate system, here we use local coordinate system
ICS=ICS_list{3};

DIS_list={'pos2dist';'pos2dist_Lei';'MGreatCircle';'wgs84'};% Matlab Great Circle and Ellipsiod, here we may use sphere to approximate the distance
dis_method=DIS_list{2};

Ae_method=1;        % method for computing triangle element area, Ae
                    % 1= plane approximation; 2=spherical triangle
CTIP=1;             % 0 = no tidal potential forcing
                    % 1 = add tidal equilibrium forcing
                    % 2 = add tidal self-attraction and loading coefficient
                    % to calculate the incremental, CTIP does not matter,
                    % it only matter for forwards model
internal_inverse=1;

cent01=load('/disks/NASUSER/lshi/file_for_Ed/LTEs/unified/test_bc_in/code_2020_04_25/sfbofs/SFB_18HC_input/alpha_coops_37HC.mat','alpha');% COOPS sequence
    aHC=cent01.alpha;clear cent01; % aHC (size 4x37), include 1) name, 2) angular speed (unit: degree/hour), 3) angular speed (unit: /s), 4) nodal factor, amplitude (non-dimensional), 5) nodal factor, phase (degree)
    whichHC=[4 6 1 2];% ['K1';'O1';'M2';'S2'];
    aHC=aHC([1 3],whichHC);
TPC_list=aHC(1,:); % TPC_list={'K1';'O1';'M2';'S2'}; % 'K(1)'    'M(2)'    'N(2)'    'O(1)'    'S(2)'    'P(1)'    'Q(1)'    'K(2)'    '2MK(3)'    'M(4)'    'MK(3)'    'M(6)'    'MS(4)'    'MN(4)'    'M(3)'    '2SM(2)'    'M(8)'    'S(4)'
% load /disks/NASUSER/lshi/file_for_Ed/LTEs/unified/test_bc_in/code_2020_04_25/sfbofs/SFB_input/STNAMPPHA.mat;
load /disks/NASUSER/lshi/file_for_Ed/LTEs/unified/test_bc_in/code_2020_04_25/wcadcirc/sta_data/STNAMPPHA.mat;
whichstn=[1:155];% include all 149 station from Rachel plus 6 deep water sensor
STNAMP=STNAMP(whichstn,whichHC);STNPHA=STNPHA(whichstn,whichHC);% row is station, column is tidal constituents (COOPS sequence)
whichnode=STNLOC(whichstn,5);% STNLOC contain location info column 1) sequence, 2) COOPS stn no, 3) lon, 4) lat, 5) ADCIRC node no corropondng to the COOPS station;
which53HC=[1 2 5 6];% the sequence of HC in fort.53 file is different with 37 COOPS HC, and often in different order, so make sure you pick correct ones % check !!!!!! you may take a look of read_in_fort_51_53(infile);what the sequence of the HC
%             infile='/disks/NASUSER/lshi/file_for_Ed/LTEs/unified/test_bc_in/code_2020_04_25/regional/test06/normal/adcirc_t02_r05/fort.53';
%             infile=['/disks/NASUSER/lshi/file_for_Ed/LTEs/unified/test_bc_in/code_2020_04_25/sfbofs/test05_18HC/r00/fort.53'];
            infile=['/disks/NASUSER/lshi/file_for_Ed/LTEs/unified/test_bc_in/code_2020_04_25/wcadcirc/data/fort.53'];
            HC5153=read_in_fort_51_53(infile); % HC5153 (size(mesh node size x 2x(no of HC)), in wc case no of HC is has 13: K1,O1,P1,Q1,M2,S2,N2,K2,M4,M6,M8,S4,S6
HCAMPstn=HC5153(whichnode,(which53HC-1)*2+1);HCPHAstn=HC5153(whichnode,(which53HC-1)*2+2);
HCstnOBS=complex(STNAMP.*cos(STNPHA*pi()/180),    -STNAMP.*sin(STNPHA*pi()/180));
HCstnMOL=complex(HCAMPstn.*cos(HCPHAstn*pi()/180),-HCAMPstn.*sin(HCPHAstn*pi()/180));
HCinno=HCstnOBS-HCstnMOL; % It compute the innovation vector

dxincr=cell(length(whichHC),2);
for i=[1:length(whichHC)] %Select TIDAL POTENTIAL CONSTITUENT
    TPC=TPC_list{i},   % NAME OF TIDAL POTENTIAL CONSTITUENT
    
    HCin=cell(6,1);
    HCin{1,1}=aHC{1,i};% name 
    HCin{2,1}=aHC{2,i};% angular speed
    HCin{3,1}=whichnode;% station node
    HCin{4,1}=HCstnOBS(:,i);
    HCin{5,1}=HCstnMOL(:,i);
    HCin{6,1}=HCinno(:,i);% HCinno, innovation vector

    %---------------------------------------------------------
    tStart=tic;
    [dxb, dxf]=LTEs(HCin,fn_grd, fn_err,NLI, NLI_MFO, ERR, DAS, TPC,ICS,Ae_method,CTIP,testcase,internal_inverse,dis_method,CITE); % dxb, xdf is incremental for BC and tidal potential 
    dxincr{i,1}=dxb;dxincr{i,2}=dxf;
%    LTEs(fn_grd, fn_err,TPC,ICS,Ae_method,CTIP,testcase,internal_inverse,dis_method,CITE);
    disp(['************** total time ',num2str(toc(tStart)),' sec*******']);
end
save dxincr_dxb_dxf_r12_nscale_1_0_R0_03.mat dxincr;% save dxincr_dxb_dxf_r12_nscale_1_0_R0_0003.mat dxincr;% just a name, sometime add the value of parameters to differentiate it
 %----------------------------------------------------------------------  

 %----------------------------------------------------------------------  
function [A, Af, A6,D12,xabc, Ae] =aen2bc_A1_5(meshele,meshnod,HLBa,HLFa,w,ICS,SFacAvg,Ae_method,meshnod_orig,CTIP,internal_inverse,H,dis_method)
% 
xabc=[];
    lnodes= length(meshnod);
    leles = length(meshele);
    %A = sparse(lnodes,lnodes);
    % Triplets approach to assemble matrix A
    ntriplets=leles*9;
    I = zeros (ntriplets, 1) ;
    J = zeros (ntriplets, 1) ;
    X = zeros (ntriplets, 1) ;
    if CTIP>0
        X2 = X;
    end
    if internal_inverse
        X6=X;
        DX12=X;
        %D2=X;
        g=9.81;
        MCO=sqrt(g*max(H,1));
        MCOa=mean(MCO(meshele(:,2:4)),2);    
    
    end
    
    ntriplets = 0 ;
    %------
    aij0=[1./6. 1./12. 1./12.;...
     1./12. 1./6. 1./12.;...
     1./12. 1./12. 1./6.];
if ~strcmp(ICS,'Local2')

    x1=meshnod(meshele(:,2),2);
    x2=meshnod(meshele(:,3),2);
    x3=meshnod(meshele(:,4),2);
    %-------------------------------
    %add in v03k for global grid
    
    %loc=find( abs(x1)>(170*pi/180)  & ~(x1>0&x2>0& x3>0) & ~(x1<0 & x2<0 & x3<0));
    % v03l for Arctic Ocean
    xx0=[x1 x2 x3 ];
    xd=max(xx0,[],2)-min(xx0,[],2);
    %Cdx=180;
    Cdx=110;
    %loc=find(  xd>Cdx*pi/180  & ~(x1>0&x2>0& x3>0) & ~(x1<0 & x2<0 & x3<0));
    loc=find(  xd>Cdx*pi/180  & ~(x1>=0&x2>=0& x3>=0) & ~(x1<=0 & x2<=0 & x3<=0));

    fprintf('%d elements found along the merged boudary \nwith delt long> %d\n', length(loc),Cdx);
    if ~isempty(loc)
        loc2=find(x1(loc)<0);
        x1(loc(loc2))=x1(loc(loc2))+360*pi/180;
        loc2=find(x2(loc)<0);
        x2(loc(loc2))=x2(loc(loc2))+360*pi/180;
        loc2=find(x3(loc)<0);
        x3(loc(loc2))=x3(loc(loc2))+360*pi/180;
    end
    
    %-----------------------------
    y1=meshnod(meshele(:,2),3);
    y2=meshnod(meshele(:,3),3);
    y3=meshnod(meshele(:,4),3);
end
    Rearth=6371.0*1000.0;
    switch ICS
        case 'Cartesian'
            dxx=-[ x2-x3    x3-x1  x1-x2];
            dyy=[ y2-y3    y3-y1  y1-y2];
             %Ae0=0.5*((meshnod(meshele(:,2),2)-meshnod(meshele(:,4),2)).*(meshnod(meshele(:,3),3)-meshnod(meshele(:,2),3))...
             %-(meshnod(meshele(:,2),2)-meshnod(meshele(:,3),2)).*(meshnod(meshele(:,4),3)-meshnod(meshele(:,2),3)));
            %Ae=0.5*((x1-x3).*(y2-y1)-(x1-x2).*(y3-y1));
            %d=Ae0-Ae;
       
 
        case 'Spherical2'
            %dxx=-Rearth*[ (x2-x3)    (x3-x1)  (x1-x2)].*cos((y1+y2+y3)/3);
            %dxx=-Rearth*[ (x2-x3).*cos(0.5*(y2+y3))   (x3-x1).*cos(0.5*(y3+y1))  (x1-x2).*cos(0.5*(y1+y2))];
            d1=x2-x3;d2=x3-x1;d3=x1-x2;
            for di=1:3
                eval(['dd=d' int2str(di) ';']);
                loc=find(abs(dd)*180/pi>180);
                if ~isempty(loc)
                    fprintf('%d need to be corrected:\n from %6.1f -> ',length(loc),dd(loc)*180/pi);
                    for ii=1:length(loc);
                        iloc=loc(ii);
                        if dd(iloc)>0
                            dd(iloc)=dd(iloc)-360*pi/180;
                        else
                            dd(iloc)=dd(iloc)+360*pi/180;
                        end
                        fprintf('%6.2f\n',dd(iloc)*180/pi);
                    end
                    eval(['d' int2str(di) '=dd;']);
                end
                
            end
            
            dxx=-Rearth*[ d1.*cos(0.5*(y2+y3))   d2.*cos(0.5*(y3+y1))  d3.*cos(0.5*(y1+y2))];
            dyy= Rearth*[ y2-y3    y3-y1  y1-y2];
      
        case 'Local2'
            % Added on 03/27/2019
            %here we use local coordinate that conserve the angle between Meridian
            % passing through the
            % centroid and great circle link nodes and centroid. 
            % step 1, calculate the centroid
%             n1 = meshele(:,2); n2=meshele(:,3); n3=meshele(:,4);
             cent01=pi()/180.0;
%             x1=cent01*orilon(n1); y1=cent01*orilat(n1);
%             x2=cent01*orilon(n2); y2=cent01*orilat(n2);
%             x3=cent01*orilon(n3); y3=cent01*orilat(n3);
            
            x1=cent01*meshnod_orig(meshele(:,2),2);
            x2=cent01*meshnod_orig(meshele(:,3),2);
            x3=cent01*meshnod_orig(meshele(:,4),2);
            y1=cent01*meshnod_orig(meshele(:,2),3);
            y2=cent01*meshnod_orig(meshele(:,3),3);
            y3=cent01*meshnod_orig(meshele(:,4),3);

            % calculate the center in a unit sphere, pay attention to degree and radius
            d31=[cos(y1).*cos(x1) cos(y1).*sin(x1) sin(y1)];
            d32=[cos(y2).*cos(x2) cos(y2).*sin(x2) sin(y2)];
            d33=[cos(y3).*cos(x3) cos(y3).*sin(x3) sin(y3)];
            cent00=(d31+d32+d33)/3;
            % now convert back to logitude and latitude
            [x4,y4,d3]=cart2sph(cent00(:,1),cent00(:,2),cent00(:,3));% d3 is radius < 1, we won't use it
            dist=zeros(leles,6);
%             x1=orilon(n1); y1=orilat(n1);
%             x2=orilon(n2); y2=orilat(n2);
%             x3=orilon(n3); y3=orilat(n3);
            x1=meshnod_orig(meshele(:,2),2);
            x2=meshnod_orig(meshele(:,3),2);
            x3=meshnod_orig(meshele(:,4),2);
            y1=meshnod_orig(meshele(:,2),3);
            y2=meshnod_orig(meshele(:,3),3);
            y3=meshnod_orig(meshele(:,4),3);

            x4=x4/cent01 ; y4=y4/cent01 ; 
            switch dis_method
                case 'pos2dist_Lei'
                    dist(:,1)=pos2dist_Lei(y1,x1,y2,x2,2);
                    dist(:,2)=pos2dist_Lei(y2,x2,y3,x3,2);
                    dist(:,3)=pos2dist_Lei(y3,x3,y1,x1,2);
                    dist(:,4)=pos2dist_Lei(y1,x1,y4,x4,2);
                    dist(:,5)=pos2dist_Lei(y2,x2,y4,x4,2);
                    dist(:,6)=pos2dist_Lei(y3,x3,y4,x4,2);% whats the mehtod???
                case 'MGreatCircle'
                    dist(:,1)=distance(y1,x1,y2,x2,'degree');
                    dist(:,2)=distance(y2,x2,y3,x3,'degree');
                    dist(:,3)=distance(y3,x3,y1,x1,'degree');% whats the mehtod???
                    dist(:,4)=distance(y1,x1,y4,x4,'degree');
                    dist(:,5)=distance(y2,x2,y4,x4,'degree');
                    dist(:,6)=distance(y3,x3,y4,x4,'degree');% whats the mehtod???

                    dist=dist*pi/180*Rearth;
                case 'wgs84'
                    ellipsoid=referenceEllipsoid('wgs84');
                    dist(:,1)=distance(y1,x1,y2,x2,ellipsoid);
                    dist(:,2)=distance(y2,x2,y3,x3,ellipsoid);
                    dist(:,3)=distance(y3,x3,y1,x1,ellipsoid);% whats the mehtod???
                    dist(:,4)=distance(y1,x1,y4,x4,ellipsoid);
                    dist(:,5)=distance(y2,x2,y4,x4,ellipsoid);
                    dist(:,6)=distance(y3,x3,y4,x4,ellipsoid);% whats the mehtod???

            end
            %    angl=zeros(leles,3)
            xa=pos2angl_Lei(y1,x1,y4,x4,dist(:,4));
            xb=pos2angl_Lei(y2,x2,y4,x4,dist(:,5));
            xc=pos2angl_Lei(y3,x3,y4,x4,dist(:,6));% whats the mehtod???
            xabc.xa=xa;
            xabc.xb=xb;
            xabc.xc=xc;
            
            %------------------
            x1=xa(:,1); y1=xa(:,2);
            x2=xb(:,1); y2=xb(:,2);
            x3=xc(:,1); y3=xc(:,2);
            
            dxx=-[ x2-x3    x3-x1  x1-x2];
            dyy=[ y2-y3    y3-y1  y1-y2];

  
    end
    %----------------------------
    % Computing Ae
    switch Ae_method
        case 1 % use plane approximation
            % x and y are in meters
            fprintf('Computing element area Ae using plane triangle\n')
            %Ae=0.5*((x1-x3).*(y2-y1)-(x1-x2).*(y3-y1));
                Ae=0.5*(-dxx(:,2).*dyy(:,3)+dxx(:,3).*dyy(:,2));
        case 2
            fprintf('Computing element area Ae using spherical triangle\n');
            fprintf('base on L''Huilier''s Theorem\n\n');
            lon1=meshnod_orig(meshele(:,2),2); % lon in deg
            lon2=meshnod_orig(meshele(:,3),2);
            lon3=meshnod_orig(meshele(:,4),2);
            
            lat1=meshnod_orig(meshele(:,2),3); % lat in deg
            lat2=meshnod_orig(meshele(:,3),3);
            lat3=meshnod_orig(meshele(:,4),3);
            
            a = pos2dist_rad(lat1,lon1,lat2,lon2,2); % arc length in rad
            b = pos2dist_rad(lat2,lon2,lat3,lon3,2); % arc length in rad
            c = pos2dist_rad(lat3,lon3,lat1,lon1,2); % arc length in rad

            s=(a+b+c)/2;
            E=sqrt(tan(s/2).*tan(0.5*(s-a)).*tan(0.5*(s-b)).*tan(0.5*(s-c)));
            Ae=atan(E)*4*Rearth^2;
            %-----------------------

    end
    fprintf(1,'Building matix A...\n');
    if isfield(HLBa,'BETA1')
        HLFa1=HLFa.f1;
        HLFa2=HLFa.f2;
        HLBa1=HLBa.BETA1;
        HLBa2=HLBa.BETA2;
    else
       HLFa1=HLFa;
       HLFa2=HLFa;
       HLBa1=HLBa;
       HLBa2=HLBa;
    end
    tic   
    %----------------------
    % Add dissipation matrix    
    for e=1:leles
        n=meshele(e,2:4);
        %-----------------------------
        %dx -> dphil/dx; dy -> dphil/dy
        dx=dyy(e,:);
        dy=dxx(e,:);
        %---------(-A1-A2)----------
        if strcmp(ICS,'Spherical1')
            aij=- 0.25 * HLBa(e) * ((dx')*(dx)/SFacAvg(e)+(dy')*(dy)*SFacAvg(e))/Ae(e);
        else
            % old
            %aij=- 0.25 * HLBa(e) * ((dx')*(dx)+(dy')*(dy))/Ae(e);% aij= HLBa(e)*Ae(e)* (dx')*(dx) /(2*Ae(e)*2*Ae(e));
            % Matrix
            aij=- 0.25 *( HLBa1(e)*(dx')*(dx)+ HLBa2(e)*(dy')*(dy))/Ae(e);% aij= HLBa(e)*Ae(e)* (dx')*(dx) /(2*Ae(e)*2*Ae(e));
        end
        %--------(-A3 + A4)-------------
        %aij=aij+ 0.25 * HLFa(e) * (-(dx')*(dy)+(dy')*(dx))/Ae(e);
        aij=aij+ 0.25  * (- HLFa1(e)*(dx')*(dy)+ HLFa2(e)*(dy')*(dx))/Ae(e);
        %---------(A5)------------------
        aij5=aij+complex(0.0,w)*Ae(e)*aij0;
        %-----------------------------
        if internal_inverse
            daij=- 0.25 * MCOa(e) * ((dx')*(dx)+(dy')*(dy))/Ae(e);% aij= HLBa(e)*Ae(e)* (dx')*(dx) /(2*Ae(e)*2*Ae(e));
            aij6=Ae(e)*aij0;
        end
        %A(n,n) = A(n,n)+aij;
        %------------------
        % Triplets approach
        for krow=1:3
            for kcol=1:3
                ntriplets = ntriplets + 1 ;
                I (ntriplets) = n (krow) ;
                J (ntriplets) = n (kcol) ;
                X (ntriplets) = aij5 (krow,kcol) ;
                if CTIP>0
                    X2 (ntriplets) = aij (krow,kcol) ;
                end
                if internal_inverse
                    DX12 (ntriplets) = daij (krow,kcol) ;
                    X6 (ntriplets) = aij6 (krow,kcol) ;
                end
            end
        end
    end
    % Add self attraction loading coefficient 
    if CTIP==2   
        X=X*(1-0.085);
    end
    A = sparse (I,J,X,lnodes,lnodes) ;
    if CTIP>0
        Af = sparse (I,J,X2,lnodes,lnodes) ;
    else
        Af=0;
    end
    if internal_inverse
        A6=sparse (I,J,X6,lnodes,lnodes) ;
        D12=sparse (I,J,DX12,lnodes,lnodes) ;
    else
        A6=0;D12=0;
    end

toc
fprintf(1,'%d ntriplet',ntriplets);
  
function [PSI,PSIs,A11,nA12,Af1,nA61,AAf,B2,BFib,F6]=getWandPSI(A,Af,A6,bcs,opnod,lnodes,Feq,CE,CTIP,NLI_MFO)
%
% one thing you have to decide, the size of W, it is size of 
%       (n-m)x m
% or     n   x m
% you have to decide. 
%
% beware the W is a sparse matrix, even though every element is non-zero
%
%
% calculate WL with a size of n x m
% beware if length of opnod and nbc=bcs(:,3) is in same length
%
% Theoretically AEN*PSI = BE
% X = A\B   is the solution to the equation A*X = B
%PSI = full(A\B);clear A B;
%
% does A need to reassemble 
%A=A1+A2;%
%T A=A1+A2+A3+A4+A5;
nbc = bcs(:,3);ksta=length(nbc);
kki=[1:lnodes]';
kki=setdiff(kki,nbc);

B2=bcs(:,2);
if ~CTIP & ~strcmp(NLI_MFO,'LTE_MFO_BC_Potential_CE');
Fib=zeros(lnodes,1);
F6=zeros(lnodes,1);
elseif CTIP & ~strcmp(NLI_MFO,'LTE_MFO_BC_Potential_CE');
Fib=Feq;
F6=zeros(lnodes,1);
elseif CTIP & strcmp(NLI_MFO,'LTE_MFO_BC_Potential_CE');
Fib=Feq;
F6=CE;
else;
    error('error: we have not reach the point yet !');
end;
BFib=[B2;Fib;F6];

A11=A(kki,kki);nA12=-A(kki,nbc);
Af1=Af(kki,:);
nA61=-A6(kki,:);
AAf=[nA12 Af1 nA61];
    PSIs=A11\(AAf*BFib);
    
%-----
PSI=zeros(lnodes,1);
PSI(kki)=PSIs;
PSI(nbc)=bcs(:,2);
        
%----------------------------------------------------------------------
  
function [A11,nA12,Af1,nA61]=partition_matrix(A,Af,A6,opnod,lnodes)
%  [PSI,PSIs,A11,nA12,Af1,nA61,AAf,B2,BFib,F6]=getWandPSI(A,Af,A6,bcs,opnod,lnodes,Feq,CE,CTIP,NLI_MFO)
%
kki=[1:lnodes]';kki=setdiff(kki,opnod);

A11=A(kki,kki);
nA12=-A(kki,opnod);
Af1=Af(kki,:);
nA61=-A6(kki,:);

 %----------------------------------------------------------------------

function     [dxb, dxf]=LTEs(HCin,fn_grd, fn_err,NLI, NLI_MFO, ERR, DAS, TPC,ICS,Ae_method,CTIP,testcase,internal_inverse,dis_method,CITE);
%%% LTEs(HCin,fn_grd, fn_err,NLI, NLI_MFO, ERR, DAS, TPC,ICS,Ae_method,CTIP,testcase,internal_inverse,dis_method,CITE)
	
% [lnodes,leles,f,Ae,meshele,meshnod,bathy,opnod]=Bohai_input_domain_preparation(testcase,fn_grd,ICS)
% halt the large scale modification for now.
dxb=[];dxf=[];

  fprintf(1,'*************LTEs PROGRAM***************\n%s TIDAL POTENTIAL CONSTITUENT\n ',TPC);
  
  disp(['grid file  : ' fn_grd]);
  disp(['error file : ' fn_err]);

%  read grid, convert into Quoddy-type format
disp([' ']); 
disp(['read grid file']);

[meshele, meshnod,bathy, opnod] = xmgrid2quoddy02(fn_grd); % for broken_square case, opnod == [], no open boundary
bathy(:,2)=max(1,bathy(:,2));
orilat=meshnod(:,3);
deg2rad=pi/180;
meshnod_orig=meshnod;

Rearth=6371.0*1000.0;
fprintf(1,'%s coordinates \n',ICS);

switch ICS
    case 'Cartesian' % all the confusion about the if the input file is in lon lat or something else. 
          disp(['since you define Cartesian, you have to provide a conversion to meter for meshnod']);
          disp(['  in an idealized domain, the coordinate is not in longitude and latitude']);
          disp(['  so in the particular ICS case Cartesian, you have to explicitly provide a conversion']);
          disp(['  take a look of the code in an idealized domain, the coordinate is not in longitude and latitude']);
        % convert the lon and lat to meters.
% % %         centlon=meshnod(:,2)-mean(meshnod(:,2));
% % %         centlat=meshnod(:,3)-mean(meshnod(:,3));
% % %             cent01=1000.0;
% % % 
% % %         centx=cent01*centlon;% centx=(2*pi()*6371.0*1000.0*cos(meshnod(:,3)*pi()/180.0)/360.0).*centlon;
% % %         centy=cent01*centlat;% centy=(2*pi()*6371.0*1000.0/360.0)*centlat;
% % %         meshnod(:,2)=centx;meshnod(:,3)=centy;
        SFacAvg=1;          
          
    case 'Local2'
        slam=meshnod(:,2)*deg2rad;
        sfea=meshnod(:,3)*deg2rad;
        SFacAvg=sfea*0+1;
end

lnodes= length(meshnod);
leles = length(meshele);

  disp([' ']);
disp(['read boundary file']);

% depth
    H=bathy(:,2);% H=74.0*ones(length(bathy(:,2)),1);
%    H=mean(H)*ones(lnodes,1);% to see if there is any different
%------------------------------------------------
% Friction fomular
CF=5;
switch CF
    case 5
        C=0.003./max(H,1);  % this is too big for O1 in Bohai Sea
         C=0.00060./max(H,1); % I use this for O1 in Bohai Sea, 
         C=0.00050./max(H,1); % I use this for O1 in Bohai Sea, LTE 
         C=0.00055./max(H,1); % I use this for O1 in Bohai Sea, LTE 
                              % the value should between 0.0005-0.0006. it
                              % for shallow water 0.003 juts too big,
                              % that's why it rely on extra enery
                              % dissipation for global model.
    case 1
        C=0.5*(0.0000375+0.00005);
        C=2.75*C./sqrt(H); % Bohai Case, looks better with this
    
    case 3
        C=0.5*(0.0000375+0.00005);
        C=2.75*C./sqrt(H)*2.5;
end

fprintf(1,'Friction coefficient CF=%d\n',CF)
%--------------------------------------------------
% tidal angular speed
switch TPC
    case 'K1';%    case 'K(1)'
        w=HCin{2,1};% 2.0*pi()/(3600*360/15.041069); % pay extra attention to it.
        TPK=0.141565;                     % Astronomical forcing amplitude (m)
        ETRF=0.736;                       % earth tide potential reduction factor
        %ETRF=0.693;                       % earth tide potential reduction factor

    case 'M2';%    case 'M(2)'
        w=HCin{2,1};% 2.0*pi()/(3600*360/28.9841042); % pay extra attention to it.
        TPK=0.242334;                     % Astronomical forcing amplitude (m)
        ETRF=0.693;                       % earth tide potential reduction factor

    case 'N2';%    case 'N(2)'
        w=HCin{2,1}; % pay extra attention to it.
        TPK=0.000000;                     % Astronomical forcing amplitude (m)
        ETRF=0.693;                       % earth tide potential reduction factor

    case 'O1';%    case 'O(1)'
        w=HCin{2,1};% 2.0*pi()/(3600*360/13.943035); % pay extra attention to it.
        TPK=0.100661;                     % Astronomical forcing amplitude (m)
        ETRF=0.695;                       % earth tide potential reduction factor

    case 'S2'%    case 'S(2)'
        w=HCin{2,1};% 2.0*pi()/(3600*360/30.0); % pay extra attention to it.
        TPK=0.112743;                     % Astronomical forcing amplitude (m)
        ETRF=0.693;                       % earth tide potential reduction factor

    case 'P1';%    case 'P(1)'
        w=HCin{2,1}; % pay extra attention to it.
        TPK=0.000000;                     % Astronomical forcing amplitude (m)
        ETRF=0.693;                       % earth tide potential reduction factor

    case 'Q1';%    case 'Q(1)'
        w=HCin{2,1}; % pay extra attention to it.
        TPK=0.000000;                     % Astronomical forcing amplitude (m)
        ETRF=0.693;                       % earth tide potential reduction factor

    case 'K2';%    case 'K(2)'
        w=HCin{2,1}; % pay extra attention to it.
        TPK=0.000000;                     % Astronomical forcing amplitude (m)
        ETRF=0.693;                       % earth tide potential reduction factor
end
fprintf(1,'w=%e\n',w);
%--------------------------
% f
g=9.81;
f=2.0*(2.0*pi/(24*3600))*sin(orilat*pi/180.0);
%---------------------------------------------
% Add internal wave dissipation option 03/05/2019

if CITE==0
% BETA
    BETA=complex(C,w);
% Lambda
    lambda=-g./(BETA.*BETA+f.*f);
% HLB
    HLB=H.*lambda.*BETA;
        HLBa=mean(HLB(meshele(:,2:4)),2);
% HLF
    HLF=H.*lambda.*f;
        HLFa=mean(HLF(meshele(:,2:4)),2);
else
    % Compute C matrix C11, C12,C22
    pathin=['../input_data_' testcase '/'];
    eval(['load ' pathin 'HxyNb_' testcase]) 
    Cdir=2.76;
    Cc=Nb*0;
    loc=find( H>10 );    
    Cc(loc)=Cdir/w/4/pi*sqrt(abs((Nb(loc).^2-w*w).*(Nm(loc).^2-w*w)));
    %-------------------
    % incread scale factor for Pacific
    %loc=find(meshnod_orig(:,2)>120 | meshnod_orig(:,2)<-90);
    loc=find(meshnod_orig(:,2)>120 | meshnod_orig(:,2)<-64);
        Cc(loc)=Cc(loc)*2;
%     figure(1)
%     clf
%     plot(meshnod_orig(:,2),meshnod_orig(:,3),'k.')
%     hold on
%     plot(meshnod_orig(loc,2),meshnod_orig(loc,3),'ro')
    %------------------
    Cc=real(Cc);
    C11=Cc.*Hx.^2;
    C22=Cc.*Hy.^2;
    C12=Cc.*(Hx.*Hy);
    % f
    f1=f-C12;f2=f+C12;
    
% BETA
     BETA1=complex(C+C11,w);BETA2=complex(C+C22,w);
% Lambda
    lambda=-g./(BETA1.*BETA2+f1.*f2);
% HLB
    HLB1=H.*lambda.*BETA1;
    HLB2=H.*lambda.*BETA2;
    HLBa.BETA1=mean(HLB1(meshele(:,2:4)),2);
    HLBa.BETA2=mean(HLB2(meshele(:,2:4)),2);
% HLF
    HLF1=H.*lambda.*f1;
    HLF2=H.*lambda.*f2;
    HLFa.f1=mean(HLF1(meshele(:,2:4)),2);
    HLFa.f2=mean(HLF2(meshele(:,2:4)),2);
end    
        g=9.81;
        %MCO=g*max(H,1);% MCO=g*max(mean(H)*ones(lnodes,1),1);%  
        MCO=sqrt(g*max(H,1));
        MCOa=mean(MCO(meshele(:,2:4)),2);    


%------------------------------------------------
%Assemble Matrix
[A, Af, A6,D12,xabc, Ae] =aen2bc_A1_5(meshele,meshnod,HLBa,HLFa,w,ICS,SFacAvg,Ae_method,meshnod_orig,CTIP,internal_inverse,MCOa,dis_method);
%     A1=aen2bc_A1(meshele,meshnod,HLBa,Ae);
%     A2=aen2bc_A2(meshele,meshnod,HLBa,Ae);
%     A3=aen2bc_A3(meshele,meshnod,HLFa,Ae);
%     A4=aen2bc_A4(meshele,meshnod,HLFa,Ae);
%     A5=aen2bc_A5(meshele,meshnod,Ae,w);
%----------------------------------------------------

%
%   further parition matrices A, Af, A6 (see notes)
%
%     HCin=cell(6,1);
%     HCin{1,1}=aHC{4,i};% name 
%     HCin{2,1}=aHC{1,i};% angular speed
%     HCin{3,1}=whichnode;% station node
%     HCin{4,1}=HCstnOBS(:,i);% HCinno
%     HCin{5,1}=HCstnMOL(:,i);% HCinno
%     HCin{6,1}=HCinno(:,i);% HCinno
	[A11,nA12,Af1,nA61]=partition_matrix(A,Af,A6,opnod,lnodes); % pay attention that it is only associate with open boundary.
    
%
%
% ok here we go
% 
%
BC_preprocess=0;
if BC_preprocess;
%   originally BC_preprocess is to optimize the BC before optamized BC and
%   tidal potential together. now it is obsolete, and removed for now
end;
       

%     HCin=cell(6,1);
%     HCin{1,1}=aHC{4,i};% name 
%     HCin{2,1}=aHC{1,i};% angular speed
%     HCin{3,1}=whichnode;% station node
%     HCin{4,1}=HCstnOBS(:,i);% HCinno
%     HCin{5,1}=HCstnMOL(:,i);% HCinno
%     HCin{6,1}=HCinno(:,i);% HCinno, innovation vector

    whichnode=HCin{3,1};% station node
    nsta=length(whichnode);

for nscale=[3];% for nscale=3:3;
    
        relaxation_scale=nscale*20000; % relaxation_scale=100000; % meter
        reference_diff_co=mean(MCOa); % m2/s % reference_diff_co=100; % m2/s
        total_time=(relaxation_scale*relaxation_scale/(2*reference_diff_co));
        half_dt_iter=8;%100;
        dt=total_time/(2*half_dt_iter);% dt=4243 sec?;
        
        [norm01]=std_random_method(dt,half_dt_iter,A6,D12,lnodes,whichnode,nsta,meshele,meshnod);

jackknifing=1;

if jackknifing; % it only works for strcmp(DAS,'DA_in') and strcmp(DAS,'DA_BC_in') 

    % jackknfing only fit for linear model run, so we still leave it here,
    % but will not further cite the array jdiff
    jdiff=zeros(nsta,nsta+2);
    jdiff(:,1)=HCin{4,1};

    for njack=0:0;njack,
%    for njack=0:length(inpsiLLfinal(:,2));njack,
%    for njack=0:21:length(inpsiLLfinal(:,2));njack,
mj=[1:njack-1 njack+1:nsta];
Jwhichnode=whichnode(mj);
Jnsta=length(Jwhichnode);
opnod;% related to bcs
% JinpsiLLfinal=inpsiLLfinal(mj,:);
% jFeq=Feq;
% jbcs=bcs;
    
    
if strcmp(DAS,'DA_NO');
    return;
elseif strcmp(DAS,'DA_BC');
    % interatively adjust the BC, now we all assume only iterate once only
    % do you think it is better to use a function instead of write here
    % here we still follow the old formulation without update for now.
    
    % now you can concentrate the on the boundary condition part with the
    % similar approach to the more broad approach.
    
        L0=2000; % 2 km
        B0=0.1;
        %L0=200; % 2 km

        fprintf(1,'L0=%6.2f;\n B0=%6.2f\n',L0,B0);
%%%         B=BCeCovariance(meshnod,bcs,L0,B0);
            for i=1:length(opnod);uncerB(i)=B0;end;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        R0=0.03;
        %R0=0.001;
        %R0=1000;

        fprintf(1,'R0=%8.4f\n',R0);
        R=OBeCovariance(meshnod,Jwhichnode,R0);%R0(25,25)=0.06;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % ok here we start the iteration
        % inpsi=JinpsiLLfinal;
%%%        check_bc_algorithm=0;
%%%        nbcs=jbcs;
        for kkk=1:1; % loop for k really not necessary since for nom-linear model, we only interation once. 
            
% % %             mpsi=PSI(inpsi(:,3));
% % %             dd=inpsi(:,2)-mpsi; % this two lines (including ^)  provide the dd that is innovation vector

            dd=HCin{6,1}; % ???

            %%% n1=inpsi(:,3); % node number for the stations
            HT=sparse(lnodes,Jnsta);
            H1=eye(Jnsta);
            HT(Jwhichnode,:)=H1;
            
            %%% nbc = nbcs(:,3);% that is the same as opnod that is, nbc=opnod
            ksta=length(opnod);kki=[1:lnodes]';kki=setdiff(kki,opnod);
            HT=HT(kki,:);HH=HT';
                       
            % this is to calculate the adjustment applyed to boundary condition;
%%%            if ~check_bc_algorithm;
[dxb]=DA_diff_incr_BC(dd,dt,half_dt_iter,A6,D12,norm01,lnodes,opnod,kki,HH,HT,A11,nA12,uncerB,R);
%%%            else;
%  %%%              [dxb,x,NPSI]=IO_dd(WL,inpsi,bcpsi,opnod,lnodes,B,R,W0L,CTIP);
%%%                [dxb]=IO_dd(dd,dt,half_dt_iter,A6,D12,norm01,lnodes,nbc,kki,HH,HT,A11,nA12,uncerB,B,R);
%%%            end;

% % %         if strcmp(NLI, 'NL_Import');
% % %             % import the non-linear solution, boundary condition, mesh grid etc
% % %             % import
% % %         elseif strcmp(NLI,'LTE_Local');
% % %             if strcmp(NLI_MFO, 'LTE_MFO_BC');
% % %             	[NPSI]=LTE_fwd_from_bc_potential(nbcs(:,2)+dxb,[],nbc,kki,A11,nA12,[]);
% % %             elseif strcmp(NLI_MFO, 'LTE_MFO_BC_Potential');
% % %          		[NPSI]=LTE_fwd_from_bc_potential(nbcs(:,2)+dxb,jFeq,nbc,kki,A11,nA12,Af1);
% % %             elseif strcmp(NLI_MFO, 'LTE_MFO_BC_Potential_CE');
% % % %         		[NPSI]=LTE_fwd_from_bc_potential(nbcs(:,2)+dxb,jFeq,nbc,kki,A11,nA12,Af1);
% % % %                [NPSI]=LTE_fwd_from_bc_potential_CE(nbcs(:,2)+dxb,Feq,nCE,nbc,kki,A11,nA12,Af1,nA61);
% % %                 [NPSI]=LTE_fwd_from_bc_potential_CE(nbcs(:,2)+dxb,Feq,CE,nbc,kki,A11,nA12,Af1,nA61);
% % %             end;
% % %         end;
% % %             
% % %     %----------------------------
% % %     % check error in Pacific, Atlantic and Indian Oceans
% % %     PSI=NPSI;
% % %     nbcs(:,2)=nbcs(:,2)+dxb;
% % % jdiff(:,njack+2)=PSI(inpsiLLfinal(:,3));
% % % % mean(abs(inpsiLLfinal(:,2)-PSI(inpsiLLfinal(:,3)))),
% % % % std(a1),
% % % % std(abs(JinpsiLLfinal(:,2))-abs(PSI(JinpsiLLfinal(:,3)))),
        end;
% % % 
% % % %        filename_end=['_' testcase '_1it_' TPC '_ITE' num2str(CITE) '_' ICS '_Ae' int2str(Ae_method) '_bs' int2str(bs) '_T' int2str(CTIP) '_inv' int2str(internal_inverse) '_' dis_method '_CF' int2str(CF) '.14'];
% % %         out14ok=1;
% % %         if out14ok==1;
% % % %             ikout14=['hre' filename_end];
% % % %             fprintf(1,'Write to %s\n',ikout14);
% % % %             writegrd(meshele,meshnod,real(NPSI),ikout14);
% % % % 
% % % %             figure(5);clf;plotpsis(meshele,meshnod,real(NPSI));colorbar;
% % % %             ikout14=['him' filename_end];fprintf(1,'Write to %s\n',ikout14);
% % % %             writegrd(meshele,meshnod,imag(NPSI),ikout14);
% % % %             figure(6);clf;plotpsis(meshele,meshnod,imag(NPSI));colorbar;
% % % 
% % % %             ikout14=['ha' filename_end];fprintf(1,'Write to %s\n',ikout14);
% % % %             writegrd(meshele,meshnod,abs(NPSI),ikout14);
% % %             figure(13);clf;plotpsis(meshele,meshnod,abs(NPSI));colorbar;
% % % 
% % % %             ikout14=['hp' filename_end];fprintf(1,'Write to %s\n',ikout14);
% % % % 
% % % %             % % % writegrd(meshele,meshnod,180*angle(NPSI)/pi(),ikout14);
% % % %             % % %   figure(2);clf;plotpsis(meshele,meshnod,180*angle(NPSI)/pi());colorbar;
% % % %             writegrd(meshele,meshnod,mod(180*(-1)*angle(NPSI)/pi(),360),ikout14);
% % %             figure(14);clf;plotpsis(meshele,meshnod,mod(180*(-1)*angle(NPSI)/pi(),360));colorbar;
% % %         end
% % % %         cent10=zeros(length(inpsiLLfinal(:,1)),9);
% % % %         cent10(:,1)=inpsiLLfinal(:,1);
% % % %         cent10(:,2)=inpsiLLfinal(:,3);
% % % %         cent10(:,3)=inpsiLLfinal(:,5);
% % % %         cent10(:,4)=inpsiLLfinal(:,6);
% % % %         cent10(:,5)=abs(inpsiLLfinal(:,2));
% % % %         cent10(:,6)=mod(180*(-1)*angle(inpsiLLfinal(:,2))/pi(),360);
% % % %         cent10(:,7)=abs(NPSI(cent10(:,2)));
% % % %         cent10(:,8)=mod(180*(-1)*angle(NPSI(cent10(:,2)))/pi(),360);
% % % %         cent10(:,9)=abs(inpsiLLfinal(:,2)-NPSI(cent10(:,2)));
% % % % 
% % % %         %save('endresults.dat','cent10','-ASCII');
% % % %         if out14ok==1
% % % %             filename=['endresults' strrep(filename_end(1:end-3),'.', 'p') '.dat'];
% % % %             save(filename,'cent10','-ASCII');
% % % %         end
% % % %         % Output in .mat format for Rachel's plotting
% % % %         if internal_inverse==0;
% % % %             filename=['output/endresults' strrep(filename_end(1:end-3),'.', 'p')];
% % % % 
% % % %             eval(['save ' filename ' cent10 NPSI  '])
% % % %         end
% % % %         toc
% % % %         disp('******************end function LTEs*****************************');
% % % % 

elseif strcmp(DAS,'DA_In');
    return;
elseif strcmp(DAS,'DA_CE_dyn1');
    return;    
elseif strcmp(DAS,'DA_BC_In');
    % the most important decision is which approach, 
    % in this particular application, we decide to use
    % solution 4. 
    % change both the tidal potential and boundary forcing
    
        L0=2000; % 2 km
        B0=0.1;
        %L0=200; % 2 km
        fprintf(1,'L0=%6.2f;\n B0=%6.2f\n',L0,B0);
%         B=BCeCovariance(meshnod,bcs,L0,B0);
            for i=1:length(opnod);uncerB(i)=B0;end;

%        uncer=1.*ones(lnodes,1); % 3 amphidromes !!!
%        uncer=0.1*ones(lnodes,1);
        uncer=0.02*ones(lnodes,1); % results total error std: 0.17m, Amplitude: 0.074 from 0.38m, and 0.10m; relaxation_scale=20000
%        uncer=0.2*ones(lnodes,1); % results total error std: 0.0157m, Amplitude: 0.0048; relaxation_scale=20000
                                   % M2 0.3782==0.1043 => 0.0883==0.0460; relaxation_scale=100000
                                   % K1 0.0652==0.0322 => 0.0400==0.0257; relaxation_scale=100000
        uncer=0.1*ones(lnodes,1); % results total error std: 0.0341m, Amplitude: 0.0158; relaxation_scale=20000
                                   % K1 0.0652==0.0322 => 0.0452==0.0272; relaxation_scale=100000
%        uncer=0.05*ones(lnodes,1); % results total error std: 0.0880m, Amplitude: 0.0395; relaxation_scale=20000
                                   % O1 0.0611==0.0337 => 0.0286==0.0250; relaxation_scale=20000
                                   % S2 0.0980==0.0372 => 0.0341==0.0164; relaxation_scale=20000
                                   % K1 0.0652==0.0322 => 0.0401==0.0229; relaxation_scale=20000
                                   % K1 0.0652==0.0322 => 0.0494==0.0285; relaxation_scale=100000
%        uncer=0.03*ones(lnodes,1); % results total error std: 0.1351m, Amplitude: 0.0596

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        R0=0.03;        %R0=0.001;        %R0=1000; %        R0=0.0001;   
%        R0=0.003;
%        R0=0.0003;
        fprintf(1,'R0=%8.4f\n',R0);
        R=OBeCovariance(meshnod,Jwhichnode,R0);%R0(25,25)=0.06;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % ok here we start the iteration
        %%% inpsi=JinpsiLLfinal;
        %%% nbcs=jbcs;
            for kkk=1:1
%         mpsi=PSI(inpsi(:,3));
%             dd=inpsi(:,2)-mpsi;

            dd=HCin{6,1}; % innovation vector

            %%% n1=inpsi(:,3); % node number for the stations
            HT=sparse(lnodes,Jnsta);
            H1=eye(Jnsta);
            HT(Jwhichnode,:)=H1;
            
            %%% nbc = nbcs(:,3);% that is the same as opnod that is, nbc=opnod
            ksta=length(opnod);kki=[1:lnodes]';kki=setdiff(kki,opnod);
            HT=HT(kki,:);HH=HT';
                       
            % this is to calculate the adjustment applyed to boundary condition;
        	[dxb,dxf]=DA_diff_incr_BC_In(dd,dt,half_dt_iter,A6,D12,norm01,lnodes,opnod,kki,HH,HT,A11,nA12,Af1,uncerB,uncer,R);
%        	[dxb,dxf,nFeq,Pa1,Pa2]=DA_diff_incr_BC_In_uncer(jFeq,dd,dt,half_dt_iter,A6,D12,norm01,lnodes,nbc,kki,HH,HT,A11,nA12,Af1,uncerB,uncer,R);
%        	[dxb,dxf,nFeq,Pa1,Pa2,Pa3,Pa4]=DA_diff_incr_BC_In_uncer_01(meshele,meshnod,jFeq,dd,dt,half_dt_iter,A6,D12,norm01,lnodes,nbc,kki,HH,HT,A11,nA12,Af1,uncerB,uncer,R);
% save('SFB_k1_dianosis.mat');
aaa=0;

% % %         if strcmp(NLI, 'NL_Import');
% % %             % import the non-linear solution, boundary condition, mesh grid etc
% % %             % import
% % %         elseif strcmp(NLI,'LTE_Local');
% % % %            if strcmp(NLI_MFO, 'LTE_MFO_BC_Potential');
% % %          		[NPSI]=LTE_fwd_from_bc_potential(nbcs(:,2)+dxb,nFeq,nbc,kki,A11,nA12,Af1);
% % % %            else;
% % % %                fprintf(1,'you do not use the results you compute');
% % % %            end;
% % %         end;
% % %         PSI=NPSI;
% % %         nbcs(:,2)=nbcs(:,2)+dxb;
% % %         jFeq=nFeq;
% % % jdiff(:,njack+2)=PSI(inpsiLLfinal(:,3));
% % % kkk,a1=inpsiLLfinal(:,2)-PSI(inpsiLLfinal(:,3));
% % % std(a1),
% % % std(inpsiLLfinal(:,2)-PSI(inpsiLLfinal(:,3))),
% % % % [[1:32]' abs(inpsiLLfinal(:,2)) abs(PSI(inpsiLLfinal(:,3))) abs(inpsiLLfinal(:,2))-abs(PSI(inpsiLLfinal(:,3))) ],
% % % % [[1:32]' mod(180*(-1)*angle(inpsiLLfinal(:,2))/pi(),360) mod(180*(-1)*angle(PSI(inpsiLLfinal(:,3)))/pi(),360) mod(180*(-1)*angle(inpsiLLfinal(:,2))/pi(),360)-mod(180*(-1)*angle(PSI(inpsiLLfinal(:,3)))/pi(),360)]
            end; % kkk
% % % %                     filename_end=['_' testcase '_1it_' TPC '_ITE' num2str(CITE) '_' ICS '_Ae' int2str(Ae_method) '_bs' int2str(bs) '_T' int2str(CTIP) '_inv' int2str(internal_inverse) '_' dis_method '_CF' int2str(CF) '.14'];
% % %             out14ok=1;
% % %         if out14ok==1
% % % %             ikout14=['hre' filename_end];
% % % %             fprintf(1,'Write to %s\n',ikout14);
% % % %             writegrd(meshele,meshnod,real(NPSI),ikout14);
% % % % 
% % % %             figure(5);clf;plotpsis(meshele,meshnod,real(NPSI));colorbar;
% % % %             ikout14=['him' filename_end];fprintf(1,'Write to %s\n',ikout14);
% % % %             writegrd(meshele,meshnod,imag(NPSI),ikout14);
% % % %             figure(6);clf;plotpsis(meshele,meshnod,imag(NPSI));colorbar;
% % % % 
% % % %             ikout14=['ha' filename_end];fprintf(1,'Write to %s\n',ikout14);
% % % %             writegrd(meshele,meshnod,abs(NPSI),ikout14);
% % %             figure(5);clf;plotpsis(meshele,meshnod,abs(NPSI));colorbar;
% % % 
% % % %             ikout14=['hp' filename_end];fprintf(1,'Write to %s\n',ikout14);
% % % % 
% % % %             % % % writegrd(meshele,meshnod,180*angle(NPSI)/pi(),ikout14);
% % % %             % % %   figure(2);clf;plotpsis(meshele,meshnod,180*angle(NPSI)/pi());colorbar;
% % % %             writegrd(meshele,meshnod,mod(180*(-1)*angle(NPSI)/pi(),360),ikout14);
% % %             figure(6);clf;plotpsis(meshele,meshnod,mod(180*(-1)*angle(NPSI)/pi(),360));colorbar;
% % %             
% % % %             figure(5);clf;plotpsis(meshele,meshnod,abs(Pa1));colorbar;
% % % %             % ikout14=['Pa1.14'];writegrd(meshele,meshnod,Pa1,ikout14);
% % % %             
% % % %              figure(6);clf;plotpsis(meshele,meshnod,abs(Pa2));colorbar;
% % % % %             ikout14=['Pa2.14'];writegrd(meshele,meshnod,Pa2,ikout14);
% % % %             
% % % %              figure(7);clf;plotpsis(meshele,meshnod,abs(Pa3));colorbar;
% % % % %             ikout14=['Pa2.14'];writegrd(meshele,meshnod,Pa2,ikout14);
% % % %             
% % % %              figure(8);clf;plotpsis(meshele,meshnod,abs(Pa4));colorbar;
% % % % %
% ikout14=['Pa2.14'];writegrd(meshele,meshnod,Pa2,ikout14); 
% % %         end    
    
elseif strcmp(DAS,'DA_BC_CE_dyn1');
    return;
end; % if strcmp(DAS,'DA_NO')

   end;% for njack=1:length(inpsiLLfinal(:,2));

end; % if jackknifing

% % % jdiff01=0*jdiff;for i=1:length(inpsiLLfinal(:,2))+2;jdiff01(:,i)=jdiff(:,i)-jdiff(:,1);end;
% % % jdiff02=abs(jdiff);
% % % jdiff03=abs(jdiff01);
% % % jdiff04=mod(180*(-1)*angle(jdiff)/pi(),360);
% % % 
% % % % save jdiff.mat jdiff jdiff01 jdiff02 jdiff03 jdiff04;
% % % % save(['jdiff_' num2str(nscale) '.mat'], 'jdiff', 'jdiff01', 'jdiff02', 'jdiff03', 'jdiff04');

end;

    return;
    %
    % here gives the final tidal field.
    %
corrPSI=NPSI; % note here 
% corrPSI=NPSI+repfw*beta; % note here 
                         % let abb=repfw*beta;
                         % then abb(inpsiLLfinal(:,3)) == dd
                         % also 
                         % corrPSI(inpsiLLfinal(:,3)) == inpsiLLfinal(:,2)
                         % that is even not so sure if representor repfw
                         % calculation is correct, the beta is correct
                         % (assuming sigmae is zero in Eq 7 at
                         % EE_JAOT_2002)
                         % inpsiLLfinal(:,2) is observation
                         % from the results, in general
                         % for O1
                         % in liaodong bay model is +, odd pts 8
                         % in bohai bay model -, odd pts 21
                         % in laizhou bay, odd pts 25
% latest developemnt: since here we are not identify the

%   figure();clf;plotpsis(meshele,meshnod,real(corrPSI));colorbar;
%   figure();clf;plotpsis(meshele,meshnod,imag(corrPSI));colorbar;

%T    figure(125);clf;plotpsis(meshele,meshnod,abs(corrPSI));colorbar;
 %T   figure(126);clf;plotpsis(meshele,meshnod,mod(180*(-1)*angle(corrPSI)/pi(),360));colorbar;
% Output in .mat format for Rachel's plotting 
errA=abs(abs(inpsiLLfinal(:,2))-abs(corrPSI(inpsiLLfinal(:,3))));
errA2=errA./abs(inpsiLLfinal(:,2))*100;
mean(abs(inpsiLLfinal(:,2)-corrPSI(inpsiLLfinal(:,3)))),
fprintf('--------------------------------------------------------\n')
fprintf('Inversed error at %d stations in amp: %6.3f m; \navg %6.3fm; %2.0f %%; rmse=%6.3f m\n',...
        length(errA),sum(errA),mean(errA),mean(errA2),sqrt(mean((errA.^2))));
%==============================================
 % here we start to calculate the velocity, 
% here we would follow equation_2019_02_22.doc notation
% the bottom friction should be 
%   tau_bottom=c(x,y)*V0*(u,v)
%       c(x,y)=0.003/max(h,1.0)
%       V0=max_ V*pi()/4
% if only consider M2, K1, then 
%       max_V=max((am2*am2+bm2*bm2)/(am2+bm2),(ak1*ak1+bk1*bk1)/(ak1+bk1),0.5)
%       ???
% where D is the water depth, 
%       a and b is the major and minor axis of the tida ellipse.
% and cxx, cxy and cyy would be the dissipation due to the internal wave
% 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % function try_to_calculate_some_matrix_and_vectors

% re-assign the depth
% % %     bathy(:,2)=74.00*ones(size(bathy(:,2)));
   
    
% % % x1=meshnod(meshele(:,2),2);x2=meshnod(meshele(:,3),2);x3=meshnod(meshele(:,4),2);
% % % y1=meshnod(meshele(:,2),3);y2=meshnod(meshele(:,3),3);y3=meshnod(meshele(:,4),3);

h1=corrPSI(meshele(:,2));  h2=corrPSI(meshele(:,3));  h3=corrPSI(meshele(:,4));
    % here is the calculation of slopes.
xa=xabc.xa;
xb=xabc.xb;
xc=xabc.xc;

[Fx,Fy,Aee,Arr]= element03(xa(:,1),xb(:,1),xc(:,1),xa(:,2),xb(:,2),xc(:,2),h1,h2,h3);

% here we calculate the coefficients

Hele=H(meshele(:,2:4));
cele=0.003./max(Hele, 1.0);
fele=f(meshele(:,2:4));
cxx=zeros(size(meshele(:,2:4)));cxy=cxx;cyy=cxx;
betaxele=complex(cele+cxx,w);
betayele=complex(cele+cyy,w);
f1ele=-fele+cxy;
f2ele=fele+cxy;
lambdaele=-9.81./(betaxele.*betayele-f1ele.*f2ele);
c11=betayele.*lambdaele;
c12=-f1ele.*lambdaele;
c21=-f2ele.*lambdaele;
c22=betaxele.*lambdaele;

uele=mean(c11,2).*Fx+mean(c12,2).*Fy;
vele=mean(c21,2).*Fx+mean(c22,2).*Fy;
Aele=elementspernode01a(meshele); 
fprintf(1,'Computing u and v\n')

cent03=Aele*Aee;
PSIu=Aele*(uele.*Aee)./cent03;
PSIv=Aele*(vele.*Aee)./cent03;


figure(11);clf;plotpsis(meshele,meshnod,abs(PSIu));colorbar;
figure(12);clf;plotpsis(meshele,meshnod,abs(PSIv));colorbar;
figure(13);clf;plotpsis(meshele,meshnod,sqrt(abs(PSIu).*abs(PSIu)+abs(PSIv).*abs(PSIv)));colorbar;
%-------------------------------------------
% filename=['output/endresults' strrep(filename_end(1:end-3),'.', 'p')];
% filename=[filename '_C' int2str(-log10(Csigma)) '_uv'];
% eval(['save ' filename ' cent10 NPSI corrPSI dt_iter PSIu PSIv '])
% fprintf(1,'Save to file\n%s\n',filename)

%-------------------------------------------


 
function    xa=pos2angl_Lei(y1,x1,y4,x4,dist);
% function    [xa,eab]=pos2angl_Lei(y1,x1,x4,y4,dist);
% assume that the coordinate is only one point
a1=[x4 y4];a2=[x4 y4+0.02];
b1=[x4 y4];b2=[x1 y1];
cent01=pi()/180.0;

a1=cent01*a1;a2=cent01*a2;
b1=cent01*b1;b2=cent01*b2;
d3a1=[cos(a1(:,2)).*cos(a1(:,1)) cos(a1(:,2)).*sin(a1(:,1)) sin(a1(:,2))];
d3a2=[cos(a2(:,2)).*cos(a2(:,1)) cos(a2(:,2)).*sin(a2(:,1)) sin(a2(:,2))];
d3b1=d3a1;
d3b2=[cos(b2(:,2)).*cos(b2(:,1)) cos(b2(:,2)).*sin(b2(:,1)) sin(b2(:,2))];
uva=cross(d3a1,d3a2);
uvb=cross(d3b1,d3b2);
uvc=cross(uva,uvb);
normuvc=sqrt(uvc(:,1).*uvc(:,1)+uvc(:,2).*uvc(:,2)+uvc(:,3).*uvc(:,3));
uvd=d3a1;
% % % eab=sign(dot(uvc,uvd,2))*atan2(norm(uvc),dot(uva,uvb,2));
eab=sign(dot(uvc,uvd,2)).*atan2(normuvc,dot(uva,uvb,2));
xa=[dist.*cos(eab+0.5*pi()) dist.*sin(eab+0.5*pi())];
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------

function [Fx,Fy,Ae,Ar]= element03(x1,x2,x3,y1,y2,y3,h1,h2,h3)
% [ a,b,c,A]= element(x,y,h)
% [ a,b,c,A]= element(meshnod(meshele(m,2:4),2:3),bathy(meshele(m,2:4),2))
%  Given the triangular element triplets of 
% xy(3,2), h(3) 
% h(x,y) =a + bx + cy
%
%[ a,dx,dy,A]= element(xy,h)
% h(x,y) =a + dx*x + dy*y
% A is area of the triangle

%x = xy(:,1);
%y = xy(:,2);

% From A simple guide to finite elements
% A = .5*det([ones(3,1) x y]);
% Ar = 1./(2.*A);
%  A=  0.5*( x(2).*y(3)-x(3).*y(2)   +   x(3).*y(1)-x(1).*y(3)   +   x(1).*y(2)-x(2).*y(1) );
%  A=  x(2).*y(3)-x(3).*y(2)   +   x(3).*y(1)-x(1).*y(3)   +   x(1).*y(2)-x(2).*y(1) );
%  x1=2,y1=0,h1=0,x2=0,y2=1,h2=0,x3=0,y3=0,h3=1,
 cent00 =  x1.*(y2-y3)+x2.*(y3-y1)+x3.*(y1-y2);
 cent01 =  h1.*(y2-y3)+h2.*(y3-y1)+h3.*(y1-y2);
 cent02 =  x1.*(h2-h3)+x2.*(h3-h1)+x3.*(h1-h2);
 cent00r=1./cent00;
 Fx=cent00r.*cent01;
 Fy=cent00r.*cent02;
 Ae=0.5*cent00;
 Ar=2.0*cent00r;
 
%-----------------------------------
function [rep]=deltaf(A,ins,bcs,opnod,lnodes);
% calculation of representator function
% you need to under =satnd it is about the forcing not for the boundary
% condition
    nsti=length(ins(:,1)); %nubmer of stations for observatopm;
    rep=zeros(lnodes,nsti);
    %rep=sparse(lnodes,nsti);
    %rep=sparse([],[],[],lnodes,nsti,lnodes*nsti);
    %-------------------------
    cent1=A';
    %nsta=length(bcs(:,1));
    if ~isempty(bcs)
    indices = sub2ind(size(cent1), bcs(:,3), bcs(:,3));
    cent1(bcs(:,3),:)=0;
    cent1(indices)=1;
    end
%     i=1;
%     indices1 = sub2ind(size(cent1), 1:lnodes, [1:lnodes]*0+i);
    %-----
   % tic
    for i=1:nsti
        %fprintf('%d\n',i);
%         dlt=dlt0;
%         dlt(ins(i,3))=1;
        dlt=sparse(ins(i,3),1,1,lnodes,1);
%    dlt=A6*dlt;
    %cent1=A1+A2+A3+A4+A5;
   %T cent1=A';
    %cent1=cent1';
%T    cent1(bcs(:,3),:)=0;
%          for j=1:nsta
%              cent1(bcs(j,3),bcs(j,3))=1;
%          end
%T    cent1(indices)=1;
%         toc
%    cent2=A6';
%    rep(:,i)=cent2*(cent1\dlt);
    %ii=indices1+(i-1)*lnodes;
   
  %  fprintf(1,'%d',i)
  
    rep(:,i)=cent1\dlt;
    end;
  rep=sparse(rep);
    
%     tic
%     dlt=zeros(lnodes,nsti);
%     dlt(ins(:,3),ins(:,3))=1;
%     rep2=cent1\dlt;
%     dd=rep-rep2;
%     max(max(dd))
%     min(min(dd))
%     toc
%function [repfw]=deltaf_bc_zero_forward(A1,A2,A3,A4,A5,A6,rep,ins,bcs,opnod,lnodes);
function [repfw]=deltaf_bc_zero_forward(A,A6,rep,ins,bcs,opnod,lnodes);
% calculation of representator function
% you need to under =satnd it is about the forcing not for the boundary
% condition
    nsti=length(ins(:,1));
    nsta=length(bcs(:,1));
    repfw=zeros(lnodes,nsti);
    
    %---------------
    cent1=A;
    if ~isempty(bcs)
    indices = sub2ind(size(cent1), bcs(:,3), bcs(:,3));
    cent1(bcs(:,3),:)=0;
    cent1(indices)=1;
    end
    for i=1:nsti
%     dlt=zeros(lnodes,1);
%     dlt(bcs(i,3))=1;
    dlt=rep(:,i);
    dlt=A6*dlt;
    dlt(bcs(:,3))=0;
    %cent1=A1+A2+A3+A4+A5;
    %Tcent1=A;
    %Tcent1(bcs(:,3),:)=0;
    
    %for j=1:nsta;cent1(bcs(j,3),bcs(j,3))=1;end;
    %Tcent1(indices)=1;
    repfw(:,i)=cent1\dlt;
% % %     cent2=speye(lnodes,lnodes)';
% % %     rep(:,i)=cent2*(cent1\dlt);
    end;

function [beta]=coin(repfw,ins,dd,lnodes,Csigma);
% calculation of representator function
% you need to under =satnd it is about the forcing not for the boundary
% condition
    nsti=length(ins(:,1));
    R=repfw(ins(:,3),:);
%    sigmae=0.025*eye(nsti,nsti);
 %   sigmae=0.0*eye(nsti,nsti);
  %  sigmae=0.001*eye(nsti,nsti);
    sigmae=Csigma*eye(nsti,nsti);
%    discrep=bcs(:,2); it is dd
    beta=(R+sigmae)\dd;
%--------------------end of anding inverse functions    
    
function dist = pos2dist_Lei(lag1,lon1,lag2,lon2,method)
% function dist = pos2dist(lag1,lon1,lag2,lon2,method)
% calculate distance between two points on earth's surface
% given by their latitude-longitude pair.
% Input lag1,lon1,lag2,lon2 are in degrees, without 'NSWE' indicators.
% Input method is 1 or 2. Default is 1.
% Method 1 uses plane approximation,
% only for points within several tens of kilometers (angles in rads):
% d =
% sqrt(R_equator^2*(lag1-lag2)^2 + R_polar^2*(lon1-lon2)^2*cos((lag1+lag2)/2)^2)
% Method 2 calculates sphereic geodesic distance for points farther apart,
% but ignores flattening of the earth:
% d =
% R_aver * acos(cos(lag1)cos(lag2)cos(lon1-lon2)+sin(lag1)sin(lag2))
% Output dist is in km.
% Returns -99999 if input argument(s) is/are incorrect.
% Flora Sun, University of Toronto, Jun 12, 2004.
% % % if nargin < 4
% % %     dist = -99999;
% % %     disp('Number of input arguments error! distance = -99999');
% % %     return;
% % % end
% % % if abs(lag1)>90 | abs(lag2)>90 | abs(lon1)>360 | abs(lon2)>360
% % %     dist = -99999;
% % %     disp('Degree(s) illegal! distance = -99999');
% % %     return;
% % % end
% % % if lon1 < 0
% % %     lon1 = lon1 + 360;
% % % end
% % % if lon2 < 0
% % %     lon2 = lon2 + 360;
% % % end
% % % % Default method is 1.
% % % if nargin == 4
% % %     method == 1;
% % % end
% % % if method == 1
% % %     km_per_deg_la = 111.3237;
% % %     km_per_deg_lo = 111.1350;
% % %     km_la = km_per_deg_la * (lag1-lag2);
% % %     % Always calculate the shorter arc.
% % %     if abs(lon1-lon2) > 180
% % %         dif_lo = abs(lon1-lon2)-180;
% % %     else
% % %         dif_lo = abs(lon1-lon2);
% % %     end
% % %     km_lo = km_per_deg_lo * dif_lo.* cos((lag1+lag2)*pi/360);
% % %     dist = sqrt(km_la.^2 + km_lo.^2);
% % % else
    R_aver = 6371000.0; % R_aver = 6374;
    deg2rad = pi/180;
    lag1 = lag1 * deg2rad;
    lon1 = lon1 * deg2rad;
    lag2 = lag2 * deg2rad;
    lon2 = lon2 * deg2rad;
    dist = R_aver * acos(cos(lag1).*cos(lag2).*cos(lon1-lon2) + sin(lag1).*sin(lag2));
    
   
function dist = pos2dist_rad(lag1,lon1,lag2,lon2,method)
% function dist = pos2dist(lag1,lon1,lag2,lon2,method)
% calculate distance between two points on earth's surface
% given by their latitude-longitude pair.
% Input lag1,lon1,lag2,lon2 are in degrees, without 'NSWE' indicators.
% Input method is 1 or 2. Default is 1.
% Method 1 uses plane approximation,
% only for points within several tens of kilometers (angles in rads):
% d =
% sqrt(R_equator^2*(lag1-lag2)^2 + R_polar^2*(lon1-lon2)^2*cos((lag1+lag2)/2)^2)
% Method 2 calculates sphereic geodesic distance for points farther apart,
% but ignores flattening of the earth:
% d =
% R_aver * acos(cos(lag1)cos(lag2)cos(lon1-lon2)+sin(lag1)sin(lag2))
% Output dist is in km.
% Returns -99999 if input argument(s) is/are incorrect.
% Flora Sun, University of Toronto, Jun 12, 2004.
if nargin < 4
    dist = -99999;
    disp('Number of input arguments error! distance = -99999');
    return;
end
if abs(lag1)>90 | abs(lag2)>90 | abs(lon1)>360 | abs(lon2)>360
    dist = -99999;
    disp('Degree(s) illegal! distance = -99999');
    return;
end
if lon1 < 0
    lon1 = lon1 + 360;
end
if lon2 < 0
    lon2 = lon2 + 360;
end
% Default method is 1.
if nargin == 4
    method == 1;
end
if method == 1
    km_per_deg_la = 111.3237;
    km_per_deg_lo = 111.1350;
    km_la = km_per_deg_la * (lag1-lag2);
    % Always calculate the shorter arc.
    if abs(lon1-lon2) > 180
        dif_lo = abs(lon1-lon2)-180;
    else
        dif_lo = abs(lon1-lon2);
    end
    km_lo = km_per_deg_lo * dif_lo.* cos((lag1+lag2)*pi/360);
    dist = sqrt(km_la.^2 + km_lo.^2);
else
    %R_aver = 6374;
    deg2rad = pi/180;
    lag1 = lag1 * deg2rad;
    lon1 = lon1 * deg2rad;
    lag2 = lag2 * deg2rad;
    lon2 = lon2 * deg2rad;
    %dist = R_aver * acos(cos(lag1).*cos(lag2).*cos(lon1-lon2) + sin(lag1).*sin(lag2));
    
    % Out put arc length in radians --Liujuan.tang
    dist = acos(cos(lag1).*cos(lag2).*cos(lon1-lon2) + sin(lag1).*sin(lag2));
end
    
  
%-----------------------------------------------------------------------

function B1 = elementspernode(A1,A2)
% Aele=elementspernode(meshele);
% le = elementspernode(Aele, node );
% Returns a list of the elements which share node 
% Aele [ nodes, list of elements] 

if nargin==1
 % initial call  create Aele from meshele
 meshele = A1;
 nnode= max(max(meshele(:,2:4)));
 A = sparse(nnode,length(meshele(:,1)));
 for e = 1:length(meshele(:,1))
   A(meshele(e,2:4),e) = e;
 end
 Aele=zeros(nnode,10);
 for n=1:nnode
    l = find(A(n,:)>0);
    nl = length(l);
    Aele(n,1:nl) = A(n,l);
 end
 
 
 
 B1 = Aele;
else
 % Second call with two parameters;
 % Aele, n
 Aele = A1;
 n = A2;
 
 le = Aele(n,find(Aele(n,:)>0));

 B1 = le;
end

%-----------------------------------------------------------------------

function B1 = elementspernode01a(meshele)
% Aele=elementspernode(meshele);
% le = elementspernode(Aele, node );
% Returns a list of the elements which share node 
% Aele [ nodes, list of elements] 

 % initial call  create Aele from meshele
 %meshele = A1;
 nnode= max(max(meshele(:,2:4)));
 ne = length(meshele);

 ntriplets=ne*3;
 I = zeros (ntriplets, 1) ;
 J = zeros (ntriplets, 1) ;
 X = zeros (ntriplets, 1) ;

fprintf('Computing A...\n')
ntriplets=0;
for e = 1:ne
    for k=1:3
        ntriplets = ntriplets + 1 ;
        I (ntriplets) = meshele(e,1+k) ;
        J (ntriplets) = e ;
        %X (ntriplets) = e ;
        X (ntriplets) = 1 ;
    end

end
B1 = sparse (I,J,X,nnode,ne) ;
%-----------------------------------------------------------------------

function B1 = elementspernode01(A1)
% Aele=elementspernode(meshele);
% le = elementspernode(Aele, node );
% Returns a list of the elements which share node 
% Aele [ nodes, list of elements] 

 % initial call  create Aele from meshele
 meshele = A1;
 nnode= max(max(meshele(:,2:4)));
 A = sparse(nnode,length(meshele(:,1)));
 for e = 1:length(meshele(:,1))
   A(meshele(e,2:4),e) = e;
 end
 Aele=zeros(nnode,10);
 for n=1:nnode
    l = find(A(n,:)>0);
    nl = length(l);
    Aele(n,1:nl) = A(n,l);
 end
 B1 = Aele;
%-----------------------------------------------------------------------

function B1 = elementspernode02(A1,A2)
% Aele=elementspernode(meshele);
% le = elementspernode(Aele, node );
% Returns a list of the elements which share node 
% Aele [ nodes, list of elements] 
 % Second call with two parameters;
 % Aele, n
 Aele = A1;
 n = A2;
 
 le0 = Aele(n,find(Aele(n,:)>0));
 le=nonzeros(Aele(n,:));
 B1 = le;

%-----------------------------------------------------------------------

function [ a,Fx,Fy,A,dFxdhi,dFydhi]= element(xy,h)
% [ a,b,c,A]= element(x,y,h)
% [ a,b,c,A]= element(meshnod(meshele(m,2:4),2:3),bathy(meshele(m,2:4),2))
%  Given the triangular element triplets of 
% xy(3,2), h(3) 
% h(x,y) =a + bx + cy
%
%[ a,dx,dy,A]= element(xy,h)
% h(x,y) =a + dx*x + dy*y
% A is area of the triangle
x = xy(:,1);
y = xy(:,2);

x = xy(:,1);
y = xy(:,2);

% From A simple guide to finite elements
 A = .5*det([ones(3,1) x y]);
 b= [ x(2)*y(3)-x(3)*y(2) 
      x(3)*y(1)-x(1)*y(3)
      x(1)*y(2)-x(2)*y(1)];
 c = [ y(2)-y(3)
       y(3)-y(1)
       y(1)-y(2)];
 d = [ x(3)-x(2)
       x(1)-x(3)
       x(2)-x(1)];
 
 a         = 1.0/(2.*A) *(b'*h);
 Fx        = 1.0/(2.*A) *(c'*h);
 Fy        = 1.0/(2.*A) *(d'*h);


% Fx = dx
dFxdhi= 1.0/(2.*A) *c;
dFydhi= 1.0/(2.*A) *d;

% Matlab technique:
%D = [ones(3,1) x y];
%A = D\h;
%a=A(1);dx=A(2); dy=A(3);

% area
%A = .5*det(D);
%-----------------------------------------------------------------------

%-----------------------------------------------------------------------

function AE= elementarea(xy)

% AE is area of the triangle
%disp('within element01');
%xy,h,
x = xy(:,1);
y = xy(:,2);
%disp('within element01');

% From AE simple guide to finite elements
 AE = .5*det([ones(3,1) x y]);

% Matlab technique:
%D = [ones(3,1) x y];
%AE = D\h;
%a=AE(1);dx=AE(2); dy=AE(3);

% area
%AE = .5*det(D);
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------

function [Fx,Fy,A]= element01(xy,h)
% [ a,b,c,A]= element(x,y,h)
% [ a,b,c,A]= element(meshnod(meshele(m,2:4),2:3),bathy(meshele(m,2:4),2))
%  Given the triangular element triplets of 
% xy(3,2), h(3) 
% h(x,y) =a + bx + cy
%
%[ a,dx,dy,A]= element(xy,h)
% h(x,y) =a + dx*x + dy*y
% A is area of the triangle
%disp('within element01');
%xy,h,
x = xy(:,1);
y = xy(:,2);
%disp('within element01');

% From A simple guide to finite elements
 A = .5*det([ones(3,1) x y]);
 Ar = 1./(2.*A);
% b= [ x(2)*y(3)-x(3)*y(2) 
%      x(3)*y(1)-x(1)*y(3)
%      x(1)*y(2)-x(2)*y(1)];
 c = [ y(2)-y(3)
       y(3)-y(1)
       y(1)-y(2)];
 d = [ x(3)-x(2)
       x(1)-x(3)
       x(2)-x(1)];
%   c,d,h,
% size(c),size(h),
% a         = Ar *(b'*h);
 Fx        = Ar *(c'*h);
 Fy        = Ar *(d'*h);



% Matlab technique:
%D = [ones(3,1) x y];
%A = D\h;
%a=A(1);dx=A(2); dy=A(3);

% area
%A = .5*det(D);
%-----------------------------------------------------------------------

function [Fx,Fy,A]= element02(xy,h)
% [ a,b,c,A]= element(x,y,h)
% [ a,b,c,A]= element(meshnod(meshele(m,2:4),2:3),bathy(meshele(m,2:4),2))
%  Given the triangular element triplets of 
% xy(3,2), h(3) 
% h(x,y) =a + bx + cy
%
%[ a,dx,dy,A]= element(xy,h)
% h(x,y) =a + dx*x + dy*y
% A is area of the triangle

%x = xy(:,1);
%y = xy(:,2);

% From A simple guide to finite elements
 A = .5*det([ones(3,1) x y]);
 Ar = 1./(2.*A);
% b= [ x(2)*y(3)-x(3)*y(2) 
%      x(3)*y(1)-x(1)*y(3)
%      x(1)*y(2)-x(2)*y(1)];
 c = [ xy(2,2)-xy(3,2)       xy(3,2)-xy(1,2)       xy(1,2)-xy(2,2)];
 d = [ xy(3,1)-xy(2,1)       xy(1,1)-xy(3,1)       xy(2,1)-xy(1,1)];
 
% a         = Ar *(b'*h);
 Fx        = Ar *(c'*h);
 Fy        = Ar *(d'*h);



% Matlab technique:
%D = [ones(3,1) x y];
%A = D\h;
%a=A(1);dx=A(2); dy=A(3);

% area
%A = .5*det(D);
%-----------------------------------------------------------------------

function [bndeleLei,polys,S] = findbndy(meshele,meshnod,plottingon)
%function [bndelenew,polys,S] = findbndy(meshele,meshnod,plottingon)
% [bndelenew] = findbndy(meshele,meshnod);
% [bndelenew] = findbndy(meshele,meshnod,plottingon);
% This finds the boundary elements of meshele
%  I.E.  all segments with only one element
disp(['mesh/meshtools/findbndy:  Finding the boundary, Watch for red triangles'])

if nargin==2
plottingon = 0;
else
plottingon = 0; %1;
end

%  Make sure the smallest index comes first for all meshelements
for i = 1:length(meshele)
if (meshele(i,2)<meshele(i,3) & meshele(i,2)<meshele(i,4) )
%  Correct order
elseif (meshele(i,3)<meshele(i,2) & meshele(i,3)<meshele(i,4) )
me2 = meshele(i,2);
meshele(i,2) = meshele(i,3);
meshele(i,3) = meshele(i,4);
meshele(i,4) = me2;
elseif (meshele(i,4)<meshele(i,3) & meshele(i,4)<meshele(i,2) )
me2 = meshele(i,2);
meshele(i,2) = meshele(i,4);
meshele(i,4) = meshele(i,3);
meshele(i,3) = me2;
end
end

polys = meshele;    
polys(:,1)=meshele(:,4);   %   triangles close on themselves

[m,n]=size(polys);
c = ones(4*m,1);
c(:) = polys;
l = min(c):max(c);    % this will be harder    
mc = min(c)-1 ;

%size(meshnod),l,this is l not 1.
xy=meshnod(l,2:3);
i1 = polys(:,1)-mc;
j1 = polys(:,2)-mc;
i2 = polys(:,2)-mc;
j2 = polys(:,3)-mc;
i3 = polys(:,3)-mc;
j3 = polys(:,4)-mc;
i4 = polys(:,4)-mc;
j4 = polys(:,1)-mc;
S = sparse([i1;i2;i3],[j1;j2;j3],1);
%spy(S,'g')

St = S';
Sd = S-St;
%hold on
%spy(Sd,'r')
[m,n]=size(Sd);

disp(['start loop on non zeros'])

[i,j,junk]=find(Sd==1);
bndele=[ [1:length(i)]' i j];bndeleLei=[ [1:length(i)]' i j,zeros(length(i),1) ones(length(i),1)];

if plottingon == 1
figure(2)
clf
gplot(S,xy,'w')
hold on
plot([meshnod(i,2) meshnod(j,2)]',[meshnod(i,3) meshnod(j,3)]','g')
hold off
axis('equal')
orient landscape
grid
end



bndelenew = bndele;
n = length(bndele);
% ld = find(meshnod(:,3)<10000);
l =  find(meshnod(:,2)==min(meshnod(:,2))) ;
%start with the smallest x 

if length(l)>1
l = l(1);
end

ilast = find(bndele(:,2)==l);
%ilast = 1
bndelenew(1,:)= bndele(ilast,:);

N = n;
for j = 2:n
l=find(bndele(:,2)==bndele(ilast,3));

if length(l)==1
bndelenew(j,:)= bndele(l,:);
bndele(ilast,:) = -bndele(ilast,:);
ilast = l;
else
N = N-1;
    if plottingon == 1
    disp(['length(l)~=1  ',int2str(ilast),'  ', int2str(bndele(ilast,3))...
     ,'  N =',int2str(N)])
    end
end

end
bndelenew = bndelenew(1:N,:);
bndelenew(:,1) = [1:N]';
n = N;

if plottingon == 1


%disp('size bndelenew')
%size(bndelenew)

figure(3)
clf
ax = [ -2000 2000 21500 23500 ; 6000 8000 19500 21500; ...
  -2000 2000 0 4000;  18000 22000 20000 23000];
plot([meshnod(bndelenew(:,2),2) meshnod(bndelenew(:,3),2)]',...
  [meshnod(bndelenew(:,2),3) meshnod(bndelenew(:,3),3)]','g')

for i = 1:n
       text(meshnod(bndelenew(i,2),2),meshnod(bndelenew(i,2),3),[int2str(i)])
end
title('bndele numbers     Use zoom to examine')
axis('equal')
grid

figure(4)
clf
plot([meshnod(bndelenew(:,2),2) meshnod(bndelenew(:,3),2)]',...
  [meshnod(bndelenew(:,2),3) meshnod(bndelenew(:,3),3)]','g')

for i = 1:n
        text(meshnod(bndelenew(i,2),2),meshnod(bndelenew(i,2),3),...
         [int2str(bndelenew(i,2))])
end
title('meshnod numbers     Use zoom to examine')
axis('equal')
grid
end    % if plottingon == 1



bndelenew = [bndelenew zeros(n,1) ones(n,1)];
% new one needs to be changed
% Set shoreline to value of  one
% Set open water to value of five
%  Do clever find statements to get indexs of ocean in lopen
l1 = find(meshnod(bndelenew(:,2),2)>80000);
l2 = find(meshnod(bndelenew(:,2),3) <-194000);
lopen = [l1 ; l2];
%disp(' ocean ')
%size(lopen)
bndelenew(lopen,5)=ones(length(lopen),1)*5;    % open sea


lplus = find(bndele(:,1)>0);
 if plottingon == 1 | (length(lplus) > 1)
%figure(5)
l = find(bndelenew(:,5)==5);
nl = find(bndelenew(:,5)~=5);
size(bndelenew);
size(meshnod);
if ~isempty(l)
  plot([meshnod(bndelenew(l,2),2) meshnod(bndelenew(l,3),2)]',...
  [meshnod(bndelenew(l,2),3) meshnod(bndelenew(l,3),3)]','b')
hold on
end
if ~isempty(nl)
  plot([meshnod(bndelenew(nl,2),2) meshnod(bndelenew(nl,3),2)]',...
  [meshnod(bndelenew(nl,2),3) meshnod(bndelenew(nl,3),3)]','g')
hold off
end

% Bad ones from bndele
lplus = find(bndele(:,1)>0);
if length(lplus) > 1
  bndele(lplus,:);
hold on
  plot([meshnod(bndele(lplus,2),2) meshnod(bndele(lplus,3),2)]',...
  [meshnod(bndele(lplus,2),3) meshnod(bndele(lplus,3),3)]','r')
  plot([meshnod(bndele(lplus,2),2) meshnod(bndele(lplus,3),2)]',...
  [meshnod(bndele(lplus,2),3) meshnod(bndele(lplus,3),3)]','r*')
title(' Interior triangles in Red','Color',[1 0 0])
%disp(' interior triangles, probably not meshed correctly:')
%[bndele(lplus,2) bndele(lplus,3)]
disp('Make a micro meshele and feed through rghthand')
disp('meshelemicro = [ 1 1070 981 1071 ; 2 2210 2211 2123 ]')
disp('meshelemicro = rghthand(meshelemicro,meshnod)') 
disp('meshele = [ meshele ; meshelemicro];')

hold off
end

end %  if plottingon == 1 figure 5

%-----------------------------------------------------------------------


%-----------------------------------------------------------------------

function plotpsi(meshele,meshnod,PSI,bcgrad)

tic
plot(meshnod(:,2),meshnod(:,3))
L = axis;
%L=[ -76.3285  -76.2121   37.2962   37.4823]; axis(L)
dc=(max(PSI)-min(PSI))/12;
contours=[min(PSI):dc:max(PSI)];
%contours=[0:5:120]
contmesh(L,meshele,meshnod,PSI,contours);
% hold on; ploteles(meshele,meshnod,2,'k'); hold off
%colorbar vert
%axis([ -76.5352  -76.1891   37.9611   38.6021])
% L=[ -77.5  -74.5   33.5   37.0]; axis(L)
title(['aen2bc '])
zoom on
disp(['plotting took ',num2str(toc),' sec'])
if nargin==4
% plot the gradient vectors

x = (meshnod(bcgrad(:,1),2)+meshnod(bcgrad(:,2),2))/2;
y = (meshnod(bcgrad(:,1),3)+meshnod(bcgrad(:,2),3))/2;
u=real(bcgrad(:,3)); v = imag(bcgrad(:,3));

hold on
%axis([-10 70 -10 70]);
axis equal
%vect(x,y,u,v,max(sqrt(u.*u+v.*v))/10,'k')
hold off

end

%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------

%-----------------------------------------------------------------------

function [A, meshele, meshnod,bathy] = xmgrid2quoddy(filename_xmg)

% [A, meshele, meshnod,bathy] = xmgrid2quoddy(filename_xmg);
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
A.ele=meshele(:,2:4);
A.nod=meshnod(:,2:3);
A.bathy=bathy(:,2);
A.filename=filename_xmg;
%-----------------------------------------------------------------------

function [meshele, meshnod,bathy] = xmgrid2quoddy01(filename_xmg)

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

function [meshele, meshnod, bathy, opnod] = xmgrid2quoddy02(filename_xmg)
% [meshele, meshnod,bathy] = xmgrid2quoddy(filename_xmg);
% Convert the xmgredit5 mesh file into quoddy arrays.
%
fid = fopen(filename_xmg,'rt');
fprintf(1,'%s\n',filename_xmg);
line = fgetl(fid);

[nums] = fscanf(fid,'%i %i',[2]);
numele=nums(1); 
numnode=nums(2);
A = fscanf(fid,'%f',[4,numnode])';

meshnod= A(:,1:3);
bathy=A(:,[1 4]);

A = fscanf(fid,'%f',[5,numele])';

meshele=A(:,[1 3 4 5]);

 tline = fgetl(fid);
 fprintf(1,'%s\n',tline);
 tline = fgetl(fid);
 fprintf(1,'%s\n',tline);

 c=strsplit(tline,' ');
 nb=str2num(c{1});
 tline = fgetl(fid); fprintf(1,'%s\n',tline);
    
 c=strsplit(tline,' ');
   nt=str2num(c{1});
TL=0;opnod=zeros(nt,1);

for n=1:nb;
 tline = fgetl(fid);
 fprintf(1,'%s\n',tline);
 c=strsplit(tline,' ');
   cent01=str2num(c{1});
    for i=1:cent01;
        tline = fgetl(fid);
        c=strsplit(tline,' ');
        opnod(TL+i)=str2num(c{1});
    end;
 TL=TL+cent01;
 fprintf('TL=%d\n',TL);
end;

fclose(fid);
clear A;

%---------------------------------------------------------------------

%-----------------------------------------------------------------------% -----------------------------------------------------

% -----------

function [meshele, meshnod,bathy,mesh2nesh,nesh2mesh]=reorderthesequence(meshele, meshnod,bathy);
for i=1:length(meshnod);
    nnn=meshnod(i,1);
    mesh2nesh(nnn)=i;
    nesh2mesh(i)=nnn;
    meshnod(i,1)=i;
    bathy(i,1)=i;
end;

for i=1:length(meshele);
    meshele(i,1)=i;
    meshele(i,2:4)=mesh2nesh(meshele(i,2:4));
end;


%-----------------------------------------------------------------------
        function B=BCeCovariance(meshnod,bcpsi,L0,B0);
% this is to calculate the error covariance matrix for boundary condition.
% as usual, the the error would be used decrease exponentially with a
% coefficient BCd
        cent01=meshnod(bcpsi(:,3),2);
        cent02=repmat(cent01,1,length(cent01));
        cent03=cent02-cent02';
        
        dent01=meshnod(bcpsi(:,3),3);
        dent02=repmat(dent01,1,length(dent01));
        dent03=dent02-dent02';
        
        eent01=sqrt(cent03.*cent03+dent03.*dent03);
        
        eent01=exp(-eent01/L0);
        eent01(eent01<exp(-4))=0.0;
        B=B0*eent01;
%-----------------------------------------------------------------------
        
        function R=OBeCovariance(meshnod,whichnode,R0);
        R=R0*eye(length(whichnode));
%-----------------------------------------------------------------------
% 
%         function [dxb,x,NPSI]=IO_dd(WL,inpsi,bcpsi,opnod,lnodes,B,R,W0L,CTIP);
%     n1=inpsi(:,3); % node number for the stations
%     H=sparse(lnodes,length(inpsi(:,3)));
%     H1=eye(length(inpsi(:,3)));
%     H(n1,:)=H1;H=H';
%     
%     HW=H*WL;
%     xb=bcpsi(:,2);
%     d=inpsi(:,2); % station data in complex form
%     dd=d-HW*xb;
%     %fprintf(1,'sum(abs(dd0)=%f; %d stations at\n',sum(abs(dd)),length(dd));
%     
%     % Add tidal forcing term
%     if CTIP
%         dd=dd-W0L(inpsi(:,3));
%     end
%     
%     fprintf(1,'sum(abs(dd)=%f; %d stations at\n',sum(abs(dd)),length(dd));
%     BHW=B*HW';
%     RB=R+HW*BHW;
%     dxb=BHW*(RB\dd);
% %    sum(abs(dxb)),
%     x=xb+dxb;
%     NPSI=WL*x;
%     % Add tidal forcing term
%     if CTIP
%         NPSI=NPSI+W0L;
%     end
%     %--------------
%     % Check station amplitude error
%     errA=abs(inpsi(:,2))-abs(NPSI(n1));
%     sum(errA);
%     errA1=abs(errA);
% 
%     errA2=errA1./abs(inpsi(:,2))*100;
% 
% %     cent05=mod(180*(+1)*angle(inpsi(:,2))/pi(),360)-60;
% %     errP=cent05-mod(180*(-1)*angle(NPSI(n1))/pi(),360);
%     fprintf('******************\nTotal error at %d stations in amp: %6.3f m;\n ',...
%         length(errA),sum(errA1));
%     fprintf('avg %6.3fm; %4.1f %% ;rmse=%6.3f m\n',...
%         mean(errA1),mean(errA2),sqrt(mean((errA.^2))));
%     %----------------------------
%     % check error in Pacific, Atlantic and Indian Oceans
%     slon=inpsi(:,5);
%     slat=inpsi(:,6);
%     list_ocean={'Pacific';'Atlantic';'Indian'};
%     col=['ro';'ks';'b*'];
% %     figure(2)
% %     clf
%     for i=1:3
%         switch i
%             case 1 % Pacific Ocean
%                 tocean='Pacific';
%                 loc=find(slon>120 | slon<-100 | (slon>-100 &slon<-70& slat<10) | (slon<-90 & slat<15));
%             case 2 % Atlantic
%                 tocean='Atlantic';
%                 loc=find(slon>-70 & slon<30  );
%             case 3 % Indian Ocean
%                 
%                 loc=find( slon>30  & (slon<120));
%               
%         end
% %         hold on
% %         plot(slon(loc),slat(loc),col(i,:))
%         
%     fprintf('---------------------------\n%s Ocean: %d stations in amp: %6.3f m;\n ',...
%         list_ocean{i},length(errA(loc)),sum(errA1(loc)));
%     fprintf('avg %6.3fm; %4.1f %% ;rmse=%6.3f m\n',...
%         mean(errA1(loc)),mean(errA2(loc)),sqrt(mean((errA(loc).^2))));
%         
%     end
%    % p
%     %     fprintf('Total error at %d stations in phase: %6.3f m; \navg %6.3fm; rmse=%6.3f m\n',...
% %         length(errP),sum(errP),mean(errP),sqrt(mean((errP.^2))));
% 
%     
        %-----------------------------------------------------------------------
    

% % % 
% % %         function [dxb,nxbc]=DA_incr_bc(xbc,dd,HH,HT,A11,nA12,B,R);
% % %             % what is original boundary condition
% % %         
% % %     fprintf(1,'sum(abs(dd)=%f; %d stations at\n',sum(abs(dd)),length(dd));
% % %     
% % %     HWT=nA12'*((A11')\HT);
% % %     HW=HWT';
% % %     
% % %     BHWT=B*HWT;
% % %     RB=R+HW*BHWT;
% % %     dxb=BHWT*(RB\dd);
% % %     nxbc=xbc+dxb;
% % %         function [dxf,nFeq]=DA_diff_incr_In(Feq,dd,nbc,kki,HH,HT,A11,Af11,Af12,B,R);
% % %             % what is original boundary condition
% % %         
% % %     fprintf(1,'sum(abs(dd)=%f; %d stations at\n',sum(abs(dd)),length(dd));
% % %     
% % %     Af1=[Af11 Af12];
% % %     HWT=Af1'*((A11')\HT);
% % %     HW=HWT';
% % %     
% % %     BHWT=B*HWT;
% % %     RB=R+HW*BHWT;
% % %     dxf=BHWT*(RB\dd);
% % %     nFeq=Feq;
% % %     nFeq(kki)=Feq(kki)+dxf(1:length(kki));
% % %     nFeq(nbc)=Feq(nbc)+dxf(length(kki)+1:length(kki)+length(nbc));
    

 
