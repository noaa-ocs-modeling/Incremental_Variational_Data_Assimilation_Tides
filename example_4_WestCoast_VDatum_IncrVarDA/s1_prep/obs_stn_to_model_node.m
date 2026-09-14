% ok here estabilish obs stn to model 1-to-1 connection
% /disks/NASUSER/lshi/file_for_Neil
% ec2012_v3d_chk.grd	 read_input.m	     sea-level-rise.txt
% Mesh.grd		 sea_level_rise.dat
% obs_stn_to_model_node.m  sea_level_rise.txt

function obs_stn_to_model_node();
% here is the deal, all data read in would be in its original data sequence
% without deletion. and the corrosponding name and 
% whichstn select the stn (in the sequence), whichstn=Nstn
% whichHC select the HC in the STNHC,        whcihHC=Nhc
% the dimension of of the STNAMP =[Nstn, Nhc],STNPHA =[Nstn, Nhc],
% STNLOC=[Nstn, 3], column 1 = sequence (not necessary), can be = 0
%                   column 2 = coops station number
%                   column 3 = Longitude
%                   column 4 = Latitude
%                   column 5 = meshgrid node number (=0, if not
%                              corrosponding with mesh directly
% STNLOC=zeros(length(aaa(:,1),5);STNLCO(:,2)=aaa(:,4);STNLCO(:,3)=aaa(:,5);STNLCO(:,4)=aaa(:,6);STNLCO(:,5)=aaa(:,3);
% STNHC = cell[5,Nhc], column 1 = name
%                      column 2 = angular speed (degree/hour)
%                      column 3 = angular speed (/s)
%                      column 4 = node factor, amplitude (non-dimensional)
%                      column 5 = node factor, phase (degree)
% default order for STNHC (CO-OPS)
% co-ops order
% 1	M2	28.984104
% 2	S2	30
% 3	N2	28.43973
% 4	K1	15.041069
% 5	M4	57.96821
% 6	O1	13.943035
% 7	M6	86.95232
% 8	MK3	44.025173
% 9	S4	60
% 10	MN4	57.423832
% 11	NU2	28.512583
% 12	S6	90
% 13	MU2	27.968208
% 14	2N2	27.895355
% 15	OO1	16.139101
% 16	LAM2	29.455626
% 17	S1	15
% 18	M1	14.496694
% {'M2','S2','N2','K1,'M4','O1','M6','MK3','S4','MN4','NU2','S6','MU2','2N2','OO1','LAM2','S1','M1','J1',
%  'MM','SSA','SA','MSF','MF','RHO','Q1','T2','R2','2Q1','P1','2SM2','M3','L2','2MK3','K2','M8','MS4']
% 19	J1	15.5854435
% 20	MM	0.5443747
% 21	SSA	0.0821373
% 22	SA	0.0410686
% 23	MSF	1.0158958
% 24	MF	1.0980331
% 25	RHO	13.471515
% 26	Q1	13.398661
% 27	T2	29.958933
% 28	R2	30.041067
% 29	2Q1	12.854286
% 30	P1	14.958931
% 31	2SM2	31.015896
% 32	M3	43.47616
% 33	L2	29.528479
% 34	2MK3	42.92714
% 35	K2	30.082138
% 36	M8	115.93642
% 37	MS4	58.984104
% set the STNHC 
cent01=load('/disks/NASUSER/lshi/file_for_Ed/LTEs/unified/test_bc_in/code_2020_04_25/sfbofs/SFB_18HC_input/alpha_37HC.mat','alpha');% the sequence is diff with COOPS
% 1     {[1.4052e-04]}    {[0.9990]}    {[118.5000]}    {'M(2)'     }
% 2     {[1.3788e-04]}    {[0.9990]}    {[233.9000]}    {'N(2)'     }
% 3     {[1.4544e-04]}    {[     1]}    {[       0]}    {'S(2)'     }
% 4     {[6.7598e-05]}    {[1.0280]}    {[140.3000]}    {'O(1)'     }
% 5     {[7.2921e-05]}    {[1.0170]}    {[334.3000]}    {'K(1)'     }
% 6     {[1.4584e-04]}    {[1.0220]}    {[128.8000]}    {'K(2)'     }
% 7     {[1.4316e-04]}    {[0.9210]}    {[198.8000]}    {'L(2)'     }
% 8     {[1.3524e-04]}    {[0.9990]}    {[349.3000]}    {'2N(2)'    }
% 9    {[1.4564e-04]}    {[     1]}    {[132.8000]}    {'R(2)'     }
%10     {[1.4525e-04]}    {[     1]}    {[ 47.2000]}    {'T(2)'     }
%11     {[1.4280e-04]}    {[0.9990]}    {[297.5000]}    {'Lambda(2)'}
%12     {[1.3559e-04]}    {[0.9990]}    {[234.9000]}    {'Mu(2)'    }
%13     {[1.3823e-04]}    {[0.9990]}    {[119.5000]}    {'Nu(2)'    }
%14     {[7.5560e-05]}    {[1.0330]}    {[222.8000]}    {'J(1)'     }
%15     {[7.0282e-05]}    {[1.7510]}    {[131.9000]}    {'M(1)'     }
%16     {[7.8245e-05]}    {[1.0870]}    {[356.2000]}    {'OO(1)'    }
%17     {[7.2523e-05]}    {[     1]}    {[ 34.4000]}    {'P(1)'     }
%18     {[6.4959e-05]}    {[1.0280]}    {[255.7000]}    {'Q(1)'     }
%19     {[6.2319e-05]}    {[1.0280]}    {[      11]}    {'2Q(1)'    }
%20     {[6.5312e-05]}    {[1.0280]}    {[141.3000]}    {'Rho(1)'   }
%21     {[2.8104e-04]}    {[0.9990]}    {[     237]}    {'M(4)'     }
%22     {[4.2156e-04]}    {[0.9980]}    {[355.5000]}    {'M(6)'     }
%23     {[5.6208e-04]}    {[0.9970]}    {[     114]}    {'M(8)'     }
%24     {[2.9089e-04]}    {[     1]}    {[       0]}    {'S(4)'     }
%25     {[4.3633e-04]}    {[     1]}    {[       0]}    {'S(6)'     }
%26     {[2.1078e-04]}    {[0.9990]}    {[177.7000]}    {'M(3)'     }
%27     {[7.2722e-05]}    {[     1]}    {[     180]}    {'S(1)'     }
%28     {[2.1344e-04]}    {[1.0170]}    {[ 92.8000]}    {'MK(3)'    }
%29     {[2.0812e-04]}    {[1.0160]}    {[262.7000]}    {'2MK(3)'   }
%30     {[2.7840e-04]}    {[0.9990]}    {[352.4000]}    {'MN(4)'    }
%31     {[2.8596e-04]}    {[0.9990]}    {[118.5000]}    {'MS(4)'    }
%32     {[1.5037e-04]}    {[0.9990]}    {[241.5000]}    {'2SM(2)'   }
%33     {[5.3234e-06]}    {[1.0570]}    {[      18]}    {'Mf'       }
%34     {[4.9252e-06]}    {[0.9990]}    {[241.5000]}    {'Msf'      }
%35     {[2.6392e-06]}    {[0.9960]}    {[244.6000]}    {'Mm'       }
%36     {[1.9911e-07]}    {[     1]}    {[235.6000]}    {'Sa'       }
%37     {[3.9821e-07]}    {[     1]}    {[111.1000]}    {'Ssa'      }
aHC=cent01.alpha;clear cent01;cent01=zeros(4,37);
    for i=1:37;cent01(1,i)=(180/pi())*3600*aHC{1,i};cent01(2,i)=aHC{1,i};cent01(3,i)=aHC{2,i};cent01(4,i)=aHC{3,i};end;
    cent01=cent01(:,[1 3 2 5 21 4 22 28 24 30 13 25 12 8 16 11 27 15 14 35 37 36 34 33 20 18 10 9 19 17 32 26 7 29 6 23 31]);
    Nhc=37;
    alpha=cell(5,Nhc);
    alpha(1,:)={'M2','S2','N2','K1','M4','O1','M6','MK3','S4','MN4','NU2','S6','MU2','2N2','OO1','LAM2','S1','M1','J1',...
                'MM','SSA','SA','MSF','MF','RHO','Q1','T2','R2','2Q1','P1','2SM2','M3','L2','2MK3','K2','M8','MS4'}
    for i=1:Nhc;alpha{2,i}=cent01(1,i);alpha{3,i}=cent01(2,i);alpha{4,i}=cent01(3,i);alpha{5,i}=cent01(4,i);end;
    save('/disks/NASUSER/lshi/file_for_Ed/LTEs/unified/test_bc_in/code_2020_04_25/sfbofs/SFB_18HC_input/alpha_coops_37HC.mat','alpha');
    load('/disks/NASUSER/lshi/file_for_Ed/LTEs/unified/test_bc_in/code_2020_04_25/sfbofs/SFB_18HC_input/alpha_coops_37HC.mat','alpha');
%     {'M2'  }    {[ 28.9841]}    {[1.4052e-04]}    {[0.9990]}    {[118.5000]}
%     {'S2'  }    {[ 30.0000]}    {[1.4544e-04]}    {[     1]}    {[       0]}
%     {'N2'  }    {[ 28.4397]}    {[1.3788e-04]}    {[0.9990]}    {[233.9000]}
%     {'K1'  }    {[ 15.0411]}    {[7.2921e-05]}    {[1.0170]}    {[334.3000]}
%     {'M4'  }    {[ 57.9682]}    {[2.8104e-04]}    {[0.9990]}    {[     237]}
%     {'O1'  }    {[ 13.9430]}    {[6.7598e-05]}    {[1.0280]}    {[140.3000]}
%     {'M6'  }    {[ 86.9523]}    {[4.2156e-04]}    {[0.9980]}    {[355.5000]}
%     {'MK3' }    {[ 44.0252]}    {[2.1344e-04]}    {[1.0170]}    {[ 92.8000]}
%     {'S4'  }    {[ 60.0000]}    {[2.9089e-04]}    {[     1]}    {[       0]}
%     {'MN4' }    {[ 57.4238]}    {[2.7840e-04]}    {[0.9990]}    {[352.4000]}
%     {'NU2' }    {[ 28.5126]}    {[1.3823e-04]}    {[0.9990]}    {[119.5000]}
%     {'S6'  }    {[ 90.0000]}    {[4.3633e-04]}    {[     1]}    {[       0]}
%     {'MU2' }    {[ 27.9682]}    {[1.3559e-04]}    {[0.9990]}    {[234.9000]}
%     {'2N2' }    {[ 27.8954]}    {[1.3524e-04]}    {[0.9990]}    {[349.3000]}
%     {'OO1' }    {[ 16.1391]}    {[7.8245e-05]}    {[1.0870]}    {[356.2000]}
%     {'LAM2'}    {[ 29.4556]}    {[1.4280e-04]}    {[0.9990]}    {[297.5000]}
%     {'S1'  }    {[ 15.0000]}    {[7.2722e-05]}    {[     1]}    {[     180]}
%     {'M1'  }    {[ 14.4967]}    {[7.0282e-05]}    {[1.7510]}    {[131.9000]}
%     {'J1'  }    {[ 15.5854]}    {[7.5560e-05]}    {[1.0330]}    {[222.8000]}
%     {'MM'  }    {[  0.5444]}    {[2.6392e-06]}    {[0.9960]}    {[244.6000]}
%     {'SSA' }    {[  0.0821]}    {[3.9821e-07]}    {[     1]}    {[111.1000]}
%     {'SA'  }    {[  0.0411]}    {[1.9911e-07]}    {[     1]}    {[235.6000]}
%     {'MSF' }    {[  1.0159]}    {[4.9252e-06]}    {[0.9990]}    {[241.5000]}
%     {'MF'  }    {[  1.0980]}    {[5.3234e-06]}    {[1.0570]}    {[      18]}
%     {'RHO' }    {[ 13.4715]}    {[6.5312e-05]}    {[1.0280]}    {[141.3000]}
%     {'Q1'  }    {[ 13.3987]}    {[6.4959e-05]}    {[1.0280]}    {[255.7000]}
%     {'T2'  }    {[ 29.9589]}    {[1.4525e-04]}    {[     1]}    {[ 47.2000]}
%     {'R2'  }    {[ 30.0411]}    {[1.4564e-04]}    {[     1]}    {[132.8000]}
%     {'2Q1' }    {[ 12.8543]}    {[6.2319e-05]}    {[1.0280]}    {[      11]}
%     {'P1'  }    {[ 14.9589]}    {[7.2523e-05]}    {[     1]}    {[ 34.4000]}
%     {'2SM2'}    {[ 31.0159]}    {[1.5037e-04]}    {[0.9990]}    {[241.5000]}
%     {'M3'  }    {[ 43.4762]}    {[2.1078e-04]}    {[0.9990]}    {[177.7000]}
%     {'L2'  }    {[ 29.5285]}    {[1.4316e-04]}    {[0.9210]}    {[198.8000]}
%     {'2MK3'}    {[ 42.9271]}    {[2.0812e-04]}    {[1.0160]}    {[262.7000]}
%     {'K2'  }    {[ 30.0821]}    {[1.4584e-04]}    {[1.0220]}    {[128.8000]}
%     {'M8'  }    {[115.9364]}    {[5.6208e-04]}    {[0.9970]}    {[     114]}
%     {'MS4' }    {[ 58.9841]}    {[2.8596e-04]}    {[0.9990]}    {[118.5000]}


a=1
% 
% read xlsx file
infile='/disks/NASUSER/lshi/NOAA_COOPS_tidal_datum/Copy_of_USA_WestCoast_Complete_List_xGEOID17B_copy.xlsx';% dataset is identical %load('/disks/NASUSER/lshi/file_for_Ed/LTEs/unified/test_bc_in/code_2020_04_25/wcadcirc/data/coops149HC_WestCoast_20180309.mat','xyp');

[num, txt, raw]= xlsread(infile,'Harmonics');
aaaa=0;
cent01=num(:,[1 3 4 5 7 8]);
num=zeros(length(cent01(:,1)),7);num(:,1)=[1:length(cent01(:,1))];num(:,[2:7])=cent01;clear cent01;
allhc=num;

% load /disks/NASUSER/lshi/NOAA_COOPS_tidal_datum/Harmonics_all_correct.txt;
% load /disks/NASUSER/lshi/NOAA_COOPS_tidal_datum/Harmonics_all_16July2013_tidyup.txt;allhc=Harmonics_all_16July2013_tidyup;
lon=allhc(:,3);
lat=allhc(:,4);
stnnam=allhc(:,2);
%%% allhc=allhc(find(lon>-123.1874 & lon < -121.434 & lat > 37.3793 & lat < 38.3199 & stnnam ~= 9414958),:);this is to limit the station only to San Francisco Bay
allhc1L=allhc(1:37:length(allhc(:,1)),:);
% this for database Harmonics_all_16July2013_tidyup.txt; stnseq=[16 1 2 3 4 12 5 11 6 7 9 8 10 14 24 21 19 20 22 13 18 17 23];% stnseq=[16 15 1 2 3 4 12 5 11 6 7 9 8 10 14 24 21 19 20 22 13 18 17 23];
% this is for database Copy_of_USA_WestCoast_Complete_List_xGEOID17B_copy.xlsx
% % % cent01=zeros(37*32,7);iii=0;for i=[2:24 26 27 28 30 32:36];iii=iii+1;cent01(37*(iii-1)+1:37*iii,:)=allhc(37*(i-1)+1:37*i,:);end;
% % % allhc=cent01;
% % % allhc1L=allhc(1:37:length(allhc(:,1)),:);
save allhc1L_new.dat allhc1L -ASCII;
% allhc1L=allhc(1:37:length(allhc(:,1)),[4 1 3 2 6]);

% % %             % load /disks/NASUSER/lshi/NOAA_COOPS_tidal_datum/merged_20130127_tidal_DATUM_Chesapeake_Delaware_Bay.txt;
% % %             % load /disks/NASUSER/lshi/NOAA_COOPS_tidal_datum/merr_20181106_num01.txt;
% % %             % load /disks/NASUSER/lshi/file_for_Neil/sea_level_rise.dat;
% % %             load /disks/NASUSER/lshi/file_for_Neil/estofs_domain_interp/sea_level_rise.dat;
% % %             merged=sea_level_rise;% merged=merged_20130127_tidal_DATUM_Chesapeake_Delaware_Bay;
% % %             ax=merged(:,1);ay=merged(:,2);
ax=allhc1L(:,3);ay=allhc1L(:,4);
% % % bbb=merged(ax > -74.686 & ax < -70.981 & ay > 39.221 & ay < 42.841,:);
% % % size(bbb),size(merged),
% % % merged=bbb;

% % %     load /disks/NASUSER/lshi/VDATUM/datum_stat_interp/observed_error_COOPS/uncer.txt;
% % %     merged(:,2:3)=uncer(:,2:3);
m=length(ax); % m=length(merged(:,1));
% co-ops order
% 1	M2	28.984104
% 2	S2	30
% 3	N2	28.43973
% 4	K1	15.041069
% 5	M4	57.96821
% 6	O1	13.943035
% 7	M6	86.95232
% 8	MK3	44.025173
% 9	S4	60
% 10	MN4	57.423832
% 11	NU2	28.512583
% 12	S6	90
% 13	MU2	27.968208
% 14	2N2	27.895355
% 15	OO1	16.139101
% 16	LAM2	29.455626
% 17	S1	15
% 18	M1	14.496694
% 19	J1	15.5854435
% 20	MM	0.5443747
% 21	SSA	0.0821373
% 22	SA	0.0410686
% 23	MSF	1.0158958
% 24	MF	1.0980331
% 25	RHO	13.471515
% 26	Q1	13.398661
% 27	T2	29.958933
% 28	R2	30.041067
% 29	2Q1	12.854286
% 30	P1	14.958931
% 31	2SM2	31.015896
% 32	M3	43.47616
% 33	L2	29.528479
% 34	2MK3	42.92714
% 35	K2	30.082138
% 36	M8	115.93642
% 37	MS4	58.984104
% old order of STNAMP and STNPHA, corresponding to 37 co-ops order % [4 1 3 6 2 30 26 35 5 7 34 8 37 10 32 31 36 9]; 18 total; in sequence
% 'K(1)'    'M(2)'    'N(2)'    'O(1)'    'S(2)'    'P(1)'    'Q(1)'    'K(2)'    '2MK(3)'    'M(4)'    'MK(3)'    'M(6)'    'MS(4)'    'MN(4)'    'M(3)'    '2SM(2)'    'M(8)'    'S(4)'
% new order of STNAMP and STNPHA, corresponding to 37 co-ops order % [4 6 30 26 1 2 3 35 5 7 34 8 37 10 32 31 36 9]; 18 total; in sequence
% 'K(1)'    'O(1)'    'P(1)'    'Q(1)'    'M(2)'    'S(2)'    'N(2)'    'K(2)'    '2MK(3)'    'M(4)'    'MK(3)'    'M(6)'    'MS(4)'    'MN(4)'    'M(3)'    '2SM(2)'    'M(8)'    'S(4)'

NTIF=37;STNAMP=zeros(m,NTIF);STNPHA=zeros(m,NTIF);
for TPC=1:NTIF ;
NHC=TPC;
STNAMP(:,TPC)=allhc(NHC:37:length(allhc(:,1)),6);
STNPHA(:,TPC)=allhc(NHC:37:length(allhc(:,1)),7);
end;
% TPC=2 ;NHC=6;STNAMP(:,TPC)=allhc(NHC:37:length(allhc(:,1)),6);STNPHA(:,TPC)=allhc(NHC:37:length(allhc(:,1)),7);
% TPC=3 ;NHC=30;STNAMP(:,TPC)=allhc(NHC:37:length(allhc(:,1)),6);STNPHA(:,TPC)=allhc(NHC:37:length(allhc(:,1)),7);
% TPC=4 ;NHC=26;STNAMP(:,TPC)=allhc(NHC:37:length(allhc(:,1)),6);STNPHA(:,TPC)=allhc(NHC:37:length(allhc(:,1)),7);
% TPC=5 ;NHC=1;STNAMP(:,TPC)=allhc(NHC:37:length(allhc(:,1)),6);STNPHA(:,TPC)=allhc(NHC:37:length(allhc(:,1)),7);
% TPC=6 ;NHC=2;STNAMP(:,TPC)=allhc(NHC:37:length(allhc(:,1)),6);STNPHA(:,TPC)=allhc(NHC:37:length(allhc(:,1)),7);
% TPC=7 ;NHC=3;STNAMP(:,TPC)=allhc(NHC:37:length(allhc(:,1)),6);STNPHA(:,TPC)=allhc(NHC:37:length(allhc(:,1)),7);
% TPC=8 ;NHC=35;STNAMP(:,TPC)=allhc(NHC:37:length(allhc(:,1)),6);STNPHA(:,TPC)=allhc(NHC:37:length(allhc(:,1)),7);
% TPC=9 ;NHC=5;STNAMP(:,TPC)=allhc(NHC:37:length(allhc(:,1)),6);STNPHA(:,TPC)=allhc(NHC:37:length(allhc(:,1)),7);
% TPC=10;NHC=7;STNAMP(:,TPC)=allhc(NHC:37:length(allhc(:,1)),6);STNPHA(:,TPC)=allhc(NHC:37:length(allhc(:,1)),7);
% TPC=11;NHC=34;STNAMP(:,TPC)=allhc(NHC:37:length(allhc(:,1)),6);STNPHA(:,TPC)=allhc(NHC:37:length(allhc(:,1)),7);
% TPC=12;NHC=8;STNAMP(:,TPC)=allhc(NHC:37:length(allhc(:,1)),6);STNPHA(:,TPC)=allhc(NHC:37:length(allhc(:,1)),7);
% TPC=13;NHC=37;STNAMP(:,TPC)=allhc(NHC:37:length(allhc(:,1)),6);STNPHA(:,TPC)=allhc(NHC:37:length(allhc(:,1)),7);
% TPC=14;NHC=10;STNAMP(:,TPC)=allhc(NHC:37:length(allhc(:,1)),6);STNPHA(:,TPC)=allhc(NHC:37:length(allhc(:,1)),7);
% TPC=15;NHC=32;STNAMP(:,TPC)=allhc(NHC:37:length(allhc(:,1)),6);STNPHA(:,TPC)=allhc(NHC:37:length(allhc(:,1)),7);
% TPC=16;NHC=31;STNAMP(:,TPC)=allhc(NHC:37:length(allhc(:,1)),6);STNPHA(:,TPC)=allhc(NHC:37:length(allhc(:,1)),7);
% TPC=17;NHC=36;STNAMP(:,TPC)=allhc(NHC:37:length(allhc(:,1)),6);STNPHA(:,TPC)=allhc(NHC:37:length(allhc(:,1)),7);
% TPC=18;NHC=9;STNAMP(:,TPC)=allhc(NHC:37:length(allhc(:,1)),6);STNPHA(:,TPC)=allhc(NHC:37:length(allhc(:,1)),7);
% 1 K1	4 ; [4 1 3 6 2 30 26 35 5 7 34 8 37 10 32 31 36 9]; 18 total; in sequence
% 4 O1	6
% 6 P1	30
% 7 Q1	26
% 2 M2	1
% 5 S2	2
% 3 N2	3
% 8 K2	35
% 9 M4	5
% 10 M6	7
% 11 2MK3	34
% 12 MK3	8
% 13 MS4	37
% 14 MN4	10
% 15 M3	32
% 16 2SM2	31
% 17 M8	36
% 18 S4	9

merged=allhc1L;

% fn_grd   = '/disks/NASUSER/lshi/Chesapeake_Delaware_Bay/fort.14';% chesapeake bay and Delaware bay
% fn_grd   = '/disks/NASUSER/lshi/file_for_Wei/testNY/code/ny_simple.14';% new york
fn_grd   = '/disks/NASUSER/lshi/file_for_Neil/Mesh.grd';% new york
fn_grd   = '/disks/NASUSER/lshi/file_for_Neil/ec2012_v3d_chk.grd';% new york
fn_grd='/disks/NASUSER/lshi/file_for_Ed/LTEs/unified/test_bc_in/code_2020_04_25/sfbofs/SFB_input/fort.14_land_bc_change_to_0_1';
fn_grd   ='/disks/NASUSER/lshi/file_for_Ed/LTEs/unified/test_bc_in/code_2020_04_25/wcadcirc/data/fort.14';

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
    cent01=(meshnod(:,2)-merged(j,3))*cos(merged(j,4)*pi()/180.0)*1000.0*6400.0*2.0*pi()/360.0;% cent01=(meshnod(:,2)-merged(j,3))*cos(38.0*pi()/180.0)*1000.0*6400.0*2.0*pi()/360.0;
    cent02=(meshnod(:,3)-merged(j,4))*1000.0*6400.0*2.0*pi()/360.0;
    cent03=sqrt(cent01.*cent01+cent02.*cent02);
%    cent03=distance(meshnod(:,3),meshnod(:,2),merged(j,4),merged(j,3))*6400*1000;
    cent04=min(cent03);indist(j,1)=cent04;
    ind(j,7)=cent04;
    ind(j,3)=find(cent03==cent04,1,'first');
end;
save indist.txt indist -ASCII;

stnloc=zeros(m,3);
stnloc(:,1)=360.0+ax;stnloc(:,2)=ay;stnloc(:,3)=[1:m]';
save stnloc.dat stnloc -ASCII;

% % % that is all point
% % scatter(ind(:,5),ind(:,6),5,ind(:,2));
% % scatter(meshnod(ind(:,3),2),meshnod(ind(:,3),3),5,ind(:,2));
% % for j=1:m;
% %     plot([ind(j,5),meshnod(ind(j,3),2)],[ind(j,6),meshnod(ind(j,3),3)],'g');
% % end;
% % triplot(meshele(:,2:4),meshnod(:,2),meshnod(:,3));
% % axis([-78 -73 35.5 40.5]);
% % 
% 
% % that is all point
% triplot(meshele(:,2:4),meshnod(:,2),meshnod(:,3));
% % scatter(ind(indist>250.0,5),ind(indist>250.0,6),50,ind(indist>250.0,2));
% % scatter(meshnod(ind(indist>250.0,3),2),meshnod(ind(indist>250.0,3),3),50,ind(indist>250.0,2));
% onandoff=[49 83 81 45 19 62 63 64 78 82 83 90 95 97 99 93 94 100 102 114 115 119 124 127 126 128 132 136 139 140 137 141 147 146 148 149 155 156];
% onandoff=sort(unique(onandoff));
% onaoff=ones(m,1);
% onaoff(onandoff,:)=0;onaoff=logical(onaoff);onaoff=onaoff & indist>250.0;
% 
% scatter(ind(onaoff,5),ind(onaoff,6),100,'r','fill');
% scatter(meshnod(ind(onaoff,3),2),meshnod(ind(onaoff,3),3),100,'g','fill');
% for j=1:m;
%     if onaoff(j);
%     plot([ind(j,5),meshnod(ind(j,3),2)],[ind(j,6),meshnod(ind(j,3),3)],'r','linewidth',2);
% %    text(meshnod(ind(j,3),2),meshnod(ind(j,3),3),num2str(j),'color',[1 0 0],'fontsize',16);
%     text(ind(j,5),ind(j,6),num2str(j),'color',[1 0 0],'fontsize',16);
%     end;
% end;
% axis([-78 -73 35.5 40.5]);set(gca,'dataaspectratio',[1,cos(38.0*pi()/180) 1]);
% % scatter(ind(indist>250.0,5),ind(indist>250.0,6),100,'r','fill');
% % scatter(meshnod(ind(indist>250.0,3),2),meshnod(ind(indist>250.0,3),3),100,'g','fill');
% % for j=1:m;
% %     if indist(j)>250.0;
% %     plot([ind(j,5),meshnod(ind(j,3),2)],[ind(j,6),meshnod(ind(j,3),3)],'r','linewidth',2);
% % %    text(meshnod(ind(j,3),2),meshnod(ind(j,3),3),num2str(j),'color',[1 0 0],'fontsize',16);
% %     text(ind(j,5),ind(j,6),num2str(j),'color',[1 0 0],'fontsize',16);
% %     end;
% % end;
% % axis([-78 -73 35.5 40.5]);set(gca,'dataaspectratio',[1,cos(38.0*pi()/180) 1]);
% % % onandoff=[49 83 81 45 19 62 63 64 78 82 83 90 95 97 99 93 94 100 102 114 115 119 124 127 126 128 132 136 139 140 137 141 147 146 148 149 155 156];
% % % onandoff=sort(unique(onandoff));
% % % onaoff=ones(m,1);
% % % onaoff(onandoff,:)=0;onaoff=logical(onaoff);%onaoff=onaoff & indist>250.0;
% % % aaa=ind(onaoff,:);

cent01=ind(:,7);
aaa=ind(cent01<2000,[1:6]);
STNAMP=STNAMP(cent01<2000,:);STNPHA=STNPHA(cent01<2000,:);
% expand aaa and STNAMP TSNPHA, M2, that is to put control point from the
% non-COOPS data, for example deep sea presure (bpr, bottom pressure recorder) gauge
bbb=[-120.699	32.247	24	102	585627;...
     -127.013	39.331	24	77	695286;...
     -128.900	42.604	24	75	722170;...
     -128.778	45.859	24	74	739372;...
     -130.813	46.773	8	107	716991;...
     -129.617	48.762	24	78	726222];% there is 6 bpr within the domain from south to north
bbb1=0; 
aaa(m+1:m+6,1:6)=zeros(6,6);
    aaa(m+1:m+6,1)=[m+1:m+6]';aaa(m+1:m+6,3)=bbb(:,5);aaa(m+1:m+6,4)=bbb(:,4);aaa(m+1:m+6,5)=bbb(:,1);aaa(m+1:m+6,6)=bbb(:,2);
STNAMP(m+1:m+6,1:NTIF)=nan*zeros(6,NTIF);STNPHA(m+1:m+6,1:NTIF)=nan*zeros(6,NTIF);
    fesAP=load('/disks/NASUSER/lshi/file_for_Ed/LTEs/unified/test_bc_in/code_2020_04_25/wcadcirc/data/bpr_v130121_cm_4_20180928.mat','K1');
i=4;    STNAMP(m+1:m+6,i)=fesAP.K1(bbb(:,4),1)/100;STNPHA(m+1:m+6,i)=fesAP.K1(bbb(:,4),2);
    fesAP=load('/disks/NASUSER/lshi/file_for_Ed/LTEs/unified/test_bc_in/code_2020_04_25/wcadcirc/data/bpr_v130121_cm_4_20180928.mat','O1');
i=6;    STNAMP(m+1:m+6,i)=fesAP.O1(bbb(:,4),1)/100;STNPHA(m+1:m+6,i)=fesAP.O1(bbb(:,4),2);
    fesAP=load('/disks/NASUSER/lshi/file_for_Ed/LTEs/unified/test_bc_in/code_2020_04_25/wcadcirc/data/bpr_v130121_cm_4_20180928.mat','M2');
i=1;    STNAMP(m+1:m+6,i)=fesAP.M2(bbb(:,4),1)/100;STNPHA(m+1:m+6,i)=fesAP.M2(bbb(:,4),2);
    fesAP=load('/disks/NASUSER/lshi/file_for_Ed/LTEs/unified/test_bc_in/code_2020_04_25/wcadcirc/data/bpr_v130121_cm_4_20180928.mat','S2');
i=2;    STNAMP(m+1:m+6,i)=fesAP.S2(bbb(:,4),1)/100;STNPHA(m+1:m+6,i)=fesAP.S2(bbb(:,4),2);
abcdef=0
% % % % expand aaa and STNAMP TSNPHA, K1
% % % % expand aaa and STNAMP TSNPHA, M2, that is to put control point from the
% % % % non-COOPS data, for example fes2014 database (interpolate to two points,
% % % % node number 689 and 1110)
% % % load Datasets_fes2104_pha.dat; fespha=Datasets_fes2104_pha;
% % % load Datasets_fes2104_amp.dat; fesamp=Datasets_fes2104_amp;
% % % aaa(m+1,1:6)=zeros(1,6);aaa(m+2,1:6)=zeros(1,6);
% % %     aaa(m+1:m+2,1)=[m+1 m+2]';aaa(m+1:m+2,3)=[689 1110]';aaa(m+1:m+2,4)=[2014 2014]';aaa(m+1:m+2,5)=meshnod([689 1110],2);aaa(m+1:m+2,6)=meshnod([689 1110],3);
% % % STNAMP(m+1,1:5)=zeros(1,5);STNAMP(m+2,1:5)=zeros(1,5);STNAMP(m+1:m+2,2)=fesamp([689 1110]);
% % % STNPHA(m+1,1:5)=zeros(1,5);STNPHA(m+2,1:5)=zeros(1,5);STNPHA(m+1:m+2,2)=360+fespha([689 1110]);
% % % % expand aaa and STNAMP TSNPHA, K1
% % % load Datasets_fes2104_k1_pha.dat; fespha=Datasets_fes2104_k1_pha;
% % % load Datasets_fes2104_k1_amp.dat; fesamp=Datasets_fes2104_k1_amp;
% % % STNAMP(m+1:m+2,1)=fesamp([689 1110]);
% % % STNPHA(m+1:m+2,1)=360+fespha([689 1110]);
% % % % expand aaa and STNAMP TSNPHA, N2
% % % infile='n2';fesamp=load(['Datasets_fes2104_' infile '_amp.dat']);fespha=load(['Datasets_fes2104_' infile '_pha.dat']);
% % % i=3;STNAMP(m+1:m+2,i)=fesamp([689 1110]);STNPHA(m+1:m+2,i)=360+fespha([689 1110]);
% % % % expand aaa and STNAMP TSNPHA, O1
% % % infile='o1';fesamp=load(['Datasets_fes2104_' infile '_amp.dat']);fespha=load(['Datasets_fes2104_' infile '_pha.dat']);
% % % i=4;STNAMP(m+1:m+2,i)=fesamp([689 1110]);STNPHA(m+1:m+2,i)=360+fespha([689 1110]);
% % % % expand aaa and STNAMP TSNPHA, S2
% % % infile='s2';fesamp=load(['Datasets_fes2104_' infile '_amp.dat']);fespha=load(['Datasets_fes2104_' infile '_pha.dat']);
% % % i=5;STNAMP(m+1:m+2,i)=fesamp([689 1110]);STNPHA(m+1:m+2,i)=360+fespha([689 1110]);
% % % % expand aaa and STNAMP TSNPHA, P1
% % % infile='p1';fesamp=load(['Datasets_fes2104_' infile '_amp.dat']);fespha=load(['Datasets_fes2104_' infile '_pha.dat']);
% % % i=6;STNAMP(m+1:m+2,i)=fesamp([689 1110]);STNPHA(m+1:m+2,i)=360+fespha([689 1110]);
% % % % expand aaa and STNAMP TSNPHA, Q1
% % % infile='q1';fesamp=load(['Datasets_fes2104_' infile '_amp.dat']);fespha=load(['Datasets_fes2104_' infile '_pha.dat']);
% % % i=7;STNAMP(m+1:m+2,i)=fesamp([689 1110]);STNPHA(m+1:m+2,i)=360+fespha([689 1110]);
% % % % expand aaa and STNAMP TSNPHA, K2
% % % infile='k2';fesamp=load(['Datasets_fes2104_' infile '_amp.dat']);fespha=load(['Datasets_fes2104_' infile '_pha.dat']);
% % % i=8;STNAMP(m+1:m+2,i)=fesamp([689 1110]);STNPHA(m+1:m+2,i)=360+fespha([689 1110]);

STNLOC=zeros(length(aaa(:,1)),5);STNLOC(:,1)=aaa(:,1);STNLOC(:,2)=aaa(:,4);STNLOC(:,3)=aaa(:,5);STNLOC(:,4)=aaa(:,6);STNLOC(:,5)=aaa(:,3);
STNCLX=complex(STNAMP.*cos(STNPHA*pi()/180),    -STNAMP.*sin(STNPHA*pi()/180));
save('STNAMPPHA.mat','STNAMP','STNPHA','STNCLX','STNLOC');
        cent01=STNAMP(:,5).*cos(STNPHA(:,5)*pi()/180);
        cent02=STNAMP(:,5).*sin(STNPHA(:,5)*pi()/180);
        cent03=complex(cent01,cent02);
        %%% stnseq1=[25 24 15 1 2 3 4 12 5 11 6 7 9 8 10 14 23 20 18 19 21 13 17 16 22];
%         stnseq1=[34 33 20 1 2 3 4 13 5 11 6 7 9 8 10 19 31 27 25 26 29 17 24 22 30];
%         stnseq2=[34 33 20 1 2 3 4 13 5 11 6 7 9 8 10 19 31 27 25 26 29 17 24 22 30 16 18 15 14 12 21 28 32 23];
%         stnseq2=[16 18 15 14 12 21 28 32 23];
%         stnseq3=[34 33 20 1 13 8 19 23 26 22 30];
        stnseq1=[1:149];
        stnseq2=[150:155];
        stnseq3=[150];

%        plot(cent03(stnseq1),'r');hold on;
       plot(cent03(stnseq1),'r*');hold on;grid on;
%       plot(cent03(stnseq2),'b');
       plot(cent03(stnseq2),'b^')
       %plot(cent03(stnseq2),'b');
       plot(cent03(stnseq3),'g^')
%        plot(cent03(1,1),'b*');hold on;
        grid on;% axis([-1 1 -0 1]);
        axis equal;

% % % % save CH_datum.err aaa -ASCII;
% % % fid = fopen('sea_level_rise_OB.err','w+'); % fid = fopen('sea_level_rise_OB.err','w+');
% % % for j=1:length(aaa);
% % %     fprintf(fid,'%i %12.6f %i %i %12.6f %12.6f\n',[aaa(j,:)]);
% % % % % %     fprintf(fid,'%i %12.6f %i %i %12.6f %12.6f %12.2f\n',[aaa(j,:)]);
% % % end;
% % % fclose(fid);



% load /disks/NASUSER/lshi/file_for_Ed/LTEs/unified/test_bc_in/code_adcirc/results_02/STNHC_zero_nodal.txt;STNHC=STNHC_zero_nodal;
% read in observation data, the sequence is K1 M2 N2 O1 S2

cent01=load('/disks/NASUSER/lshi/file_for_Ed/LTEs/unified/test_bc_in/code_2020_04_25/sfbofs/SFB_18HC_input/alpha_coops_37HC.mat','alpha');% the sequence is diff with COOPS
    aHC=cent01.alpha;clear cent01;
% 1     {'M2'  }    {[ 28.9841]}    {[1.4052e-04]}    {[0.9990]}    {[118.5000]}
% 2     {'S2'  }    {[ 30.0000]}    {[1.4544e-04]}    {[     1]}    {[       0]}
% 3     {'N2'  }    {[ 28.4397]}    {[1.3788e-04]}    {[0.9990]}    {[233.9000]}
% 4     {'K1'  }    {[ 15.0411]}    {[7.2921e-05]}    {[1.0170]}    {[334.3000]}
% 5     {'M4'  }    {[ 57.9682]}    {[2.8104e-04]}    {[0.9990]}    {[     237]}
% 6     {'O1'  }    {[ 13.9430]}    {[6.7598e-05]}    {[1.0280]}    {[140.3000]}
% 7     {'M6'  }    {[ 86.9523]}    {[4.2156e-04]}    {[0.9980]}    {[355.5000]}
% 8     {'MK3' }    {[ 44.0252]}    {[2.1344e-04]}    {[1.0170]}    {[ 92.8000]}
% 9     {'S4'  }    {[ 60.0000]}    {[2.9089e-04]}    {[     1]}    {[       0]}
% 10    {'MN4' }    {[ 57.4238]}    {[2.7840e-04]}    {[0.9990]}    {[352.4000]}
% 11    {'NU2' }    {[ 28.5126]}    {[1.3823e-04]}    {[0.9990]}    {[119.5000]}
% 12    {'S6'  }    {[ 90.0000]}    {[4.3633e-04]}    {[     1]}    {[       0]}
% 13    {'MU2' }    {[ 27.9682]}    {[1.3559e-04]}    {[0.9990]}    {[234.9000]}
% 14    {'2N2' }    {[ 27.8954]}    {[1.3524e-04]}    {[0.9990]}    {[349.3000]}
% 15    {'OO1' }    {[ 16.1391]}    {[7.8245e-05]}    {[1.0870]}    {[356.2000]}
% 16    {'LAM2'}    {[ 29.4556]}    {[1.4280e-04]}    {[0.9990]}    {[297.5000]}
% 17    {'S1'  }    {[ 15.0000]}    {[7.2722e-05]}    {[     1]}    {[     180]}
% 18    {'M1'  }    {[ 14.4967]}    {[7.0282e-05]}    {[1.7510]}    {[131.9000]}
% 19    {'J1'  }    {[ 15.5854]}    {[7.5560e-05]}    {[1.0330]}    {[222.8000]}
% 20    {'MM'  }    {[  0.5444]}    {[2.6392e-06]}    {[0.9960]}    {[244.6000]}
% 21    {'SSA' }    {[  0.0821]}    {[3.9821e-07]}    {[     1]}    {[111.1000]}
% 22    {'SA'  }    {[  0.0411]}    {[1.9911e-07]}    {[     1]}    {[235.6000]}
% 23    {'MSF' }    {[  1.0159]}    {[4.9252e-06]}    {[0.9990]}    {[241.5000]}
% 24    {'MF'  }    {[  1.0980]}    {[5.3234e-06]}    {[1.0570]}    {[      18]}
% 25    {'RHO' }    {[ 13.4715]}    {[6.5312e-05]}    {[1.0280]}    {[141.3000]}
% 26    {'Q1'  }    {[ 13.3987]}    {[6.4959e-05]}    {[1.0280]}    {[255.7000]}
% 27    {'T2'  }    {[ 29.9589]}    {[1.4525e-04]}    {[     1]}    {[ 47.2000]}
% 28    {'R2'  }    {[ 30.0411]}    {[1.4564e-04]}    {[     1]}    {[132.8000]}
% 29    {'2Q1' }    {[ 12.8543]}    {[6.2319e-05]}    {[1.0280]}    {[      11]}
% 30    {'P1'  }    {[ 14.9589]}    {[7.2523e-05]}    {[     1]}    {[ 34.4000]}
% 31    {'2SM2'}    {[ 31.0159]}    {[1.5037e-04]}    {[0.9990]}    {[241.5000]}
% 32    {'M3'  }    {[ 43.4762]}    {[2.1078e-04]}    {[0.9990]}    {[177.7000]}
% 33    {'L2'  }    {[ 29.5285]}    {[1.4316e-04]}    {[0.9210]}    {[198.8000]}
% 34    {'2MK3'}    {[ 42.9271]}    {[2.0812e-04]}    {[1.0160]}    {[262.7000]}
% 35    {'K2'  }    {[ 30.0821]}    {[1.4584e-04]}    {[1.0220]}    {[128.8000]}
% 36    {'M8'  }    {[115.9364]}    {[5.6208e-04]}    {[0.9970]}    {[     114]}
% 37    {'MS4' }    {[ 58.9841]}    {[2.8596e-04]}    {[0.9990]}    {[118.5000]}
whichHC=[4 6 30 26 1 2 3 35];%[5 4 17 18 1 3 2 6];%[5 1 2 4 3 17 18 6];% 18 total old whichHC=[5 1 2 4 3 17 18 6 29 21 28 22 31 30 26 32 23 24];
whichstn=[1:155];
% whichHC represent the order 'K(1)'    'O(1)'    'P(1)'    'Q(1)'    'M(2)'    'S(2)'    'N(2)'    'K(2)' 
% if whichHC=[5 4 17 18 1 3 2 6 29 21 28 22 31 30 26 32 23 24];%[5 1 2 4 3 17 18 6];% 18 total old whichHC=[5 1 2 4 3 17 18 6 29 21 28 22 31 30 26 32 23 24];
% whichHC represent the order 'K(1)'    'O(1)'    'P(1)'    'Q(1)'    'M(2)'    'S(2)'    'N(2)'    'K(2)'    '2MK(3)'    'M(4)'    'MK(3)'    'M(6)'    'MS(4)'    'MN(4)'    'M(3)'    '2SM(2)'    'M(8)'    'S(4)'
% do we need input .obs file?


% 
% 
% NTIF=8;
% alpha=cell(4,NTIF);
% alpha{1,1}='K1 SAEL';alpha{1,2}='M2 SAEL';alpha{1,3}='N2 SAEL';alpha{1,4}='O1 SAEL';alpha{1,5}='S2 SAEL';
% alpha{2,1}=0.000072921158358;alpha{2,2}=0.000140518902509;alpha{2,3}=0.000137879699487;alpha{2,4}=0.000067597744151;alpha{2,5}=0.000145444104333;
% alpha{4,1}='K1';alpha{4,2}='M2';alpha{4,3}='N2';alpha{4,4}='O1';alpha{4,5}='S2';
% mc=0;
% for k=[1 2 3 4 5];mc=mc+1;
%     aaa(:,2)=STNAMP(:,k);
%     fn_err_ha=['fort_14' '_B_' alpha{4,k} '_ha.obs'];
%         fid = fopen(fn_err_ha,'w+');
%         fprintf(fid,'%s\n',[alpha{4,k} ' ha']);
%         fprintf(fid,'%i\n',length(STNAMP(:,1)));
%         for i=1:length(STNAMP(:,1));
%            fprintf(fid,'%i\t %8.3f\t %i\t %i\t %12.6f\t %12.6f\n',[aaa(i,:)]);
% %            fprintf(fid,'%i\t %8.3f\t %i\t %i\t %12.6f\t %12.6f\n',i,STNHC(i,6*(mc-1)+3),b2a(inLL(i,3),2),0,inLL(i,5),inLL(i,6));
%         end;
%         fclose(fid);
%     aaa(:,2)=STNPHA(:,k);
%     fn_err_hp=['fort_14' '_B_' alpha{4,k} '_hp.obs'];
%         fid = fopen(fn_err_hp,'w+');
%         fprintf(fid,'%s\n',[alpha{4,k} ' hp']);
%         fprintf(fid,'%i\n',length(STNAMP(:,1)));
%         for i=1:length(STNAMP(:,1));
%            fprintf(fid,'%i\t %8.3f\t %i\t %i\t %12.6f\t %12.6f\n',[aaa(i,:)]);
% %            fprintf(fid,'%i\t %8.3f\t %i\t %i\t %12.6f\t %12.6f\n',i,STNHC(i,6*(mc-1)+4),b2a(inLL(i,3),2),0,inLL(i,5),inLL(i,6));
%         end;
%         fclose(fid);
% end;

% % % 
% % % fn_grd   = '/disks/NASUSER/lshi/file_for_Wei/testNY/code/ny_mh.14';% new york
% % % [meshele, meshnod,bathy] = xmgrid2quoddy01(fn_grd);
% % % 
% % % aaa01=aaa;aaa01(:,2)=aaa(:,2)-bathy(aaa(:,3),2);
% % % fid = fopen('CH_datum_MHW_01.err','w+');
% % % for j=1:length(aaa);
% % %     fprintf(fid,'%i %12.6f %i %i %12.6f %12.6f\n',[aaa01(j,:)]);
% % % % % %     fprintf(fid,'%i %12.6f %i %i %12.6f %12.6f %12.2f\n',[aaa(j,:)]);
% % % end;
% % % fclose(fid);
% % % 
% % % aaa01=aaa;aaa01(:,2)=bathy(aaa(:,3),2);
% % % fid = fopen('CH_datum_MHW_MD.err','w+');
% % % for j=1:length(aaa);
% % %     fprintf(fid,'%i %12.6f %i %i %12.6f %12.6f\n',[aaa01(j,:)]);
% % % % % %     fprintf(fid,'%i %12.6f %i %i %12.6f %12.6f %12.2f\n',[aaa(j,:)]);
% % % end;
% % % fclose(fid);






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
