function write_fort_13();
%
fn_grd='/disks/NASUSER/lshi/file_for_Ed/LTEs/unified/test_bc_in/code_2020_04_25/sfbofs/SFB_input/fort.14_land_bc_change_to_0_1';
[meshele, meshnod, bathy] = xmgrid2quoddy01(fn_grd); % for broken_square case, opnod == [], no open boundary
lnodes= length(meshnod);
leles = length(meshele);
% 
% % here we set a uniform FC at beginning
% svbf=0.5*(0.0025+0.00375)*ones(lnodes,1);
% % for i=1:lnodes;
% %     svbf(i)=0.5*(0.0025+0.00375);
% % end;
% 
% load Datasets01.dat;% svbf=Datasets01;
% nodd=[105:119 204:219 314 326:327 1172:1179];
% svbf(nodd)=Datasets01(nodd);
% 
% % specify particular point of high bottom friction
% % svbf([166 164 269 267 266 372 373])=0.02;

load Datasets_e_i_9.dat;svbf=Datasets_e_i_9;
% svbf([166 164 163 269 267 266 372 373 485])=0.02;
% svbf([166 164 163 269 267 266 372 373 485 272 273 274 381 382])=0.02;
svbf([166 164 163 269 267 266 372 373 485 272 273 274 276 381 382 158 260 259])=0.03;


infile='fort.13_e_i_9';
writefort13(lnodes,svbf,infile);

%-----------------------------------------------------------------------
function writefort13(lnodes,svbf,infile);

% lnodes= length(meshnod);
% leles = length(meshele);

fid = fopen(infile,'w+');

fprintf(fid,'%s\n','pnw');
fprintf(fid,'%i\n',lnodes);
fprintf(fid,'%i\n',1);
fprintf(fid,'%s\n','quadratic_friction_coefficient_at_sea_floor');
fprintf(fid,'%i\n',1);
fprintf(fid,'%i\n',1);

fprintf(fid,'%20.16f\n',0.0025);
fprintf(fid,'%s\n','quadratic_friction_coefficient_at_sea_floor');
fprintf(fid,'%i\n',lnodes);
for i=1:lnodes;
   fprintf(fid,'%i\t %12.6f\n',i,svbf(i));
end;
fclose(fid);

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

%-----------------------------------------------------------------------
function writegrd(meshele,meshnod,PSI,infile)
%function writegrd(meshele,meshnod,PSI)
% writes a finite element grid in xmgredit format

%fid = fopen('new.grd','w+');
fid = fopen(infile,'w+');
fprintf(fid,'%s\n','new.grd');

[nnod,junk]=size(meshnod);
for i=1:3
   nnmat(:,i) = meshnod(:,i);
end
nnmat(:,4) = PSI(:,1);
[nele,junk]=size(meshele);
nemat(:,1) = meshele(:,1);
nemat(:,2) = 3;
nemat(:,3) = meshele(:,2);
nemat(:,4) = meshele(:,3);
nemat(:,5) = meshele(:,4);
fprintf(fid,'%i %i\n',[nele nnod]);
for i=1:nnod
   fprintf(fid,'%i\t %12.6f\t %12.6f\t %12.6f\n',[nnmat(i,:)]);
end
for i=1:nele
   fprintf(fid,'%i %i %i %i %i\n', [nemat(i,:)]);
end
fclose(fid);

%-----------------------------------------------------------------------
function writefort24(meshele,meshnod,NTIF,alpha,SALTAMP,SALTPHA,infile);

% lnodes= length(meshnod);
% leles = length(meshele);
[nnod,junk]=size(meshnod);
[nele,junk]=size(meshele);

fid = fopen(infile,'w+');
for k=1:NTIF;
fprintf(fid,'%s\n',alpha{1,k});
fprintf(fid,'%20.16f\n',alpha{2,k});
fprintf(fid,'%i\n',1);
fprintf(fid,'%s\n',alpha{4,k});

for i=1:nnod;
   fprintf(fid,'%i\t %12.6f\t %12.6f\n',i,SALTAMP(i,k),SALTPHA(i,k));
end;
end;
fclose(fid);

%-----------------------------------------------------------------------



