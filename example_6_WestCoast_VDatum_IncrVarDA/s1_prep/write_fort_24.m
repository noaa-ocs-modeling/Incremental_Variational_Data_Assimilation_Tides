function write_fort_24();

fn_grd='/disks/NASUSER/lshi/file_for_Ed/LTEs/unified/test_bc_in/code_2020_04_25/sfbofs/SFB_input/fort.14_land_bc_change_to_0_1';
[meshele, meshnod, bathy] = xmgrid2quoddy01(fn_grd); % for broken_square case, opnod == [], no open boundary
lnodes= length(meshnod);
leles = length(meshele);

NTIF=5;
alpha=cell(4,NTIF);
alpha{1,1}='K1 SAEL';alpha{1,2}='M2 SAEL';alpha{1,3}='N2 SAEL';alpha{1,4}='O1 SAEL';alpha{1,5}='S2 SAEL';
alpha{2,1}=0.000072921158358;
alpha{2,2}=0.000140518902509;
alpha{2,3}=0.000137879699487;
alpha{2,4}=0.000067597744151;
alpha{2,5}=0.000145444104333;
alpha{4,1}='K1';alpha{4,2}='M2';alpha{4,3}='N2';alpha{4,4}='O1';alpha{4,5}='S2';

SALTAMP=zeros(lnodes,NTIF);SALTPHA=zeros(lnodes,NTIF);

infile='fort.24_start_zero';

writefort24(meshele,meshnod,NTIF,alpha,SALTAMP,SALTPHA,infile);

% a2b=zeros(lnodes,2);a2b(:,1)=[1:lnodes]';
% for i=1:lnodes;
%     cx=meshnod(i,2);cy=meshnod(i,3);
%     cent01=meshnod1(:,2)-cx;cent02=meshnod1(:,3)-cy;
%     dist=sqrt(cent01.*cent01+cent02.*cent02);
%     [M,I]=min(dist);
%     a2b(i,2)=I;
% end;
% bath_diff=bathy(:,2)-bathy1(a2b(:,2),2);
% 
% b2a=a2b;
% [M,I]=sort(a2b(:,2));
% b2a(:,2)=I;
% 
% save a2b.dat a2b -ASCII;
% save b2a.dat b2a -ASCII;
% save a2b2a.mat a2b b2a;
% 
% ikout14='bath_diff.14';
%          writegrd(meshele,meshnod,bath_diff,ikout14);
%              figure(6);clf;plotpsis(meshele,meshnod,bath_diff);colorbar;
%              figure(7);clf;plotpsis(meshele,meshnod,bathy(:,2));colorbar;
%              figure(8);clf;plotpsis(meshele,meshnod,bathy1(a2b(:,2),2));colorbar;
% print -dtiff bath_diff.tiff;

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



