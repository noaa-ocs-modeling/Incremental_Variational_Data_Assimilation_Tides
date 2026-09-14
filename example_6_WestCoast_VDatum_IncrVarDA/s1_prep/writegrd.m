
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

