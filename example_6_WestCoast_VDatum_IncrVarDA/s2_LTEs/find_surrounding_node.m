function Jnsn=find_surrounding_node(meshele,Jwhichnode);
% this is to find out the surrounding node of nodes (Jwhichnode). and it will double
% count the node, if this node in two elements surrounding the node. 
% the output 
% fn_grd   ='/disks/NASUSER/lshi/file_for_Ed/LTEs/unified/test_bc_in/code_2020_04_25/wcadcirc/data/fort.14';
% [meshele, meshnod, bathy, opnod] = xmgrid2quoddy02(fn_grd);% [meshele, meshnod, bathy] = xmgrid2quoddy01(fn_grd); % for broken_square case, opnod == [], no open boundary
% lnodes= length(meshnod);
% leles = length(meshele);

Jnsta=length(Jwhichnode);

Jnsn=cell(Jnsta);
k=1;
    for j=1:Jnsta
        jj=Jwhichnode(j);
        cent01=find(meshele(:,2)==jj);
        cent02=meshele(cent01,[3 4]);cent02=reshape(cent02,[],1);
        cent01=find(meshele(:,3)==jj);
        cent03=meshele(cent01,[2 4]);cent03=reshape(cent03,[],1);
        cent01=find(meshele(:,4)==jj);
        cent04=meshele(cent01,[2 3]);cent04=reshape(cent04,[],1);
        
        Jnsn{j}=[cent02;cent03;cent04];
    end
    