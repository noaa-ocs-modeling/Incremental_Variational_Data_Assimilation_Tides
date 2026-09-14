function all2=f_split_global(all,z,lon0)
if ~isempty(z)
    all.nodes(:,end)=z;
end
        meshele=all.elements;
        x1=all.nodes(meshele(:,3),2);
        x2=all.nodes(meshele(:,4),2);
        x3=all.nodes(meshele(:,5),2);
        y1=all.nodes(meshele(:,3),3);
        y2=all.nodes(meshele(:,4),3);
        y3=all.nodes(meshele(:,5),3);
        yy=[y1 y2 y3];
        xx=[x1 x2 x3 ];
        xd=max(xx,[],2)-min(xx,[],2);
        %x2-x3
        %dxx=-Rearth*[ (x2-x3).*cos(0.5*(y2+y3))   (x3-x1).*cos(0.5*(y3+y1))  (x1-x2).*cos(0.5*(y1+y2))];
        %dxx=[x2-x3 x3-x1 x1-x2];
        
        %loc=find(xd<300 );
        %-------------------------
        %loc2=find( y1>80 & xd>100 & xd<180 & ~(x1>0&x2>0& x3>0) & ~(x1<0 & x2<0 & x3<0));
        loc=find(  xd>80  & ~(x1>0&x2>0& x3>0) & ~(x1<0 & x2<0 & x3<0));
        %loc=find(  xd>180  & ~(x1>0&x2>0& x3>0) & ~(x1<0 & x2<0 & x3<0));
        %loc2=find( abs(x1)>(170)  & ~(x1>0&x2>0& x3>0) & ~(x1<0 & x2<0 & x3<0));

        nloc=length(loc);
        nn=all.line2(2);
        %xx(loc,:)
        fprintf('%d elements found along the boudary\n', nloc);
        newxyz=[];
        newelement=[];
        newnode=[];
        k=0;
        for i=1:nloc
            iloc=loc(i);
            ixx=xx(iloc,:);
            ielement=meshele(iloc,:);
            loc2=find(ixx<0);
            newelement(end+1,:)=all.elements(iloc,:);
            for j=1:length(loc2)
                k=k+1;
                inode=ielement(loc2(j)+2);
                newxyz(end+1,1:3)=[all.nodes(inode,2)+360 all.nodes(inode,3) all.nodes(inode,4)];
                newnode(end+1)=k+nn;
                newelement(end,2+loc2(j))=newnode(end);
            end
            %------------------
            loc2=find(ixx>0);
            %newelement(end+1,:)=all.elements(iloc,:);
            for j=1:length(loc2)
                k=k+1;
                inode=ielement(loc2(j)+2);
                newxyz(end+1,1:3)=[all.nodes(inode,2)-360 all.nodes(inode,3) all.nodes(inode,4)];
                newnode(end+1)=k+nn;
                %newelement(end,2+loc2(j))=newnode(end);
                all.elements(iloc,2+loc2(j))=newnode(end);
            end
        end
        
        newx=[all.nodes(:,2);newxyz(:,1)];
        newy=[all.nodes(:,3);newxyz(:,2)];
        newz=[all.nodes(:,4);newxyz(:,3)];
        
        %loc=find(xd<360 );
        allnewelements=all.elements;
        %allnewelements(loc,:)=[];
        allnewelements=[allnewelements;newelement];
 %       allnewelements=[all.elements;newelement];
%         trimesh(all.elements(loc,3:5),all.nodes(:,2),all.nodes(:,3),all.nodes(:,4)*0);
%         hold on
%         trimesh(newelement(:,3:5),newx,newy,newz*0+10);
%         hold on
%         trimesh(allnewelements(:,3:5),newx,newy,newz*0+20);
%         hold on
%         trimesh(allnewelements(:,3:5),newx+360,newy,newz*0+20);
        nn=length(newx);
        newx2=[newx ;newx+360];
        newy2=[newy ;newy];
        newz2=[newz ;newz];
        newxyz2=[newx2 newy2 newz2];
        
        newelements2=allnewelements;
        newelements2(:,3:end)=newelements2(:,3:end)+nn;
        newe2=[allnewelements;newelements2];
        %-------------------
        [C,IA,IC]=unique(newxyz2,'rows','stable');
        nC=length(C);
        
        fprintf('%d nodes unique; %d nodes removed \n',nC,length(IC)-nC);
        all2.nodes=(1:nC)';
        all2.nodes(:,2:4)=C;
        
        newe2=IC(newe2(:,3:5));
        [newe2]=unique(newe2,'rows','stable');
        if ~isempty(lon0);
            x=all2.nodes(newe2,2);
            x=reshape(x,[],3);
            xmin=min(x,[],2);
            xmax=max(x,[],2);
            loc=find(xmax<(lon0-5) | xmin>(360+lon0+5));
            newe2(loc,:)=[];
        end
            
        all2.elements=(1:length(newe2))';
        all2.elements(:,2)=all2.elements(:,1)*0+3;
        all2.elements(:,3:5)=newe2;
        all2.line2=[length(newe2) nC];
        all2.line1='Splitted grid';
        
%         all2.x=newx2;
%         all2.y=newy2;
%         all2.z=newz2;
%         all2.elements=newe2(:,3:end);
        
