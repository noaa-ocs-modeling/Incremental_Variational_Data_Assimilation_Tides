clear,
reloadok=1;
addcok=1;
plotgageok=1;
testcase='SFWC300s3'

splitok=0;
if reloadok==1
    pathin=['input_data_' testcase '/'];
    temp=dir([pathin '*14*'])
    gridfile=temp.name;
    all=f_read_bathy_adcirc_grd(pathin,gridfile,0);
    temp=dir([pathin '*inpsiLLfinal*']);
    dartfile= temp.name;
    dartdata=load([pathin dartfile]);
end
%---------------------------------
% plotting 

TPC_list={'O1';'M2';'S2';'K1'};
CSlist={'Cartesian';'ADCIRC';'Spherical';'Spherical2';'Local';'Local2'};
Aelist=[1 2];
data=[];
%load cm
load cm_jet_blue
colormap(cm)
tlist={'Amplitude';'Phase'};
%cm=flipud(cm);
dx=0.04;
xy_lim=[min(all.nodes(:,2))-dx max(all.nodes(:,2))+dx  min(all.nodes(:,3))-dx max(all.nodes(:,3))+dx];
plotcase=1;
switch plotcase
    case 1
        xy_lim=[ -134-1  -115    30-1    52+1];
        ampc_lim=[0.2 1];
        pcv=[30:10:330];
        dxtick=5;
        ax_po=[0.08 0.08 0.7 0.85];
        gcf_po=[10 100 [500 500]/1000*1600];
                uvc_lim=[0 0.5];

    case 2
        xy_lim=[ -123.15  -121.2    37.3    38.9];
        ampc_lim=[0 1];
        pcv=[30:10:330];
        dxtick=0.5;
        ax_po=[0.08 0.08 0.8 0.85];
        gcf_po=[10 100 [800 500]/1000*1600];
        uvc_lim=[0 1];


end
tcok=1;
ta=12;tt=14;
bkc=[0 0 0]+1;
ttc=[0 0 0];
x=[all.nodes(:,2)-360;all.nodes(:,2); all.nodes(:,2)+360];
y=[all.nodes(:,3) ;all.nodes(:,3) ; all.nodes(:,3)];
inv1='inv1';
corrok=1;
for ii=2
    ITPC=TPC_list{ii};
    
for i=6
    for j=2
        %ifile=['endresults_bohai_1it_' ITPC '_4p375e-05_' CSlist{i} '_Ae' int2str(j) '.mat'];
        %temp=dir(['output/endresults_' testcase '*' ITPC '*' CSlist{i} '_Ae' int2str(j) '_*bs*_T2*' inv1 '*C*.mat']);
        %temp=dir(['output/endresults_' testcase '*' ITPC '*' CSlist{i} '_Ae' int2str(j) '_*bs*_T2*' inv1 '*wgs*C*.mat']);
        %temp=dir(['output/endresults_' testcase '*' ITPC '*' CSlist{i} '_Ae' int2str(j) '_*bs*_T2*' inv1 '*.mat']);
        
        iname=['output/endresults_' testcase '*' ITPC '*ITE1*' CSlist{i} '_Ae' int2str(j) '_*bs*_T2_*' inv1 '*'];
        iname=['output/endresults_' testcase '*1000*' ITPC '*ITE1*' CSlist{i} '_Ae' int2str(j) '_*bs*_T2_*' inv1 '*CF7*'];
        iname=['out_' testcase '/endresults_' testcase '*' ITPC '*ITE1*' CSlist{i} '_Ae' int2str(j) '_*bs*_T2_*' inv1 '*CF8*'];
        iname=['out_' testcase '/endresults_' testcase '*1000*' ITPC '*ITE1*' CSlist{i} '_Ae' int2str(j) '_*bs*_T2_*' inv1 '*CF7*'];

        temp=dir(iname);
        ifile=temp.name;
        %eval(['load output/' ifile]);
        eval(['load out_'  testcase '/' ifile]);
        fprintf(1,'load %s\n',ifile)
        if corrok
            NPSI=corrPSI;
            addt='corrNPSI';
        else
            addt='NPSI';
        end
        
        %----------------------
        for k=1:2%:3
            if k~=2
                figure(k)
                set(gcf,'position',gcf_po)
                %set(gcf,'paperpositionmode','auto','InvertHardcopy','off')
                set(gcf,'paperpositionmode','auto','Scrollable','off')
                clf
            end
            hax=axes('position',ax_po);
            set(hax,'color',bkc)
            switch k
                case 1
                    z=abs(NPSI);
                    %all2=f_split_global(all,z);
                    %z2=[z1; z1;z1];
                    %z=griddata(x,y,z2,all0.nodes(:,2),all0.nodes(:,3));
                    cv=[ 0:0.05:1.6];
                     c_lim=[0 1.5];
                     c_lim=ampc_lim;
                switch ITPC
                    case 'N2'
                        cv=[ 0:0.02:1.6];
                        c_lim=[0 .5];
                    case 'O1'
                        c_lim=[0 .5];
                    case 'M2'
                        cv=[ 0.1:0.1:1.6];
                    case 'S2'
                        cv=[ 0:0.04:1.6];
                        c_lim=[0 .5];
                    case 'K1'
                        c_lim=[0 .5];
                end
               
                %c_lim=[0 0.7];
                sf=sprintf('%.2f',max(z));
               
               
                case 2
                z=mod(180*(-1)*angle(NPSI)/pi(),360);

                c_lim=[0 360];
                cv=[30:10:330];
                cv=pcv;
                case 3
                    z=sqrt(abs(PSIu).^2+abs(PSIv).^2);
                    c_lim=uvc_lim;
                %cv=[60:30:330];
        end
        %all2=f_split_global(all,z,lon0);
        all2=all;
        all2.nodes(:,end)=z;
        %colormap(jet)
        colormap(cm)
        %colormap(flipud(jet))
        if k~=2
        trisurf(all2.elements(:,3:5),all2.nodes(:,2),all2.nodes(:,3),all2.nodes(:,end));
        %trisurf(all2.elements,all2.x,all2.y,all2.z);
        
        hold on
        %trisurf(all.elements(loc,3:5),all.nodes(:,2)+360,all.nodes(:,3),z);
        view([0 90])
        hc=colorbar;
        kk=1;
        set(hc,'position',[ax_po(kk,1)+ax_po(kk,3)+0.06, ax_po(kk,2),0.03,ax_po(kk,4)*.8],'fontsize',tt,'AxisLocation','in')
        caxis(c_lim)
        %axis tight
        %axis equal
        shading interp
        set(gca,'fontsize',ta,'color',bkc)
        %load C:\Users\Liujuan.Tang\Documents\VDatum\WestCoast\tide_stations\coops149HC_WestCoast_20180309 xyp
        eval(['load ' pathin 'coops149HC_WestCoast_20180309 xyp'])
        if k==1 & plotgageok==1
            %if plotcase==1
                mt='o';
            %else
            %    mt
            %---
            %-----Plot DART location
            hold on
            %plot3(dartdata(:,5), dartdata(:,6),dartdata(:,6)*0+10,'k+','markersize',8,'linewidth',1 );
            if k==1 & plotcase ==2
            %-----Plott 259 stations with datum data
            load C:\Users\Liujuan.Tang\Documents\VDatum\WestCoast\tide_stations\obs259_wc_20180307 
            hold on
            plot3(obs.x,obs.y,obs.x*10,mt,'color',[0 0 0]+.5)
            for ista=1:length(obs.x);
                text(obs.x(ista),obs.y(ista),int2str(obs.id(ista)));
            end
            %--------plotting 149 stations have HC data
            
            %load C:\Users\Liujuan.Tang\Documents\VDatum\WestCoast\tide_stations\coops149HC_WestCoast_20180309 xyp
            eval(['load ' pathin 'coops149HC_WestCoast_20180309 xyp'])
            hold on
            plot3(xyp(:,2),xyp(:,1),xyp(:,1)*10,mt,'color','r')
            for ista=1:length(xyp);
                %text(xyp(ista,2),xyp(ista,1),int2str(xyp(ista,3)));
            end
            end
            
            %-------Plot 38 open coast tide stations
             eval(['load ' pathin '/' testcase '_snode locmin '])
            plot3(xyp(locmin,2),xyp(locmin,1),xyp(locmin,1)*10,mt,'color','k')
            
        end
        %title([strrep(testcase, '_','-') ' ' ITPC '  '  tlist{k} ' ' tlist{k+1} '  ' strrep(ifile(16:end-4),'_','-') ' ' addt],'fontsize',12)
        title([strrep(testcase, '_','-') ' ' ITPC '  '    strrep(ifile(28:end-4),'_','-') ' ' addt],'fontsize',12)
        box on
        grid off
        axis(xy_lim)
        axis on
        end
        %-----------------
        plotdartok=1;
        if k==1 & plotdartok==1
            axis off
            hax=axes('position',ax_po);
            set(hax,'color','none')
            set(gca,'color','none')
            %text(118, 40.5,['max=' sf ' m'])
            outfile=[pathin '/' testcase '_inpsiLLfinal_' ITPC '.dat' ];
            data=load(outfile);
            xyp=data(:,[5 6]);
            loc=find(xyp(:,1)>180);
            xyp(loc,1)=xyp(loc,1)-360;
            hold on
            %plot3(xyp(:,1),xyp(:,2),data(:,2),'ro','markersize',5)
            scatter(xyp(:,1),xyp(:,2),8,data(:,2)/100,'markerfacecolor','flat','MarkerEdgeColor','black')
            loc=find(xyp(:,1)<0);
            xyp(loc,1)=xyp(loc,1)+360;


            %plot3(xyp(:,1),xyp(:,2),data(:,2),'ro','markersize',5)
            %scatter(xyp(:,1),xyp(:,2),30,data(:,2)/100,'markerfacecolor','flat','MarkerEdgeColor','black')
            caxis(c_lim)
            axis(xy_lim)

        end
        %pp
        %-----------------
        if tcok==1 & k>1
            axis off
            hax=axes('position',ax_po);
            set(hax,'color','none')
            set(gca,'color','none')
             nodes=all2.nodes(:,2:3);
             elements=[all2.elements(:,3:5)];
            hold on
            [cout,hout] = tricontour(nodes,elements,all2.nodes(:,end),cv);
            hc=colorbar;
            kk=1;
            set(hc,'position',[ax_po(kk,1)+ax_po(kk,3)+0.06, ax_po(kk,2),0.03,ax_po(kk,4)*.8],'fontsize',tt,'ytick',0:30:360,'xcolor',ttc,'ycolor',ttc)
            if k==1
                title(hc,'(m)     (deg)','fontsize',tt)
            elseif k==3
                title(hc,'(m/s)    ','fontsize',tt)
            else
                title(hc,'(m)     (deg)','fontsize',tt)

            end
            caxis(c_lim)
            axis(xy_lim)
        end
        figname=['f' testcase '_' ifile(18:end-4) '_' int2str(k) '_' addt '.png'];
        grid off
        %-------------------------------------
        xtick=-180:dxtick:200;
        %xtt=[];
        for ix=1:length(xtick)
            if (xtick(ix)==0 | xtick(ix)==180)
                xtt(ix)={[num2str(xtick(ix)) '^o']};
            elseif  xtick(ix)<180 & xtick(ix)>0
                xtt(ix)={[num2str(xtick(ix)) '^oE']};
            else
                xtt(ix)={[num2str(-xtick(ix)) '^oW']};
                if 360-xtick(ix)<0;
                    xtt(ix)={[num2str(-(360-xtick(ix))) '^oE']};
                end
                    
            end
        end
        ytick=-90:dxtick:90;
        for ix=1:length(ytick)
            if (ytick(ix)==0 )
                ytt(ix)={[num2str(ytick(ix)) '^o']};
            elseif  ytick(ix)<0
                ytt(ix)={[num2str(-ytick(ix)) '^oS']};
            else
                ytt(ix)={[num2str(ytick(ix)) '^oN']};
            end
        end
        %-----------------------------------
         set(gca,'xtick',xtick,'xticklabel',xtt,'ytick',ytick,'yticklabel',ytt,'xcolor',ttc,'ycolor',ttc,'fontsize',ta)
        box on
        %xlim([0 360]+30)
        if k==2
            boxok=1;
if boxok
    x_lim=xlim;
    y_lim=ylim;
    hold on
    plot(x_lim([1 2 2 1 1]),y_lim([1 1 2 2 1]),'k','linewidth',1)
end

            %eval(['print -dpng -r300 outfigs/' figname]) 
           %eval(['print -dtiff -r300 outfigs/' figname])
        end
        %pause
        
        end

    end

end
end
pp
figure(1)
%-----------------
% find station
list_notinsta=[
    9415584 
    9415447 
    9415438 
    9415338 
    9415623
    9415266
    9414551
    9415498
    %------------
%     9415126
%     9415052
%     9415379
%     9414874
%     9414525
    ];
outfile='SF_TideStations_nocovered.txt';
fid=fopen(outfile,'wt')
for i=1:length( list_notinsta)
    iid=list_notinsta(i);
    loc=find(obs.id==iid);
    fprintf(fid,'%d %10.4f %10.5f\n',iid,obs.x(loc),obs.y(loc));
    hold on
    plot(obs.x(loc),obs.y(loc),'r*')
end
fclose(fid);
    
pp

sfbp=[-122.7396   38.2047
 -122.5092   37.8165
 -122.3969   37.5012
 -121.9036   37.3741
 -121.2181   37.6518
 -121.3761   38.8259
  -121.7532   38.8706
  
  -122.7396   38.2047
 ];

hold on
plot(sfbp(:,1),sfbp(:,2))