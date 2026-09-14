clear,clf
reloadok=0;
addcok=1;
%testcase='G5k5s';
testcase='G5km6';
testcase='GA5k151';
%testcase='GA2p5km151';
% testcase='GA1p8km151';
% testcase='GN2p5km_p9_151';
% testcase='G300m_2p5km_151';
testcase='GFEN2p5km300m_151';
 %testcase='GFWE2p5km300m_151';

%testcase='GW4km_151';
%testcase='GFTB5km_151';

%testcase='GF5km_151';
%testcase='G5kR1';
splitok=0;
if reloadok==1
    pathin=['../input_data_' testcase '/'];
    temp=dir([pathin '*14*'])
    gridfile=temp.name;
    all=f_read_bathy_adcirc_grd(pathin,gridfile,0);
else
    pathin=['input_data_' testcase '/'];
    temp=dir([pathin '*14*']);
    gridfile=temp.name;
    gridfile=[gridfile(1:end-3) '_all.mat' ]
    eval(['load ' pathin gridfile ' all']);5
end
%---------------------------------
% plotting 

TPC_list={ 'O1';'M2';'S2';'K1';...
    'K2';'N2';'P1';'Q1';...
    'M4';'MS4';%'MSF';M6
    'MU2';'NU2';'L2'};
CSlist={'Cartesian';'ADCIRC';'Spherical';'Spherical2';'Local';'Local2'};
Aelist=[1 2];
data=[];
%load cm
load cm_jet_blue
colormap(cm)
colormap jet
tlist={'Amplitude';'Phase'};
%cm=flipud(cm);
dx=0.04;
xy_lim=[min(all.nodes(:,2))-dx max(all.nodes(:,2))+dx  min(all.nodes(:,3))-dx max(all.nodes(:,3))+dx];
lon0=30;
xy_lim(1:2)=[0 360]+lon0;
xy_lim(3:4)=[-1 1]*90;

%xy_lim=[114 125 37 41];
tcok=1;
ta=10;tt=12;
bkc=[0 0 0]+1;
ttc=[0 0 0];
ax_po=[0.05 0.08 0.93 0.85];
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
        % if ii<9
        % iname=['output/endresults_' testcase '*' ITPC '*ITE1*' CSlist{i} '_Ae' int2str(j) '_*bs*_T2_*' inv1 '*'];
        % else
        % iname=['output/endresults_' testcase '*' ITPC '*ITE*' CSlist{i} '_Ae' int2str(j) '_*bs*_T0_*' inv1 '*CF8*'];
        % end
         iname=['out_' testcase '/endresults_' testcase '*' ITPC '*ITE1*' CSlist{i} '_Ae' int2str(j) '_*bs*_T2_*' inv1 '*'];

        temp=dir(iname);
        ifile=temp.name;
        eval(['load out_' testcase '/' ifile]);
        fprintf(1,'load %s\n',ifile)
        if corrok
            NPSI=corrPSI;
            addt='corrNPSI';
        else
            addt='NPSI';
        end
        %----------------------
        for k=1:2
            if k~=2
                figure(k)
                clf
                %set(gcf,'visible','off')
                %set(gcf)
                %set(gcf,'position',[10 100 [1000 500]/1000*1600],'color','w')
                set(gcf,'position',[10 10 [900 500]/1000*1500],'paperpositionmode','auto','Scrollable','off')
               

                
            end
            hax=axes('position',ax_po);
            set(hax,'color',bkc)
            if k==1
                    z=abs(NPSI);
                    %all2=f_split_global(all,z);
                    %z2=[z1; z1;z1];
                    %z=griddata(x,y,z2,all0.nodes(:,2),all0.nodes(:,3));
                    cv=[ 0:0.05:1.6];
                     c_lim=[0 1.5];
                     c_lim=[0 0.05];
                     c_lim=[0 0.2];
                     s5=0.2;
                switch ITPC
                    case 'MU2'
                        c_lim=[0 .05];
                    case 'N2'
                        cv=[ 0:0.02:1.6];
                        c_lim=[0 .2];
                    case 'O1'
                        c_lim=[0 .5];
                    case 'M2'
                        c_lim=[ 0 1];
                    case 'S2'
                        cv=[ 0:0.04:1.6];
                        c_lim=[0 .5];
                    case 'K1'
                        c_lim=[0 .5];
                end
               
                %c_lim=[0 0.7];
                sf=sprintf('%.2f',max(z));
                %all2=f_split_global(all,z);
               
            else
                z=mod(180*(-1)*angle(NPSI)/pi(),360);

                c_lim=[0 360];
                cv=[60:60:360];
                cv=[30:30:330];
                %cv=[60:30:330];
        end
        all2=f_split_global(all,z,lon0);
        %colormap(jet)
        colormap(cm)
        %colormap(flipud(jet))
        if k~=2
        trisurf(all2.elements(:,3:5),all2.nodes(:,2),all2.nodes(:,3),all2.nodes(:,end));
        %trisurf(all2.elements,all2.x,all2.y,all2.z);
        
        hold on
        %trisurf(all.elements(loc,3:5),all.nodes(:,2)+360,all.nodes(:,3),z);
        view([0 90])
        caxis(c_lim)
        shading  interp
        %------------------------
        % This is to double check dart location
        plotdartok=1;
        if k==1 & plotdartok==0;
            
            %text(118, 40.5,['max=' sf ' m'])
            outfile=[pathin '/' testcase '_inpsiLLfinal_' ITPC '.dat' ];
            data=load(outfile);
            xyp=data(:,[5 6]);
            loc=find(xyp(:,1)>180);
            xyp(loc,1)=xyp(loc,1)-360;
            hold on
            plot3(xyp(:,1),xyp(:,2),data(:,2),'ro','markersize',5)
            %scatter(xyp(:,1),xyp(:,2),8,data(:,2)/100,'markerfacecolor','flat','MarkerEdgeColor','black')
            loc=find(xyp(:,1)<0);
            xyp(loc,1)=xyp(loc,1)+360;
            %axis(xy_lim)
            
            plot3(xyp(:,1),xyp(:,2),data(:,2),'ro','markersize',5)
            %scatter(xyp(:,1),xyp(:,2),4,data(:,2))
        end
        %---------------------------
         set(gca,'fontsize',ta,'color',bkc)

        hc=colorbar;
        kk=1;
        hc_po=[.21, .64,0.01,.18];
        set(hc,'position',hc_po,'ytick',[ 0:.25:1],'fontsize',ta,'AxisLocation','in')
        y1=ylabel(hc,'Amplitude (m)');
        %yl.Position(1)=min(xlim(hc))-100;
        %hc.Label.VerticalAlignment = 'bottom';

        %set(hc,'position',[ax_po(kk,1)+ax_po(kk,3)+0.03, ax_po(kk,2)+.1,0.03,ax_po(kk,4)*.7],'fontsize',tt,'AxisLocation','in')
        
        
       
        
        %title([strrep(testcase, '_','-') ' ' ITPC '  '  tlist{k} ' ' tlist{k+1} '  ' strrep(ifile(16:end-4),'_','-') ' ' addt],'fontsize',16)
        title( [ITPC '  '  tlist{k} ' ' tlist{k+1} '  '  addt ' ' strrep(testcase, '_','-') ],'fontsize',16)
        box on
        grid off
        axis(xy_lim)
        axis on
        %---------------
        
        if k==1 & plotdartok==1;
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
            scatter(xyp(:,1),xyp(:,2),30,data(:,2)/100,'markerfacecolor','flat','MarkerEdgeColor','black')
            axis(xy_lim)
            
        end
        %------------------
        end
        if tcok==1 & k>1
            axis off
            hax=axes('position',ax_po);
            set(hax,'color','none')
            set(gca,'color','none')
             nodes=all2.nodes(:,2:3);
             elements=[all2.elements(:,3:5)];
            hold on
            [cout,hout] = tricontour(nodes,elements,all2.nodes(:,end),cv);
            %tricontour(elements,nodes(:,1),nodes(:,2),all2.nodes(:,end),cv);         
            %tricontour(elements,nodes(:,1),nodes(:,2),all2.nodes(:,end),cv);     
            caxis(c_lim)
            axis(xy_lim)

            hc=colorbar;
            ylabel(hc,'Phase (deg)')
            %hc.Title.String='(m)     (deg)';
            kk=1;
            %set(hc,'position',[ax_po(kk,1)+ax_po(kk,3)+0.03, ax_po(kk,2),0.03,ax_po(kk,4)/1.5],'fontsize',ta,'ytick',0:30:360,'xcolor',ttc,'ycolor',ttc)
            set(hc,'position',hc_po,'fontsize',ta,'ytick',0:90:360,'xcolor',ttc,'ycolor',ttc)
            %title(hc,'(m)     (deg)','fontsize',ta)
            %text(261.3802  -15.5834,'(m)     (deg)')
            
        end
        figname=['outfigs/f' testcase '_' ifile(18:end-4) '_' int2str(k) '_' addt '.png'];
        grid off
        %-------------------------------------
        xtick=-180:30:200+360;
        %xtt=[];
        for ix=1:length(xtick)
            if (xtick(ix)==0 | xtick(ix)==180)
                xtt(ix)={[int2str(xtick(ix)) '^o']};
            elseif  xtick(ix)<180
                xtt(ix)={[int2str(xtick(ix)) '^oE']};
            else
                xtt(ix)={[int2str(360-xtick(ix)) '^oW']};
                if 360-xtick(ix)<0;
                    xtt(ix)={[int2str(-(360-xtick(ix))) '^oE']};
                end
                    
            end
        end
        ytick=-90:30:90;
        for ix=1:length(ytick)
            if (ytick(ix)==0 )
                ytt(ix)={[int2str(ytick(ix)) '^o']};
            elseif  ytick(ix)<0
                ytt(ix)={[int2str(-ytick(ix)) '^oS']};
            else
                ytt(ix)={[int2str(ytick(ix)) '^oN']};
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
                plot(x_lim([1 2 2 1 1]),y_lim([1 1 2 2 1]),'k','linewidth',2)
            end

            eval(['print -dpng -r150 ' figname]) 
            %exportgraphics(gcf, figname, 'Resolution', 150);
            %eval(['print -dpng -r300 outfigs/' figname])
           %eval(['print -dtiff -r300 outfigs/' figname])
        end
       
        
        end

    end

end
end

