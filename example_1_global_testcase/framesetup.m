function [ax_po ,hb_po]=framesetup(gcf_po,na,nb,dfx,dfy)
%clf
%gcf_po=[4 343 1044 616];
%set(gcf,'Position',  gcf_po);
%hb_po=0;
%na=3;nb=3;
d1=1/nb;   d2=1/na;
d11=.85*d1;
%;
%d22=.8*d2;ddy=.15*d2;
d22=.9*d2;
ddy=d22*dfy;
ddx=d11*dfx;
d10=d1*.12;
d20=d10*gcf_po(3)/gcf_po(4)/1.5;
%ddy=.2*d22;

xt=0:0.2:1;
xtl=NaN*xt;
%------------------------------------------------------------------
iplot=0;
barw=20/gcf_po(3);
for iframe=1:na

    for jframe=1:nb
    iplot=iplot+1;
        %ax_po(iframe,:)=[d10*jframe+d11*(jframe-1)  1-d22*(iframe)-d20 d11 d22 ];
        %ax_po(iplot,:)=[d10+d11*(jframe-1)  1-d22*(iframe)-d20+0.05-((iframe-1)*ddy) d11 d22 ];
        ax_po(iplot,:)=[d10+d11*(jframe-1)+(jframe-1)*ddx  1-d22*(iframe)-d20-((iframe-1)*ddy) d11 d22 ];
        hb_po(iplot,:)=[d10+d11*(jframe-1)+(jframe-1)*ddx+d11+ddx/3  1-d22*(iframe)-d20-((iframe-1)*ddy) barw d22*.8 ];
    end
end
%save ax_hb ax_po hb_po