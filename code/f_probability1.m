function [umult,SIG_umult,f] = f_probability1(THX,MASK)

%Parameterization from Huss et al 2007 (unteraargletscher model) and
%distributed based on pixel thicknesses

% THXdist=cur.F.thx(cur.mask);


%% set up param distributions

N=30;
rs1=rand(N,1);
% 
% u_b = rs1(:,1);
% thx=50+400.*rs1(:,2);
% n=2.8+0.4.*rs1(:,3)

rs2=randn(N,1);

% rs3=lognrnd(-2,0.75,N,1);
u_b = (rs1(:,1));
u_b(u_b>1)=1; u_b(u_b<0)=0;
n=3+(0.2/3).*rs2(:,1);
thx=THX(MASK);

uS=50;

f=zeros(length(thx),N);
%%
for ip=1:N
    disp(ip)
    for iTH=1:length(thx)
        cTHX=thx(iTH);
        zs=linspace(1,cTHX,100);
        for z=1:length(zs);
           cSt(z) = u_b(ip)*uS+uS*(1-u_b(ip)).*(1-((cTHX-zs(z))/cTHX)^(n(ip)+1));
        end
        cSm(iTH)=trapz(cSt.*diff([0 zs]))./(uS.*cTHX);
        clear cSt
    end
    f(:,ip)=cSm;
    clear cSm
end

umult=NaN.*THX;
SIG_umult=f;
umult(MASK)=nanmean(f,2);
SIG_umult(MASK)=nanstd(f,[],2);



