% This program converts frequency-direction distributions to
% wavenumber-direction distributions.
clear
close all
nu = 4;% factor in drift vel., u. nu = 4 works best and gives surface drift U0 of 3.2% of U10. 
color = ['b';'g';'r';'c';'m';'k';'b';'g';'r';'c';'m';'k';'b';'g';'r';'c';'m';'k';'b';'g';'r';'c';'m';'k';'b';'g';'r';'c';'m';'k';];
sym = ['.';'+';'o';'*';'s';'d';'.';'+';'o';];
name = 'waves';
dr=pi/180;dth = dr;% 1 degree res.
gamma = 0.074; % Surface tension.
rho = 1000; % Water density.
D = 12;%water depth in m.
D3 = D/3;% zero current depth
pim4 = (2*pi)^(-4);
rhorat = 1.25/1000.0;
kappa = 0.4;

kgmolwt = 0.029; %gm molecular wt. for air/1000 in Kgm/mol.  Earth
RRR = 8.314; % Universal gas constant in J/K/mol 
mua = 1.837e-5;%Dynamic viscosity of air at 20C [Pa s];
g=9.80665;

% compute cosine of angles.
th = ([1:360]-180)*pi/180;
cth = cos(th').^2*ones(1,240);
cth1 = cos(th')*ones(1,240);
sth = sin(th').^2*ones(1,240);

load wind
nr = 0;% number of spectra in set.
% runs = [23 24 28 34 36 49 50 55 87 89 91 92 112 114 115 116 127 128 132 139 140 142 143 144 145 146 156 157 158 160 161 172 173 175 185 186];
runs = 5;%[23 24 28 50 55 87 89 91 92 112 114 115 116 127 128 132 139 143 144 145 146 156 157 158 160 161 172 173 175 185 186];
% for run = runs;
% for run=1:200;%[189 82 62]
dr=pi/180;

for run=runs
    fnsave=['kspect\ksp_',int2str(run)];
eval(['load ',fnsave,' K wn dwn run runl Hs corr wsstats'])

f=wn;
beta=[];smm=[];dir=[];
F=K;%Normalize by number of points and directional resolution.
fm=find(sum(F)==max(sum(F(:,:)))); % Index of peak

for fr = 1:35
%fr = Frequency number.
fa=round(5+2*abs(f(fr)-f(fm))/f(fm));% Fit +- this angle.
if fa > 20 ; fa = 20; end

% Test pattern to check direction finding of peak.
if 1==2
F(:,1)=sech(5*dr*([1:360]'-111)).^2;
F(:,2)=sech(5*dr*([1:360]'-222)).^2;
end


% Move peak of DS to center.
s1=F(:,fr);
j0=0; j1=0; j2=0; pp(1)=0; smk=999;
fs1=find(s1>0);
if length(fs1) > 1;% Ignore spectral lines with less than 2 entries.
s2=[F(271:360,fr); F(:,fr); F(1:90,fr)];
smooth=RUNAVF(s2,1,10);
smooth=smooth(91:450);
smax=find(smooth == max(smooth));
smove=smax; % First location of peak.
if smove < 180;
% sss(181:360)=s1(smove:smove+179);
sss(180:360)=s1(smove:smove+180);
% sss(1:181-smove)=s1(smove+180:360);
sss(1:180-smove)=s1(smove+181:360);
% sss(180-smove:180)=s1(1:smove+1);
sss(180-smove+1:179)=s1(1:smove-1);
else
sss(541-smove:360)=s1(1:smove-180);
sss(1:180)=s1(smove-179:smove);
sss(181:540-smove)=s1(smove+1:360);
end

j0=smove-180;

sss=sss;
% Spectrum has been shifted by j0.


ss=dr*cumsum(sss);
sm=ss(360);smk=sm;
ss=ss/sm;
sh=max(ss)/2;
j1=find(abs(ss-sh)==min(abs(ss-sh)));j1=j1(1);j=j1+90;
figure(1);clf;plot([1:360],[ss; sss/max(sss)],[1 360],[sh sh],'--')

s2=[sss(271:360) sss sss(1:90)];
s=zeros(240,1);
s=s2(j-120:j+119);
ss=dr*cumsum(s);
sm=ss(240);
ss=ss/sm;
sh=max(ss)/2;
j2=find(abs(ss-sh)==min(abs(ss-sh)));j2=j2(1);j2=j2-120;
ss=(ss-sh)/sh;
pp=polyfit([120+j2-fa:120+j2+fa]*dr,ss([120+j2-fa:120+j2+fa]),1);
se=zeros(240,1);
se(1:240)=sech(pp(1)*dr*([-120-j2:119-j2]));se=se.*se;
figure(2);plot([1:240],RUNAVF(s,1,15),[1:240],se*sm*pp(1)/2,'.g')
xlabel(['wavelength = ',num2str(2*pi./wn(fr)),' m'])





% ss=dr*cumsum(sss);
% sm=ss(360);smk=sm;
% ss=ss/sm;
% sh=max(ss)/2;
% j1=find(abs(ss-sh)==min(abs(ss-sh)));j1=j1(1);j=j1+90;
% figure(1);clf;plot([1:360],[ss; sss/max(sss)],[1 360],[sh sh],'--')
% 
% s2=[sss(271:360) sss sss(1:90)];
% s=zeros(240,1);
% s=s2(j-120:j+119);
% ss=dr*cumsum(s);
% sm=ss(240);
% ss=ss/sm;
% sh=max(ss)/2;
% j2=find(abs(ss-sh)==min(abs(ss-sh)));j2=j2(1);j2=j2-120;
% ss=(ss-sh)/sh;
% pp=polyfit([120+j2-fa:120+j2+fa]*dr,ss([120+j2-fa:120+j2+fa]),1);
% se=zeros(240,1);
% se(1:240)=sech(pp(1)*dr*([-120-j2:119-j2]));se=se.*se;
% figure(2);plot([1:240],s,[1:240],se*sm*pp(1)/2,'.g')
% pause(1)
pause
end; % end of "if mean(s1"
beta=[beta pp(1)];
smm=[smm smk];
dir=[dir j0+j1+j2];
end


eval(['save betak\betak',int2str(run),' beta smm dir wn fm'])

end
