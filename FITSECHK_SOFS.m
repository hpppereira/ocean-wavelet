% This program reads in wavenumber-direction distributions from WDM in folder "kspect" and determines 
% direction (dir) and beta by fitting to (sech(beta*theta))^2 at each
% wavenumber. The directionally averaged (running over 5 degrees),F, is
% stored in the folder "betak" along with beta, dir, wn and dwn.
clear

dr=pi/180;
for runstart = 6;
nsum = 1;% Number of 10 min files to sum. Integer > 0.

F = [];nf = 0;
for run=runstart:runstart+nsum-1;%

        fnsave=['kspect\ksp_',int2str(run)];
eval(['load ',fnsave,' K wn dwn o run runl Hs run corr wsstats'])

m=runl;
res = 5;

% Smooth E over 5 degrees.
% F = E;
flen = length(wn);
for fn = 1:flen
     FF(:,fn) = RUNAVF(K(:,fn),1,res);% Smooth over res degrees.
end
if nf == 0; F = 0*K;end
F = F + FF;nf=nf+1;
end

F = F/nf;







f=wn;
beta=[];smm=[];dir=[];
fm=find(sum(F)==max(sum(F(:,10:end)))); % Index of peak

for fr = 1:o
%fr = Frequency number.
fa=round(5+2*abs(f(fr)-f(fm))/f(fm));% Fit +- this angle.
if fa > 20 ; fa = 20; end


% Move peak of DS to center.
s1=F(:,fr);
j0=0; j1=0; j2=0; pp(1)=0; smk=999;
fs1=find(s1>0);
if length(fs1) > 1;% Ignore spectral lines with less than 2 entries.
s2=[F(271:360,fr); F(:,fr); F(1:90,fr)];
smooth=RUNAVF(s2,1,10);
smooth=smooth(91:450);
% Find mean direction rather than max direction.
sdir = atan2(sum(smooth'.*sin((90-[1:360])*dr)),sum(smooth'.*cos((90-[1:360])*dr)));
% smax=find(smooth == max(smooth));
% smove=smax; % First location of peak.
smove = MOD360(round(sdir./dr));
smove = MOD360(90 - smove);% return to geographic coords.

if smove < 180;
sss(180:360)=s1(smove:smove+180);
sss(1:180-smove)=s1(smove+181:360);
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
figure(2);clf;plot([1:240],RUNAVF(s,1,15),[1:240],se*sm*pp(1)/2,'.g')
xlabel(['wavelength = ',num2str(2*pi./wn(fr)),' m'])






pause(1)
% pause
end; % end of "if mean(s1"
beta=[beta pp(1)];
smm=[smm smk];
dir=[dir j0+j1+j2];
end


eval(['save betak\betak_',int2str(run),' beta smm nsum dir wn dwn fm F'])
end
