% Plot polar frequency spectra in polar coordinates.
% Last edited:
% 2010.02.24 run on given days
% 2010.02.22 akm & MAD
% 2016.07.07 MAD
clear
close all
%% 

 

%% 
th = [1:360];
dth = pi/180;dr = pi/180;
cth = cos(th*dr); sth = sin(th*dr);
power=2;%fspectra are multiplied by f^'power'
maxf = 0.3 ; % resolution of polar plot
fh = 0.5;% High freq cut-off for period calculations.
fl = 0.02;% Low freq cut-off for period calculations (but: f(1)=0.0313).
%

runs=[6];
HSa = [];HM0a = [];TMa = [];TPa = [];WDMa = [];WDPa = [];WDPMa = []; 
ET = zeros(360,20); %Mean spectra.
HS = [];TM = [];WDM = [];TP = [];WDPM = [];WSST = [];TM02 = [];
HM0 = [];

for run = runs
figure(400);clf;
eval(['load fspectd\fsp_',int2str(run),' FCC FSECH corrFCC corrFSECH f df dir bbeta wnm SCDf Hs'])

% E = FCC;% Corrected spread spectrum
E = FSECH;% Corrected spread spectrum - sech
% Determine dfs and integrate
fnv=f/f(1);
fnvf=find(fnv>1.9999 & fnv < 2.0001);
nv=fnvf-1;
df=f*log(2)/nv;
lenf = length(f);
Hs=4*sqrt(sum(sum(E).*df.*f)*pi/180);

% Compute mean period.(Hs=4*sqrt(sum(sum(E).*df.*f)*pi/180);)
ifh = max(find(f < fh))  
ifl = min(find(f > fl))  
M0 = sum(sum(E(:,ifl:ifh),1)*dth.*f(ifl:ifh).*df(ifl:ifh));
M1 = sum(sum(E(:,ifl:ifh),1)*dth.*f(ifl:ifh).^2.*df(ifl:ifh));
Tm = M0/M1;
TM = [TM Tm];
Hm0 = 4*sqrt(M0);
HM0 = [HM0 Hm0];

% Compute mean second moment period, "TM02".(Hs=4*sqrt(sum(sum(E).*df.*f)*pi/180);)
ifh = max(find(f < fh))  
ifl = min(find(f > fl))  
M0 = sum(sum(E(:,ifl:ifh),1)*dth.*f(ifl:ifh).*df(ifl:ifh));
M1 = sum(sum(E(:,ifl:ifh),1)*dth.*f(ifl:ifh).^2.*df(ifl:ifh));
M2 = sum(sum(E(:,ifl:ifh),1)*dth.*f(ifl:ifh).^3.*df(ifl:ifh));
Tm2 = sqrt(M0/M2);
TM02 = [TM02 Tm2];

% Compute  mean direction.
Mc = sum(sum(cth'*ones(1,ifh-ifl+1).*E(:,ifl:ifh)).*dth.*f(ifl:ifh).*df(ifl:ifh));
Ms = sum(sum(sth'*ones(1,ifh-ifl+1).*E(:,ifl:ifh)).*dth.*f(ifl:ifh).*df(ifl:ifh));
Dir = atan2(Ms,Mc)/dr;
Dir = MOD360(Dir);

wdm = Dir;
WDM = [WDM Dir];


% Compute peak period.
Ef = sum(E,1)*dth.*f;
[b,a] = max(Ef);
if a < 2 | a > lenf - 1; a = 2;end
a1 = a-1;a2 = a; a3 = a+1;
sfE = f(a-1).*Ef(a-1)+ f(a).*Ef(a) + f(a+1).*Ef(a+1);
sE = Ef(a-1)+ Ef(a) + Ef(a+1);
fP = sfE./sE;
TP = [TP 1/fP];

% Compute  mean peak direction.
Mc1 = sum(cth'*ones(1,1).*E(:,a1)).*dth.*f(a1).*df(a1);
Ms1 = sum(sth'*ones(1,1).*E(:,a1)).*dth.*f(a1).*df(a1);
Mc2 = sum(cth'*ones(1,1).*E(:,a2)).*dth.*f(a2).*df(a2);
Ms2 = sum(sth'*ones(1,1).*E(:,a2)).*dth.*f(a2).*df(a2);
Mc3 = sum(cth'*ones(1,1).*E(:,a3)).*dth.*f(a3).*df(a3);
Ms3 = sum(sth'*ones(1,1).*E(:,a3)).*dth.*f(a3).*df(a3);
Mc = (Mc1 + Mc2 + Mc3)/3;
Ms = (Ms1 + Ms2 + Ms3)/3;
Dir = atan2(Ms,Mc)/dr;
Dir = MOD360(Dir);

wdpm = Dir;
WDPM = [WDPM Dir];

% Smooth spectrum in direction
for jk = 1:lenf
    E(:,jk) = RUNAVF(E(:,jk),1,10);
end

dr = pi/180;
th = ([1:360])*dr ;%
cth = cos(th); sth = sin(th);
k = f(1:lenf);
kx = k'*cth;ky = k'*sth;
kx = kx';ky = ky';kk = abs(kx + i*ky);
F = E(:,1:lenf);
file = 1;
polar_geo(th,maxf*ones(size(th)),'-k');hold on
FF(:,:) = F(MOD360(90 - [1:360]),:);% Return to math convention for contour.
hold on;contour(kx,ky,kk.^power.*FF);colorbar
ths =sprintf(' Hs=%4.1fm',Hs)
ttm =sprintf(' Tm=%4.1fs',Tm)
ttp =sprintf(' Tp=%4.1fs',1/fP)
twdm=sprintf(' WDM=%3.0f',wdm)
twdpm=sprintf(' WDPM=%3.0f',wdpm)
% tpower=sprintf(' power, n=%3.0f',power)
title([ths,ttm,ttp,twdm,twdpm])
xlabel(['frequency^{',int2str(power),'} \times polar spectrum.  Run ',int2str(run)],'fontsize',12)
pause

print('-dpng',['polarplot_',int2str(run),'_freq'])

end ; % end of run loop




