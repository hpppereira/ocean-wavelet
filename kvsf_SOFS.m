%To plot output from WDM: weighted mean wavenumber in each frequency band.
%That is to determine the observed dispersion relation -- the relationship
%between frequency and observed wavenumber.

clear
g = 9.81;
gamma = 0.074; % Surface tension.
rho = 1025; % Water density.
dr = pi/180;

runstart = 3;
nruns = 3;

for run = runstart:runstart+nruns-1;
% % maxwn=0.16;
% % maxwn=0.32;
% maxwn= 1;
% kfac=120/maxwn;
% np = 3;
% dr=pi/180;
% set maximum wavenumber: usually to fix(2*pi/diameter of array)
maxwn=0.6;
kfac=120/maxwn;
figure(10*run);clf;
%
eval(['ex2=exist([''wdms/yw'',int2str(run),''.mat'']);'])
if ex2==2
run
eval(['load wdms/yw',int2str(run),' AAp ddd kkmp f d']);% ' AAp ddd kkmp f d'
D = d;%water depth at buoy in m.
  eval(['load wdms\data_',int2str(run),' ws WE NS heave'])

% nv = Number of voices in Wavelet transformation.
fnv=f/f(1);
fnvf=find(fnv>1.9999 & fnv < 2.0001);
nv=fnvf-1;


stdh=std(detrend(heave));
varh=stdh*stdh;
wsstats = STATSM(ws);
% 




jjj = run

ii = find(kkmp > 60.0);AAp(ii) = 0;
AA=AAp; 

fnv=f/f(1);
fnvf=find(fnv>1.9999 & fnv < 2.0001);
nv=fnvf-1;
ddd=MOD360(ddd);
ii=find(ddd == 0); ddd(ii)=ddd(ii)+360;
[m,n]=size(ddd);
flen = n;
runl=m;
df=f*log(2)/nv;

% Find spectral peak, fm.
   eval(['load fspect\fsp_',int2str(run)]);
    fm=find(sum(E).*f*dr==max(sum(E).*f*dr)); % Index of peak
fm;
fst = fm-1;fe = fst + 7;if fe > flen;fe = flen;end;%Start and stop frequencies.
for j=1:n
j;
% Compute weighted mean k = wnm = wavenumber measured.
wnm(j) = sum(AA(:,j).*kkmp(:,j))/sum(AA(:,j));
end
wn = WAVEKT(f,D,gamma,rho);wn = wn';% Theoretical wavenumber

% Plot dispersion relation -- measured and theoretical.
figure(28);clf;
 hold on;h3(1)= plot(f,wn,'.-k'); tr1 = ['theoretical dispersion.'];hold on
 hold on;h3(2)=plot(f,wnm,'.--r');tr2 = ['observed dispersion.'];hold on

set(gca,'fontsize',15)
% axis([0 0.7 0 2])
title(['k vs f ---> dispersion  Run ',int2str(run)],'fontsize',14)
ylabel('wavenumber [m^{-1}]','fontsize',14)
xlabel('frequency [Hz]','fontsize',14)
hleg = legend([h3],tr1,tr2,'location','SouthEast');

eval(['print -djpeg99 dispersion',int2str(run)]) 


% Plot contours of slope amplitude on k and f axes.

waven=[1:kfac*maxwn]/kfac;dk=1/kfac;wavenn=[1:kfac*maxwn+1]/kfac;
k=round(kkmp*kfac);
ik=find(k > maxwn*kfac);
kkmp(ik)=((maxwn*kfac+1)/kfac)*ones(size(ik));
AMP = [];AMP = zeros(121,n);
KAMP = [];KAMP = zeros(121,n);
for j = 1:n
    k=round(kkmp(:,j)*kfac);
  for kk = 1:maxwn*kfac+1
 ik=[];ik=find(k==kk);
  if ~isempty(ik);
         AMP(kk,j)=AMP(kk,j) + sum(AA(ik,j).^2)/runl;
         KAAMP(kk,j)=sqrt(AMP(kk,j)).*kk./kfac;
  end
  
  end
end

 figure(30);clf;
 [m1,n1]=size(KAAMP);
 m1 = m1 - 10;
 n1 = n1 - 2;


 pcolor(wavenn(1:m1),f(1:n1),(KAAMP(1:m1,1:n1)'));shading interp;colorbar
 hold on;h3(1)= plot(wn,f,'.-k'); tr1 = ['theoretical dispersion.'];hold on
 hold on;h3(2)=plot(2*wn,2*f,'.--k');tr2 = ['2nd bound harmonics.'];hold on

set(gca,'fontsize',15)
axis([0 0.6 0 0.7])
title(['slope contours  Run ',int2str(run)],'fontsize',14)
xlabel('wavenumber [m^{-1}]','fontsize',14)
ylabel('frequency [Hz]','fontsize',14)
hleg = legend([h3],tr1,tr2,'location','NorthWest');

eval(['print -djpeg99 dispersion_contour',int2str(run)]) 



end
pause(2)
end

