% This program takes output of WDM and summarizes in wavenumber polar space.
% and stores the wavenumber-direction spectrum ("ksp") in the folder
% "kspect". The wavenumbers are logarithmically spaced from k1 to k2 in o
% bins. These are set in lines 9 & 10.

dr=pi/180;
% set maximum wavenumber: usually to (2*pi/diameter of array), k2
% set minimum wavenumber: usually to (2*pi/Length of longest waves to be analysed), k1
o = 35;% number of wavenumbers.
k1 = 0.03;k2=6;dlnk=(log(k2)-log(k1))/(o-1);%Set wavenumbers logarithmically spaced.
k = exp(log(k1)+[0:o-1]*dlnk);
dk = k.*dlnk;

figure(10*run);clf;
%
eval(['ex2=exist([''yw'',int2str(run),''.mat'']);'])
if ex2==2
run
eval(['load yw',int2str(run)])

stdh=mean(std(detrend(ws)));
varh=stdh*stdh;
wsstats = stats(ws);

clear ws data

[m,n]=size(AAp);
runl=m;
K=zeros(360,o);
ddd=mod360(ddd);
ii=find(ddd == 0); ddd(ii)=ddd(ii)+360;

%nv = Number of voices in Wavelet transformation.
fnv=f/f(1);
fnvf=find(fnv>1.9999 & fnv < 2.0001);
nv=fnvf-1;
df=f*log(2)/nv;
for kk = 1:o
 run
 kk
 ik=[];ik=find(kkmp >= k(kk)-dk(kk)/2 & kkmp < k(kk)+dk(kk)/2);
  if ~isempty(ik);
   for di = 1:360
    id=[];id=find(ddd(ik)==di);
    K(di,kk)=K(di,kk) + sum(AAp(ik(id)).^2)/nv*1.03565/runl;
   end
  end
end
K=(K./(ones(360,1)*(k.*dk)))*180/pi;
corr=varh/(sum(sum(K).*k.*dk)*pi/180);
corr

Hs=4*sqrt(sum(sum(K).*k.*dk)*pi/180);

fnsave=['kspect\ksp_',int2str(run)];
eval(['save ',fnsave,' K k dk run runl Hs varh corr wsstats'])

subplot(221)
loglog(k,sum(K(:,1:o)).*k*pi/180,'.-'); grid on
xlabel('wavenumber [m^{-1}]')
ylabel('spectral density [m^3]')
subplot(222)
contour(1:360,k(1:o),K(:,1:o)'.*(k(1:o)'.^1*ones(1,360)));grid
title('Energy spectrum')
ylabel('wavenumber [m^{-1}]')
xlabel('direction [degrees]')
subplot(223)
contour(1:360,k(1:o),K(:,1:o)'.*(k(1:o)'.^3*ones(1,360)));grid
title('k^2 \times Energy spectrum')
ylabel('wavenumber [m^{-1}]')
xlabel('direction [degrees]')
subplot(224)
contour(1:360,k(1:o),K(:,1:o)'.*(k(1:o)'.^4*ones(1,360)));grid
title('k^3 \times Energy spectrum')
ylabel('wavenumber [m^{-1}]')
xlabel('direction [degrees]')
pause(1)

clear K k run runl Hs corr AAp ddd kk df di ik k dk
else
end
