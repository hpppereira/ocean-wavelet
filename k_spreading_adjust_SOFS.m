% This program adjusts the spreading of the wind-sea to agree with the cross/down wave 
% slope ratio obtained for the frequency spectrum and plots normalized direction 
% distributions (FNn) from WDM(n=1) & WDM_corr(n=4)
% from one wavenumber below the peak to 6 wavenumbers above it.
% The swell wavenumbers are not adjusted. The wind-sea is adjusted from 3 wavenumbers below the peak. 
% The corrected spectrum FCC(k,phi) is stored in the folder "kspectd" as ksp_run#.
% Also stored:  FCC wn dwn dir bbeta wnm SCDk Hs.
% where wnm is energy-averaged wavenumber in each frequency band, f. SCDk
% is std(cross slope)/std(down slope) vs wavenumber. wnm and SCDf are read
% from the folder: "fspectd" created by the program: "f_spreading_adjust_SWIFT". 
% They are used here to determine the spreading with wavenumber of the wind-sea.


clear
color = ['b';'g';'r';'c';'m';'k';'b';'g';'r';'c';'m';'k';'b';'g';'r';'c';'m';'k';'b';'g';'r';'c';'m';'k';'b';'g';'r';'c';'m';'k';];
sym = ['.';'+';'o';'*';'s';'d';'.';'+';'o';];

clear
g = 9.81;
gamma = 0.074; % Surface tension.
dr = pi/180;
sf = 1/0.14/3;% Sampling frequency
runstart = 6;%First run
nsum = 1;% Number of runs

for run = runstart:runstart+nsum-1;
dr=pi/180;
% set maximum wavenumber: usually to fix(2*pi/diameter of array)
%
eval(['ex2=exist([''wdms/yw'',int2str(run),''.mat'']);'])
if ex2==2
run
fnsave=['fspectd\fsp_',int2str(run)];
eval(['load ',fnsave,' wnm SCDf'])
eval(['load wdms/yw',int2str(run),' AAp ddd kkmp f d']);% ' AAp ddd kkmp f d'
  eval(['load wdms\data_',int2str(run),' ws WE NS heave'])
eval(['load betak\betak_',int2str(run),' beta smm nsum dir wn dwn fm F'])
betakeep = beta;
% Find start of wind-sea from diff(tilt) > 0.01; max(heave-spectrum)
% after that yields f_p of wind-sea.
sp = SPECTF(WE,NS,.42,2);
sph = SPECTF(heave,.42,2);
sphh=sph(:,2);
ft=sp(:,1);
tilt = sp(:,2)+sp(:,3);
dtilt=[0; diff(tilt)];
jt = find(dtilt>0.01);
[dum,jt2]=max(sphh(jt(1):end));
fpk = ft(jt(1)+jt2);
fpk = WAVEK(fpk,d);
f = wn;df = dwn;

[dum,fpm]=min(abs(fpk-f));% fpm is index of peak wavenumber, fpk.

stdh=std(detrend(heave));
varh=stdh*stdh;

    flen = length(f);
    dirw = dir;% To locate legend
    fm=find(sum(F).*f*dr==max(sum(F).*f*dr)); % Index of peak
    swx = []; swy = [];
     for k = 1:flen
        DIR=(90-dir(k))*dr;% return to math coords to calculate CD

        F1([1:360],k) = F(MOD360([1:360]),k);
     end
%      % Add pulse to check directions
%      F1([1:360],50) = 0;F1(50,50)=1;dir(50)=50;
     eta=heave;
     Hs = 4*std(eta)
spw=f.*sum(F)*dr;
figure(3);clf;
h3(1)=loglog(f,spw,'*r');
grid on;hold on;
tr2 = ['WDM k-spectrum. Run = ',int2str(run)];hold on
set(gca,'fontsize',15)
title('Wavenumber spectra from WDM','fontsize',14)
xlabel('wavenumber [rad/m]','fontsize',14)
ylabel('Spectrum [m^3 / rad]','fontsize',14)
hleg = legend([h3],tr2,'location','SouthWest');

eval(['print -dpng WDM_k_spectra',int2str(run)]) 
pause(2)

% Convert CD(f) from folder "fspectd" to CD(k).
SCDk = interp1(wnm,SCDf,wn,'linear','extrap');

    % Compute corrected WDM based on CD = crosswave/downwave slope variances.
    CD = SCDk.*SCDk;
%     CD1 = swy./swx;% CD1 = crosswave/downwave slope variances for WDM spreading.

 load invbetavsSCD;% Polynomial of 1/beta vs SQRT(CD) created by "betavsSCD.m".
Beta = 1.0./(interp1(SCD,invbeta,sqrt(CD),'linear','extrap'));
f
fl = find(Beta > 40);Beta(fl)=40;
fl = find(Beta < 0.9);Beta(fl)=0.9;
fws = fpm - 3;%Start of wind-sea 3 wavenumbers below peak.
for k = fws:flen;% Loop over frequencies and move peak to 180 degrees.


F1([1:360],k) = F1(MOD360([1:360]+180+dir(k)),k);% Shift mean direction to 180';
end


 for j = 1:flen
        FN1(:,j) = F1(:,j)./(max(F1(:,j)));
        CD1(j) = sum(F1(90:270,j).*sin(([90:270]-180)*dr)'.^2)/sum(F1(90:270,j).*cos(([90:270]-180)*dr)'.^2);% cross/down mss.

 end
     
 
beta = 1.0./(interp1(SCD,invbeta,sqrt(abs(CD1)),'linear','extrap'));

fl = find(beta > 40);beta(fl)=40;
fl = find(beta < 0.9);beta(fl)=0.9;

    for k = fws:flen;% Loop over frequencies

ab = 0.3;% SQRT of Lower limit for correcting SECH^2.
secB = sech(Beta(k).*(([1:360]-180)*dr)');secBB = secB;fl=find(secB < ab);secB(fl)=ab;SecB=sum(secB.^2);
secb = sech(beta(k).*(([1:360]-180)*dr)');fl=find(secb < ab);secb(fl)=ab;Secb=sum(secb.^2);
Fc(:,k) = Secb/SecB.*F1(:,k).*(secB./secb).^2;% Directional spread corrected frequency spectrum Fc.
FSECH(:,k) = Fc(180,k).*secBB.^2;% sech^2 directional spread corrected frequency spectrum FSECH.
CDc(k) = sum(Fc(90:270,k).*sin(([90:270]-180)*dr)'.^2)/sum(Fc(90:270,k).*cos(([90:270]-180)*dr)'.^2);% Diagnostic check.
end
% %Iterate to find beta to yield observed CD.
for iter = 1:5 
beta = 1.0./(interp1(SCD,invbeta,sqrt(CDc),'linear','extrap'));
fl = find(beta > 40);beta(fl)=40;
fl = find(beta < 0.9);beta(fl)=0.9;

for k = fws:flen;% Loop over frequencies

secB = sech(Beta(k).*(([1:360]-180)*dr)');fl=find(secB < ab);secB(fl)=ab;SecB=sum(secB.^2);
secb = sech(beta(k).*(([1:360]-180)*dr)');fl=find(secb < ab);secb(fl)=ab;Secb=sum(secb.^2);
Fc(:,k) = Secb/SecB.*Fc(:,k).*(secB./secb).^2;% Directional spread corrected frequency spectrum Fc.
CDc(k) = sum(Fc(90:270,k).*sin(([90:270]-180)*dr)'.^2)/sum(Fc(90:270,k).*cos(([90:270]-180)*dr)'.^2);% Diagnostic check.
bbeta(k) = beta(k);
end
end
% 
% Renormalize directional spreads
FCC = F1;
for k = fws:flen;% only windsea
    Fc(:,k)  = Fc(:,k)*sum(F1(:,k))./sum(Fc(:,k));
    FSECH(:,k)  = FSECH(:,k)*sum(F1(:,k))./sum(FSECH(:,k));
FCC(1:360,k) = Fc(MOD360([1:360]+180-dir(k)),k);% Shift back to mean direction
FN5(:,k) = FSECH(:,k)./(max(FSECH(:,k)));
FSECH(1:360,k) = FSECH(MOD360([1:360]+180-dir(k)),k);% Shift back to mean direction
end
% 
% Check variance conservation.
Var = sum(sum(FCC(:,:))*dr.*f.*df);
corrFCC = 4*sqrt(Var)/Hs
Var = sum(sum(FSECH(:,:))*dr.*f.*df);
corrFSECH = 4*sqrt(Var)/Hs
    
% Save spread corrected WDM for windsea. Swell not corrected.
fnsave=['kspectd\ksp_',int2str(run)];wn = f;dwn = df;
eval(['save ',fnsave,' FCC FSECH corrFCC corrFSECH wn dwn dir bbeta wnm SCDk Hs'])

       for j = 1:flen
        F4([1:360],j) = Fc(MOD360([1:360]),j);
        FN4(:,j) = F4(:,j)./(max(F4(:,j)));
     end

    figure(1);clf;plot(f,dir,'*--k','linewidth',2);grid on
    set(gca,'fontsize',15)
    xlabel('wavenumber [rad/m]','fontsize',20)
    ylabel('Direction from [deg.]','fontsize',20)
    title(['mean direction of waves. Run ',int2str(run)],'fontsize',20)
    eval(['print -dpng k_Directions',int2str(run)]) 

    figure(2);clf;subplot(211)
    figure(2);hold on;h1=plot(f,betakeep,'.-b');grid on;t1 = ['WDM'];
    figure(2);hold on;h3=plot(f,bbeta,'o-r');grid on;t3 = ['WDM_{corr}'];
    figure(2);hold on;h4=plot(f,Beta,'.-k');grid on;t4 = ['(\sigma_x/\sigma_y)^2'];

      set(gca,'fontsize',15)
    ylabel('\beta','fontsize',20)
    title(['spreading (\beta) of waves. Run ',int2str(run)],'fontsize',20)
        hleg = legend([h1,h3,h4],t1,t3,t4,'location','NorthEast');


    figure(2);subplot(212)
    beta = 1./(sum(FN1(:,:))*dr);
    hold on;h1=plot(f,beta,'.-b');grid on;t1 = ['WDM'];
    
    Bbeta = (Beta/2).*coth(pi.*Beta/2);
    figure(2);hold on;h4=plot(f,Bbeta,'.-k');grid on;t4 = ['(\sigma_x/\sigma_y)^2'];

      set(gca,'fontsize',15)
    xlabel('wavenumber [rad/m]','fontsize',20)
    ylabel('A','fontsize',20)
    title(['spreading (A) of waves. Run ',int2str(run)],'fontsize',20)
    hleg = legend([h1,h4],t1,t4,'location','NorthEast');
    eval(['print -dpng k_A_beta',int2str(run)]) 

   
    % Plot 8 panels of FNn from fm-1 to fm+6
    figure(3);clf;
    fe = fpm + 6;if fe > 56;fe = 56;end
    n = 0;
    for j = fpm-1:fe
        n=n+1;
        if n < 9;
        subplot(4,2,n);
         hold on;h3 = plot([1:360],[FN1(:,j)],'.-b');grid on;t3 = ['WDM'];
        hold on;h4 = plot([1:360],[FN4(:,j)],'.-k');grid on;t4 = ['WDM_{corr}'];
        hold on;h5 = plot([1:360],[FN5(:,j)],'.-r');grid on;t5 = ['WDM_{sech}'];

    axis([0 360 0 1.2])
        xlabel('direction (from) [deg.]','fontsize',10)
    title(['k/k_{p} = ',num2str(f(j)/f(fpm),3),',       \sigma_{x} / \sigma_{y} = ',...
        num2str(1./SCDk(j),3),',   \lambda = ',num2str(2*pi./f(j)),' m'],'fontsize',10)
       if n == 1 & dirw(fpm-1)<=180;hleg = legend([h4,h3,h5],t4,t3,t5,'location','NorthEast');end
   if n == 1 & dirw(fpm-1)>180;hleg = legend([h4,h3,h5],t4,t3,t5,'location','NorthWest');end
 
        end
    end
    
        eval(['print -dpng k_dirsp',int2str(run)]) 
        
        


end
pause(2)
end
    
  
    

