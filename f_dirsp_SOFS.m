% This program plots normalized direction distributions (FNn) from WDM(n=1), MLM(n=2), MEM(n=3) & WDM_corr(n=4) 
% from 1 frequency below the peak to 7 above the peak.

clear
color = ['b';'g';'r';'c';'m';'k';'b';'g';'r';'c';'m';'k';'b';'g';'r';'c';'m';'k';'b';'g';'r';'c';'m';'k';'b';'g';'r';'c';'m';'k';];
sym = ['.';'+';'o';'*';'s';'d';'.';'+';'o';];
tic
dr=pi/180;
g = 9.81;
gamma = 0.074; % Surface tension.
rho = 1000; % Water density.
D = 12;%water depth in m.

% load wind

for run=[3]
% for run=[62]
    figure(1);clf;
    figure(2);clf;
    figure(3);clf;

    
  eval(['load \SOFS\yw',int2str(run)])
  eval(['load \SOFS\fspect\fsp_',int2str(run)])
  eval(['load \SOFS\beta\beta_',int2str(run)])
  
      eval(['load \SOFS\wdms\data_',int2str(run),' ws WE NS heave'])

    % Compare Fourier and Wavelet 1-D frequency Spectra.



%sf=sample_rate;
sf = 1/0.14;
% if length(data) > 10244095
%     
% %data=binavg(data,2);   % binavg;  Modify as necessary to subsample time series.
% 
% % Reduce data to a multiple of 4096 points: wavelet analysis done 4096
% % point chunks; i.e. 2^12 points
% tlen = fix(length(data)/4096);
%    data=data(1:4096*tlen,:);
% 
% % Apply calibrations [to metres] to the wave staffs and set in desired order.
% 
% ws=[.58467*data(:,3) .57756*data(:,1) .59762*data(:,4) ...
%    .58998*data(:,6) .57155*data(:,5) .57363*data(:,2) ];
% 
% stdh=mean(std(detrend(ws)));
% varh=stdh*stdh;
% Hs = 4*stdh;
% wsstats = STATSM(ws)

% % Get wind data from "wind" and from Met87/R87nnn.
% WD = 99; U1 = 99; WD2 = 99; U2 = 99; Cd2 = 99;
%     j8 = find(RR==run);
%     if length(j8) > 0;
%     WD2 = MEANANG(WinD(j8).*dr)./dr;WD2 = MOD360(WD2); 
%     U2 = mean(Ud(j8));
%     Cd2 = mean(Cd(j8));
%     end
%nv = Number of voices in Wavelet transformation.
fnv=f/f(1);
fnvf=find(fnv>1.9999 & fnv < 2.0001);
nv=fnvf-1;
df=f*log(2)/nv;
dlnf = log(2)./nv;

dwn = dlnf.*f;
wnm = exp(log(f)-dlnf/2);% lower limit at each wn
wnp = exp(log(f)+dlnf/2);% upper limit at each wn
wn = f;
% 


% Find mss(f) from timeseries of slopes rotated into wave propagation directions
dws = [];cws = [];rws = [];
dirm = 90 - dir;% Geog to math directions.
for k = 1:20;% Loop over frequencies
DIR=dirm(k)*dr;% For numerical simulations DIR = pi.
sx=-WE*cos(DIR)+NS*sin(DIR);% Downwave slope time series.
sy=WE*sin(DIR)+NS*cos(DIR);% Crosswave slope time series.
sq=SPECTF(sx,sy,1/sf,2,1);
fq = sq(:,1);
sqx = sq(:,2);
sqy = sq(:,3);
% Sum linear freqs into logarithmic bins.
ii = find(fq>=wnm(k) & fq<wnp(k));
sxx(k) = sum(sqx(ii));
syy(k) = sum(sqy(ii));
ft = sq(:,1);
ftlen = length(ft);
% dws(k,:) = sq(:,2);
% cws(k,:) = sq(:,3);
% rws(k) = dws(k)./cws(k);
end
DC = sxx./syy;
% for dif = 1:ftlen
% %     ji = find(rws(:,dif)>0.995*max(rws(:,dif)));
% %     dir(dif) = MOD360(90-ji(1));
%     ji = find(diff(diff(rws(180:360,dif)))==min(diff(diff(rws(180:360,dif)))));
%     dir(dif) = MOD360(90-ji(1)+1);
% end
%     j1 = find(dir > 180);dir(j1) = dir(j1) - 360;
% % figure(31);clf;plot([1:360],dws(:,40)./cws(:,40),'.-')
% % figure(31);clf;plot([1:360],dws(:,40)./cws(:,40),'.-')
% figure(32);clf;plot([1:360],[rws(:,20)],'.-');grid on
% figure(37);clf;plot([1:358],[diff(diff(rws(:,20)))],'.-');grid on
% figure(34);clf;plot(ft,dir,'.-')
% 
%     if run < 10;eval(['load Met87/R8700',int2str(run),'av']);end
%     if run < 100 & run > 9;eval(['load Met87/R870',int2str(run),'av']);end
%     if run < 1000 & run > 99;eval(['load Met87/R87',int2str(run),'av']);end
%     eval(['load wdms/fspect/fsp_',int2str(run)]);% for "runl"
%     eval(['load betaf\beta',int2str(run)]);
%     lenR = fix(runl/4/30);% Number of 30 second bits in Met87 data.
%        WD = MEANANG(av(1:lenR,1).*dr)./dr; WD = MOD360(WD);
%     U1 = mean(av(1:lenR,2));
    
    nv=4;%Number of voices in Wavelet transformation.
    df=f*log(2)/nv;
    flen = length(f);
    dirw = dir;% To locate legend
    fm=find(sum(F).*f*dr==max(sum(F).*f*dr)); % Index of peak
    figure(1);hold on;plot(f/f(fm),MOD360(dir),'.-b');grid on
     for j = 1:flen
        F1([1:360],j) = F(MOD360([0:359]),j);
        FN1(:,j) = F1(:,j)./(max(F1(:,j)));
     end
     eta=ws(:,1);
sp=SPECTF(eta,1/sf,32,1);
spw=f.*sum(F)*dr;
figure(3);clf;
h3(1)=loglog(sp(:,1),sp(:,2),'.-b');hold on
h3(2)=loglog(f,spw,'*r');
grid on;hold on;
tr1 = ['FFT spectrum. Run = ',int2str(run)];hold on
tr2 = ['WDM spectrum. Run = ',int2str(run)];hold on
set(gca,'fontsize',15)
title('comparison of spectra vs frequency','fontsize',14)
xlabel('frequency [Hz]','fontsize',14)
ylabel('Spectrum [m^2 / Hz]','fontsize',14)
hleg = legend([h3],tr1,tr2,'location','SouthWest');

eval(['print -djpeg99 \SOFS/WDM_FFT_spectra',int2str(run)]) 
%pause

    % Compute corrected WDM based on CD = crosswave/downwave slope variances.
    CD = syy./sxx;
%     CD1 = swy./swx;% CD1 = crosswave/downwave slope variances for WDM spreading.
%     load invbetavsCD;% Polynomial of 1/beta vs CD.
% Beta = 1.0./(polyval(pbc,CD));
% nb = length(pbc);
% fl = find(Beta > 10);Beta(fl)=10;
% fl = find(Beta < 0.1);Beta(fl)=0.1;
 load invbetavsSCD;% Polynomial of 1/beta vs SQRT(CD) created by "betavsSCD.m".
% Beta = 1.0./(polyval(pbc,sqrt(CD)));
Beta = 1.0./(interp1(SCD,invbeta,sqrt(CD),'linear','extrap'));

nb = length(pbc);
fl = find(Beta > 40);Beta(fl)=40;
fl = find(Beta < 0.1);Beta(fl)=0.1;
for k = 1:flen;% Loop over frequencies and move peak to 180 degrees.
    smove=dir(k); % First location of peak.
    s1 = F(:,k);
if smove < 180;
    smove(smove<=0)=1;
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

F1([1:360],k) = sss';


end

 for j = 1:flen
        FN1(:,j) = F1(:,j)./(max(F1(:,j)));
        CD1(j) = sum(F1(90:270,j).*sin(([90:270]-180)*dr)'.^2)/sum(F1(90:270,j).*cos(([90:270]-180)*dr)'.^2);% cross/down mss.

 end
     

 % beta = 1.0./(polyval(pbc,sqrt(CD1)));
beta = 1.0./(interp1(SCD,invbeta,sqrt(CD1),'linear','extrap'));

fl = find(beta > 40);beta(fl)=40;
fl = find(beta < 0.1);beta(fl)=0.1;

for k = 1:flen;% Loop over frequencies

ab = 0.0001;% SQRT of Lower limit for correcting SECH^2.
secB = sech(Beta(k).*(([1:360]-180)*dr)');fl=find(secB < ab);secB(fl)=ab;SecB=sum(secB.^2);
secb = sech(beta(k).*(([1:360]-180)*dr)');fl=find(secb < ab);secb(fl)=ab;Secb=sum(secb.^2);
Fc(:,k) = Secb/SecB.*F1(:,k).*(secB./secb).^2;% Directional spread corrected frequency spectrum Fc.
CDc(k) = sum(Fc(90:270,k).*sin(([90:270]-180)*dr)'.^2)/sum(Fc(90:270,k).*cos(([90:270]-180)*dr)'.^2);% Diagnostic check.
end
% %Iterate to find beta to yield observed CD.
for iter = 1:1
% beta = 1.0./(polyval(pbc,sqrt(CDc)));
beta = 1.0./(interp1(SCD,invbeta,sqrt(CDc),'linear','extrap'));
fl = find(beta > 40);beta(fl)=40;
fl = find(beta < 0.1);beta(fl)=0.1;

for k = 1:flen;% Loop over frequencies

secB = sech(Beta(k).*(([1:360]-180)*dr)');fl=find(secB < ab);secB(fl)=ab;SecB=sum(secB.^2);
secb = sech(beta(k).*(([1:360]-180)*dr)');fl=find(secb < ab);secb(fl)=ab;Secb=sum(secb.^2);
Fc(:,k) = Secb/SecB.*Fc(:,k).*(secB./secb).^2;% Directional spread corrected frequency spectrum Fc.
CDc(k) = sum(Fc(90:270,k).*sin(([90:270]-180)*dr)'.^2)/sum(Fc(90:270,k).*cos(([90:270]-180)*dr)'.^2);% Diagnostic check.
bbeta(k) = beta(k);
end
end

% 
% Renormalize directional spreads
for k = 1:flen
    Fc(:,k)  = Fc(:,k)*sum(F1(:,k))/sum(Fc(:,k));
end
% 
% Check variance conservation.
Var = sum(sum(Fc(:,:))*dr.*f.*df);
corr = 4*sqrt(Var)/Hs
% Save spread corrected WDM all frequencies moved to 180 deg.
fnsave=['fspectc\fsp_',int2str(run)];
eval(['save ',fnsave,' Fc f df dir bbeta Hs CD1 CDc'])

for k = 1:flen;% Loop over frequencies and move peak back to dir degrees.
    smove=dir(k); % First location of peak.
    s1 = Fc(:,k);
if smove < 181;
    smove(smove<=0)=1;
sss(smove:smove+180)=s1(180:360);
sss(smove+180:360)=s1(1:181-smove);
sss(1:smove-1)=s1(180-smove+1:179);
else
sss(1:smove-180)=s1(541-smove:360);
sss(smove-179:smove)=s1(1:180);
sss(smove+1:360)=s1(181:540-smove);
end

Fcd([1:360],k) = sss';


end
fnsave=['fspectd\fsp_',int2str(run)];
eval(['save ',fnsave,' Fcd f df dir bbeta Hs CD1 CDc'])


       for j = 1:flen
        F4([1:360],j) = Fc(:,j);
        FN4(:,j) = F4(:,j)./((F4(180,j)));
     end
  % Plot CDc vs CD as check
  figure(5);clf;plot(CD,CDc,'*r',[0 1],[0 1],'--b');grid on;
%   pause

% % MLM
%     eval(['load betafmlm\beta',int2str(run)]);
%     figure(1);hold on;plot(f/f(fm),MOD360(dir-180),'.-g');grid on
%     % Convert to geographic from directions.
%         for k = 1:flen;% Loop over frequencies and move peak to 180 degrees.
%     smove=MOD360(dir(k))+1; % First location of peak.
%     s1 = F(:,k);
% if smove < 180;
% sss(181:360)=s1(smove:smove+179);
% sss(1:181-smove)=s1(smove+180:360);
% sss(180-smove:180)=s1(1:smove+1);
% else
% sss(541-smove:360)=s1(1:smove-180);
% sss(1:180)=s1(smove-179:smove);
% sss(181:540-smove)=s1(smove+1:360);
% end
% F2([1:360],k) = sss';
%         end
% for j = 1:flen
%         FN2(:,j) = F2(:,j)./(max(F2(:,j)));
%     end
%         CD2 = swy./swx;% CD2 = crosswave/downwave slope variances for MLM spreading.
% 
% % MEM        
%     eval(['load betafmem\beta',int2str(run)]);
%     figure(1);hold on;plot(f/f(fm),MOD360(90-dir),'.-r');grid on
%      % Convert to geographic from directions.
%      for k = 1:flen;% Loop over frequencies and move peak to 180 degrees.
%     smove=MOD360(dir(k))+1; % First location of peak.
%     s1 = F(:,k);
% if smove < 180;
% sss(181:360)=s1(smove:smove+179);
% sss(1:181-smove)=s1(smove+180:360);
% sss(180-smove:180)=s1(1:smove+1);
% else
% sss(541-smove:360)=s1(1:smove-180);
% sss(1:180)=s1(smove-179:smove);
% sss(181:540-smove)=s1(smove+1:360);
% end
% F3([1:360],k) = sss';
%      end
% 
%     for j = 1:flen
%         FN3(:,j) = F3(:,j)./(max(F3(:,j)));
%     end
%         CD3 = swy./swx;% CD3 = crosswave/downwave slope variances for MEM spreading.
WD = dir(18);
    figure(1);hold on;plot([2 f(end)/f(fm)],[WD WD],'--k','linewidth',2);grid on
    set(gca,'fontsize',15)
    text(2,WD+10,'Wind')
    xlabel('f/f_p','fontsize',20)
    ylabel('Direction from [deg.]','fontsize',20)
    axis([0 15 0 360])
    title(['mean direction of waves. Run ',int2str(run),'.  U = ',num2str(99,2)],'fontsize',20)
    eval(['print -djpeg99 \SOFS/Directions',int2str(run)]) 
    pause(5)
    figure(2);subplot(211)
     eval(['load beta\beta_',int2str(run)]);
    fm=find(sum(F).*f*dr==max(sum(F).*f*dr)); % Index of peak
    figure(2);hold on;h1=plot(f/f(fm),beta,'.-b');grid on;t1 = ['WDM'];
%     eval(['load betafmlm\beta',int2str(run)]);
%     figure(2);hold on;h2=plot(f/f(fm),beta,'.-g');grid on;t2 = ['MLM'];
%     eval(['load betafmem\beta',int2str(run)]);
%     figure(2);hold on;h3=plot(f/f(fm),beta,'.-r');grid on;t3 = ['MEM'];
    figure(2);hold on;h4=plot(f/f(fm),Beta,'.-k');grid on;t4 = ['(\sigma_x/\sigma_y)^2'];

      set(gca,'fontsize',15)
    xlabel('f/f_p','fontsize',20)
    ylabel('\beta','fontsize',20)
    title(['spreading (\beta) of waves. Run ',int2str(run)],'fontsize',20)
        hleg = legend([h1,h4],t1,t4,'location','NorthEast');


    figure(2);subplot(212)
    beta = 1./(sum(FN1(:,:))*dr);
    hold on;h1=plot(f/f(fm),beta,'.-b');grid on;t1 = ['WDM'];
%         beta = 1./(sum(FN2(:,:))*dr);
%     hold on;h2=plot(f/f(fm),beta,'.-g');grid on;t2 = ['MLM'];
%         beta = 1./(sum(FN3(:,:))*dr);
%     hold on;h3=plot(f/f(fm),beta,'.-r');grid on;t3 = ['MEM'];
%     figure(2);hold on;plot(f/f(fm),beta,'.-g');grid on
%     figure(2);hold on;plot(f/f(fm),beta,'.-r');grid on
    Bbeta = (Beta/2).*coth(pi.*Beta/2);
    figure(2);hold on;h4=plot(f/f(fm),Bbeta,'.-k');grid on;t4 = ['(\sigma_x/\sigma_y)^2'];

      set(gca,'fontsize',15)
    xlabel('f/f_p','fontsize',20)
    ylabel('A','fontsize',20)
    title(['spreading (A) of waves. Run ',int2str(run)],'fontsize',20)
    hleg = legend([h1,h4],t1,t4,'location','NorthEast');
    eval(['print -djpeg99 \SOFS/A_beta',int2str(run)]) 

    
    % Plot 8 panels of FNn from fm-1 to fm+7
    figure(3);clf;
    fe = fm + 7;if fe > 20;fe = 20;end
    n = 0;
    for j = fm -1:fe
        n=n+1;
        if n < 9;
        subplot(4,2,n);
%         h1 = plot([1:360],[FN3(:,j)],'.-r');grid on;t1 = ['MEM'];
%         hold on;h2 = plot([1:360],[FN2(:,j)],'.-g');grid on;t2 = ['MLM'];
        hold on;h3 = plot([1:360],[F1(:,j)],'.-b');grid on;t3 = ['WDM'];
        hold on;h4 = plot([1:360],[F4(:,j)],'.-k');grid on;t4 = ['WDM_{corr}'];
%     axis([0 360 0 1.5])
        xlabel('direction (from) [deg.]','fontsize',10)
    title(['f/f_p = ',num2str(f(j)/f(fm),2),',       \sigma_{x} / \sigma_{y} = ',num2str(sqrt(sxx(j)./syy(j)),2)],'fontsize',10)
   if n == 1 & dirw(fm-1)<=180;hleg = legend([h4,h3],t4,t3,'location','NorthEast');end
   if n == 1 & dirw(fm-1)>180;hleg = legend([h4,h3],t4,t3,'location','NorthWest');end
        end
    end
    
        eval(['print -djpeg99 \SOFS/dirsp',int2str(run)]) 
        %pause
    toc    

% Plot 9 panels of FNn from fm-1 to fm+7
    figure(4);clf;
    load betamod;% PP and mdif PPC mCD
    fe = fm + 7;if fe > 20;fe = 20;end
    n = 0;
    for j1 = fm -1:fe
        n=n+1;
        subplot(3,3,n);
        for j = 1:9
    dcd = mCD(j);
 if CD(j1) > dcd
     pp = PPC(:,j);
%      bet = polyval(pp,CD(j1));
     bet= 1.0./(interp1(SCD,invbeta,sqrt(CD(j1)),'linear','extrap'));

     dth = (j-1)*5;
% Construct bimodal spectrum from sum of 2 sech^2, separated by 2*dth [degrees].
F(:,1)=sech(bet*dr*([1:360]'-180-dth)).^2 + 0.0*rand(1,360)';
F(:,2)=sech(bet*dr*([1:360]'-180+dth)).^2 + 0.0*rand(1,360)';
F2 = (F(:,1)+F(:,2))/2;
figure(4);plot([1:360],F2,[color(j),'.-']);grid on;hold on;
 end
        end
         xlabel('direction (from) [deg.]','fontsize',10)
    title(['f/f_p = ',num2str(f(j1)/f(fm),2)],'fontsize',10)
    end
   
    
        eval(['print -djpeg99 \SOFS/dirsp_pos',int2str(run)]) 
        %pause
        
               
          % Plot 16 panels of FNn from fm-2 to 18
    figure(5);clf;
        eval(['load betaf\beta',int2str(run)]);

    load betamod;% PP and mdif PPC mCD
        load invcosnvsCD;% pbc
        pp2 = pbc;

    fe = fm + 7;if fe > 20;fe = 20;end
    n = 0;n2 = 1.0;
%     for j1 = fm -1:fe
    for j1 = fm-2:1:18
        n=n+1;
%         if n == 16;n2 = 1.2;end
        subplot(4,4,n);
           pp = PPC(:,1);
%      bet = polyval(pp,CD1(j1));
bet= 1.0./(interp1(SCD,invbeta,sqrt(CD1(j1)),'linear','extrap'));

  %   bet = 1.1*beta(j1);
     N = 1.2/polyval(pp2,CD1(j1));
% Construct bimodal spectrum from sum of 2 sech^2, separated by 2*dth [degrees].
Fs=sech(bet*dr*([1:360]'-180)).^2 ;
Fcos=cos(dr*([1:360]'-180)/2).^N ;

figure(5);plot([1:360],n2*FN1(1:360,j1),[color(3),'*-']);grid on;hold on;
 figure(5);plot([1:360],Fcos,[color(1),'.']);grid on;hold on;
figure(5);plot([1:360],Fs(1:360),[color(2),'.']);grid on;hold on;

 axis([0 360 0 1.5])
         xlabel('direction (from) [deg.]','fontsize',10)
    title(['f/f_p = ',num2str(f(j1)/f(fm),3),'   \beta = ',num2str(bet,3)],'fontsize',10)
    end
   
    
        eval(['print -djpeg99 \SOFS/dirsp_WDM',int2str(run)]) 
        %pause
        
              % Plot 16 panels of FNn from fm-2 to 18
    figure(5);clf;
        eval(['load betafmlm\beta',int2str(run)]);

    load betamod;% PP and mdif PPC mCD
        load invcosnvsCD;% pbc
        pp2 = pbc;

    fe = fm + 7;if fe > 20;fe = 20;end
    n = 0;
%     for j1 = fm -1:fe
    for j1 = fm-2:1:18
        n=n+1;
        subplot(4,4,n);
           pp = PPC(:,1);
%      bet = polyval(pp,CD1(j1));
bet = 1.0./(interp1(SCD,invbeta,sqrt(CD1(j1)),'linear','extrap'));

 %    bet = 1.1*beta(j1);
     N = 1.2/polyval(pp2,CD1(j1));
% Construct bimodal spectrum from sum of 2 sech^2, separated by 2*dth [degrees].
Fs=sech(bet*dr*([1:360]'-180)).^2 ;
Fcos=cos(dr*([1:360]'-180)/2).^N ;

figure(5);plot([1:360],FN2(1:360,j1),[color(6),'*-']);grid on;hold on;
% figure(5);plot([1:360],Fs(1:360),[color(2),'.-']);grid on;hold on;
% figure(5);plot([1:360],Fcos,[color(3),'.-']);grid on;hold on;
axis([0 360 0 1])
         xlabel('direction (from) [deg.]','fontsize',10)
    title(['f/f_p = ',num2str(f(j1)/f(fm),3),'   \beta = ',num2str(bet,3)],'fontsize',10)
    end
   
    
        eval(['print -djpeg99 \SOFS/dirsp_MLM',int2str(run)]) 
        %pause
        
              % Plot 16 panels of FNn from fm-2 to 18
    figure(5);clf;
        eval(['load betafmem\beta',int2str(run)]);

    load betamod;% PP and mdif PPC mCD
        load invcosnvsCD;% pbc
        pp2 = pbc;

    fe = fm + 7;if fe > 20;fe = 20;end
    n = 0;
%     for j1 = fm -1:fe
    for j1 = fm-2:1:18
        n=n+1;
        subplot(4,4,n);
           pp = PPC(:,1);
%      bet = polyval(pp,CD1(j1));
bet = 1.0./(interp1(SCD,invbeta,sqrt(CD1(j1)),'linear','extrap'));

%      bet = 1.1*beta(j1);
     N = 1.2/polyval(pp2,CD1(j1));
% Construct bimodal spectrum from sum of 2 sech^2, separated by 2*dth [degrees].
Fs=sech(bet*dr*([1:360]'-180)).^2 ;
Fcos=cos(dr*([1:360]'-180)/2).^N ;

figure(5);plot([1:360],FN3(1:360,j1),[color(6),'*-']);grid on;hold on;
% figure(5);plot([1:360],Fs(1:360),[color(2),'.-']);grid on;hold on;
% figure(5);plot([1:360],Fcos,[color(3),'.-']);grid on;hold on;
axis([0 360 0 1])
         xlabel('direction (from) [deg.]','fontsize',10)
    title(['f/f_p = ',num2str(f(j1)/f(fm),3),'   \beta = ',num2str(bet,3)],'fontsize',10)
    end
   
    
        eval(['print -djpeg99 \SOFS/dirsp_MEM',int2str(run)]) 
        %pause
        
                  
          % Plot 16 panels of FNn from fm-2 to 18
    figure(5);clf;
        eval(['load betaf\beta',int2str(run)]);

    load betamod;% PP and mdif PPC mCD
        load invcosnvsCD;% pbc
        pp2 = pbc;

    fe = fm + 12;if fe > 20;fe = 20;end

    n = 0;
%     for j1 = fm -1:fe
       fep = fm + 12;if fep > 20;fep = 20;end
    for j1 = fm-2:1:17
%     for j1 = fm:1:20
        n=n+1;
        subplot(4,4,n);
           pp = PPC(:,1);
       bet = bbeta(j1);
     N = 1.2/polyval(pp2,CD1(j1));


     Fs=sech(bbeta(j1)*dr*([1:360]'-180)).^2 ;
Fcos=cos(dr*([1:360]'-180)/2).^N ;

figure(5);plot([1:360],FN4(1:360,j1),[color(6),'*-']);grid on;hold on;
figure(5);plot([1:360],Fs(1:360),[color(5),'-'],'linewidth',2);grid on;hold on;
% figure(5);plot([1:360],Fcos,[color(3),'.-']);grid on;hold on;
axis([0 360 0 1.5])
         xlabel('direction (from) [deg.]','fontsize',10)
    title(['f/f_p = ',num2str(f(j1)/f(fm),3),'   \beta = ',num2str(bbeta(j1))],'fontsize',10)
    end
   
    
        eval(['print -djpeg99 \SOFS/dirsp_WDM_corr',int2str(run)]) 
        %pause
       
 % Plot crosswave/downwave ratio for different spreads CDn vs 9 frequencies from fm-1 to fm+7
    figure(7);clf;
           h1 = plot(f(fm-2:fe)/f(fm),1./CD3(fm-2:fe),'*-r');grid on;t1 = ['MEM'];hold on;
           h2 = plot(f(fm-2:fe)/f(fm),1./CD2(fm-2:fe),'*-g');grid on;t2 = ['MLM'];hold on;
           h3 = plot(f(fm-2:fe)/f(fm),1./CD1(fm-2:fe),'*-b');grid on;t3 = ['WDM'];hold on;
           h4 = plot(f(fm-2:fe)/f(fm),1./CD(fm-2:fe),'s-k');grid on;t4 = ['(\sigma_x/\sigma_y)^2'];hold on;
           h5 = plot(f(fm-2:fe)/f(fm),1./CDc(fm-2:fe),'*-k');grid on;t5 = ['WDM_{corr'];hold on;
      
    ylabel('downwave/crosswave variance ratio','fontsize',10)
    xlabel('f/f_p ','fontsize',10)
   hleg = legend([h5,h4,h3,h2,h1],t5,t4,t3,t2,t1,'location','NorthEast');
    
    
        eval(['print -djpeg99 \SOFS/CDratio',int2str(run)]) 
        %pause
end;%if length > 16384
toc
