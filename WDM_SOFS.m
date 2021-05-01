% WDM.m  Wavelet Directional Method 
% Driver for producing WDM directional wave spectra plots
%           from a specified array of 3 or more wave staffs.
%
%  
% M.Donelan / Bergen 1994.11.26; Modified by A.Babanin and M.Donelan
% 2012.05 and earlier
clf
clear
d = 400;%depth of mooring in meters.
edepth = 0.05;% effective depth of drag ---> velocities.

np = 4;% Set number of wave staffs.
npp = (np*(np-1))/2;
tic
%  load data array with each wave staff time series as a column.
for run = [6]


    
    
   % To read SOFS TRI wave raw data.
%Largest waves on 30th Sep, 2013 are at 0600 hrs.
fid = fopen('Sep/201309030600.UVH','rt');
%fid = fopen('Sep/201309030800.UVH','rt');
% fid = fopen('Jul/201307300600.UVH','rt');
for j = 1:11
    fgetl(fid);
end
x = zeros(10000,4);aa = 0;
for j = 1:10000;
a = fgetl(fid);
if a ~= -1;
b = str2num(a);
x(j,:) = b(1,:);
end
if a == -1;aa = aa + 1;end 
if aa == 1;je = j-1;end
end
data = x(1:je,:); 
 
    
data=BINAVG(data,3);   % binavg;  Modify as necessary to subsample time series.

% Reduce data to a multiple of 1024 points: wavelet analysis done 1024
% point chunks; i.e. 2^10 points
tlen = fix(length(data)/1024);
   data=data(1:1024*tlen,:);

% Apply calibrations [to metres] to the wave staffs and set in desired order.
% Radius taken to be 0.5 m set up/down & left/right wave staffs: E,N,W,S.
R1 = 0.5;
WE = BuoyVelocityToSlope(data(:,4),0.14*3,d,edepth);% 400 m depth assumed. To the West
NS = BuoyVelocityToSlope(data(:,3),0.14*3,d,edepth);% 400 m depth assumed. To the North

ws=[data(:,2)-R1.*WE  data(:,2)-R1.*NS  data(:,2)+R1.*WE  data(:,2)+R1.*NS ];
heave = data(:,2);
  
      ws = detrend(ws);% Remove mean and trend.
      heave = detrend(heave);
  eval(['save wdms\data_',int2str(run),' ws WE NS heave'])

%   a and R are the angle [rad] and position [m] of the wave staffs 
%       from the array centre; i.e. Polar coordinates of the staffs   
%   np is number of probes
%
A = [90:90:360]; dr=pi/180;  A = dr*A; % Geographic coords: N = 0; E = 90.
R = sqrt(2).*[R1 R1 R1 R1] ;   % Radius of array (m). i.e. Polar coordinates of the staffs

X = R.*cos(A);
Y = R.*sin(A);

l=0;x=[];y=[];
for j=1:np-1
for k=(j+1):np
   l=l+1;
   x(l)=X(k)-X(j);
   y(l)=Y(k)-Y(j);
end
end

i=sqrt(-1);
r=abs(x+i*y);
a=atan2(y,x);

l=0;rr=[];
for j=1:npp-1
for k=(j+1):npp
   l=l+1;
   rr(l)=a(j)-a(k);
   csj(l)=cos(a(j));
   csk(l)=cos(a(k));
   snj(l)=sin(a(j));
   snk(l)=sin(a(k));
   rk(l)=r(k);
   rj(l)=r(j);
end
end
rr=rr*180/pi;
ii=find(rr<0);
rr(ii)=rr(ii)+360;
ii=find(rr>70 & rr<110);
jj=find(rr>250 & rr<290);
ij=sort([ii jj]);

rr
clear ii jj X Y x y rr a r
ns=1/0.14/3;n=1024;lf=.0625/2;hf=0.5;nv=4;% Set parameters required by Wavelet.m 

for i1=1:n:length(ws(:,1))-n+1
i1
  lp=log(lf)/log(2);lp=floor(lp);
  hp=log(hf)/log(2);hp=ceil(hp);

 
  AMP = [];
  for jj=1:np
     eval(['wx',int2str(jj),' = [];'])
     eval(['B',int2str(jj),' = [];'])
     eval(['AMP',int2str(jj),' = [];'])
     x=ws(i1:i1+n-1,jj);x=x(:);
     eval(['[wx',int2str(jj),',f]=WAVELET(x,lp,hp,nv,ns);'])
     eval(['B',int2str(jj),'=angle(wx',int2str(jj),');'])
     eval(['AMP',int2str(jj),'=abs(wx',int2str(jj),');'])
  end
  AMP = AMP1;
  for jj = 2:np
      eval(['AMP = AMP + AMP',int2str(jj),';'])
      eval(['clear AMP',int2str(jj),' wx',int2str(jj)])
  end
  AMP = AMP/np;
  clear AMP1 wx1 x
%Compute AMP from heave only to avoid pitch and roll errors.
  eval(['wx0 = [];'])
     eval(['AMP0 = [];'])
     x=heave(i1:i1+n-1);x=x(:);
     eval(['[wx0,f] = WAVELET(x,lp,hp,nv,ns);'])
     eval(['AMP0 = abs(wx0);'])
     AMP = AMP0;

mxmf=max(find(f < ns/2));
kkm=[];thm=[];kks=[];ths=[];
%********************************************
for mf=1:mxmf
%********************************************

l=0;b=[];
for j=1:np-1
for k=(j+1):np
   l=l+1;
   eval(['b(:,l)=B',int2str(k),'(:,mf)-B',int2str(j),'(:,mf);'])
end
end
ill=find(b(:) > pi); b(ill)=b(ill)-2*pi;
ill=find(b(:) < -pi); b(ill)=b(ill)+2*pi;
lb=length(b);
l=0;AA=[];
for j=1:npp-1
for k=(j+1):npp
   l=l+1;
   if length(find(ij == l)) ==0
      AA(:,l)=zeros(lb,1); bk(:,l)=zeros(lb,1);bj(:,l)=zeros(lb,1);
   else
      AA(:,l)=rk(l)/rj(l)*b(:,j)./b(:,k);
      bk(:,l)=b(:,k);bj(:,l)=b(:,j);
   end
end
end

th=[];ll=0;

for l=ij
   ll=ll+1;
   th(:,ll)=atan2((AA(:,l)*csk(l)-csj(l)),(snj(l)-AA(:,l)*snk(l)));
   kk(:,ll)=(snk(l)*bj(:,l)/rj(l)-snj(l)*bk(:,l)/rk(l)) ./ ...
           (csj(l)*snk(l) - csk(l)*snj(l))./cos(th(:,ll));
end
th = th+(kk<0)*pi;
kk = abs(kk);
thm = [thm MEANANG(th')'];
ths = [ths std(th')'];
kkm = [kkm mean(kk')'];
kks = [kks std(kk')'];

end        %%%  end of 'for mf=1:mxmf'  loop (frequency bin loop)
clear AA b bj bk ill th kk j k l

if i1 == 1 
   AAp=AMP;
   ddd=round(thm*180/pi);
   kkmp=kkm;
   eval(['save wdms/yw',int2str(run),' AAp ddd kkmp f d'])
else
  eval(['load wdms/yw',int2str(run)])
  AAp=[AAp;AMP];
  kkmp=[kkmp;kkm];
  ddd=[ddd;round(thm*180/pi)];
  eval(['save wdms/yw',int2str(run),' AAp ddd kkmp f d'])
end
clear AAp ddd kkmp f
end

clear B1
 for jj = 2:np
      eval(['clear B',int2str(jj)])
  end

clear csj csk dr hf hp kkm kks lb ll lf mf mxmf mw A R
clear n rj rk snj snk thm ths wn t i1 ij jj lp AMP


WAVEPLOTS
WAVENUM_LOG_SOFS
toc
end   % end of 'if length(data) > 4096' loop


