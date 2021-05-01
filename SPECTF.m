function S = spectf(x,y,dt,Nfa,a0)
%        S = SPECTf(x,dt,Nfa)
%        S = SPECTf(x,y,dt,Nfa)   cross spectrum
% 
% Frequency averaged power spectrum estimate,  GEOPHYSICAL NORMALIZATION
% Trend is removed, Blackman-Harris window is used. K.K.Kahma 1990-05-19
%
%     x , y  = data vectors
%     dt = sampling interval in seconds
%     Nfa = number of elementary frequency bands which are averaged
%
%     S(:,1) = f      (1/second == Hz)
%     S(:,2) = Sxx    (unit*unit*second)
%
% If cross spectrum is calculated 
%     S(:,3) = Syy
%     S(:,4) = Sxy
%     S(:,5) = phase angle = 180/pi*atan2(-imag(Sxy),real(Sxy))
%     S(:,6) = coherence   = abs(Sxy./sqrt(Sxx.*Syy))
%
%     positive phase means x leads y

%        S = SPECTF(x,y,dt,Nfa,a0)
% Elementary frequency bands 0:a0-1 (matlab index 1:a0) are ignored. 
% Default a0 is 0, i.e. all bands including zero (mean value) are incuded.
% 

x = x(:).';      	% Make sure x is a row vector
N  = max(size(x));      % Number of data points
window=blackhar(N).';

if max(size(y)) ~= N,
   if (max(size(y)) == 1) | (nargin < 5)

% ***************
   % Spectrum 
% ***************

   if (nargin < 4), Nfa = 0; end    % default a0
   if (nargin < 3), dt = 31; end    % default Nfa
   a0 = Nfa; Nfa = dt; dt = y; 

   Nfft=0; maxb=0; C=0; df=0;         % To define these variables before Xx
   Xx = fft(window.*detrend(x));
   Nfft = length(Xx);                 % Number of points in FFT
   maxb = Nfft/2+1;
   Xx(maxb+1:Nfft)=[];
   Xx(maxb) = Xx(maxb)/2;

   C = dt/(Nfa*pi*norm(window)^2);    % Scaling coefficient
   df = 2*pi/(dt*Nfft);

   if Nfa==1
      f = [a0:maxb-1]*df;
      Pxx = (abs(Xx(a0+1:maxb)).^2)*C;
   else
      if Nfa > 20
%       When Nfa is large enough this is as fast as vectorized
%       averaging and it requires far less memory    
        m=0; a=a0+1; b=a0+Nfa;
        while b <= maxb
           m=m+1;
           Pxx(m) = sum(abs(Xx(a:b)).^2)*C;
           f(m) = df*((a+b-2)/2);
           a=a+Nfa; b=b+Nfa;
        end
      else
        m=fix((maxb-a0) / Nfa);
        f=([1:m]*Nfa+(a0-0.5-Nfa/2))*df;
        b=a0+m*Nfa;

%        Old bin averaging loop
%        sx=zeros(m,Nfa);
%        for i=1:Nfa 
%           sx(:,i) = abs(Xx(a0+i:Nfa:b)).^2;
%        end
%        Pxx=(sum(sx.')*C);

        sx=zeros(Nfa,m);
        sx(:) = abs(Xx(a0+1:b)).^2;  
        Pxx=(sum(sx)*C);
      end
      a=a0+1+m*Nfa;
      if a <= maxb
         m=m+1;
         c = maxb+1-a; 
         Pxx(m) = sum(abs(Xx(a:maxb)).^2)*C*Nfa/c;
         f(m) = df*(a+maxb-2)/2;
      end
    end 
    clear Xx window
    S = [f/2/pi;2*pi*Pxx].';
 
  else

  error('x and y are not of same size'); end

else

% **********************
   % Cross spectrum
% **********************

   if (nargin < 5), a0 = 0; end    % default a0
   if (nargin < 4), Nfa = 31; end  % default Nfa

   y = y(:).';
   Nfft=0; maxb=0; C=0; df=0;
   Xx = fft(window.*detrend(x));
   Nfft = length(Xx);                 % Number of points in FFT
   maxb = Nfft/2+1;
   Xx(maxb+1:Nfft)=[];
   Xx(maxb) = Xx(maxb)/2;

   C = dt/(Nfa*pi*norm(window)^2);    % Scaling coefficient
   df = 2*pi/(dt*Nfft);

   Yy = fft(window.*detrend(y));
   Yy(maxb) = Yy(maxb)/2;
   Yy(maxb+1:Nfft)=[];

   if Nfa==1
      f = [a0:maxb-1]*df;
      Pxx = (abs(Xx(a0+1:maxb)).^2)*C;
      Pyy = (abs(Yy(a0+1:maxb)).^2)*C;
      Pxy = (conj(Xx(a0+1:maxb)).*Yy(a0+1:maxb))*C;
   else
      if Nfa > 20
         m=0; a=a0+1; b=a0+Nfa;
         while b <= maxb
            m=m+1;
            Pxx(m) = sum(abs(Xx(a:b)).^2)*C;
            Pyy(m) = sum(abs(Yy(a:b)).^2)*C;
            Pxy(m) = sum(conj(Xx(a:b)).*Yy(a:b))*C;
            f(m) = df*((a+b-2)/2);
            a=a+Nfa; b=b+Nfa;
         end
      else
         m=fix((maxb-a0) / Nfa);
         f=([1:m]*Nfa+(a0-0.5-Nfa/2))*df;
         b=a0+m*Nfa;
%         sx=zeros(m,Nfa);
%         for i=1:Nfa 
%            sx(:,i)  = abs(Xx(a0+i:Nfa:b)).^2;
%            sy(:,i)  = abs(Yy(a0+i:Nfa:b)).^2;
%            sxy(:,i) = conj(Xx(a0+i:Nfa:b)).*Yy(a0+i:Nfa:b);
%         end
         sx=zeros(Nfa,m);
         sx(:) = abs(Xx(a0+1:b)).^2;  
         Pxx=(sum(sx)*C);
         sx(:) = abs(Yy(a0+1:b)).^2;  
         Pyy=(sum(sx)*C);
         sx(:) = conj(Xx(a0+1:b)).*Yy(a0+1:b);  
         Pxy=(sum(sx)*C);
         a=a0+1+m*Nfa;
      end
   
      if a <= maxb
         m=m+1; 
         c = maxb+1-a; 
         Pxx(m) = sum(abs(Xx(a:maxb)).^2)*C*Nfa/c;
         Pyy(m) = sum(abs(Yy(a:maxb)).^2)*C*Nfa/c;
         Pxy(m) = sum(conj(Xx(a:maxb)).*Yy(a:maxb))*C*Nfa/c;
         f(m) = df*(a+maxb-2)/2;
      end
   end
   phase = 180/pi*atan2(-imag(Pxy),real(Pxy));
   coh   = abs(Pxy./sqrt(Pxx.*Pyy));
   clear Xx Yy window sx sy sxy
   S = [f/2/pi;2*pi*Pxx;2*pi*Pyy;2*pi*Pxy;phase;coh].';
end
