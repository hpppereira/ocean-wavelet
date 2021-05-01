function x=tdelay(x,T_s,T_d,T_c);
%
% tdelay.m;   1995-03.  M. Donelan, CCIW.
%
% function x=tdelay(x,T_s,T_d,T_c);
%
%   This routine corrects for a time delay in the signal x.
%
%  x, the time series, must be an even number of points long.
%  T_s is sampling interval in sec.
%  T_d is the time delay in sec.
%  T_c is the cutoff period [sec] if low pass filtering desired [optional]


x=fft(x);
[n,m] = size(x);

w=2*pi*(1:n/2)*(1/T_s)/n;

% Compute phase shift corresponding to the time delay, T_d.
phdelay=T_d*w;
phdelay=[0 phdelay -phdelay(n/2-1:-1:1)]';
phdelay=phdelay*ones(1,m);
phase=angle(x)+phdelay;mag=abs(x);
x=mag.*cos(phase)+sqrt(-1)*mag.*sin(phase);

if nargin ==4
   ww = [0 w w(n/2-1:-1:1)]';
   i  = find(ww >= 2*pi/T_c); 
   x(i) = zeros(length(i),1);
end



x=real(ifft(x));
x(1,:) = x(2,:);
x(n,:) = x(n-1,:);
