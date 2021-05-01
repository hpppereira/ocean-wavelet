
function slope=BuoyVelocityToSlope(Vel,T_s,depth,edepth);
%
% VelocityToSlope.m;   2016-09.  M. Donelan, Duncan, BC.
%
% function slope=VelocityToSlope(Vel,T_s,depth);
%
%   This routine computes the slopes of the water surfave from the time series of horizontal velocity components
%
%  Vel, the velocity component time series, must be an even number of points long.
%  T_s is sampling interval in sec.
%  depth is that at buoy in m.
%  edepth = effective depth of buoy velocities in m.
dr = pi/180;% degrees to radians.
x=fft(Vel);
[n,m] = size(x);

w=2*pi*(1:n/2)*(1/T_s)/n;

% First convert buoy velocity to orbital velocity using transfer function
% derived through Morison's equation.
load buoy_orbital_velocity_transfer_fcn Uf URatio Uphase
Uphdelay = -dr*(interp1(2*pi*Uf(3:end),Uphase(3:end),w,'linear','extrap'));% Phase correction
Uphdelay=[0 Uphdelay -Uphdelay(n/2-1:-1:1)]';
Uphdelay=Uphdelay*ones(1,m);
w=[0 w w(n/2-1:-1:1)]';
w(w<0.1) = 0.1;

Beta = 1.0./(interp1(2*pi*Uf(3:end),URatio(3:end),w,'linear','extrap'));% Magnitude correction
phase=angle(x) + Uphdelay;mag=abs(x.*Beta); % Shift phase and modify magnitude.
v=mag.*cos(phase)+sqrt(-1)*mag.*sin(phase);


% x=real(ifft(v));
% x(1,:) = x(2,:);
% x(n,:) = x(n-1,:);

% x is now the orbital velocity (v). Convert to slope.
x = v;

w=2*pi*(1:n/2)*(1/T_s)/n;
w=[0 w w(n/2-1:-1:1)]';
w(w<0.1) = 0.1;

k = WAVEK(w/2/pi,depth);
c = w./k;
c(c<0.2) = 0.2;

v = exp(k.*edepth).*x./c.*tanh(k.*depth); % Converts velocity (adjusted to surface) to slope.
phdelay = -pi./2.*ones(1,n/2);
phdelay=[0 phdelay -phdelay(n/2-1:-1:1)]';
phdelay=phdelay*ones(1,m);
phase=angle(v) + phdelay;mag=abs(v); % Shift phase by 90 deg.
v=mag.*cos(phase)+sqrt(-1)*mag.*sin(phase);


x=real(ifft(v));
x(1,:) = x(2,:);
x(n,:) = x(n-1,:);
slope = x;
