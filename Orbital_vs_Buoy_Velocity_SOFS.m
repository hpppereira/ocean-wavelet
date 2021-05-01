% Calculates buoy x-axis velocity, UB for a given orbital velocity, U; by integrating acceleration = Drag Force/mass.
clear
close all
R = 0.5;% radius of buoy, m -- TriAxys in SOFS
M = 230;% Kgm mass of buoy
rho = 1025;% Kgm/m^3 -- density of sea water
dt = 0.42;%seconds
A = 0.5*pi*R^2;% Drag Area, m^2, projected bottom half of sphere.
V = 2/3*pi*R^3;% Half sphere volume for inertial calcs.
% load UvsUB
load wdms/data_6
[ u, w ] = ORBITAL( heave, 1/0.14/3, 400, heave-0.1 );
UB = 0;
ax = 0;
U = u;
Ub(1) = UB;
au = DIFF1BYF(U,0.42,0.3);% orbital accel for inertial force calc.
ax = au(1);
for j = 1:2048
    for k = 1:20
    ax = 1.5*(0.2*A*rho*abs(U(j)-UB)*(U(j)-UB) + 0.4*V*rho*(au(j)-ax))/M;
    end
    UB = UB + dt*ax;
    Ub(j) = UB;
end
sp = SPECTF(U(65:2048),Ub(65:2048),0.42,16);
Uf = sp(:,1);
SU = sp(:,2);
SUb = sp(:,3);
URatio = sqrt(SUb./SU);
Uphase = sp(:,5);
coh = sp(:,6);
save buoy_orbital_velocity_transfer_fcn Uf URatio Uphase
figure(1);clf;semilogx(Uf,URatio,'.-');grid on
pause
figure(2);clf;semilogx(Uf,Uphase,'.-');grid on
pause
figure(3);clf;semilogx(Uf,coh,'.-');grid on
