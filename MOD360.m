function x = mod360(x);% Changed to make mod360(360) = 360 not 0.
%
%  function x = mod360(x)
%
%    Given array of angles x [deg], returns x in range 1-360. i.e. mod(x,360).
%
x = mod(x,360);
x(x==0) = 360;

