function w = blackhar(n)
%
%  BLACKHAR(N) returns the N-point Blackman-Harris window as a column vector
%  K.Kahma 1989-07-20

m = (0:n-1)' * ( 2*pi/(n-1) ) ;

w = (.35875 - .48829*cos(m) + .14128*cos(2*m) - 0.01168*cos(3*m) ) ;
