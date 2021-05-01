function [sireal,siimag] = diff1byF( si, dt, cut_off )
%        [sireal,siimag] = diff1byF( si, dt, cut_off )
%
% This function differentiates a signal in the time domain by an 
% equivalent product in the frequency program, and cut off energy
% of periods less than cut-off.
%
%   si      = input signal;
%   dt      = sampling rate [s];
%   cut_off = cut_off period [s].

% Frequency range estimation.
%
si = si(:) ; lsi = length(si) ; clsi = ceil( lsi/2 ) ;
f = ( -clsi : 1 : clsi-1 )' / lsi / dt ;
i = sqrt( -1 ) ;

% Signal transformation to the frequency domain.
%
SI = fft( detrend( si ) );
SIshift = fftshift( SI ) ;

% Differentiation.
%
SI(1:clsi,1)    = SIshift(1:clsi)     .* ( i * 2 * pi * f(1:clsi)     ) ;
SI(clsi+1,1)    = 0 ;
SI(clsi+2:lsi,1) = SIshift(clsi+2:lsi) .* ( i * 2 * pi * f(clsi+2:lsi) ) ;

% Energy elimination at frequencies greater than 1/cut-off.
%
fc = 1 / cut_off ; ind = find( f <= -fc | f >= fc ) ; 
if length(ind) == 1, SI(ind) = 0; else, SI(ind) = zeros(size(ind)) ; end

% Signal transformation back to the temporal domain.
%
si = ifft( fftshift( SI ) ) ;
sireal = real( si ) ;
siimag = imag( si ) ;
