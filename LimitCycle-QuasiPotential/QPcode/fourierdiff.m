function w = fourierdiff(v,m,interval,N)
%FOURIERDIFF   Fourier spectral differentiation of a periodic real function.
%   W = FOURIERDIFF(V,M,INTERVAL,N) returns the M-th derivative of function
%   V sampled at N (even) equidistant points on INTERVAL.
%   W = FOURIERDIFF(V,M,INTERVAL) returns the M-th derivative of function V
%   sampled at 128 equidistant points on INTERVAL.
%   W = FOURIERDIFF(V,M) returns the M-th derivative of function V sampled 
%   at 128 equidistant points on the default interval [0 2*pi].
%   W = FOURIERDIFF(V) returns the first derivative of function V sampled 
%   at 128 equidistant points on the default interval [0 2*pi].
%
%   M must be a non-negative integer and N must be even.
%   Note that the approximation is poor if V is not periodic.
%
%   See also   FFT

%   The algorithm is based on the book
%       Trefethen, L.N.: Spectral Methods in MATLAB, SIAM, 2000
%
%   Zoltn Csti
%   2015/04/13


% Check input
if nargin < 1
    error('MATLAB:fourierdiff:fewInputs', 'At least one input argument is needed.');
elseif nargin < 2
    m = 1;
    interval = [0 2*pi];
    N = 128;
elseif nargin < 3
    interval = [0 2*pi];
    N = 128;
elseif nargin < 4
    N = 128;
end
if nargin > 4
    error('MATLAB:fourierdiff:manyInputs', 'Too many input arguments.');
end
if mod(N,2) ~= 0 % N is odd
    error('MATLAB:fourierdiff:oddN', 'N must be even.');
end
if floor(m) ~= ceil(m) ... % m is not an integer
    || m < 0
    error('MATLAB:fourierdiff:notInteger', 'm must be a non-negative integer.');
end

v = v(:);
% Wavenumbers (with ordering suitable for MATLAB's FFT)
k = [0:N/2 -N/2+1:-1]';
% Discrete Fourier-transform of v via the FFT
v_hat = fft(v);
% Perform differentiation in the Fourier space
a = interval(1);
b = interval(2);
normalization = (2*pi/(b-a))^m; % mapping from [a b] to [0 2*pi]
w_hat = normalization*(1i*k).^m.*v_hat;
% For odd derivatives there is a loss of symmetry
if mod(m,2) == 1
     w_hat(N/2) = 0;
end
% Return to the physical space
w = real(ifft(w_hat));