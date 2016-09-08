function [Omega, x, y, idx] = make_domain_2d(delta, shape)
% MAKE_DOMAIN_2D  Creates a simple 2d domain (compact subset 
%                 of R^2) for interpolation.
%
%  Parameters:
%    delta : (scalar) spacing in each dimension (assumed isotropic)
%    shape : shape of the domain; see below for valid options.

% mjp, sept 2016

if nargin < 2, shape = 'square'; end
if nargin < 1, delta = 0.02; end

[X,Y] = meshgrid(-1:delta:1, -1:delta:1);

x = X(:);  y = Y(:);  Omega = [x y];

switch(shape)
  case 'square'
    % pass
    idx = 1:length(x);
    
  case 'triangle'
    idx = find(sum((Omega+1)/2, 2) <= 1);
    Omega = Omega(idx,:);

  otherwise
    error('unsupported shape');
end    
