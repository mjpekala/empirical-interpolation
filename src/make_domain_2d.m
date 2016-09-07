function Omega = make_domain_2d(varargin)
% MAKE_DOMAIN_2D  Creates a simple 2d domain (compact subset 
%                 of R^2) for interpolation.
%

shapes = {'square', 'triangle'};

p = inputParser;
addOptional(p, 'delta', .02, @(x) x > 0);
addParameter(p, 'shape', 'square', @(x) ismember(x, shapes));
parse(p, varargin{:});
p = p.Results;


x = -1:p.delta:1;
[X,Y] = meshgrid(x, x);
Omega = [X(:) Y(:)];

switch(p.shape)
  case 'square'
    % pass
    
  case 'triangle'
    bits = sum( (Omega+1)/2, 2) > 1;
    Omega(bits,:) = [];

  otherwise
    error('unsupported shape');
end    
