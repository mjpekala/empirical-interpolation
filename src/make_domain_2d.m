function [Omega, x, y, idx] = make_domain_2d(varargin)
% MAKE_DOMAIN_2D  Creates a simple 2d domain (compact subset 
%                 of R^2) for interpolation.
%

shapes = {'square', 'triangle'};

p = inputParser;
addOptional(p, 'delta', .02, @(x) x > 0);
addParameter(p, 'shape', 'square', @(x) ismember(x, shapes));
parse(p, varargin{:});
p = p.Results;

[X,Y] = meshgrid(-1:p.delta:1, -1:p.delta:1);

x = X(:);  y = Y(:);  Omega = [x y];

switch(p.shape)
  case 'square'
    % pass
    idx = 1:length(x);
    
  case 'triangle'
    idx = find(sum((Omega+1)/2, 2) <= 1);
    Omega = Omega(idx,:);

  otherwise
    error('unsupported shape');
end    
