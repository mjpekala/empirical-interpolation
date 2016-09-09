function h = mesh2(z_subset, domain_info)
% MESH2  Like mesh(), but supports non-rectangular domains.

% Fill in the subset of the domain for which we have z values.
z = NaN*ones(size(domain_info.x));
z(domain_info.idx) = z_subset;

% We assume the overall domain is square (usually [-1,1]^2)
n = sqrt(numel(domain_info.x));
to_square = @(v) reshape(v, n, n);

h = mesh(to_square(domain_info.x), to_square(domain_info.y), to_square(z));
