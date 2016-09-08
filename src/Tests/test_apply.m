% TEST_APPLY  Unit tests for the apply() function.

% mjp, sept 2016

addpath('..');

% 1d case
f = {};
f{1} = @(x) ones(size(x));
f{2} = @(x) x;
f{3} = @(x) x.^2;

x = (-1:.01:1)';
Z = apply(f, x);
assert(all(Z(:,1) == 1));
assert(all(Z(:,2) == x));


% 2d case
f = {};
f{1} = @(X) X(:,1) + X(:,2);
f{2} = @(X) X(:,1).*X(:,2);
f{3} = @(X) X(:,1).^2 + X(:,2).^2;

X = [1 1 ; 2 2; 3 3];

Z = apply(f, X);
assert(size(Z,1) == size(X,1));
assert(size(Z,2) == length(f));

for ii = 1:length(f)
    assert(all(Z(:,ii) == f{ii}(X)));
end


