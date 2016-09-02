% TEST_APPLY  Unit tests for the apply() function.

addpath('..');

f = {};
f{1} = @(X) X(:,1) + X(:,2);
f{2} = @(X) X(:,1).*X(:,2);
f{3} = @(X) X(:,1).^2 + X(:,2).^2;

X = [1 1 ; 2 2; 3 3];

Z = apply(f, X);

for ii = 1:length(f)
    assert(all(Z(:,ii) == f{ii}(X)));
end


