function W_n = make_polynomial_basis(d, max_degree)
% MAKE_POLYNOMIAL_BASIS  Returns functions that implement a
%                        polynomial basis in d dimensions.
%
%   d          : The number of spatial dimensions.
%   max_degree : The highest degree polynomial to include.
%

% mjp, sept 2016.

switch(d)
  case 1
    W_n = cell(max_degree+1,1);
    for ii = 0:max_degree
        W_n{ii+1} = @(x) x.^ii;
    end
    
  case 2
    % the code below is a way to enumerate the space:
    %   ii+jj <= n      for 0 <= n <= max_degree
    [ii,jj] = ind2sub([max_degree+1, max_degree+1], ... 
                      find(tril(ones(max_degree+1,max_degree+1))));
    jj = jj - 1;
    ii = max_degree+1-ii;

    % sort by degree (for aesthetics and so the average function appears first)
    s = sum([ii jj], 2);
    [~,idx] = sort(s, 'ascend');
    ii = ii(idx);  jj = jj(idx);

    W_n = cell(length(ii),1);
    for kk = 1:length(W_n)
        W_n{kk} = @(X2d)  X2d(:,1).^ii(kk) .* X2d(:,2).^jj(kk);
    end
    
  otherwise
    error('sorry, not yet implemented');
end

