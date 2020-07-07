function [sigma] = validateCovMatrix(sig)

% [sigma] = validateCovMatrix(sig)
%
% -- INPUT --
% sig:      sample covariance matrix
%
% -- OUTPUT --
% sigma:    positive-definite covariance matrix
%

EPS = 10^-0;
ZERO = 10^-1;

sigma = sig;
[r err] = cholcov(sigma, 0);

%if (err ~= 0)
    % the covariance matrix is not positive definite!
check = 0;
cnt = 0;
while check == 0
    cnt = cnt + 1;
    [v d] = eig(sigma);
    
    % set any of the eigenvalues that are <= 0 to some small positive value
    for n = 1:size(d,1)
        if (d(n, n) <= ZERO)
            d(n, n) = EPS;
        end
    end
%     if cnt > 10
%         d = d + 1e-5 * eye(size(d,1));
%     end
    
    % recompose the covariance matrix, now it should be positive definite.
    sigma = v*d*v';
    sigma = real(sigma);
    [r err] = cholcov(sigma, 0);
    if (err ~= 0)
        check = 0;
    else
        check = 1;
    end
end
%end