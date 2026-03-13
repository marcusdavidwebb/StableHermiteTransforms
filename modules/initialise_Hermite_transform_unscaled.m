function [d, Q] = initialise_Hermite_transform_unscaled(x)
    % Builds the orthogonal matrix Q and weight vector d 
    % without using scaling to avoid overflow
    % 
    % the coeffs2vals transform is d .* (Q * cfs) and
    % the val2coeffs transform is Q' * (vals ./ d).
    % x must be a column vector containing the Gauss-Hermite nodes
    
    N = length(x);
    Q = zeros(N);
    hjm1 = ones(N,1) * pi^(-1/4);    % h_0(x) (first Hermite polynomial)
    Q(:,1) = hjm1; 
    if N > 1
        hj = sqrt(2) * x .* hjm1;    % h_1(x) (second Hermite polynomial)
        Q(:,2) = hj;
    end

    % Loop to compute higher-order Hermite polynomials
    for j = 3:N
        % Use recurrence relation for Hermite polynomials
        [hjm1, hj] = deal(hj, sqrt(2/(j-1)) * x .* hj - sqrt((j-2)/(j-1)) * hjm1);

        % Assign the current Hermite polynomial to the matrix
        Q(:,j) = hj;  
    end

    % Compute the weights
    d = sqrt(N) * abs(Q(:,N)) .* exp(- x.^2 / 2);

    % Normalize Q matrix (convert to unscaled Hermite functions)
    for j = 1:N
        Q(:,j) = Q(:,j) ./ abs(Q(:,N)) .* exp(- 0.5 * log(N));
    end
end