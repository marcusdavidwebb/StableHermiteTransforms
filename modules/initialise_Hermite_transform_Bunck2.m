function [d, Q] = initialise_Hermite_transform_Bunck2(N)
%INITIALISE_HERMITE_TRANSFORM_BUNCK2  Bunck recurrence with binary scaling.
% For positive integer N, vals = d .* (Q' * cfs) evaluates the expansion
% at the N Gauss-Hermite nodes, and cfs = Q * (vals ./ d) inverts it.
% See Bunck, "A fast algorithm for evaluation of normalized Hermite
% functions", BIT (2009), appendix.

x = hermpts(N);
degree = (1:N-1)';
alpha = sqrt(2 ./ degree);
beta = sqrt((degree-1) ./ degree);

% Split exp(-x^2/2) into a representable seed and blocks of 2^-512.
% Even at the outermost nodes the seed is at least pi^(-1/4) * 2^-512.
blockBits = 512;
blockFactor = 2^(-blockBits);
logBlock = blockBits * log(2);
halfSquare = 0.5 * x.^2;
blocks = floor(halfSquare / logBlock);
seed = pi^(-1/4) * exp(blocks * logBlock - halfSquare);

Q = zeros(N);
d = zeros(N,1);
values = zeros(N,1);
sqrtN = sqrt(N);
for k = 1:N
    node = x(k);
    remaining = blocks(k);
    scale = pow2(-blockBits * remaining);
    previous = 0;
    current = seed(k);
    values(1) = current * scale;

    for j = 1:N-1
        next = (alpha(j) * node) * current - beta(j) * previous;
        previous = current;
        current = next;
        if remaining > 0 && abs(current) >= 1
            % Apply one pending block to both terms of the recurrence.
            previous = previous * blockFactor;
            current = current * blockFactor;
            remaining = remaining - 1;
            scale = pow2(-blockBits * remaining);
        end
        % While blocks remain, |current| < 1. If scale rounds to zero,
        % so does the Hermite function, but the recurrence keeps its seed.
        values(j+1) = current * scale;
    end

    d(k) = sqrtN * abs(values(N));
    Q(:,k) = values / d(k);
end
end
