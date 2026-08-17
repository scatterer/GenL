function [result] = fastMatrixExponentiation_single(A, k)


    % Analytical attempt

    % [V, D] = eig(A); % V contains eigenvectors, D contains eigenvalues as a diagonal matrix
    % 
    % % Step 2: Identify the dominant eigenvalue and its eigenvector
    % [lambda_max, idx] = max(abs(diag(D))); % Find dominant eigenvalue
    % v1 = V(:, idx); % Dominant eigenvector
    % 
    % % Step 3: Normalize the eigenvector to stabilize computation
    % %v1 = v1 / v1(1); % Normalize such that the first component is 1 (optional)
    % 
    % % Step 4: Compute the stabilized ratio y21/y11
    % result_alt = v1(2) / v1(1); % Use components of the dominant eigenvector

    % [V, D] = eig(A);
    % d = min(k,1e2);
    % r = rem(k,d);
    % q = (k-r)/d;
    % 
    % result_alt = [1 0; 0 1];
    % invV = inv(V);
    % for i = 1:q
    %   result_alt = V * D^d * invV * result_alt;
    %   result_alt = result_alt/max(max(real(result_alt)));
    % end
    % result_alt = V * D^r * invV * result_alt;
    % result_alt = result_alt/max(max(real(result_alt)));

    %ratio_alt = result_alt(2,1)/result_alt(1,1);

   % end

    % Optimized Exponentiation by Squaring for 2x2 complex matrices
   
    % Base case for k = 0 (return identity matrix)
    if k == 0
        result = eye(2);
        return;
    end
    
    % Initialize result as identity matrix
    result = eye(2);
    % Fast exponentiation loop (iterative approach to avoid recursion)
    while k > 0
        if mod(k, 2) == 1
            % If k is odd, multiply result by current A
            result       = result * A;
        end
        % Square the matrix
        A  = A  * A;
        A  = A/max(max(real(A)));

        % Divide k by 2 (integer division)
        k = floor(k / 2);
    end
end