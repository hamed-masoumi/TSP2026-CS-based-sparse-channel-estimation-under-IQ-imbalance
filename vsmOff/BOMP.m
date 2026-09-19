function xhat = BOMP(A, y, np)
% BOMP  Block Orthogonal Matching Pursuit.
%
%   xhat = BOMP(A, y, np)
%
%   Performs block orthogonal matching pursuit using pairs of columns
%   [A(:,k), A(:,N/2+k)] as blocks.

    % Matrix dimensions
    [~, N] = size(A);
    Nblk = N / 2;

    % A_normalized = normalize(A, "norm", 2);

    % Construct blocks:
    %   D{k}     = [A(:,k), A(:,N/2+k)]
    %   Dindx{k} = [k, N/2+k]
    D = cell(Nblk, 1);
    Dindx = cell(Nblk, 1);

    for blockIdx = 1:Nblk
        D{blockIdx} = [A(:, blockIdx), A(:, Nblk + blockIdx)];
        Dindx{blockIdx} = [blockIdx, Nblk + blockIdx];
    end

    % Initialization
    AA = [];
    ind_max = [];
    y_orig = y;

    M = length(y);
    eps = M * np;
    err(1) = norm(y_orig)^2;

    itr = 1;

    % Handle non-zero input
    if norm(y_orig) ~= 0

        while (err >= eps) && (itr <= M)

            % -------------------------------------------------------------
            % Support detection
            % -------------------------------------------------------------
            tmp1 = zeros(Nblk, 1);

            for blockIdx = 1:Nblk
                tmp1(blockIdx) = norm(D{blockIdx}' * y);
            end

            tmp1_sorted = sort(tmp1, "descend");

            % Stop if all block correlations are negligible
            if norm(tmp1) <= (1e-10) * norm(y_orig)
                break;
            end

            % Select the block with maximum correlation
            maxCandidates = find(tmp1_sorted(1) == tmp1);
            ind_max(itr) = maxCandidates(1);

            % Avoid selecting the same block twice
            if sum(ind_max == ind_max(itr)) > 1
                ind_max(itr) = [];
                break;
            end

            % -------------------------------------------------------------
            % Update active matrix and estimate coefficients
            % -------------------------------------------------------------
            AA = [AA, D{ind_max(itr)}];

            xhat1 = pinv(AA) * y_orig;

            % Compute residual
            r = y_orig - AA * xhat1;

            % Stopping criterion
            err = norm(r)^2;

            % Use residual for the next support-detection iteration
            y = r;

            itr = itr + 1;
        end

        % -------------------------------------------------------------
        % Construct final output vector
        % -------------------------------------------------------------
        if isempty(ind_max)
            xhat = zeros(length(A(1, :)), 1);
        else
            xhat = zeros(length(A(1, :)), 1);

            for blockIdx = 1:length(ind_max)
                coeffIdx = (blockIdx - 1) * 2 + (1:2);

                xhat(Dindx{ind_max(blockIdx)}, 1) = xhat1(coeffIdx, 1);
            end
        end

    else
        % Zero input
        xhat = zeros(length(A(1, :)), 1);
    end

end
