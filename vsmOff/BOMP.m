function xhat = BOMP(A,y,np)

[~,N] = size(A);

A_normalized = normalize(A,"norm",2);
% Construct blocks for the CS matrix [A,A^c] as [a_n, a^c_n]
Nblk = N / 2;
for itr_blk = 1:Nblk
    % D_normalized{itr_blk} = [A_normalized(:,itr_blk), A_normalized(:,N/2+itr_blk)];
    D{itr_blk} = [A(:,itr_blk), A(:,N/2+itr_blk)];
    Dindx{itr_blk} = [itr_blk, N/2+itr_blk];
end

AA = [];
ind_max = [];
y_orig=y;
M = length(y);
eps = M*np;
err(1) = (norm(y_orig))^2;



itr = 1;
if norm(y_orig) ~= 0

    while (err >= eps)&&(itr<=(M))

        % Support detection
        for itr_blk = 1 : Nblk
            tmp1(itr_blk,1) = norm((D{itr_blk})' * y);
        end

        tmp1_sorted = sort(tmp1,'descend');
        if norm(tmp1) <= (1e-10)*norm(y_orig)
            break;
        end

        zz = find(tmp1_sorted(1)==tmp1);
        ind_max(itr) = zz(1);



        if sum(ind_max==ind_max(itr))>1  % To avoid choosing the same index twice
            ind_max(itr) = [];
            break;
        end

        % ### To account for mirror components ###
        AA = [AA,D{ind_max(itr)}]; % To account for mirror components

        xhat1 = (pinv(AA))*y_orig;

        if max(abs(xhat1))>1e4
            flagg = 1;
        end


        r = y_orig - AA*xhat1;

        err = (norm(r))^2;    % Stopping criterion based on noise level (eps2)

        y=r;

        itr = itr + 1;

    end

    if isempty(ind_max)
        xhat = zeros(length(A(1,:)),1);
    else
        xhat = zeros(length(A(1,:)),1);
        for itr_blk = 1:length(ind_max)
            xhat(Dindx{ind_max(itr_blk)},1) = xhat1( ((itr_blk-1)*2 + 1):(itr_blk*2),1);
        end
    end

else

    xhat = zeros(length(A(1,:)),1);
end

end
