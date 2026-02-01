function xhat = SAOMP(A,y,np)

[~,N] = size(A);

AA = [];
ind_max = [];
y_orig=y;
M = length(y);
eps = M*np;
err(1) = (norm(y_orig))^2;

A_normalized = normalize(A,"norm",2);



itr = 1;
col_indx = [];
if norm(y_orig) ~= 0
    while (err >= eps)&&(itr<=(M))
        tmp1 = abs(A_normalized'*y);
        tmp1(col_indx) = 0;
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
        flag_mirror = 0;
        if ind_max(itr) > N/2
            ind_max(itr+1) = ind_max(itr) - N/2;
            itr = itr + 1;
            flag_mirror = 1;
        else
            ind_max(itr+1) = ind_max(itr) + N/2;
            itr = itr + 1;
            flag_mirror = 1;
        end
        % -----------------------
        if flag_mirror == 1 % To account for mirror components
            AA = [AA,A(:,ind_max((itr-1):itr))]; % To account for mirror components
        else
            AA = [AA,A(:,ind_max(itr))];
        end

        xhat1 = (pinv(AA))*y_orig;


        r = y_orig - AA*xhat1;

        err = (norm(r))^2;    % Stopping criterion based on noise level (eps2)

        y=r;

        itr = itr + 1;

    end

    if isempty(ind_max)
        xhat = zeros(length(A(1,:)),1);
    else
        xhat = zeros(length(A(1,:)),1);
        for itr1 = 1:length(ind_max)
            xhat(ind_max(itr1),1) = xhat1(itr1);
        end
    end

else
    xhat = zeros(length(A(1,:)),1);
end
supp_detect = ind_max;
end
