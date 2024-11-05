
function c_cor = SOFT_DECODER_GROUPE1(c_ds_flip, H, P1_ds, MAX_ITER)
    numC = length(H(1,:));
    numF = length(H(:,1));

    Q1 = zeros(numC,1);
    Q0 = zeros(numC,1);
    q1 = zeros(numF, numC);
    q0 = zeros(numF, numC);
    for i = 1:length(q1)
        Q1(i) = P1_ds(i);
        Q0(i) = 1 - P1_ds(i);
        for j = 1:numF
            q1(j,i) = P1_ds(i);
            q0(j,i) = 1 - P1_ds(i);
        end
    end
    r0 = zeros(numF, numC);
    c_est = c_ds_flip;

    %% debug

    fprintf('Start:\n');
    fprintf('Q0: %s\n', mat2str(Q0));
    fprintf('Q1: %s\n', mat2str(Q1));
    fprintf('c_est: %s\n', mat2str(c_est));
    fprintf('Parity Check Fail: %d\n', parityCheckFail(c_est, H));

    %for i = 1:MAX_ITER
        % boucle principale
        %1 : q de c vers f
        %2 : calcul r
        %3 : r de f vers c
        %4 : update des q
    %end
    i = 1;
    while i <= MAX_ITER && parityCheckFail(c_est, H)
        
        % modif : chaque f -> un r par c (dont il est exclu)
        % update r
        for j = 1:numF
            for y = 1:numC
                q_temp = zeros(numC,1);
                for n = 1:numC
                    if H(j, n) == 1 && n ~= y
                        q_temp(n) = q1(j,n);
                    end
                end
                [r0(j, y),~] = cal_r(q_temp);
            end
        end

        % modif : chaque c -> un q par f (dont il est exclu)
        % update q
        for j = 1:numC
            for y = 1:numF
                r_temp = zeros(numF,1);
                for n = 1:numF
                    if H(y, n) == 1 && n ~= y
                        r_temp(n) = r0(n,j);
                    end
                end
                [q0(y,j), q1(y,j)] = cal_q(P1_ds(j), r_temp);
            end
        end

        % update Q
        for j = 1:numC
            r_temp = zeros(numF,1);
            for y = 1:numF
                if H(y, j) == 1
                    r_temp(y) = r0(y,j);
                end
            end
            [Q0(y,j), Q1(y,j)] = cal_q(P1_ds(j), r_temp);
        end
        
        c_est = estimate(Q0, Q1);

        %% debug
        fprintf('Iteration %d:\n', i);
        fprintf('Q0: %s\n', mat2str(Q0));
        fprintf('Q1: %s\n', mat2str(Q1));
        fprintf('c_est: %s\n', mat2str(c_est));
        fprintf('Parity Check Fail: %d\n', parityCheckFail(c_est, H));
        
        i = i + 1;
    end
    
    c_cor = logical(c_est);

end

function [q0, q1] = norm_q(q0i, q1i)
    norm = q0i + q1i;
    q0 = q0i / norm;
    q1 = q1i / norm;
end

function [r0, r1] = cal_r(q1)
    prod = 1;
    for i = 1:length(q1)
        if q1(i) ~= 0
            prod = prod*(1-2*q1(i));
        end
    end
    r0 = 1/2 + (1/2)*prod;
    r1 = 1 - r0;
end

function [q0, q1] = cal_q(p, r0)
    q0i = cal_q0(p, r0);
    q1i = cal_q1(p, r0);
    [q0, q1] = norm_q(q0i, q1i);
end

function q0 = cal_q0(p, r0)
    prod = 1;
    for i = 1:length(r0)
        if r0(i) ~= 0
            prod = prod*r0(i);
        end
    end
    q0 = (1 - p)*prod;
end

function q1 = cal_q1(p, r0)
    prod = 1;
    for i = 1:length(r0)
        if r0(i) ~= 0
            prod = prod*(1 - r0(i));
        end
    end
    q1 = p*prod;
end

function c_est = estimate(Q0, Q1)
    numC = length(Q0(:,1));
    c_temp = zeros(numC, 1);
    for i = 1:numC
        if Q1(i) > Q0(i)
            c_temp(i) = 1;
        else
            c_temp(i) = 0;
        end
    end
    c_est = c_temp;
end

function condition = parityCheckFail(c_est, H)

    M = mod(H*c_est, 2);
    sum = 0;
    for i = 1:length(M)
        sum = sum + M(i);
    end

    if sum == 0
        condition = 0;
    else
        condition = 1;
    end

end
