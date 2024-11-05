
function c_cor = SOFT_DECODER_GROUPE1__(c_ds_flip, H, P1_ds, MAX_ITER)
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

    %for i = 1:MAX_ITER
        % boucle principale
        %1 : q de c vers f
        %2 : calcul r
        %3 : r de f vers c
        %4 : update des q
    %end
    i = 1;
    while i <= MAX_ITER && 1 %ajouter condition
        
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

        %condition ?
        % update Q
        for j = 1:numC
            r_temp = zeros(numF,1);
            for y = 1:numF
                if H(y, n) == 1
                    r_temp(n) = r0(y,j);
                end
            end
            [Q0(y,j), Q1(y,j)] = cal_q(P1_ds(j), r_temp);
        end

        i = i + 1;
    end

    % finalisation (décision)
    %c_temp = c_ds_flip;
    %for i = 1:length(c_temp)
        %if Q(i, 1) > Q(i, 2)
            %c_temp(i) = 1;
        %else
            %c_temp(i) = 0;
        %end
    %end
    
    %c_cor = c_temp;

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

disp('ok')