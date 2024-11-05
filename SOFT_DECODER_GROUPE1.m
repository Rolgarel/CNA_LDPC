
function c_cor = SOFT_DECODER_GROUPE1__(c_ds_flip, H, P1_ds, MAX_ITER)
    % initialisation
    %Q = P1_ds;
    %for i = 1:MAX_ITER
        % boucle principale
        %1 : q de c vers f
        %2 : calcul r
        %3 : r de f vers c
        %4 : update des q
    %end

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
        prod = prod*(1-2*q1(i));
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
        prod = prod*r0(i);
    end
    q0 = (1 - p)*prod;
end

function q1 = cal_q1(p, r0)
    prod = 1;
    for i = 1:length(r0)
        prod = prod*(1 - r0(i));
    end
    q1 = p*prod;
end

