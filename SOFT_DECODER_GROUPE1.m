%% SOFT_DEC0DER_GROUPE1.m
% =========================================================================
% *Authors:* Raphaël DEMBIK, Camille ROBINSON, Pierre-Angelo PEYRIE
% *Date:* 2024, Novembre
% =========================================================================
% This script is an implementation of the LDPC soft decoder.
% It contains the function that implement the decoder and a number
% of sub-functions used in the main function.
% =========================================================================

function c_cor = SOFT_DECODER_GROUPE1(c_ds_flip, H, P1_ds, MAX_ITER)
% SOFT_DECODER_GROUPE1 - Soft-decode a codeword encoded with an LDPC code
%   
%   Parameters :
%       C_DS_FLIP  Flipped codeword to be decoded
%       H          Parity-check matrix
%       P1_DS      Initial probability of having "1" on variable-node i 
%                  knowing y(i)
%       MAX_ITER   Maximum number of iterations of the decoding process
%   Return :
%       C_COR      Decoded codeword

    %% INITIALISE VARIABLES 

    numC = length(H(1,:)); % number of variable-nodes
    numF = length(H(:,1)); % number of check-nodes

    
    Q1 = []; % variable-nodes' amount of belief in "1"
    Q0 = []; % variable-nodes' amount of belief in "0"

    % Messages matrix from variable-nodes i to check-nodes j
    q1 = zeros(numF, numC); % amount of belief in "1"
    q0 = zeros(numF, numC); % amount of belief in "0"

    for c = 1:numC
        Q1 = [Q1, P1_ds(c)];
        Q0 = [Q0, 1 - P1_ds(c)];
        for f = 1:numF
            if H(f, c) == 1
                q1(f,c) = P1_ds(c);
                q0(f,c) = 1 - P1_ds(c);
            end
        end
    end

    r0 = zeros(numF, numC); % check-nodes j responses to variable-nodes i
    
    c_est = c_ds_flip; % estimation of the decoded codeword

    %% debug

    %fprintf('Start:\n');
    %fprintf('Q0: %s\n', mat2str(Q0));
    %fprintf('Q1: %s\n', mat2str(Q1));
    %fprintf('q0: %s\n', mat2str(q0));
    %fprintf('q1: %s\n', mat2str(q1));
    %fprintf('r0: %s\n', mat2str(r0));
    %fprintf('c_est: %s\n', mat2str(c_est));
    %fprintf('Parity Check Fail: %d\n', parityCheckFail(c_est, H));

    %for i = 1:MAX_ITER
        % boucle principale
        %1 : q de c vers f
        %2 : calcul r
        %3 : r de f vers c
        %4 : update des q
    %end


    %% MAIN DECODING LOOP

    i = 1;
    while i <= MAX_ITER && parityCheckFail(c_est, H)
        
        % modif : chaque f -> un r par c (dont il est exclu)
        % Update r0 : responses from check-nodes to their variable-nodes
        for f = 1:numF
            for c = 1:numC
                q_temp = [];
                for j = 1:numC
                    if H(f, j) == 1 && j ~= c
                        q_temp = [q_temp, q1(f,j)];
                    end
                end
                [r0(f, c),~] = cal_r(q_temp);
            end
        end

        % modif : chaque c -> un q par f (dont il est exclu)
        % Update q0 and q1 : messages from variable-nodes to their check-nodes
        for c = 1:numC
            for f = 1:numF
                r_temp = [];
                for j = 1:numF
                    if H(j, c) == 1 && j ~= f
                        r_temp = [r_temp, r0(j,c)];
                    end
                end
                [q0(c,f), q1(c,f)] = cal_q(P1_ds(c), r_temp);
            end
        end

        % Update Q0 and Q1 : check-nodes' amount of belief in "0" and "1"
        for c = 1:numC
            r_temp = [];
            for f = 1:numF
                if H(f, c) == 1
                    r_temp = [r_temp, r0(f,c)];
                end
            end
            [Q0(c,1), Q1(c,1)] = cal_q(P1_ds(c), r_temp);
        end
        
        c_est = estimate(Q0, Q1); % estimation of the decoded codeword

        %% debug
        %fprintf('Iteration %d:\n', i);
        %fprintf('Q0: %s\n', mat2str(Q0));
        %fprintf('Q1: %s\n', mat2str(Q1));
        %fprintf('q0: %s\n', mat2str(q0));
        %fprintf('q1: %s\n', mat2str(q1));
        %fprintf('r0: %s\n', mat2str(r0));
        %fprintf('c_est: %s\n', mat2str(c_est));
        %fprintf('Parity Check Fail: %d\n', parityCheckFail(c_est, H));
        
        i = i + 1;
    end
    
    c_cor = c_est; % final estimation of the decoded codeword

end

function [q0, q1] = norm_q(q0i, q1i)
    norm = q0i + q1i;
    q0 = q0i / norm;
    q1 = q1i / norm;
end


function [r0, r1] = cal_r(q1)
% CAL_R - Compute r0 and r1, responses of check-nodes to variable nodes
%   
%   Parameters :
%       q1 : varible-nodes' amount of belief in "1" sent to their 
%            check-nodes
%   Return :
%       [r0, r1] : response of the check-nodes about their amount of belief 
%                  in "0" and "1"

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


function c_est = estimate(Q0, Q1)
% ESTIMATE - Compute the estimation of the decoded codeword as for
%            the current amount of belief in "0" and "1" of the variable-nodes
%
%   Parameters :
%       Q0 : amount of belief in "0" according to the check-nodes
%       Q1 : amount of belief in "1" according to the check-nodes
%   Return :
%       c_est : estimation of the decoded codeword

    numC = length(Q0(:,1));
    c_temp = zeros(numC, 1);
    for i = 1:numC
        if Q1(i) > Q0(i)
            c_temp(i) = 1;
        else
            c_temp(i) = 0;
        end
    end
    c_est = logical(c_temp);
end


function condition = parityCheckFail(c_est, H)
% parityCheckFail - Verify the parity equations represented by H for the
%                   estimated decoded codeword
%
%   Parameters :
%       H : parity-check matrix
%       c_est : estimation of the decoded codeword
%   Return :
%       condition : 1 if at least one of the equations fails, 0 otherwise

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
