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


    %% MAIN DECODING LOOP

    i = 1;
    while i <= MAX_ITER && parityCheckFail(c_est, H)
        
        % Update r0 : responses from check-nodes to their variable-nodes
        for f = 1:numF
            for c = 1:numC
                if H(f, c) == 1
                    q_temp = [];
                    for j = 1:numC
                        if H(f, j) == 1 && j ~= c
                            q_temp = [q_temp, q1(f,j)];
                        end
                    end
                    [r0(f, c),~] = cal_r(q_temp);
                end
            end
        end

        % Update q0 and q1 : messages from variable-nodes to their check-nodes
        for c = 1:numC
            for f = 1:numF
                if H(f, c) == 1
                    r_temp = [];
                    for j = 1:numF
                        if H(j, c) == 1 && j ~= f
                            r_temp = [r_temp, r0(j,c)];
                        end
                    end
                    [q0(f,c), q1(f,c)] = cal_q(P1_ds(c), r_temp);
                end
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
    r0 = 1/2 + (1/2)*prod(1-2*q1);
    r1 = 1 - r0;
end

function [q0, q1] = cal_q(p, r0)
%Cal_Q - Compute q0 and q1, queries of variable nodes to check-nodes
%
%   Parameters :
%       r0 : response of the check-nodes about their amount of belief 
%            in "0"
%   Return :
%       [q0, q1] : varible-nodes' amount of belief in "0" and "1" sent to 
%       their check-nodes          
    q0i = (1 - p)*prod(r0);
    q1i = p*prod(1 - r0); % r1 = 1 - r0
    [q0, q1] = norm_q(q0i, q1i);
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

    if sum(M) == 0
        condition = 0;
    else
        condition = 1;
    end

end
