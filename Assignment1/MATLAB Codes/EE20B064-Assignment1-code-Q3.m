% EE4140 Assignment 1 Q3 [Author: Kishore Rajendran EE20B064]

% snr varies from 0 to 16 dB in steps of 2dB
snr = 0 :2: 16;
% deltas (Decoding delays)
dels = [3 6 10 20 40];

% M-ary PAM with L tap filter
M = 4;
L = 3;
node_num = M^(L-1);
N = 100002;     % Number of symbols 
d = 1/sqrt(5);   
% 4-ary PAM symbols
I = zeros(1,N-2);
rv = rand(1, N-2);
I(rv < 0.25) = -3*d;
I(rv >= 0.25 & rv < 0.5) = -d;
I(rv >= 0.5 & rv < 0.75) = d;
I(rv >= 0.75) = 3*d;
% Last 'L-1' symbols are known (Tail symbols/ Flush bits)
I = [I 3*d 3*d];
% Array containing all possible symbols for M(=4)-ary PAM
all_d = [-3*d -d d 3*d];

% Defining the filter
f0 = 0.8/sqrt(2);
f1 = -1/sqrt(2);
f2 = 0.6/sqrt(2);
filter = [f0 f1 f2];

% for loop with dels
for delta = 1:length(dels)
    del = dels(delta);
    [nodes, r_cap] = precompute_rcaps(M, filter, all_d);
    ser_plotter(snr, I, M, nodes, all_d, filter, r_cap, N, del);
end
% Legend for the plots
legend('del=3', 'del=6', 'del=10', 'del=20', 'del=40')

function [node_vals, r_cap_vals] = precompute_rcaps(M, filter, all_d)
    % Generating nodes of the trellis and all the paths it can take
    L = length(filter);
    combs = combvec(all_d, all_d);
    nodes = combs;
    node_num = M^(L-1);
    % Repeatedly making all possible combinations to obtain all nodes of the trellis
    for i = 1:L-2
        combs = combvec(all_d, combs);
        if i == L-3
            nodes = combs;
        end
    end
    nodes = nodes';
    combs = combs';
    % Pre-computing r_caps (Since they remain the same throughout the trellis)
    r_cap = zeros('like', combs);
    for i=1:L
        r_cap = r_cap + filter(i)*combs(:,i);
    end
    r_cap_vals = reshape(r_cap, [], node_num)';
    node_vals = nodes;
end

function ser_plotter(snr, I, M, nodes, all_d, filter, r_cap, N, del)
    % Generating the list of std deviations for noise
    stdevs = 10.^-(snr/20);
    % Defining arrays to store results
    mlse_seq = zeros(length(stdevs), N);
    errorcount = zeros(length(stdevs), 1);
    Pe_ser = ones(length(stdevs), 1);
    % for loop with stdev
    for l = 1:length(stdevs)
        sd = stdevs(l);
        % Generating Gaussian noise
        v = normrnd(0, sd, size(I));
        % Calculating received sequence
        r_conv = conv(I, filter);
        r = r_conv(1:N) + v;
        % Calling the Viterbi decoder
        mlse_seq(l, :) = reshape(viterbi(M, nodes, all_d, filter, r, r_cap, N, del), 1, []);
        % Comparing the sequences and counting errors (over unknown symbols)
        L = length(filter);
        mismatches = mlse_seq(l, 1:end-(L-1)) ~= I(1, 1:end-(L-1));
        errorcount(l) = sum(mismatches, 2);
        Pe_ser(l) = errorcount(l)/(N-(L-1));
    end
    plot(snr, log10(Pe_ser))
    hold on
    title('log_{10}(P_{e}) vs SNR')
    xlabel('SNR (dB)')
    ylabel('log_10 (P_{e})')
end

function seq = viterbi(M, nodes, all_d, filter, r, r_cap, N, del)
% Viterbi Algorithm
    L = length(filter);
    node_num = M^(L-1);
    CM = zeros(node_num, 1);
    ss_final = zeros(1, N);
    % Pre-filling the first L-1 columns of survivor sequence before trellis reaches a constant size
    ss = fliplr(changem(nodes, 1:M, all_d));
    ss_new = ss;
    % Calculating initial CM
    for i=1:L-1
        rcapi = zeros(node_num, 1);
        for j=i:-1:1
            rcapi = rcapi + filter(j)*nodes(:, i-j+1);
        end
        tm = (r(i)-rcapi).^2;
        CM(:, 1) = CM(:, 1) + tm;
    end
    % Looping through the trellis
    for j=L:1:N
        if mod(j, 10000) == 0
            disp(j)
        end
        % Filling TM exiting each node and comp = tm + cm
        TM = (r(j) - r_cap).^2;
        comp = TM + CM;
        % Finding comp's entering a node
        enter_num = M^(L-2);
        mod_new_vals = repmat(1:M, 1, enter_num);
        a = 0:M-1;
        % comp_enter are all the branches entering a particular node
        comp_enter = comp(enter_num*(a)+1, :);
        for i=2:enter_num
            comp_enter = [comp_enter comp(enter_num*(a)+i, :)];
        end
        comp_enter = comp_enter';
        % Finding min metric for each node
        [tm_min, min_ind] = min(comp_enter, [], 2);
        % Updating survivor sequence for each node
        for k=1:node_num
            ss_new(k, 1:min(j,del+1)) = [ss((enter_num*(min_ind(k)-1)+ceil(k/M)), 1:min(j,del+1)-1) mod_new_vals(k)];
        end
        CM(:, 1) = tm_min;
        ss = ss_new;
        % Implementing decoding delay
        if j>del
            [~, min_buffi] = min(CM(:, 1));
            ss_final(1, j-del) = ss_new(min_buffi, 1);
            % Pruned ss
            ss = ss_new(:, 2:end);
        end
    end

    % Taking survivor sequence from last node since flush bits are 3*d, 3*d
    ss_flush = ss(node_num, :);
    % Combining it with the already estimated MLSE sequence with decoding delay
    ss_flush = [ss_final(1:N-del) ss_flush(end-del+1:end)];
    % To replace indices with the actual symbol values
    seq = changem(ss_flush, all_d, 1:M);
end