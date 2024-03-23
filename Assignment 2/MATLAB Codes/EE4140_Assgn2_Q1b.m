% EE4140 Assignment 2 Q1b [Author: Kishore Rajendran EE20B064]

% Part (b)
% snr = 10dB
snr = 10;
s_var = 1;
n_var = 10.^-(snr/10);

% M-ary PAM with L tap filter (channel)
M = 2;
L = 6;
d = 1;

% Defining the filter
C = sqrt(1^2 + 0.95^2 + 0.5^2 + 0.15^2 + 0.2^2 + 0.1^2);
f0 = 1/C;
f1 = -0.95/C;
f2 = 0.5/C;
f3 = 0.15/C;
f4 = -0.2/C;
f5 = -0.1/C;
filter = [f0 f1 f2 f3 f4 f5];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Part b1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% N tap equalizer with decoding delay del
N = 20;
del = 9;

% Generating Autocorrelation matrix (for L=6 unity gain filter, can be generalized) (NxN matrix)
t1 = filter(1)*filter(2) + filter(2)*filter(3) + filter(3)*filter(4) + filter(4)*filter(5) + filter(5)*filter(6);
t2 = filter(1)*filter(3) + filter(2)*filter(4) + filter(3)*filter(5) + filter(4)*filter(6);
t3 = filter(1)*filter(4) + filter(2)*filter(5) + filter(3)*filter(6);
t4 = filter(1)*filter(5) + filter(2)*filter(6);
t5 = filter(1)*filter(6);
R = (toeplitz([s_var+n_var, t1*s_var, t2*s_var, t3*s_var, t4*s_var, t5*s_var, zeros(1, N-L)])).';

% Generating Crosscorrelation matrix (for del) (Nx1 matrix)
p = [zeros(1, del-L+1), fliplr(filter)*s_var, zeros(1, N-del-1)].';
p = p(end-N+1:end);

% Finding optimal weights (inv(Ryy) * rIy = Ryy\rIy)
w_opt_b1 = R\p;

% Storing corresponding values of Jmin
Jmin_b1 = s_var - p'*w_opt_b1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Part b2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% N1 + N2 tap Decision Feedback Equalizer with decoding delay del
% N1 tap Feed-forward filter, N2 tap Feedback filter
N1 = 15;
N2 = 5;
del = 1;

% Generating matrix R11 (for L=6 unity gain filter, can be generalized) (N1xN1 matrix)
t1 = filter(1)*filter(2) + filter(2)*filter(3) + filter(3)*filter(4) + filter(4)*filter(5) + filter(5)*filter(6);
t2 = filter(1)*filter(3) + filter(2)*filter(4) + filter(3)*filter(5) + filter(4)*filter(6);
t3 = filter(1)*filter(4) + filter(2)*filter(5) + filter(3)*filter(6);
t4 = filter(1)*filter(5) + filter(2)*filter(6);
t5 = filter(1)*filter(6);
R_11 = (toeplitz([s_var+n_var, t1*s_var, t2*s_var, t3*s_var, t4*s_var, t5*s_var, zeros(1, N1-L)])).';

% Generating matrix R12
row = zeros(1,N2);
col = zeros(N1,1);
for i=1:N2
    if(del+i>=0 && del+i<L)
        row(i) = -s_var*filter(del+i+1);
    end
end
for i=1:N1
    if(del-i+2>=0 && del-i+2<L)
        col(i) = -s_var*filter(del-i+3);
    end
end
R_12 = toeplitz(col, row);

% Generating matrix R21 (Transpose of R12)
R_21 = R_12.';
% Generating matrix R22 (Identity Matrix * s_var)
R_22 = s_var*eye(N2);

% Generating Autocorrelation matrix R (N1+N2 x N1+N2 matrix)
R = [R_11 R_12; R_21 R_22];

% Generating Crosscorrelation matrix (for del) (N1+N2 x 1 matrix)
p = [zeros(1, del-L+1), fliplr(filter)*s_var, zeros(1, N1+N2-del-1)].';
p = p(end-(N1+N2)+1:end);

% Finding optimal weights (inv(Ryy) * rIy = Ryy\rIy)
w_opt_b2 = R\p;

% Storing corresponding values of Jmin
Jmin_b2 = s_var - p'*w_opt_b2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Part b3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% (i) SWF-LE

% num number of 4-PAM symbols
num = 50000;

% Best values as found in
N = 20;
del = 9;

% snr : 0dB to 20dB in steps of 2dB
snr = 0 :2: 20;
stdevs = 10.^-(snr/20);

% Generating 4-PAM Tx symbols
I = (2*randi([0, M-1], 1, num+N+L-2)-(M-1))*d;
constellation = (2*(0:M-1)-(M-1))*d;

% Passing Tx symbols through the filter
f_rx = conv(I, filter);
f_rx = f_rx(L:num+N+L-2);

Pe_ser_swf_le = zeros(1, length(stdevs));
decoded_swf = zeros(1, length(I));

for l = 1:length(stdevs)
    sd = stdevs(l);
    n_var = sd^2;

    % Adding noise to the received symbols
    y = f_rx + normrnd(0, sd, size(f_rx));

    % Equalized symbols
    eq_swf = conv(y, w_opt_b1.');
    eq_swf = eq_swf(N:num+N-1);
    % Nearest neighbour decoding
    for i = 1:num
        [minval,minidx] = min(abs(constellation-eq_swf(i)));
        decoded_swf(i) = constellation(minidx);
    end
    % Comparing and counting errors
    Pe_ser_swf_le(l) = sum(decoded_swf(1:num) ~= I(N+L-1-del:end-del))/num;
end

% (ii) SWF-DFE

% num number of 4-PAM symbols
num = 50000;

% Best values as found in
N1 = 15;
N2 = 5;
del = 1;

% Feedforward and Feedback weights
w_ff = w_opt_b2(1:N1);
w_fb = w_opt_b2(N1+1:end);

% snr : 0dB to 20dB in steps of 2dB
snr = 0 :2: 20;
stdevs = 10.^-(snr/20);

% Generating 4-PAM Tx symbols
I = (2*randi([0, M-1], 1, num+N1+L-2)-(M-1))*d;
constellation = (2*(0:M-1)-(M-1))*d;

% Passing Tx symbols through the filter
f_rx = conv(I, filter);
f_rx = f_rx(L:num+N1+L-2);

Pe_ser_swf_dfe = zeros(1, length(stdevs));
decoded_swf = zeros(1, length(I));

for l = 1:length(stdevs)
    sd = stdevs(l);
    n_var = sd^2;

    % Adding noise to the received symbols
    y = f_rx + normrnd(0, sd, size(f_rx));

    % Equalized symbols
    eq_swf1 = conv(y, w_ff.');
    eq_swf1 = eq_swf1(N1:num+N1-1);

    % The feedback decision that must be removed
    dec_fb = zeros(1,N2);

    % Nearest neighbour decoding
    for i = 1:num
        eq_swf2 = (w_fb.')*(dec_fb.');
        eq_swf = eq_swf1(i) - eq_swf2;
        [minval,minidx] = min(abs(constellation-eq_swf));
        decoded_swf(i) = constellation(minidx);

        dec_fb = [decoded_swf(i) dec_fb(1:end-1)];
    end
    % Comparing and counting errors
    Pe_ser_swf_dfe(l) = sum(decoded_swf(1:num) ~= I(N1+L-1-del:end-del))/num;
end

% All SER plots plotted in the following code

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Part b4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% Viterbi decoder
M = 2;          % M-PAM symbols
del = 30;       % Decoding delay
N = 50004;     % Number of symbols 
d = 1;          % d for the M-ary PAM symbols
% Generating the 2-ary PAM symbols
I = zeros(1,N-5);
rv = rand(1, N-5);
I(rv < 0.5) = -d;
I(rv >= 0.5) = d;
% Last 'L-1' symbols are known (Tail symbols/ Flush bits)
I = [I d d d d d];
r_conv = conv(I, filter);
noiseless_r = r_conv(1:N);
% Array containing all possible symbols for M(=2)-ary PAM
all_d = [-d d];

[nodes, r_cap] = precompute_rcaps(M, filter, all_d);
ser_plotter(snr, I, noiseless_r, M, nodes, all_d, filter, r_cap, N, del);
% Plotting SER plots for Equalizers LE and DFE too
plot(snr,log10(Pe_ser_swf_dfe));
plot(snr,log10(Pe_ser_swf_le));
legend('VA','DFE','LE');
title('log_{10}(P_{e}) vs SNR')
xlabel('SNR (dB)')
ylabel('log_{10}(P_{e})')

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

function ser_plotter(snr, I, noiseless_r, M, nodes, all_d, filter, r_cap, N, del)
    % Generating the list of std deviations for noise
    stdevs = 10.^-(snr/20);
    % Defining arrays to store results
    mlse_seq = zeros(length(stdevs), N);
    errorcount = zeros(length(stdevs), 1);
    Pe_ser = ones(length(stdevs), 1);
    % for loop with stdev
    for l = 1:length(stdevs)
        sd = stdevs(l);
        % Generating Gaussian noise and adding it to the noiseless received sequence
        v = normrnd(0, sd, size(I));
        r = noiseless_r + v;
        % Calling the Viterbi decoder
        mlse_seq(l, :) = reshape(viterbi(M, nodes, all_d, filter, r, r_cap, N, del), 1, []);
        % Comparing the sequences and counting errors (over unknown symbols)
        L = length(filter);
        mismatches = mlse_seq(l, 1:end-(L-1)) ~= I(1, 1:end-(L-1));
        errorcount(l) = sum(mismatches, 2);
        Pe_ser(l) = errorcount(l)/(N-(L-1));
    end
    plot(snr, log10(Pe_ser))
    disp(log10(Pe_ser))
    %xlim(0, 12)
    hold on
    title('log_{10}(P_{e}) vs SNR')
    xlabel('SNR (dB)')
    ylabel('log_{10}(P_{e})')
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

    % Taking survivor sequence from last node since flush bits are d,d,d,d
    ss_flush = ss(node_num, :);
    % Combining it with the already estimated MLSE sequence with decoding delay
    ss_flush = [ss_final(1:N-del) ss_flush(end-del+1:end)];
    % To replace indices with the actual symbol values
    seq = changem(ss_flush, all_d, 1:M);
end
