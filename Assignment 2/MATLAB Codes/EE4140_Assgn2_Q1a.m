% EE4140 Assignment 2 Q1a [Author: Kishore Rajendran EE20B064]

% Part (a)
% snr = 10dB
snr = 10;
s_var = 1;
n_var = 10.^-(snr/10);

% M-ary PAM with L tap filter (channel)
M = 4;
L = 3;
d = 1/sqrt(5);

% Defining the filter
f0 = 0.8/sqrt(2);
f1 = - 1/sqrt(2);
f2 = 0.6/sqrt(2);
filter = [f0 f1 f2];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Part a1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% N tap equalizer with decoding delay del
N = 3;
del = 0;

% Generating Autocorrelation matrix (for L=3 unity gain filter, can be generalized) (NxN matrix)
t1 = filter(1)*filter(2) + filter(2)*filter(3);
t2 = filter(1)*filter(3);
R = (toeplitz([s_var+n_var, t1*s_var, t2*s_var, zeros(1, N-L)])).';

% Generating Crosscorrelation matrix (for del) (Nx1 matrix)
p = [zeros(1, del-L+1), fliplr(filter)*s_var, zeros(1, N-del-1)].';
p = p(end-N+1:end);

% Finding optimal weights (inv(Ryy) * rIy = Ryy\rIy)
w_opt_a1 = R\p;

% Storing corresponding values of Jmin
Jmin_a1 = s_var - p'*w_opt_a1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Part a2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% N tap equalizer with decoding delay del
N = 10;
del = 0;

% Generating Autocorrelation matrix (for L=3 unity gain filter, can be generalized) (NxN matrix)
t1 = filter(1)*filter(2) + filter(2)*filter(3);
t2 = filter(1)*filter(3);
R = (toeplitz([s_var+n_var, t1*s_var, t2*s_var, zeros(1, N-L)])).';

% Generating Crosscorrelation matrix (for del) (Nx1 matrix)
p = [zeros(1, del-L+1), fliplr(filter)*s_var, zeros(1, N-del-1)].';
p = p(end-N+1:end);

% Finding optimal weights (inv(Ryy) * rIy = Ryy\rIy)
w_opt_a2 = R\p;

% Storing corresponding values of Jmin
Jmin_a2 = s_var - p'*w_opt_a2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Part a3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% N tap equalizer with decoding delay del
N = 10;
del = 5;

% Generating Autocorrelation matrix (for L=3 unity gain filter, can be generalized) (NxN matrix)
t1 = filter(1)*filter(2) + filter(2)*filter(3);
t2 = filter(1)*filter(3);
R = (toeplitz([s_var+n_var, t1*s_var, t2*s_var, zeros(1, N-L)])).';

% Generating Crosscorrelation matrix (for del) (Nx1 matrix)
p = [zeros(1, del-L+1), fliplr(filter)*s_var, zeros(1, N-del-1)].';
p = p(end-N+1:end);

% Finding optimal weights (inv(Ryy) * rIy = Ryy\rIy)
w_opt_a3 = R\p;

% Storing corresponding values of Jmin
Jmin_a3 = s_var - p'*w_opt_a3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Part a4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% We know from class that for min value of Jmin, 
% Both N and Delta hit a minimum then Jmin again starts increasing

% N tap equalizer with decoding delay del
N_best = 20;
del_best = 9;

% Generating Autocorrelation matrix (for L=3 unity gain filter, can be generalized) (NxN matrix)
t1 = filter(1)*filter(2) + filter(2)*filter(3);
t2 = filter(1)*filter(3);
R = (toeplitz([s_var+n_var, t1*s_var, t2*s_var, zeros(1, N_best-L)])).';

% Generating Crosscorrelation matrix (for del) (Nx1 matrix)
p = [zeros(1, del_best-L+1), fliplr(filter)*s_var, zeros(1, N_best-del_best-1)].';
p = p(end-N_best+1:end);

% Finding optimal weights (inv(Ryy) * rIy = Ryy\rIy)
w_opt_a4 = R\p;

% Storing corresponding values of Jmin
Jmin_a4 = s_var - p'*w_opt_a4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Part a6 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% N tap equalizer with decoding delay del
N = 10;
del = 5;

% P number of snapshots
pvals = [20 100 500];
Jmin_vals_a6 = zeros(1, length(pvals));

% for loop with different values of P
for i = 1:length(pvals)
    P = pvals(i);
    
    % Generating 'P+N+L-2' 4-PAM Tx symbols
    I = (2*randi([0, M-1], 1, P+N+L-2)-(M-1))*d;
    
    % Passing Tx symbols through the filter
    f_rx = conv(I, filter);
    f_rx = f_rx(L:P+N+L-2);
    
    % Adding noise to the received symbols
    y = f_rx + (normrnd(0, sqrt(n_var), size(f_rx)));
    dec = I(N+L-1-del:end-del);
    R = zeros(N,N);
    p = zeros(N,1);

    % Using Time Averaging to estimate R and p
    for k = 1:P
        y_k = fliplr(y(k:k+N-1)).';
        R = R + (y_k*y_k.')/P;
        p = p + (y_k*dec(k))/P;
    end

    % Finding optimal weights (inv(Ryy) * rIy = Ryy\rIy)
    w_opt_a6 = R\p;

    % Storing corresponding values of Jmin
    Jmin_vals_a6(i) = s_var - p'*w_opt_a6;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Part a5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% num number of 4-PAM symbols
num = 50000;

% Best values as found in (a4)
N = N_best;
del = del_best;

% snr : 0dB to 20dB in steps of 2dB
snr = 0 :2: 20;
stdevs = 10.^-(snr/20);

% Generating 4-PAM Tx symbols
I = (2*randi([0, M-1], 1, num+N+L-2)-(M-1))*d;
constellation = (2*(0:M-1)-(M-1))*d;

% Passing Tx symbols through the filter
f_rx = conv(I, filter);
f_rx = f_rx(L:num+N+L-2);

Pe_ser_swf = zeros(1, length(stdevs));
eq_swf = zeros(1, length(I));
decoded_swf = zeros(1, length(I));

for l = 1:length(stdevs)
    sd = stdevs(l);

    % Adding noise to the received symbols
    y = f_rx + normrnd(0, sd, size(f_rx));

    % Equalized symbols
    eq_swf = conv(y, w_opt_a4.');
    eq_swf = eq_swf(N:num+N-1);
    % Nearest neighbour decoding
    for i = 1:num
        [minval,minidx] = min(abs(constellation-eq_swf(i)));
        decoded_swf(i) = constellation(minidx);
    end
    % Comparing and counting errors
    Pe_ser_swf(l) = sum(decoded_swf(1:num) ~= I(N+L-1-del:end-del))/num;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Part a7 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% Best values as found in (a4)
N = N_best;
del = del_best;

% Calculating wopt for SNR = 10dB
snr = 10;
n_var = 10.^-(snr/10);

% Generating optimal weights using TAWF (P = 500)
P = 500;
    
% Generating 'P+N+L-2' 4-PAM Tx symbols
I = (2*randi([0, M-1], 1, P+N+L-2)-(M-1))*d;

% Passing Tx symbols through the filter
f_rx = conv(I, filter);
f_rx = f_rx(L:P+N+L-2);

% Adding noise to the received symbols
y = f_rx + (normrnd(0, sqrt(n_var), size(f_rx)));
dec = I(N+L-1-del:end-del);
R = zeros(N,N);
p = zeros(N,1);

% Using Time Averaging to estimate R and p
for k = 1:P
    y_k = fliplr(y(k:k+N-1)).';
    R = R + (y_k*y_k.')/P;
    p = p + (y_k*dec(k))/P;
end

% Finding optimal weights (inv(Ryy) * rIy = Ryy\rIy)
w_opt_a7 = R\p;

% num number of 4-PAM symbols
num = 50000;

% snr : 0dB to 20dB in steps of 2dB
snr = 0 :2: 20;
stdevs = 10.^-(snr/20);

% Generating 4-PAM Tx symbols
I = (2*randi([0, M-1], 1, num+N+L-2)-(M-1))*d;
constellation = (2*(0:M-1)-(M-1))*d;

% Passing Tx symbols through the filter
f_rx = conv(I, filter);
f_rx = f_rx(L:num+N+L-2);

Pe_ser_tawf = zeros(1, length(stdevs));
eq_tawf = zeros(1, length(I));
decoded_tawf = zeros(1, length(I));

for l = 1:length(stdevs)
    sd = stdevs(l);
    n_var = sd^2;

    % Adding noise to the received symbols
    y = f_rx + normrnd(0, sd, size(f_rx));

    % Equalized symbols
    eq_tawf = conv(y, w_opt_a7.');
    eq_tawf = eq_tawf(N:num+N-1);
    % Nearest neighbour decoding
    for i = 1:num
        [minval,minidx] = min(abs(constellation-eq_tawf(i)));
        decoded_tawf(i) = constellation(minidx);
    end
    % Comparing and counting errors
    Pe_ser_tawf(l) = sum(decoded_tawf(1:num) ~= I(N+L-1-del:end-del))/num;
end

% SER plots for (a5) and (a7)
plot(snr,log10(Pe_ser_swf));
hold on
plot(snr,log10(Pe_ser_tawf));
legend('SWF LE','TAWF LE');
title('log_{10}(P_{e}) vs SNR')
xlabel('SNR (dB)')
ylabel('log_{10}(P_{e})')
