% EE4140 Assignment 1 Q2 [Author: Kishore Rajendran EE20B064]

% snr varies from 0 to 14 dB in steps of 2dB
snr = 0 :2: 14;
stdevs = 10.^-(snr/20);

% dmin for the 3 signals s.t Eb = 1
d1 = 1;
d2 = 1;
d3 = sqrt(0.4);
d = [d1 d2 d3];
q = zeros(length(stdevs), 3);
Pe_t = zeros(length(stdevs), 3);

% Part 1
for k = 1:length(stdevs)
    sd = stdevs(k);
    % Calculating Pe (Theoretical value) 
    for j = 1:3
        q(k, j) = erfc(d(j)/(sqrt(2)*sd))/2;
        if j == 1       % 2-PAM
            Pe_t(k, j) = q(k, j);
        elseif j == 2   % 4-QAM
            Pe_t(k, j) = 2*q(k,j) - q(k, j)^2;
        elseif j == 3   % 16-QAM
            Pe_t(k, j) = 3*q(k, j) - 2.25*q(k, j)^2;
        end
    end
end

figure(1)
plot(snr, log10(Pe_t))
title('Theoretical log_{10}(P_{e}) vs SNR')
xlabel('SNR (dB)')
ylabel('log_10 (P_{e})')
legend('2-PAM', '4-QAM', '16-QAM')

% Part 2 (For 4-QAM)
Pe_ua = zeros(length(stdevs), 1);
Pe_unn = zeros(length(stdevs), 1);
Pe_cbound = zeros(length(stdevs), 1);

q1 = zeros(length(stdevs), 1);
q1_bound = zeros(length(stdevs), 1);
q2 = zeros(length(stdevs), 1);

for k = 1:length(stdevs)
    sd = stdevs(k);
    q1(k) = erfc(d2/(sqrt(2)*sd))/2;
    q1_bound(k) = exp(-(d2/(sqrt(2)*sd))^2);
    q2(k) = erfc(d2/sd)/2;
    Pe_ua(k) = 2*q1(k) + q2(k);
    Pe_unn(k) = 2*q1(k);
    Pe_cbound(k) = 2*q1_bound(k);
end

figure(2)
plot(snr, log10(Pe_ua))
hold on
plot(snr, log10(Pe_unn))
plot(snr, log10(Pe_cbound))
plot(snr, log10(Pe_t(:, 2)))
hold on
title('log_{10}(P_{e}) vs SNR for 4-QAM')
xlabel('SNR (dB)')
ylabel('log_10 (P_{e})')


% Part 3 (SER simulation for 4-QAM)
d4 = 1;
N = 1e5;
tx = zeros(length(stdevs), N);
rx = zeros(length(stdevs), N);
rx_d = zeros(length(stdevs), N);
errorcount = zeros(length(stdevs), 1);

Pe_ser4 = zeros(length(stdevs), 1);

for k = 1:length(stdevs)
    sd = stdevs(k);
    % Generating Tx symbols
    tx(k, :) = (2*randi([0, 1], N, 1)-1)*d4 + 1j*(2*randi([0, 1], N, 1)-1)*d4;
    % Generating Rx symbols (Adding noise)
    rx(k, :) = tx(k, :) + sd*(randn(size(tx(k,:))) + 1i*randn(size(tx(k,:))));
    % Decoding Rx symbols
    rx_d(k, :) = d4*(sign(real(rx(k,:)))) + 1i*d4*(sign(imag(rx(k,:))));
    % Comparing and counting errors
    mismatches = rx_d(k, :) ~= tx(k, :);
    errorcount(k) = sum(mismatches, 2);
    Pe_ser4(k) = errorcount(k)/N;
end

plot(snr, log10(Pe_ser4))
legend('All pairwise symbol errors', 'Only nearest neighbours', 'Chernoff Bound', 'Theoretical', 'SER')

% Part 4 (SER simulation for 16-QAM)
d16 = d3;
N = 1e5;
tx_16 = zeros(length(stdevs), N);
rx_16 = zeros(length(stdevs), N);
rx_d16 = zeros(length(stdevs), N);
errorcount16 = zeros(length(stdevs), 1);

Pe_ser16 = zeros(length(stdevs), 1);

for k = 1:length(stdevs)
    sd = stdevs(k);
    % Generating Tx symbols
    tx_16(k, :) = (2*randi([0, 2], N, 1)-1)*d16 + 1j*(2*randi([0, 2], N, 1)-1)*d16;
    % Generating Rx symbols (Adding noise)
    rx_16(k, :) = tx_16(k, :) + normrnd(0, sd, size(tx_16(k,:))) + 1i*normrnd(0, sd, size(tx_16(k,:)));
    % Decoding Rx symbols
    rx_d16(k, :) = d16*(sign(real(rx_16(k,:)))+2*(real(rx_16(k,:))>2*d16)+(-2)*(real(rx_16(k,:))<-2*d16)) + 1i*d16*(sign(imag(rx_16(k,:)))+2*(imag(rx_16(k,:))>2*d16)+(-2)*(imag(rx_16(k,:))<-2*d16));
    % Comparing and counting errors
    mismatches16 = rx_d16(k, :) ~= tx_16(k, :);
    errorcount16(k) = sum(mismatches16, 2);
    Pe_ser16(k) = errorcount16(k)/N;
end

figure(3)
plot(snr, log10(Pe_t(:, 3)))
hold on
plot(snr, log10(Pe_ser16))

% Part 5 (Union Bound for 16-QAM)
Pe_nn16 = zeros(length(stdevs), 1);
q_16 = zeros(length(stdevs), 1);

for k = 1:length(stdevs)
    sd = stdevs(k);
    q_16(k) = erfc(d3/(sqrt(2)*sd))/2;
    Pe_nn16(k) = 3*q_16(k);
end

plot(snr, log10(Pe_nn16))
legend('Theoretical', 'SER', 'Only nearest neighbours')
title('log_{10}(P_{e}) vs SNR for 16-QAM')
xlabel('SNR (dB)')
ylabel('log_10 (P_{e})')