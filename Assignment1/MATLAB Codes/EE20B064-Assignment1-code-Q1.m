% EE4140 Assignment 1 Q1 [Author: Kishore Rajendran EE20B064]

% Defining a random BPSK input bitstream (N = 16)
ip_bits = [1, -1, 1, 1, 1, -1, 1, -1, -1, 1, -1, -1, 1, 1, 1, -1];
% Parts 1 and 2
figure(1)
plot_output(0, 4, 2, 16, ip_bits)
plot_output(0, 4, 4, 16, ip_bits)
legend('L=2', 'L=4')

% Part 3
figure(2)
plot_output(0, 8, 4, 16, ip_bits)
plot_output(0.5, 8, 4, 16, ip_bits)
plot_output(1, 8, 4, 16, ip_bits)
legend('\beta = 0', '\beta = 0.5', '\beta = 1');

% Part 4
figure(3)
rcpulse_PSD(0.5, 8, 4, 1024, 100)
rectpulse_PSD(8, 1024, 100)
legend('RC pulse-shape', 'Rect pulse-shape')

function plot_output(beta, J, L, N, ip_bits)
    % Interleaving the input bits with 'J-1' zeroes
    index = 1;
    ip_bitstream = zeros(1,J*N);
    for i = 1 :J: N*J
        ip_bitstream(i) = ip_bits(index);
        index = index + 1;
    end

    % Defining the RC filter
    rcfilter = zeros(1,2*J*L);
    index = 1;
    for i = -L*J :1: L*J
        rcfilter(index) = (sin(i*pi/J)/(i*pi/J)) * (cos(i*pi*beta/J)/(1 - (2*i*beta/J)^2));
        % To get rid of nan at x = 0
        rcfilter(L*J+1) = 1;
        % To get rid of inf when 2nd term is undefined
        rcfilter(isinf(rcfilter)) = beta*sin(pi/(2*beta))/2;
        index = index + 1;
    end
    
    % Output sequence = input convloved with filter's impulse response
    output = conv(ip_bitstream, rcfilter);
    op_tspan = -L*J : (-L*J+length(output)-1);
    plot(op_tspan, output, '-o', 'MarkerIndices', 2*J+1:J:length(output)-2*J-1)
    grid on
    title('Output sequence x(kTs)')
    xlabel('t')
    ylabel('Amplitude')
    hold on
end

function rcpulse_PSD(beta, J, L, N, R)
    output = zeros(R, 8256);
    out_fft_mag2 = zeros(R, 8192);
    for runs = 1: 1: R
        % Generating a random bitstream
        ip_bits = zeros(1, N);
        for k = 1 :1: N
            rv = unifrnd(0, 1);
            if rv < 0.5
                ip_bits(k) = -1;
            else
                ip_bits(k) = 1;
            end
        end
        % Interleaving the input bits with 'J-1' zeroes
        index = 1;
        ip_bitstream = zeros(1,J*N);
        for i = 1 :J: N*J
            ip_bitstream(i) = ip_bits(index);
            index = index + 1;
        end

        % Defining the RC filter
        rcfilter = zeros(1,2*J*L);
        index = 1;
        for i = -L*J :1: L*J
            rcfilter(index) = (sin(i*pi/J)/(i*pi/J)) * (cos(i*pi*beta/J)/(1 - (2*i*beta/J)^2));
            % To get rid of nan at x = 0
            rcfilter(L*J+1) = 1;
            % To get rid of inf when 2nd term is undefined
            rcfilter(isinf(rcfilter)) = beta*sin(pi/(2*beta))/2;
            index = index + 1;
        end
        
        % Output sequence = input convloved with filter's impulse response
        output(runs, :) = conv(ip_bitstream, rcfilter);
        % Only use one side because signal is real-valued
        fft_val = fft(output(runs, end-8192+1:end), 8192); % To take last 8192 samples of each op seq
        fft_val_s = fftshift(fft_val);
        out_fft_mag2(runs, :) = abs(fft_val_s).^2;
    end
    % Averaging the R FFT magnitudes
    rc_psd = mean(out_fft_mag2);
    % Dividing by R to get PSD estimate
    rc_psd = rc_psd/R;
    fshift = (-8192/2:8192/2-1)*(J/8192);
    plot(fshift, rc_psd)
    title('PSD Estimate of x(kT_{s})')
    xlabel('f')
    ylabel('S_{X}(f)')
    hold on
end

function rectpulse_PSD(J, N, R)
    output = zeros(R, 8199);
    out_fft_mag2 = zeros(R, 8192);
    for runs = 1: 1: R
        % Generating a random bitstream
        ip_bits = zeros(1, N);
        for k = 1 :1: N
            rv = unifrnd(0, 1);
            if rv < 0.5
                ip_bits(k) = -1;
            else
                ip_bits(k) = 1;
            end
        end
        % Interleaving the input bits with 'J-1' zeroes
        index = 1;
        ip_bitstream = zeros(1,J*N);
        for i = 1 :J: N*J
            ip_bitstream(i) = ip_bits(index);
            index = index + 1;
        end

        % Defining the Rectangular filter
        rectfilter = ones(1,J);
        
        % Output sequence = input convloved with filter's impulse response
        output(runs, :) = conv(ip_bitstream, rectfilter);
        % Only use one side because signal is real-valued
        fft_val = fft(output(runs, end-8192+1:end), 8192);
        fft_val_s = fftshift(fft_val);
        out_fft_mag2(runs, :) = abs(fft_val_s).^2;
    end
    % Averaging the R FFT magnitudes
    rect_psd = mean(out_fft_mag2);
    % Dividing by R to get PSD estimate
    rect_psd = rect_psd/R;
    fshift = (-8192/2:8192/2-1)*(J/8192);
    plot(fshift, rect_psd)
    title('PSD Estimate of x(kT_{s})')
    xlabel('f')
    ylabel('S_{X}(f)')
    hold on
end