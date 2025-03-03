clear all;
close all;

% Define range of Eb/N0 in dB
EbN0_dB = -2:0.2:10;  % From -2 dB to 10 dB with 0.2 dB steps
EbN0 = 10.^(EbN0_dB/10);  % Convert to linear scale

% Number of Monte Carlo samples
L = 100000;

% Initialize capacity array
C = zeros(size(EbN0));
C_bsc = zeros(size(EbN0));
C_arb = zeros(size(EbN0));

% 
R = 1;
i = size(EbN0, 2);
while (i > 0)
    err = qfunc(sqrt(2*R*EbN0(i)));
    C_arb(i) = 0.5*log2(1 + 2*R*EbN0(i));
    er = C_arb(i) - R;
    % Ensure capacity stays larger than R
    if (C_arb(i) - 0.001 > R)
        R = max((C_arb(i) + R) / 2, 0);
    elseif (er < -0.0000001)
        R = max((C_arb(i) + R) / 2, 0);
    else
        R = 1;
        i = i - 1;
    end
end

R = 1;
i = size(EbN0, 2);
% Compute capacity for each Eb/N0 value
while (i > 0)
    % Convert Eb/N0 to sigma^2
    % Eb/N0 = 1/(2*R*sigma^2), we'll approximate R with C iteratively
    % Initial sigma^2 estimate
    sigma2 = 1/(2*EbN0(i)*R);  % Initial R approximation = 0.5
    
    % Generate samples
    x = sign(randn(L,1));  % Generate ±1 with equal probability
    z = sqrt(sigma2)*randn(L,1);  % AWGN noise
    y = x + z;  % Channel output
    
    % Compute p(y|x=+1) and p(y|x=-1)
    py_x1 = (1/sqrt(2*pi*sigma2)) * exp(-((y-1).^2)/(2*sigma2));
    py_xm1 = (1/sqrt(2*pi*sigma2)) * exp(-((y+1).^2)/(2*sigma2));
    
    py = 0.5 * (py_x1 + py_xm1);  % Compute p(y)
    
    % Monte Carlo estimation of E{-log2(p(y))}
    h_Y = -mean(log2(py));
    
    h_Z = 0.5 * log2(2*pi*exp(1)*sigma2); % Differential entropy of noise 
    
    C(i) = h_Y - h_Z;  % Capacity
    er = C(i) - R;
    % Ensure capacity stays larger than R
    if (C(i) - 0.0001 > R)
        R = max((C(i) + R) / 2, 0);
    elseif (er < -0.0000001)
        R = max((C(i) + R) / 2, 0);
    else
        i = i - 1;
        if(i > 0)
            R = C_arb(i);
        end
    end
    
end

err = size(EbN0);
R = 1;
i = size(EbN0, 2);
while (i > 0)
    err = qfunc(sqrt(2*R*EbN0(i)));
    C_bsc(i) = 1 + err*log2(err) +  (1 - err)*log2(1 - err);
    er = C_bsc(i) - R;
    % Ensure capacity stays larger than R
    if (C_bsc(i) - 0.001 > R)
        R = max((C_bsc(i) + R) / 2, 0);
    elseif (er < -0.0000001)
        R = max((C_bsc(i) + R) / 2, 0);
    else
        i = i - 1;
        if(i > 0)
            R = C(i);
        end
    end
end


% Plot results
figure;
plot(EbN0_dB,C,'LineWidth',2,'MarkerSize',8)
hold on
plot(EbN0_dB,C_bsc,'LineWidth',2,'MarkerSize',8)
plot(EbN0_dB,C_arb,'LineWidth',2,'MarkerSize',8)
grid on;
xlabel('E_b/N_0 (dB)');
ylabel('Capacity (bits/channel use)');
title('BI-AWGN Channel Capacity');
axis([-2 10 0 1]);
legend("Soft Decision", "Hard Decision", "Unconstrained");
