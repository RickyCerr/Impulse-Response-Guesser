
sig_len = 1000; % Sets the amount of samples 

x = randn(1,sig_len); % Creats a random signal of zero-mean pseudo-random numbers from 1 to 1000

h = [-0.57534026E-02,0.99026691E-03,0.75733471E-02,-0.65141204E-02,0.13960509E-01,0.22951644E-02,-0.19994041E-01,0.71369656E-02,-0.39657373E-01,0.11260066E-01,0.66233635E-01,-0.10497202E-01,0.85136160E-01,-0.12024988E+00,-0.29678580E+00,0.30410913E+00,0.30410913E+00,-0.29678580E+00,-0.12024988E+00,0.85136160E-01,-0.10497202E-01,0.66233635E-01,0.11260066E-01,-0.39657373E-01,0.71369656E-02,-0.19994041E-01,0.22951644E-02,0.13960509E-01,-0.65141204E-02,0.75733471E-02,0.99026691E-03,-0.57534026E-02];
% our original h signal (sinc like function) to be convolved with our input
% signal and then later retrieved after noise is added

n1 = (randn(1,sig_len)/1000); % Creating two random zero-mean pseudo-random noise signals of length 1000, one is scaled by 0.001 in amplitude, the other is scaled by 1
n2 = randn(1,sig_len);

sig_end = length(h); % creating a variable length for the impulse response signal (this should be 32 in our case)

h_transpose = h.'; % transposing the matrix h, converting it from a 32x1 to a 1x32 (for matrix math)
    
X_window = zeros(1,sig_end); % creating the sliding window matrix the same size as the impulse response signal, this will be populated with values from the input signal
    
y = zeros(1, length(x) + sig_end - 1); % creating the vector that will store the convultion result
    
for i = 1:length(x)+sig_end-1 % perfrom the sliding window opperation, it will itterate through the entire input signal (1000) and the transient response (+ 32 - 1)
     if i <= length(x)
        X_window = [x(i), X_window(1:end-1)]; % Shift the window and update with the current input value
     end
     if i > length(x) % if the window is beyond the input signals length, populate the remaining window points with zeros
         X_window = [0, X_window(1:end-1)];
     end
    y(i) = sum(X_window .* h); % perform the dot product of the window and the impulse response every itteration and sum it up, store this total in one element of the vector y
end

y_condensed = y(1:1000); % delete the transient state


y1 = n1 + y_condensed; % add the low noise to the output
y2 = n2 + y_condensed; % add the high noise to the output (two different singals)

A = zeros(length(x), length(h)); % 1000 x 32 matrix of all zeros

for i = 1:length(x) % perform two sums, one from 1 to 1000 and the other from 1 to 32, effectively reaching the entire A matrix
    for j = 1:min(i, length(h))
        A(i, j) = x(i - j + 1); % populate A with the values of the input signal, each time shifting down a row and collumn
    end
end

h1_estimate = (A' * A) \ (A' * y1.'); % derive the reconsturcted h by using the least squares derivation method 
h2_estimate = (A' * A) \ (A' * y2.');

h_original_freq = real(fftshift(fft(h))); % take the real shifted values of the fourier transforms of all the h signals (original and reconstructed) for comparison
h1_freq = real(fftshift(fft(h1_estimate)));
h2_freq = real(fftshift(fft(h2_estimate)));

figure(1);

subplot(4,1,1);
plot(x);
title('x');
xlabel('Time');
ylabel('Amplitude');
subplot(4,1,2);
plot(h);
title('h');
xlabel('Time');
ylabel('Amplitude');
subplot(4,1,3);
plot(n1);
title('noise 1 (1/1000)');
xlabel('Time');
ylabel('Amplitude');
subplot(4,1,4);
plot(n2);
title('noise 2');
xlabel('Time');
ylabel('Amplitude');

figure(2);

subplot(4,1,1);
plot(y);
title('y');
xlabel('Time');
ylabel('Amplitude');
subplot(4,1,2);
plot(y_condensed);
title('y condensed');
xlabel('Time');
ylabel('Amplitude');
subplot(4,1,3);
plot(y1);
title('y1 (added noise 1)');
xlabel('Time');
ylabel('Amplitude');
subplot(4,1,4);
plot(y2);
title('y2 (added noise 2)');
xlabel('Time');
ylabel('Amplitude');

figure(3);

subplot(3,1,1);
plot(A);
title('x matrix');
xlabel('Time');
ylabel('Amplitude');
subplot(3,1,2);
plot(h1_estimate);
title('recovered h from y1');
xlabel('Time');
ylabel('Amplitude');
subplot(3,1,3);
plot(h2_estimate);
title('recovered h from y2');
xlabel('Time');
ylabel('Amplitude');

figure(4);

subplot(3,1,1);
stem(h_original_freq);
title('original h (freq)');
xlabel('w');
ylabel('Mag');
subplot(3,1,2);
stem(h1_freq);
title('h1 low noise (freq)');
xlabel('w');
ylabel('Mag');
subplot(3,1,3);
stem(h2_freq);
title('h2 higher noise (freq)');
xlabel('w');
ylabel('Mag');
