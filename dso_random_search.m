format long
%Length of the propagated sequence
length = 5;
%Number of sequences generated
num_trials = 10000000;
b = zeros(length,num_trials);
c = zeros(3);
%The coefficients of the polynomial for the Toda Flow, with c(1)
%corresponding to the second highest degree term, c(2) to the third, and so
%on.
c(1) = 0;
c(2) = -4;
c(3) = 0;
%Upper bound for when the program evaluates whether the sequence b has
%remained bounded
stat_up_check = 5;
%The radius of the interval centered at zero where the first 3 values of
%the seequence b are sampled
bound = 2;
a_check = zeros(length - 4, num_trials);
%Generate a random sequence and test its last element for having a low magnitude
b(1:3,:) = bound * (rand(3, num_trials) - 0.5);

%The non-b_{n + 2} terms in the General Quartic a_n difference equation
C = 2 * b(1,:) .* b(2,:) + b(1,:) * c(1) + 3 * b(2,:) * c(1) - 3 * b(3,:) * c(1) + b(2,:) * c(3) - b(3,:) * c(3) + b(2,:).^2 * c(2) + b(2,:).^3 * c(1) - b(3,:).^2 * c(2) - b(3,:).^3 * c(1) + b(1,:).^2 + 5 * b(2,:).^2 - 5 * b(3,:).^2 + b(2,:).^4 - b(3,:).^4;
%Generate the last initial condition using the a_n difference equation
b(4,:) = (-1/2) * (2 * b(3,:) + c(1) + sqrt((2 * b(3,:) + c(1)).^2 + 4 * C));
%Generate the rest of the sequence iteratively using the General Quartic b_n difference equation
for i = 5:length
    C = 2 * b(i - 4 + 1,:) .* b(i - 4 + 2,:) + b(i - 4 + 1,:) * c(1) + 3 * b(i - 4 + 2,:) * c(1) - 3 * b(i - 4 + 3,:) * c(1) + b(i - 4 + 2,:) * c(3) - b(i - 4 + 3,:) * c(3) + b(i - 4 + 2,:).^2 * c(2) + b(i - 4 + 2,:).^3 * c(1) - b(i - 4 + 3,:).^2 * c(2) - b(i - 4 + 3,:).^3 * c(1) + b(i - 4 + 1,:).^2 + 5 * b(i - 4 + 2,:).^2 - 5 * b(i - 4 + 3,:).^2 + b(i - 4 + 2,:).^4 - b(i - 4 + 3,:).^4;
    b(i,:) = (2 * b(i - 4,:) + 8 * b(i - 3,:) - 8 * b(i - 1,:) + 2 * b(i - 3,:) * c(2) - 2 * b(i - 1,:) * c(2) + 2 * b(i - 3,:) .* b(i - 2,:).^2 + 2 * b(i - 3,:).^2 .* b(i - 2,:) - 2 * b(i - 2,:) .* b(i - 1,:).^2 - 2 * b(i - 2,:).^2 .* b(i - 1,:) + 2 * b(i - 3,:).^2 * c(1) - 2 * b(i - 1,:).^2 * c(1) + 2 * b(i - 3,:).^3 - 2 * b(i - 1,:).^3 + 2 * b(i - 3,:) .* b(i - 2,:) * c(1) - 2 * b(i - 2,:) .* b(i - 1,:) * c(1)) / 2;
    a_check(i - 4, :) = (b(i - 4 + 4,:) == (-1/2) * (2 * b(i - 4 + 3,:) + c(1) + sqrt((2 * b(i - 4 + 3,:) + c(1)).^2 + 4 * C)) | b(i - 4 + 4,:) == (-1/2) * (2 * b(i - 4 + 3,:) + c(1) - sqrt((2 * b(i - 4 + 3,:) + c(1)).^2 + 4 * C))); 
end
is_stationary = (abs(b(length,:)) < stat_up_check & abs(b(length - 1,:)) < stat_up_check & abs(b(length - 2,:)) < stat_up_check & sum(a_check,1) == (length - 4));
stationary = b(1:3,is_stationary)
