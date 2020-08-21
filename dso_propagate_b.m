format long
length = 20;
b = zeros(length,1);

%First 3 values of the main diagonal sequence which will be propagated.
b(1) = 0;
b(2) = -1;
b(3) = -1;
c = zeros(3);

%The coefficients of the polynomial for the Toda Flow, with c(1)
%corresponding to the second highest degree term, c(2) to the third, and so
%on.
c(1) = 1;
c(2) = -5;
c(3) = -1;
%The non-b_{n + 2} terms in the General Quartic a_n difference equation
C = 2 * b(1) * b(2) + b(1) * c(1) + 3 * b(2) * c(1) - 3 * b(3) * c(1) + b(2) * c(3) - b(3) * c(3) + b(2)^2 * c(2) + b(2)^3 * c(1) - b(3)^2 * c(2) - b(3)^3 * c(1) + b(1)^2 + 5 * b(2)^2 - 5 * b(3)^2 + b(2)^4 - b(3)^4;

%Generate the last initial condition using the a_n difference equation
%The sqrt can either be added or subtracted.
b(4) = (-1/2) * (2 * b(3) + c(1) + sqrt((2 * b(3) + c(1))^2 + 4 * C));


%Generate the rest of the sequence iteratively using the General Quartic b_n difference equation
for i = 5:length
    b(i) = (2 * b(i - 4) + 8 * b(i - 3) - 8 * b(i - 1) + 2 * b(i - 3) * c(2) - 2 * b(i - 1) * c(2) + 2 * b(i - 3) * b(i - 2)^2 + 2 * b(i - 3)^2 * b(i - 2) - 2 * b(i - 2) * b(i - 1)^2 - 2 * b(i - 2)^2 * b(i - 1) + 2 * b(i - 3)^2 * c(1) - 2 * b(i - 1)^2 * c(1) + 2 * b(i - 3)^3 - 2 * b(i - 1)^3 + 2 * b(i - 3) * b(i - 2) * c(1) - 2 * b(i - 2) * b(i - 1) * c(1)) / 2;     
end

b