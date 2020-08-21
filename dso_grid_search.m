%Generates every possible length 3 sequence with elements from a lattice subset of
%the real line. Uses the quartic difference equations to propagate the
%sequence forward. Checks whether this sequence remains bounded (indicates
%a valid Jacobi Operator stationary for the quartic determined by the
%coefficient vector (c)

format long
%The final length of the main diagonal vector after the propagations.
length = 30;
%One less than the number of elements in the lattice subset of the real
%line
jumps = 8;
%The distance of the endpoints of the lattice from the origin.
radius = 2;
%The bound on the final few elements of the sequence that determines
%whether we accept the sequence as remaining bounded.
stat_up_check = 5;

a_check = zeros(length - 4, (jumps + 1)^3);
c = zeros(3);
%The coefficients of the quartic which we are checking if a DSO is
%stationary for. c(1) is the cubic term coefficent, c(2), is for the
%quadratic term, and c(3) is for the linear term.
c(1) = 1;
c(2) = -11;
c(3) = -1;

b = zeros(length, (jumps + 1)^3);
%Generates the first 3 terms in the sequence by taking every possible
%length 3 sequence using elements from the lattice subset of the real line.
for i = 1:(jumps + 1)^3
    b(1,i) = (mod(floor((i - 1) / (jumps + 1)^2), jumps + 1) - jumps / 2) * 2 * radius / jumps;
    b(2,i) = (mod(floor((i - 1) / (jumps + 1)), jumps + 1) - jumps / 2) * 2 * radius / jumps;
    b(3,i) = (mod(i - 1, jumps + 1) - jumps / 2) * 2 * radius / jumps;
end
b2 = b;
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
stationary = b(1:8,is_stationary);

%Add the stationary initial conditions for the negative of the square root

%The non-b_{n + 2} terms in the General Quartic a_n difference equation
C = 2 * b2(1,:) .* b2(2,:) + b2(1,:) * c(1) + 3 * b2(2,:) * c(1) - 3 * b2(3,:) * c(1) + b2(2,:) * c(3) - b2(3,:) * c(3) + b2(2,:).^2 * c(2) + b2(2,:).^3 * c(1) - b2(3,:).^2 * c(2) - b2(3,:).^3 * c(1) + b2(1,:).^2 + 5 * b2(2,:).^2 - 5 * b2(3,:).^2 + b2(2,:).^4 - b2(3,:).^4;
%Generate the last initial condition using the a_n difference equation
b2(4,:) = (-1/2) * (2 * b2(3,:) + c(1) - sqrt((2 * b2(3,:) + c(1)).^2 + 4 * C));
%Generate the rest of the sequence iteratively using the General Quartic b_n difference equation
for i = 5:length
    C = 2 * b2(i - 4 + 1,:) .* b2(i - 4 + 2,:) + b2(i - 4 + 1,:) * c(1) + 3 * b2(i - 4 + 2,:) * c(1) - 3 * b2(i - 4 + 3,:) * c(1) + b2(i - 4 + 2,:) * c(3) - b2(i - 4 + 3,:) * c(3) + b2(i - 4 + 2,:).^2 * c(2) + b2(i - 4 + 2,:).^3 * c(1) - b2(i - 4 + 3,:).^2 * c(2) - b2(i - 4 + 3,:).^3 * c(1) + b2(i - 4 + 1,:).^2 + 5 * b2(i - 4 + 2,:).^2 - 5 * b2(i - 4 + 3,:).^2 + b2(i - 4 + 2,:).^4 - b2(i - 4 + 3,:).^4;
    b2(i,:) = (2 * b2(i - 4,:) + 8 * b2(i - 3,:) - 8 * b2(i - 1,:) + 2 * b2(i - 3,:) * c(2) - 2 * b2(i - 1,:) * c(2) + 2 * b2(i - 3,:) .* b2(i - 2,:).^2 + 2 * b2(i - 3,:).^2 .* b2(i - 2,:) - 2 * b2(i - 2,:) .* b2(i - 1,:).^2 - 2 * b2(i - 2,:).^2 .* b2(i - 1,:) + 2 * b2(i - 3,:).^2 * c(1) - 2 * b2(i - 1,:).^2 * c(1) + 2 * b2(i - 3,:).^3 - 2 * b2(i - 1,:).^3 + 2 * b2(i - 3,:) .* b2(i - 2,:) * c(1) - 2 * b2(i - 2,:) .* b2(i - 1,:) * c(1)) / 2;
    a_check(i - 4, :) = (b2(i - 4 + 4,:) == (-1/2) * (2 * b2(i - 4 + 3,:) + c(1) + sqrt((2 * b2(i - 4 + 3,:) + c(1)).^2 + 4 * C)) | b2(i - 4 + 4,:) == (-1/2) * (2 * b2(i - 4 + 3,:) + c(1) - sqrt((2 * b2(i - 4 + 3,:) + c(1)).^2 + 4 * C))); 
end
is_stationary = (abs(b2(length,:)) < stat_up_check & abs(b2(length - 1,:)) < stat_up_check & abs(b2(length - 2,:)) < stat_up_check & sum(a_check,1) == (length - 4));
stationary = [stationary b2(1:8,is_stationary)];
%Optional lines that take out constant sequences
not_constant = stationary(1,:) ~= stationary(2,:) | stationary(2,:) ~= stationary(3,:);
stationary = stationary(1:8, not_constant)