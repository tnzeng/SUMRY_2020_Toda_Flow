%Generates a 51 x 51 matrix representation of the operator. The operator
%operations of the Toda Flow result in a polynomial of values of b_n and
%a_n in a small radius around any given value for n. As a result, this
%matrix perfectly computes difference equations for DSOs stationary for
%quartics of degree at most 10.

%The size of the finite matrix which is used as a stand-in for a Jacobi operator. Larger size is needed for higher degree polynomials.
size = 51;
%Creates a vector of indices centered at n
nminus = strcat('n-', string(floor(size/2):-1:1));
n = 'n';
nplus = strcat('n+', string(1:floor(size/2)));
nvec = [nminus n nplus];

b = str2sym(strcat('b_{', nvec, '}'));

%Only include 1 of the following 2 lines: the first for general Jacobi
%operator, the second for the DSO
a = str2sym(strcat('a_(', nvec(1,1:(size)), ')'));
%a = ones(size - 1, 1);


J = sym('J', [size size]);
for i = 1:size
    for j = 1:size
        if i == j
            J(i,j) = b(i);
        elseif abs(i - j) == 1
                J(i,j) = a(min(i,j));
        else
            J(i,j) = 0;
        end
    end
end
c = sym('c_', [1 5]);
J2 = J * J;
J3 = J2 * J;
J4 = J2 * J2;
%C = commutator(skew(J2), J)
%General Cubic
%C = commutator(skew(J3 + 0 * c(1) * J2 + c(2) * J), J);
%Monic Quartic
%C = commutator(skew(J4), J);
%Quartic + Quadratic
%C = commutator(skew(J4 + c(2) * J2), J);
%General Quartic
C = commutator(skew(J4 + c(1) * J3 + c(2) * J2 + c(3) * J), J);


%Outputs the main and off-diagonal entries of the Toda Flow
simplify(C((size + 1)/2, (size + 1)/2))
simplify(C((size + 1)/2, (size + 1)/2 + 1))
