function S = skew(J)
    %UNTITLED Summary of this function goes here
    %   Detailed explanation goes here
    S = triu(J) - tril(J);
end