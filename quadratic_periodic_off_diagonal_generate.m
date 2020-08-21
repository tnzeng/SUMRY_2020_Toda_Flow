format;
length = 15;
radius = 32;
interval = 0.015625;
vals = floor(radius / interval);
B = zeros(length, 2 * vals + 1);
for i = (-vals):vals
    B(1,i + vals + 1) = i * interval;
end
a_interval = 0.125;
a_bound = 2;
for j = 1:round(a_bound / a_interval)
    for k = 1:round(a_bound / a_interval)
        for l = 1:round(a_bound / a_interval)
            a = set_a([j * a_interval, k * a_interval, l * a_interval, 0.5], length);
            if (j == l && k == 0.5 / a_interval)
                continue;
            end
            A = a' * ones(1, 2 * vals + 1);
            C1 = a(1)^4 - 2 * a(2)^2 * a(1)^2;
            C2 = 2 * a(1)^4 * B(1,:) - 2 * a(2)^2 * a(1)^2 * B(1,:);
            C3 = a(1)^4 * B(1,:).^2  - a(2)^4 * (a(1)^2 - a(3)^2);
            B(2,:) = (-C2 - sqrt(C2.^2 - 4 * C1 * C3)) / (2 * C1);
            for i = 3:length
                B(i,:) = (a(i-2)^2 * (B(i - 1,:) + B(i - 2,:)) - a(i - 1)^2 * B(i - 1,:)) / (a(i - 1)^2);
            end
            A_check = zeros(length - 2, 2 * vals + 1);
            for i = 2:(length - 1)
                A_check(i - 1,:) = a(i - 1)^2 - a(i + 1)^2 + B(i,:).^2 - B(i + 1,:).^2;
            end
            
            A_check2 = sum((A_check == 0), 1);
            stationary = B(:,A_check2 == length - 2);
            if (size(stationary, 2) > 0)
                a
                stationary'
            end
        end
    end
end