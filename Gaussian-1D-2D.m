%% Part 1
%% Legendre Polynomial root finding

addpath('myFunctions')
syms x y z
syms t

digits(16);

l_1 = legendreP(1,x);
l_2 = legendreP(2,x);
l_3 = legendreP(3,x);
l_4 = legendreP(4,x);
l_5 = legendreP(5,x);

disp("Legendre Polynomials")
disp(l_1)
disp(l_2)
disp(l_3)
disp(l_4)
disp(l_5)
disp("")

roots_1 = vpa(solve(l_1),16);
roots_2 = vpa(solve(l_2),16);
roots_3 = vpa(solve(l_3),16);
roots_4 = vpa(solve(l_4),16);
roots_5 = vpa(solve(l_5),16);

legendrePolynomial = [string(l_1); string(l_2); string(l_3); string(l_4); string(l_5)];
first_root = [string(roots_1(1)); string(roots_2(1)); string(roots_3(1)); string(roots_4(1)); string(roots_5(1))];
second_root = ["NA";string(roots_2(2)); string(roots_3(2)); string(roots_4(2)); string(roots_5(2))];
third_root = ["NA"; "NA"; string(roots_3(3)); string(roots_4(3)); string(roots_5(3))];
fourth_root = ["NA"; "NA"; "NA"; string(roots_4(4)); string(roots_5(4))];
fifth_root = ["NA"; "NA"; "NA"; "NA"; string(roots_5(5))];

table(legendrePolynomial, first_root, second_root, third_root, fourth_root, fifth_root)


%% lookup table
function [Y_lookup, Y] = Table(X, Y, X_lookup)
    n = size(X, 2);

    for j = 2:1:n
        for i = n:-1:j
            Y(i) = (Y(i) - Y(i - 1)) / X(i) - X(i - j + 1));
        end
    end

    Y_lookup = Y(n);

    for j = n-1:-1:1
        Y_lookup = Y_lookup.*(X_lookup - X(j)) + Y(j);
    end
end

%% iterations of all 5 loops for Legendre polynomial weights
w2 = ones (2,1);
w3 = ones (3,1);
w4 = ones (4,1);
w5 = ones (5,1);

for i  1:1:2
    l = 1;
    for j = 1:1:2
        if j ~= i
            l = l * (t - roots_2(j)) / (roots_2(i) - roots_2(j));
        end
    end
    w2(i) = int(l, -1, 1);
end

for i  1:1:3
    l = 1;
    for j = 1:1:3
        if j ~= i
            l = l * (t - roots_3(j)) / (roots_3(i) - roots_3(j));
        end
    end
    w3(i) = int(l, -1, 1);
end

for i  1:1:4
    l = 1;
    for j = 1:1:4
        if j ~= i
            l = l * (t - roots_4(j)) / (roots_4(i) - roots_4(j));
        end
    end
    w4(i) = int(l, -1, 1);
end

for i 1:1:5
    l = 1;
    for j = 1:1:5
        if j ~= i
            l = l * (t - roots_5(j)) / (roots_5(i) - roots_5(j));
        end
    end
    w5(i) = int(l, -1, 1);
end

disp("Weights for P_2")
disp(w2)
disp("Weights for P_3")
disp(w3)
disp("Weights for P_4")
disp(w4)
disp("Weights for P_5")
disp(w5)


%% Part 2
%% single integral coding method
function [Y] = my_single_integral(f, a, b)
    Y = int(f, a, b);
    syms x
    syms t

    l_5 = legendreP(5,x);
    roots_5 = vpa(solve(l_5),16);
    w5 = ones(5,1);

    for i = 1:1:5
        l = 1;
        for j = 1:1:5
            if j ~= i
                l = l * ((t - roots_5(j)) / (roots_5(i) - roots_5(j)));
            end
        end
        w5(i) = int(l, -1, 1);
    end

    Y = 0;
    for i = 1:1:size(w5,1)
        x = a + ((b-a)/2 * (roots_5(i) + 1));
        Y = Y + (w5(i) * subs(f));
    end

    Y = Y * ((b-a)/2);
end

%% single integral calculation
a = 0;
b = 1;
f = x;
single_int_ex = my_single_integral(f,a,b);
fprintf("Example of the single integral function for %s from %d to %d: %s\n\n", string(f), a, b, string(single_int_ex))

%% double integral coding method
function [Y] = my_double_integral(u, a, b, g, h)
    addpath('myfunctions')
    l_5 = legendreP(5,x);
    roots_5 = vpa(solve(l_5),16);
    w5 = ones(5,1);
    for i = 1:1:5
        l = 1;
        for j = 1:1:5
            if j ~= i
                l = l * ((t - roots_5(j)) / (roots_5(i) - roots_5(j)));
            end
        end
        w5(i) = int(l, -1, 1);
    end
    Y = 0;
    for i = 1:1:size(w5,1)
        aot = a + ((b-a)/2 * (roots_5(i) + 1));

        z = aot;
        y = x;

        G = g;
        H = h;

        if class(g) == "sym"
            G = subs(G);
        end
        if class(h) == "sym"
            H = subs(H);
        end
        Y = Y + w5(i) * my_single_integral(subs(u), G, H);
    end
    Y = Y * ((b-a)/2);
end

%% Part 3
%% double integral calculation
u = (2 * y * sin(z) + (cos(z)).^2) / (sqrt(1 - y.^2));
a = 0;
b = pi/4;
g = 0;
h = sin(z);

doubl_int_ans = vpa(my_double_integral(u, a, b, g, h), 16);

fprintf("Double Integral Approximation: %s\n\n", string(doubl_int_ans))

sum = 0;

for i=0:1:9
    for j=0:1:9
        a_temp = (a+b)/10*i;
        b_temp = (a+b)/10*(i+1);
        g_temp = (h+g)/10*i;
        h_temp = (h+g)/10*(i+1);

        sum = sum + my_double_integral(u, a_temp, b_temp, g_temp, h_temp);
    end
end

doubl_int_ans_100 = vpa(sum, 16);

fprintf("Double Integral Approximation(100 subintervals): %s\n\n", string(doubl_int_ans_100))
