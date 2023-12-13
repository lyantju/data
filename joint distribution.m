%% 
yangben = xlsread('12H.xls');
X = yangben(:,1);
Y = yangben(:,2);
%% 
params_X = gevfit(X); 
k_X = params_X(1) 
sigma_X = params_X(2)
miu_X = params_X(3) 

U = gevcdf(X, k_X, sigma_X, miu_X);

Weibull_Y = wblfit(Y);
scale_Weibull_Y = Weibull_Y(1)
shape_Weibull_Y = Weibull_Y(2) 
V = wblcdf(Y, scale_Weibull_Y, shape_Weibull_Y);

%% 
rho_norm = copulafit('Gaussian',[U(:), V(:)])
[Udata,Vdata] = meshgrid(linspace(0,1,31));  
Cpdf_norm = copulapdf('Gaussian',[Udata(:), Vdata(:)],rho_norm);
Ccdf_norm = copulacdf('Gaussian',[Udata(:), Vdata(:)],rho_norm);

%% 
j = 100;
i = 1-1/j;

xp = miu_X-sigma_X/k_X*(1-(-log(i))^-k_X)
yp = wblinv(i, scale_Weibull_Y, shape_Weibull_Y)

Fp = gevcdf(xp, k_X, sigma_X, miu_X)
Fz = wblcdf(yp, scale_Weibull_Y, shape_Weibull_Y)

Tp = 1/(1-Fp);
Tz = 1/(1-Fz);

Ccdf = copulacdf('Gaussian',[Fp, Fz],rho_norm);
Tlianhe = 1/(1-Ccdf)
Ttongxian = 1/(1-Fp-Fz+Ccdf)

N = 10000000; 
XY = copularnd('Gaussian', rho_norm, N);  
C = copulacdf('Gaussian', XY, rho_norm);  
C1=copulacdf( 'Gaussian', [i, i], rho_norm);
area_points = sum(C <= C1);
area = (area_points / numel(C));
kendallT=1/(1-area)
%% 
j_values = [2, 3, 5, 10, 20, 50, 100, 200];

results = cell(length(j_values), 6);
results{1, 1} = 'j';
results{1, 2} = 'xp';
results{1, 3} = 'yp';
results{1, 4} = 'Tlianhe';
results{1, 5} = 'Ttongxian';
results{1, 6} = 'kendallT';

for idx = 1:length(j_values)
    j = j_values(idx);
    i = 1 - 1/j;

    xp = miu_X - sigma_X/k_X * (1 - (-log(i))^(-k_X));
    yp = wblinv(i, scale_Weibull_Y, shape_Weibull_Y);

    Fp = gevcdf(xp, k_X, sigma_X, miu_X);
    Fz = wblcdf(yp, scale_Weibull_Y, shape_Weibull_Y);

    Tp = 1/(1 - Fp);
    Tz = 1/(1 - Fz);

    Ccdf = copulacdf('Gaussian', [Fp, Fz], rho_norm);

    Tlianhe = 1/(1 - Ccdf);

    Ttongxian = 1/(1 - Fp - Fz + Ccdf);

    N = 10000000;
    XY = copularnd('Gaussian', rho_norm, N);
    C = copulacdf('Gaussian', XY, rho_norm);
    C1 = copulacdf('Gaussian', [i, i], rho_norm);
    area_points = sum(C <= C1);
    area = area_points / numel(C);
    kendallT = 1/(1 - area);

    results{idx+1, 1} = j;
    results{idx+1, 2} = xp;
    results{idx+1, 3} = yp;
    results{idx+1, 4} = Tlianhe;
    results{idx+1, 5} = Ttongxian;
    results{idx+1, 6} = kendallT;
end

xlswrite('result.xlsx', results);

%% 
lianheT = 20;
Ccdf_lianheT = 1-1/lianheT;
tolerance = 1e-6;
lower_bound = 0;
upper_bound = 1;

while (upper_bound - lower_bound) > tolerance
    midpoint = (lower_bound + upper_bound) / 2;
    Ccdf_midpoint = copulacdf('Gaussian', [midpoint, midpoint], rho_norm);
    
    if Ccdf_midpoint < Ccdf_lianheT
        lower_bound = midpoint;
    else
        upper_bound = midpoint;
    end
end

u1 = lower_bound
v1 = lower_bound

p1_lianheT = miu_X-sigma_X/k_X*(1-(-log(u1))^-k_X)
z1_lianheT = wblinv(v1, scale_Weibull_Y, shape_Weibull_Y)

copula_1 = @(u, v) copulacdf('Gaussian', [u v], rho_norm);

u_values1 = linspace(0, 1, 500);
v_values1 = linspace(0, 1, 500);

valid_combinations1 = [];

for i = 1:length(u_values1)
    for j = 1:length(v_values1)
        u = u_values1(i);
        v = v_values1(j);
        
        cdf = copula_1(u, v);
        
        if abs(cdf - Ccdf_lianheT )< 0.01 
            valid_combinations1 = [valid_combinations1; u v];
        end
    end
end

Ccdf_lianheT_all = copulapdf('Gaussian',[valid_combinations1(:, 1), valid_combinations1(:, 2)],rho_norm);
p1_all = miu_X-sigma_X/k_X*(1-(-log(valid_combinations1(:, 1))).^-k_X);
z1_all = wblinv(valid_combinations1(:, 2), scale_Weibull_Y, shape_Weibull_Y);
u1_all = gevpdf(p1_all, k_X, sigma_X, miu_X);
v1_all = wblpdf(z1_all, scale_Weibull_Y, shape_Weibull_Y);


product1 = Ccdf_lianheT_all.* u1_all.* v1_all;
[max_product, idx1] = max(product1);

max_u1 = valid_combinations1(idx1, 1)
max_v1 = valid_combinations1(idx1, 2)

p2_lianheT = miu_X-sigma_X/k_X*(1-(-log(max_u1))^-k_X)
z2_lianheT = wblinv(max_v1, scale_Weibull_Y, shape_Weibull_Y)
%% 
lianheT_values = [2, 3, 5, 10, 20, 50, 100, 200];

results = cell(length(lianheT_values), 8);
results{1, 1} = 'lianheT';
results{1, 2} = 'u1';
results{1, 3} = 'v1';
results{1, 4} = 'p1_lianheT';
results{1, 5} = 'z1_lianheT';
results{1, 6} = 'max_u1';
results{1, 7} = 'max_v1';
results{1, 8} = 'p2_lianheT';
results{1, 9} = 'z2_lianheT';

for idx = 1:length(lianheT_values)
    lianheT = lianheT_values(idx);
    
    Ccdf_lianheT = 1 - 1/lianheT;

    tolerance = 1e-6;
    lower_bound = 0;
    upper_bound = 1;

    while (upper_bound - lower_bound) > tolerance
        midpoint = (lower_bound + upper_bound) / 2;
        Ccdf_midpoint = copulacdf('Gaussian', [midpoint, midpoint], rho_norm);

        if Ccdf_midpoint < Ccdf_lianheT
            lower_bound = midpoint;
        else
            upper_bound = midpoint;
        end
    end

    u1 = lower_bound;
    v1 = lower_bound;

    p1_lianheT = miu_X - sigma_X/k_X * (1 - (-log(u1))^(-k_X));
    z1_lianheT = wblinv(v1, scale_Weibull_Y, shape_Weibull_Y);

    copula_1 = @(u, v) copulacdf('Gaussian', [u v], rho_norm);

    u_values1 = linspace(0, 1, 500);
    v_values1 = linspace(0, 1, 500);

    valid_combinations1 = [];

    for i = 1:length(u_values1)
        for j = 1:length(v_values1)
            u = u_values1(i);
            v = v_values1(j);

            cdf = copula_1(u, v);

            if abs(cdf - Ccdf_lianheT) < 0.01
                valid_combinations1 = [valid_combinations1; u v];
            end
        end
    end

    Ccdf_lianheT_all = copulapdf('Gaussian', [valid_combinations1(:, 1), valid_combinations1(:, 2)], rho_norm);
    p1_all = miu_X - sigma_X/k_X * (1 - (-log(valid_combinations1(:, 1))).^(-k_X));
    z1_all = wblinv(valid_combinations1(:, 2), scale_Weibull_Y, shape_Weibull_Y);
    u1_all = gevpdf(p1_all, k_X, sigma_X, miu_X);
    v1_all = wblpdf(z1_all, scale_Weibull_Y, shape_Weibull_Y);

    product1 = Ccdf_lianheT_all .* u1_all .* v1_all;
    [max_product, idx1] = max(product1);

    max_u1 = valid_combinations1(idx1, 1);
    max_v1 = valid_combinations1(idx1, 2);

    p2_lianheT = miu_X - sigma_X/k_X * (1 - (-log(max_u1))^(-k_X));
    z2_lianheT = wblinv(max_v1, scale_Weibull_Y, shape_Weibull_Y);

    results{idx+1, 1} = lianheT;
    results{idx+1, 2} = u1;
    results{idx+1, 3} = v1;
    results{idx+1, 4} = p1_lianheT;
    results{idx+1, 5} = z1_lianheT;
    results{idx+1, 6} = max_u1;
    results{idx+1, 7} = max_v1;
    results{idx+1, 8} = p2_lianheT;
    results{idx+1, 9} = z2_lianheT;
end

xlswrite('result.xlsx', results);

%% 
tongxianT = 200;

fun = @(u) 1/tongxianT - 1 + 2*u - copulacdf('Gaussian', [u, u], rho_norm);
u_solution = fsolve(fun, 0); 
v_solution = u_solution;

u2 = u_solution
v2 = v_solution

p1_tongxianT = miu_X-sigma_X/k_X*(1-(-log(u2))^-k_X)
z1_tongxianT = wblinv(v2, scale_Weibull_Y, shape_Weibull_Y)

copula_1 = @(u, v) copulacdf('Gaussian', [u v], rho_norm);

u_values2 = linspace(0, 1, 500);
v_values2 = linspace(0, 1, 500);

valid_combinations2 = [];

for i = 1:length(u_values2)
    for j = 1:length(v_values2)
        u = u_values2(i);
        v = v_values2(j);
        
        cdf = copula_1(u, v);
        Ccdf_tongxianT=(1/tongxianT)-1+u+v;
        
        if abs(cdf - Ccdf_tongxianT )< 0.001 
            valid_combinations2 = [valid_combinations2; u v];
        end
    end
end

Ccdf_tongxianT_all = copulapdf('Gaussian',[valid_combinations2(:, 1), valid_combinations2(:, 2)],rho_norm);
p2_all = miu_X-sigma_X/k_X*(1-(-log(valid_combinations2(:, 1))).^-k_X);
z2_all = wblinv(valid_combinations2(:, 2), scale_Weibull_Y, shape_Weibull_Y);
u2_all = gevpdf(p2_all, k_X, sigma_X, miu_X);
v2_all = wblpdf(z2_all, scale_Weibull_Y, shape_Weibull_Y);

product2 = Ccdf_tongxianT_all.* u2_all.* v2_all;
[max_product, idx2] = max(product2);

max_u2 = valid_combinations2(idx2, 1)
max_v2 = valid_combinations2(idx2, 2)

p2_tongxianT = miu_X-sigma_X/k_X*(1-(-log(max_u2))^-k_X)
z2_tongxianT = wblinv(max_v2, scale_Weibull_Y, shape_Weibull_Y)

%% 
tongxianT_values = [2, 3, 5, 10, 20, 50, 100, 200];

results = cell(length(tongxianT_values), 8);
results{1, 1} = 'tongxianT';
results{1, 2} = 'u2';
results{1, 3} = 'v2';
results{1, 4} = 'p1_tongxianT';
results{1, 5} = 'z1_tongxianT';
results{1, 6} = 'max_u2';
results{1, 7} = 'max_v2';
results{1, 8} = 'p2_tongxianT';
results{1, 9} = 'z2_tongxianT';

for idx = 1:length(tongxianT_values)
    tongxianT = tongxianT_values(idx);

    fun = @(u) 1/tongxianT - 1 + 2*u - copulacdf('Gaussian', [u, u], rho_norm);
    u_solution = fsolve(fun, 0); 
    v_solution = u_solution;

    u2 = u_solution;
    v2 = v_solution;

    p1_tongxianT = miu_X - sigma_X/k_X * (1 - (-log(u2))^(-k_X));
    z1_tongxianT = wblinv(v2, scale_Weibull_Y, shape_Weibull_Y);

    copula_1 = @(u, v) copulacdf('Gaussian', [u v], rho_norm);

    u_values2 = linspace(0, 1, 500);
    v_values2 = linspace(0, 1, 500);

    valid_combinations2 = [];

    for i = 1:length(u_values2)
        for j = 1:length(v_values2)
            u = u_values2(i);
            v = v_values2(j);

            cdf = copula_1(u, v);
            Ccdf_tongxianT = (1/tongxianT) - 1 + u + v;

            if abs(cdf - Ccdf_tongxianT) < 0.001
                valid_combinations2 = [valid_combinations2; u v];
            end
        end
    end

    Ccdf_tongxianT_all = copulapdf('Gaussian', [valid_combinations2(:, 1), valid_combinations2(:, 2)], rho_norm);
    p2_all = miu_X - sigma_X/k_X * (1 - (-log(valid_combinations2(:, 1))).^(-k_X));
    z2_all = wblinv(valid_combinations2(:, 2), scale_Weibull_Y, shape_Weibull_Y);
    u2_all = gevpdf(p2_all, k_X, sigma_X, miu_X);
    v2_all = wblpdf(z2_all, scale_Weibull_Y, shape_Weibull_Y);

    product2 = Ccdf_tongxianT_all .* u2_all .* v2_all;
    [max_product, idx2] = max(product2);

    max_u2 = valid_combinations2(idx2, 1);
    max_v2 = valid_combinations2(idx2, 2);

    p2_tongxianT = miu_X - sigma_X/k_X * (1 - (-log(max_u2))^(-k_X));
    z2_tongxianT = wblinv(max_v2, scale_Weibull_Y, shape_Weibull_Y);

    results{idx+1, 1} = tongxianT;
    results{idx+1, 2} = u2;
    results{idx+1, 3} = v2;
    results{idx+1, 4} = p1_tongxianT;
    results{idx+1, 5} = z1_tongxianT;
    results{idx+1, 6} = max_u2;
    results{idx+1, 7} = max_v2;
    results{idx+1, 8} = p2_tongxianT;
    results{idx+1, 9} = z2_tongxianT;
end

xlswrite('result.xlsx', results);

%% 
kendallT = 100;
target_area = 1-1/kendallT;
N = 10000000; 

XY = copularnd('Gaussian', rho_norm, N);

C = copulacdf('Gaussian', XY, rho_norm);

[C_sorted, sort_indices] = sort(C);

threshold = C_sorted(floor(target_area * N));

threshold_index = find(C == threshold, 1);

u_value_np = XY(threshold_index, 1);
v_value_np = XY(threshold_index, 2);


tolerance = 1e-6;
lower_bound = 0;
upper_bound = 1;
Ccdf_kendallT = threshold;
while (upper_bound - lower_bound) > tolerance
    midpoint = (lower_bound + upper_bound) / 2;
    Ccdf_midpoint = copulacdf('Gaussian', [midpoint, midpoint], rho_norm);
    
    if Ccdf_midpoint < Ccdf_kendallT
        lower_bound = midpoint;
    else
        upper_bound = midpoint;
    end
end

u3 = lower_bound
v3 = lower_bound

p1_kendallT = miu_X-sigma_X/k_X*(1-(-log(u3))^-k_X)
z1_kendallT = wblinv(v3, scale_Weibull_Y, shape_Weibull_Y)


copula_3 = @(u, v) copulacdf('Gaussian', [u v], rho_norm);

u_values3 = linspace(0, 1, 500);
v_values3 = linspace(0, 1, 500);

valid_combinations3 = [];

for i = 1:length(u_values3)
    for j = 1:length(v_values3)
        u = u_values3(i);
        v = v_values3(j);
        
        cdf = copula_3(u, v);
        
        if abs(cdf - Ccdf_kendallT )< 0.01 
            valid_combinations3 = [valid_combinations3; u v];
        end
    end
end

Ccdf_kendallT_all = copulapdf('Gaussian',[valid_combinations3(:, 1), valid_combinations3(:, 2)],rho_norm);
p3_all = miu_X-sigma_X/k_X*(1-(-log(valid_combinations3(:, 1))).^-k_X);
z3_all = wblinv(valid_combinations3(:, 2), scale_Weibull_Y, shape_Weibull_Y);
u3_all = gevpdf(p3_all, k_X, sigma_X, miu_X);
v3_all = wblpdf(z3_all, scale_Weibull_Y, shape_Weibull_Y);

product3 = Ccdf_kendallT_all.* u3_all.* v3_all;
[max_product, idx3] = max(product3);


max_u3 = valid_combinations3(idx3, 1)
max_v3 = valid_combinations3(idx3, 2)


p2_kendallT = miu_X-sigma_X/k_X*(1-(-log(max_u3))^-k_X)
z2_kendallT = wblinv(max_v3, scale_Weibull_Y, shape_Weibull_Y)
%% 
kendallT_values = [2, 3, 5, 10, 20, 50, 100, 200];

results = cell(length(kendallT_values), 7);
results{1, 1} = 'kendallT';
results{1, 2} = 'equal_elements';
results{1, 3} = 'p1_kendallT';
results{1, 4} = 'z1_kendallT';
results{1, 5} = 'max_u3';
results{1, 6} = 'max_v3';
results{1, 7} = 'p2_kendallT';
results{1, 8} = 'z2_kendallT';

for idx = 1:length(kendallT_values)
    kendallT = kendallT_values(idx);
    target_area = 1 - 1/kendallT;

    XY = copularnd('Gaussian', rho_norm, N);

    C = copulacdf('Gaussian', XY, rho_norm);

    [C_sorted, sort_indices] = sort(C);

    threshold = C_sorted(floor(target_area * N));

    threshold_index = find(C == threshold, 1);

    u_value_np = XY(threshold_index, 1);
    v_value_np = XY(threshold_index, 2);

    tolerance = 1e-6;
    lower_bound = 0;
    upper_bound = 1;
    Ccdf_kendallT = threshold;

    while (upper_bound - lower_bound) > tolerance
        midpoint = (lower_bound + upper_bound) / 2;
        Ccdf_midpoint = copulacdf('Gaussian', [midpoint, midpoint], rho_norm);

        if Ccdf_midpoint < Ccdf_kendallT
            lower_bound = midpoint;
        else
            upper_bound = midpoint;
        end
    end

    u3 = lower_bound;
    v3 = lower_bound;

    p1_kendallT = miu_X - sigma_X/k_X * (1 - (-log(u3))^(-k_X));
    z1_kendallT = wblinv(v3, scale_Weibull_Y, shape_Weibull_Y);

    copula_3 = @(u, v) copulacdf('Gaussian', [u v], rho_norm);

    u_values3 = linspace(0, 1, 500);
    v_values3 = linspace(0, 1, 500);

    valid_combinations3 = [];

    for i = 1:length(u_values3)
        for j = 1:length(v_values3)
            u = u_values3(i);
            v = v_values3(j);

            cdf = copula_3(u, v);

            if abs(cdf - Ccdf_kendallT) < 0.01 
                valid_combinations3 = [valid_combinations3; u v];
            end
        end
    end

    Ccdf_kendallT_all = copulapdf('Gaussian', [valid_combinations3(:, 1), valid_combinations3(:, 2)], rho_norm);
    p3_all = miu_X - sigma_X/k_X * (1 - (-log(valid_combinations3(:, 1))).^(-k_X));
    z3_all = wblinv(valid_combinations3(:, 2), scale_Weibull_Y, shape_Weibull_Y);
    u3_all = gevpdf(p3_all, k_X, sigma_X, miu_X);
    v3_all = wblpdf(z3_all, scale_Weibull_Y, shape_Weibull_Y);

    product3 = Ccdf_kendallT_all .* u3_all .* v3_all;
    [max_product, idx3] = max(product3);

    max_u3 = valid_combinations3(idx3, 1);
    max_v3 = valid_combinations3(idx3, 2);

    p2_kendallT = miu_X - sigma_X/k_X * (1 - (-log(max_u3))^(-k_X));
    z2_kendallT = wblinv(max_v3, scale_Weibull_Y, shape_Weibull_Y);

    results{idx+1, 1} = kendallT;
    results{idx+1, 2} = u3;
    results{idx+1, 3} = p1_kendallT;
    results{idx+1, 4} = z1_kendallT;
    results{idx+1, 5} = max_u3;
    results{idx+1, 6} = max_v3;
    results{idx+1, 7} = p2_kendallT;
    results{idx+1, 8} = z2_kendallT;
end

xlswrite('result.xlsx', results);
