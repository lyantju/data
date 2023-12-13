yangben = xlsread('24H.xls');
X = yangben(:,1);
Y = yangben(:,2);

[fx, Xsort] = ecdf(X);
[fy, Ysort] = ecdf(Y); 

U1 = spline(Xsort(2:end),fx(2:end),X);
V1 = spline(Ysort(2:end),fy(2:end),Y);
%% 
guangyi_X = gevfit(X); 
guangyi_Y = gevfit(Y); 

k_X = guangyi_X(1)
sigma_X = guangyi_X(2) 
miu_X = guangyi_X(3) 

k_Y = guangyi_Y(1)
sigma_Y = guangyi_Y(2) 
miu_Y = guangyi_Y(3) 

U_guangyi = gevcdf(X, k_X, sigma_X, miu_X);

V_guangyi = gevcdf(Y, k_Y, sigma_Y, miu_Y);


rmseX = sqrt(mean((U_guangyi - U1).^2))
rmseY = sqrt(mean((V_guangyi - V1).^2))

nparamX = numel(guangyi_X) ;
nobsX = numel(X);
aicX = 2 * nobsX * log(rmseX) + 2 * nparamX;
bicX = 2 * nobsX * log(rmseX) + log(nobsX) * nparamX;

nparamY = numel(guangyi_Y); 
nobsY = numel(Y); 
aicY = 2 * nobsY * log(rmseY) + 2 * nparamY;
bicY = 2 * nobsY * log(rmseY) + log(nobsY) * nparamY;

[h, p, ksstat] = kstest2(U_guangyi, U1);

fprintf('K-S: %f\n', ksstat);
fprintf('p: %f\n', p);

[h_Y, p_Y, ksstat_Y] = kstest2(V_guangyi, V1);

fprintf('K-Sͳ: %f\n', ksstat_Y);
fprintf('p: %f\n', p_Y);


%% 
norm_X = fitdist(X(:), 'Normal');
norm_Y = fitdist(Y(:), 'Normal');

U_norm = normcdf(X, norm_X.mu, norm_X.sigma);
V_norm = normcdf(Y, norm_X.mu, norm_X.sigma);


rmseX = sqrt(mean((U_norm - U1).^2))
rmseY = sqrt(mean((V_norm - V1).^2))

nparamX = 2; 
nobsX = numel(X); 
nparamY = 2; 
nobsY = numel(Y); 

aicX = 2 * nobsX * log(rmseX) + 2 * nparamX;
aicY = 2 * nobsY * log(rmseY) + 2 * nparamY;

bicX = 2 * nobsX * log(rmseX) + log(nobsX) * nparamX;
bicY = 2 * nobsY * log(rmseY) + log(nobsY) * nparamY;

[h, p, ksstat] = kstest2(U_norm, U1);
fprintf('K-S: %f\n', ksstat)
fprintf('p: %f\n', p)

[h_Y, p_Y, ksstat_Y] = kstest2(V_norm, V1);
fprintf('K-S: %f\n', ksstat_Y)
fprintf('p: %f\n', p_Y)

%% 
if any(X == 0)
    X(X == 0) = 1e-6;
end
if any(Y == 0)
    Y(Y == 0) = 1e-6;
end
Gamma_X = fitdist(X, 'Gamma');
Gamma_Y = fitdist(Y, 'Gamma');
shape_Gamma_X = Gamma_X.Params(1); 
scale_Gamma_X = Gamma_X.Params(2);
shape_Gamma_Y = Gamma_Y.Params(1) 
scale_Gamma_Y = Gamma_Y.Params(2)

U_Gamma = gamcdf(X, shape_Gamma_X, scale_Gamma_X);
V_Gamma = gamcdf(Y, shape_Gamma_Y, scale_Gamma_Y);

x = gaminv(0.95, shape_Gamma_X, scale_Gamma_X)
Fp = gamcdf(x, shape_Gamma_X, scale_Gamma_X)

rmseX = sqrt(mean((U_Gamma - U1).^2))
rmseY = sqrt(mean((V_Gamma - V1).^2))

nparamX = 2; 
nobsX = numel(X); 
nparamY = 2; 
nobsY = numel(Y); 

aicX = 2 * nobsX * log(rmseX) + 2 * nparamX;
aicY = 2 * nobsY * log(rmseY) + 2 * nparamY;
bicX = 2 * nobsX * log(rmseX) + log(nobsX) * nparamX;
bicY = 2 * nobsY * log(rmseY) + log(nobsY) * nparamY;
[h, p, ksstat] = kstest2(U_Gamma, U1);

fprintf('K-S: %f\n', ksstat)
fprintf('p: %f\n', p)

[h_Y, p_Y, ksstat_Y] = kstest2(V_Gamma, V1);

fprintf('K-S: %f\n', ksstat_Y)
fprintf('p: %f\n', p_Y)

%% 
Weibull_X = wblfit(X);
Weibull_Y = wblfit(Y);
scale_Weibull_X = Weibull_X(1); 
shape_Weibull_X = Weibull_X(2); 
scale_Weibull_Y = Weibull_Y(1)
shape_Weibull_Y = Weibull_Y(2)

U_Weibull = wblcdf(X, scale_Weibull_X, shape_Weibull_X);
V_Weibull = wblcdf(Y, scale_Weibull_Y, shape_Weibull_Y);


rmseX = sqrt(mean((U_Weibull - U1).^2))
rmseY = sqrt(mean((V_Weibull - V1).^2))
% nparamX = numel(Weibull_X);
nparamX = 2; 
nobsX = numel(X); 
nparamY = 2; 
nobsY = numel(Y); 

aicX = 2 * nobsX * log(rmseX) + 2 * nparamX;
aicY = 2 * nobsY * log(rmseY) + 2 * nparamY;
bicX = 2 * nobsX * log(rmseX) + log(nobsX) * nparamX;
bicY = 2 * nobsY * log(rmseY) + log(nobsY) * nparamY;
[h, p, ksstat] = kstest2(U_Weibull, U1);
fprintf('K-S: %f\n', ksstat)
fprintf('p: %f\n', p)

[h_Y, p_Y, ksstat_Y] = kstest2(V_Weibull, V1);
fprintf('K-S: %f\n', ksstat_Y)
fprintf('p: %f\n', p_Y)


%% 
U_hmd = ksdensity(X,X,'function','cdf');
V_hmd = ksdensity(Y,Y,'function','cdf');

rmseX = sqrt(mean((U_hmd - U1).^2))
rmseY = sqrt(mean((V_hmd - V1).^2))

[h, p, ksstat] = kstest2(U_hmd, U1);
fprintf('K-S: %f\n', ksstat)
fprintf('p: %f\n', p)

[h_Y, p_Y, ksstat_Y] = kstest2(V_hmd, V1);
fprintf('K-S: %f\n', ksstat_Y)
fprintf('p: %f\n', p_Y)
