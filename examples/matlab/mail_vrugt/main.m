% DREAM, sequential vs. 'parallel' (cf. fig. 7 of Vrugt, 2016)
N = 25;
T = 3000;
prior = @(N,d) unifrnd(-20,20,N,d);
pdf = @(x) 1/6 * normpdf(x,-8,1) + 5/6 * normpdf(x,10,1);
tar_fun = pdf(-15:0.1:15);
%[dream_x, dream_px] = dream(prior, pdf, N, T, 1);
[dream2_x, dream2_px] = par_dream(prior, pdf, N, T, 1);

% plot
% dream result
figure
% subplot(2,2,1)
% hold on
% histogram(dream_x, 'Normalization', 'pdf', 'FaceColor', 'red', 'NumBins', 50),xlim([-15 15])
% plot(-15:0.1:15, tar_fun, 'black', 'LineWidth', 2)
% hold off
% subplot(2,2,2)
% plot(reshape(dream_x(:,1,:), T, N), '.')
% % par_dream result
subplot(1,2,1)
hold on
histogram(dream2_x, 'Normalization', 'pdf', 'FaceColor', 'red', 'NumBins', 50),xlim([-15 15])
plot(-15:0.1:15, tar_fun, 'black', 'LineWidth', 2)
hold off
subplot(1,2,2)
plot(reshape(dream2_x(:,1,:), T, N), '.')


% 10-D muttimodal mixture
N = 25;
T = 10000;
pdf = @(x) 1/3 * mvnpdf(x,-5*ones(1,length(x)), diag(ones(length(x),1))) + ...
     2/3 * mvnpdf(x,5*ones(1,length(x)), diag(ones(length(x),1)));
x_tar = -10:0.1:10;
tar_fun = NaN(length(x_tar),1);
for i = 1:length(x_tar), tar_fun(i) = pdf(x_tar(i)); end
prior = @(N,d) unifrnd(-10,10,N,d);
[dream_x2, dream_px2] = dream(prior, pdf, N, T, 10);

figure
hold on
histogram(dream_x2(:,1,:), 'Normalization', 'pdf', 'FaceColor', 'red', 'NumBins', 50),xlim([-10 10])
plot(x_tar, tar_fun, 'black', 'LineWidth', 2)
hold off


% 25-D trimodal mixture
N = 25;
T = 10000;
pdf = @(x) 3/6 * mvnpdf(x, repmat(15,1,length(x)), diag(ones(length(x),1))) +...
           2/6 * mvnpdf(x, repmat(5,1,length(x)), diag(ones(length(x),1))) +...
           1/6 * mvnpdf(x, repmat(-5,1,length(x)), diag(ones(length(x),1)));
x_tar = -10:0.1:20;
tar_fun = NaN(length(x_tar),1);
for i = 1:length(x_tar), tar_fun(i) = pdf(x_tar(i)); end
prior = @(N,d) unifrnd(-10,20,N,d);
[dream_x, dream_px] = dream(prior, pdf, N, T, 15);

figure
hold on
histogram(dream_x(:,1,1), 'Normalization', 'pdf', 'FaceColor', 'red', 'NumBins', 50),xlim([-10 20])
plot(x_tar, tar_fun, 'black', 'LineWidth', 2)
hold off
