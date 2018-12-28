% examples from Paper

% RWM algorithm
prior = @(N,d) unifrnd(-10, 10, N, d);
pdf = @(x) mvnpdf(x, [0 0], [1 0.8; 0.8 1]);
[x, p_x] = rwm(prior, pdf, 50000, 2);

% plot (fig. 4)
figure
subplot(2,2,[1 3])
scatter(x(:,1),x(:,2), '.', 'black'),xlim([-4 4]),ylim([-4 4])
subplot(2,2,2)
plot(x(:,1),'DisplayName','x')
subplot(2,2,4)
plot(x(:,2),'DisplayName','x')

% AM and DE-MC comparison
prior = @(N,d) unifrnd(-20,20,N,d);
pdf = @(x) 1/6 * normpdf(x,-8,1) + 5/6 * normpdf(x,10,1);
tar_fun = pdf(-15:0.1:15);
[am_x, am_px] = am(prior, pdf, 50000, 1);
[de_x, de_px] = de_mc(prior, pdf, 10, 5000, 1);

% plot (fig. 5)
figure
subplot(1,2,1)
hold on
histogram(am_x, 'Normalization', 'pdf', 'FaceColor', 'red', 'NumBins', 50),xlim([-15 15])
plot(-15:0.1:15, tar_fun, 'black', 'LineWidth', 2)
hold off
subplot(1,2,2)
hold on
histogram(de_x, 'Normalization', 'pdf', 'FaceColor', 'red', 'NumBins', 50),xlim([-15 15])
plot(-15:0.1:15, tar_fun, 'black', 'LineWidth', 2)
hold off

% DREAM (fig. 7), sequential and parallel
[dream_x, dream_px] = dream(prior, pdf, 10, 5000, 1);
[dream2_x, dream2_px] = dream_var2(prior, pdf, 10, 5000, 1);

figure
subplot(2,2,1)
hold on
histogram(dream_x, 'Normalization', 'pdf', 'FaceColor', 'red', 'NumBins', 50),xlim([-15 15])
plot(-15:0.1:15, tar_fun, 'black', 'LineWidth', 2)
hold off
subplot(2,2,2)
plot(reshape(dream_x(:,1,:), 5000, 10), '.')

subplot(2,2,3)
hold on
histogram(dream2_x, 'Normalization', 'pdf', 'FaceColor', 'red', 'NumBins', 50),xlim([-15 15])
plot(-15:0.1:15, tar_fun, 'black', 'LineWidth', 2)
hold off
subplot(2,2,4)
plot(reshape(dream2_x(:,1,:), 5000, 10), '.')
