%% Plot data
Y_data = dlmread('cases.txt');
X_data = dlmread('X_process.txt');
time = 0:1/52:1/52*500;
time = time(1:500);

% plot Y process
figure
plot(time, Y_data)
xlabel('Time')
ylabel('Cases')

% Plot X process
figure
subplot(2,2,1)
plot(time, X_data(:,1))
ylabel('S')
xlabel('Time (years)')
subplot(2,2,2)
plot(time,X_data(:,2))
ylabel('I')
xlabel('Time (years)')
subplot(2,2,3)
plot(time,X_data(:,3))
ylabel('R')
xlabel('Time (years)')
subplot(2,2,4)
plot(time,X_data(:,4))
ylabel('H')
xlabel('Time (years)')
%% plot parameter estiamtions 
data_IF1 = dlmread('data_IF1.txt');
data_IF2 = dlmread('data_IF2.txt');
data_AIF = dlmread('data_AIF.txt');
K = 2;
start_ind = 1:3:K*3;

% [beta gamma rho]
% IF1
figure
subplot(131)
plot(data_IF1(:,start_ind))
hold on
%plot(theta_ture, 'k--')
title('\beta')
xlabel('Iteration')
ylabel('Value')
subplot(132)
plot(data_IF1(:,start_ind+1))
hold on
%plot(theta_ture, 'k--')
xlabel('Iteration')
ylabel('Value')
title('\gamma')
subplot(133)
plot(normcdf(data_IF1(:,start_ind+2)))
hold on
%plot(theta_ture, 'k--')
title('\rho')
xlabel('Iteration')
ylabel('Value')


% IF2
figure
subplot(131)
plot(data_IF2(:,start_ind))
hold on
%plot(theta_ture, 'k--')
xlabel('Iteration')
ylabel('Value')
title('\beta')
subplot(132)
plot(data_IF2(:,start_ind+1))
hold on
%plot(theta_ture, 'k--')
xlabel('Iteration')
ylabel('Value')
title('\gamma')
subplot(133)
plot(normcdf(data_IF2(:,start_ind+2)))
hold on
%plot(theta_ture, 'k--')
xlabel('Iteration')
ylabel('Value')
title('\rho')




% AIF
figure
subplot(131)
plot(data_AIF(1:50,start_ind))
hold on
%plot(theta_ture_a, 'k--')
xlabel('Iteration')
ylabel('Value')
title('\beta')
subplot(132)
plot(data_AIF(1:50,start_ind+1))
hold on
%plot(theta_ture_a, 'k--')
xlabel('Iteration')
ylabel('Value')
title('\gamma')
subplot(133)
plot(normcdf(data_AIF(1:50,start_ind+2)))
hold on
%plot(theta_ture_a, 'k--')
xlabel('Iteration')
ylabel('Value')
title('\rho')



