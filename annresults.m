clc
clear
close all

rng(4)   % fix random seed for reproducibility

% ====================================
% 1. Load Excel Data
% ====================================

data = readtable("E:\VIT\NEW ONE feb26\SLIP MODIFIED\REVISION\ANN\NUMERICAL TABLE.xlsx");

% Inputs (8 parameters)
X = table2array(data(:,1:8));

Cf = table2array(data(:,9))';
Nu = table2array(data(:,10))';  
Sh = table2array(data(:,11))';

Y = [Cf; Nu; Sh];   % 3 × samples

X = X';         % inputs must be features × samples

% ====================================
% 2. Create ANN Model
% Using Levenberg-Marquardt Algorithm
% ====================================

hiddenLayerSize = 10;   % one hidden layers

net = feedforwardnet(hiddenLayerSize,'trainlm');

% ====================================
% 3. Data Division
% ====================================

net.divideParam.trainRatio = 0.70;
net.divideParam.valRatio   = 0.15;
net.divideParam.testRatio  = 0.15;

% ====================================
% 4. Network Training Parameters
% ====================================

net.trainParam.epochs   = 200;      % maximum iterations
net.trainParam.goal     = 1e-8;     % training goal
net.trainParam.min_grad = 1e-7;

% ====================================
% 5. Train the Network
% ====================================

[net,tr] = train(net,X,Y);

% ====================================
% 6. Predict Output
% ====================================

Y_pred = net(X);

Cf_pred = Y_pred(1,:);
Nu_pred = Y_pred(2,:);
Sh_pred = Y_pred(3,:);

% ====================================
% 7. Calculate Performance
% ====================================

performance = perform(net,Y,Y_pred);
% disp('Mean Square Error (MSE)')
% disp(performance)
test_mse    = perform(net, Y(:,tr.testInd), Y_pred(:,tr.testInd));


disp('=== Performance ===')
fprintf('Overall MSE : %.6e\n', performance)
fprintf('Test-set MSE: %.6e\n', test_mse)

% ====================================
% 7b. Per-variable RMSE, MAPE and R²
% ====================================

% --- Cf ---
mse_Cf = mean((Cf - Cf_pred).^2);
rmse_Cf = sqrt(mean((Cf - Cf_pred).^2));
mae_Cf = mean(abs(Cf - Cf_pred));
mape_Cf = mean(abs((Cf - Cf_pred)./Cf))*100;
R2_Cf = 1 - sum((Cf-Cf_pred).^2)/sum((Cf-mean(Cf)).^2);

fprintf('\nCf Performance\n');
fprintf('MSE : %0.6e\n', mse_Cf);
fprintf('RMSE : %.6f\n',rmse_Cf);
fprintf('MAE  : %.6f\n',mae_Cf);
fprintf('MAPE : %.4f %%\n',mape_Cf);
fprintf('R²   : %.6f\n',R2_Cf);

% --- Nu ---
mse_Nu = mean((Nu - Nu_pred).^2);
rmse_Nu = sqrt(mean((Nu - Nu_pred).^2));
mae_Nu = mean(abs(Nu - Nu_pred));
mape_Nu = mean(abs((Nu - Nu_pred)./Nu))*100;
R2_Nu = 1 - sum((Nu-Nu_pred).^2)/sum((Nu-mean(Nu)).^2);

fprintf('\nNu Performance\n');
fprintf('MSE : %0.6e\n', mse_Nu);
fprintf('RMSE : %.6f\n',rmse_Nu);
fprintf('MAE  : %.6f\n',mae_Nu);
fprintf('MAPE : %.4f %%\n',mape_Nu);
fprintf('R²   : %.6f\n',R2_Nu);

% --- Sh ---
mse_Sh = mean((Sh - Sh_pred).^2);
rmse_Sh = sqrt(mean((Sh - Sh_pred).^2));
mae_Sh = mean(abs(Sh - Sh_pred));
mape_Sh = mean(abs((Sh - Sh_pred)./Sh))*100;
R2_Sh = 1 - sum((Sh-Sh_pred).^2)/sum((Sh-mean(Sh)).^2);

fprintf('\nSh Performance\n');
fprintf('MSE : %0.6e\n', mse_Sh);
fprintf('RMSE : %.6f\n',rmse_Sh);
fprintf('MAE  : %.6f\n',mae_Sh);
fprintf('MAPE : %.4f %%\n',mape_Sh);
fprintf('R²   : %.6f\n',R2_Sh);


disp('Dimension checking')

disp(size(Cf))
disp(size(Cf_pred))

disp(size(Nu))
disp(size(Nu_pred))

disp(size(Sh))
disp(size(Sh_pred))

% ====================================
% 8. Performance Plot
% ====================================

figure
plotperform(tr)
title('Training Performance')

% ====================================
% 9. Regression Plot
% ====================================

figure
plotregression( ...
    Y(:,tr.trainInd), Y_pred(:,tr.trainInd), 'Training', ...
    Y(:,tr.valInd),   Y_pred(:,tr.valInd),   'Validation', ...
    Y(:,tr.testInd),  Y_pred(:,tr.testInd),  'Test', ...
    Y, Y_pred, 'All')
sgtitle('ANN Regression Plot (LM Algorithm)')


% ====================================
% Training  state plot 
% ====================================
figure
plottrainstate(tr)
title('ANN Training State')


% ====================================
% 10. Error Histogram
% ====================================


% --- Figure 5: Error Histogram (with train/val/test colors) ---
figure('Name','Error Histogram','NumberTitle','off')
ploterrhist( ...
    Y(:, tr.trainInd) - Y_pred(:, tr.trainInd), 'Training', ...
    Y(:, tr.valInd)   - Y_pred(:, tr.valInd),   'Validation', ...
    Y(:, tr.testInd)  - Y_pred(:, tr.testInd),  'Test' ...
    )


% ====================================
% 11. Numerical vs Predicted Plot
% ====================================
%%%%%
figure
scatter(Cf,Cf_pred,'filled')
hold on
min_val = min([Cf Cf_pred]);
max_val = max([Cf Cf_pred]);
plot([min_val max_val],[min_val max_val],'r','LineWidth',2)
xlabel('Numerical Cf')
ylabel('Predicted CF')
title('Numerical vs Predicted Skin friction')
grid on


%%%
figure
scatter(Nu,Nu_pred,'filled')
hold on
min_val = min([Nu Nu_pred]);
max_val = max([Nu Nu_pred]);
plot([min_val max_val],[min_val max_val],'r','LineWidth',2)
xlabel('Numerical Nu')
ylabel('Predicted Nu')
title('Numerical vs Predicted Nusselt Number')
grid on

%%%
figure
scatter(Sh,Sh_pred,'filled')
hold on
min_val = min([Sh Sh_pred]);
max_val = max([Sh Sh_pred]);
plot([min_val max_val],[min_val max_val],'r','LineWidth',2)
xlabel('Numerical Sh')
ylabel('Predicted Sh')
title('Numerical vs Predicted Sherwood Number')
grid on




% ====================================
% 12. Line Comparison Plots
% ====================================

%%%
figure
plot(Cf,'ro-','LineWidth',1.5)
hold on
plot(Cf_pred,'b--','LineWidth',1.5)
xlabel('Data Index')
ylabel('Skin friction')
legend('Numerical Values','ANN Prediction')
title('Comparison of Numerical and ANN Predicted Skin friction')
grid on


figure
plot(Nu,'ro-','LineWidth',1.5)
hold on
plot(Nu_pred,'b--','LineWidth',1.5)
xlabel('Data Index')
ylabel('Nusselt Number')
legend('Numerical Values','ANN Prediction')
title('Comparison of Numerical and ANN Predicted Nusselt Number')
grid on


figure
plot(Sh,'ro-','LineWidth',1.5)
hold on
plot(Sh_pred,'b--','LineWidth',1.5)
xlabel('Data Index')
ylabel('Sherwood Number')
legend('Numerical Values','ANN Prediction')
title('Comparison of Numerical and ANN Predicted Sherwood Number')
grid on



% % =====================================
% % Performance Metrics
% % =====================================

MSE = perform(net,Y,Y_pred);

% Calculate correlation coefficient for Cf

R_Cf = corrcoef(Cf(tr.testInd),Cf_pred(tr.testInd));
R_Nu = corrcoef(Nu(tr.testInd),Nu_pred(tr.testInd));
R_Sh = corrcoef(Sh(tr.testInd),Sh_pred(tr.testInd));

disp('Correlation coefficient for Cf:')
disp(R_Cf)

disp('Correlation coefficient for Nu:')
disp(R_Nu)

disp('Correlation coefficient for Sh:')
disp(R_Sh)


% =====================================
% Display Training Results
% ====
% =================================

Cf_pred = Cf_pred';
Nu_pred = Nu_pred';
Sh_pred = Sh_pred';

Result_Table = table(Cf',Cf_pred,Nu',Nu_pred,Sh',Sh_pred);

Result_Table.Properties.VariableNames = ...
{'Cf_Numerical','Cf_Predicted','Nu_Numerical','Nu_Predicted','Sh_Numerical','Sh_Predicted'};

disp(Result_Table)

writetable(Result_Table,'ANN_Predicted_Comparison.xlsx');

height(Result_Table)

view(net)
figure('Name','ANN Architecture','NumberTitle','off')
plot(net)
title('ANN Architecture')
exportgraphics(gcf,'ANN_Architecture.png','Resolution',600);
