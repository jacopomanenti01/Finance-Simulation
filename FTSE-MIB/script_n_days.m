clear
%% Get data
% get data from FTSE MIB 
dataFTSE = hist_stock_data('12312014','09112023','BTC-EUR');
data = dataFTSE.AdjClose; % adjusted end-of-day close values {Si}
data_ok = fillmissing(data,'nearest');
N = length(data_ok);

%choose the size of the multiday return
n = 10;


% compute n-day overlapping simple returns
ri = zeros(N-n,1); %6574
for i = 1:(N-n)
    ri(i,:) = data_ok(i+n,:)/data_ok(i,:) - 1;
end

%plot empirical overlapping n-day returns vs Normal distribution
histogram(ri,100,'normalization','pdf','FaceColor','auto')
X = min(ri):.01:max(ri);
hold on
a = plot(X,normpdf(X,mean(ri),std(ri)),'k-.','LineWidth',2);
pause
axis([0.1 max(ri) 0 1])
hold off
pause

% compute n-day non-overlapping simple returns
M = floor((N - n) / n) + 1;
Rj = zeros(M,n);
for i = 1:M
    for j = 1:n
        index = (i-1)*n+j;
        if index <= N-n;
            Rj(i,j) = ri(index,:);
        end
    end
end

%% Plot against Normal Distribution
%plot empirical non-overlapping n-day returns (for j=1) vs Normal distribution
histogram(Rj(:,1),30,'normalization','pdf')
hold on
x1 = min(Rj(:,1)):.01:max(Rj(:,1));
b = plot(x1,normpdf(x1,mean(Rj(:,1)),std(Rj(:,1))),'k-.','LineWidth',2);
hold off
pause


%% split non-overlapping n-day returns into right and left tail

Rj_right = NaN(size(Rj));
right = Rj>=0;
Rj_right(right) = Rj(right);

Rj_left = NaN(size(Rj));
left = Rj<0;
Rj_left(left) = Rj(left);

%% estimating parameters of power law distribution for each j
parameters_pl_pos=zeros(n,3);
parameters_pl_neg=zeros(n,3);

%right tail
for j=1:n
    positive_returns = Rj(Rj(:, j) > 0, j);
    [alpha, xmin, L]=plfit(positive_returns);
    parameters_pl_pos(j,[1 2 3]) = [alpha, xmin, L]; %add alfa, xmin, log likelihood
    n2 = length(positive_returns);
    alpha = parameters_pl_pos(j,1);
    parameters_pl_pos(j, 4) = ((alpha - 1)/sqrt(n2)); %sd
    parameters_pl_pos(j,5) = n2; %number of positive returns
    parameters_pl_pos(j,6) = length(positive_returns(positive_returns > parameters_pl_pos(j,2))); %%number of positive returns > xmin
end


%left tail
Rj2 = -Rj;
for j=1:n
    negative_returns = Rj2(Rj2(:, j) > 0, j);
    [alpha, xmin, L]=plfit(negative_returns);
    parameters_pl_neg(j,[1 2 3]) = [alpha, xmin, L];
    n1 = length(negative_returns);
    alpha = parameters_pl_neg(j,1);
    parameters_pl_neg(j, 4) = ((alpha - 1)/sqrt(n1));
    parameters_pl_neg(j,5) = n1;
    parameters_pl_neg(j,6) = length(negative_returns(negative_returns > parameters_pl_neg(j,2)));
end

%% AVG
par_avg_pos = mean(parameters_pl_pos)
par_avg_neg = mean(parameters_pl_neg)

%% plotting the power law
% plotting the power law for negative returns against the empirical data
for j=1:n
    subplot((n/5),5,j)
    histogram(Rj2(Rj2(:, j) > parameters_pl_neg(j,2), j),20,'normalization','pdf','FaceColor','c','EdgeColor', 'black')
    hold on
    x = linspace(parameters_pl_neg(j,2), 10, 1000);
    pdf_values = (parameters_pl_neg(j,1) - 1) * (x / parameters_pl_neg(j,2)).^(-parameters_pl_neg(j,1));
    pdf_values = pdf_values / trapz(x, pdf_values);
    a = plot(x, pdf_values, 'k-', 'LineWidth', 2);
    axis([0 0.3 0 30])
    hold off
end
pause

% plotting the power law for positive returns with empirical data
for j=1:n
    subplot((n/5),5,j)
    histogram(Rj(Rj(:, j) > parameters_pl_pos(j,2), j),20,'normalization','pdf','FaceColor','c','EdgeColor', 'black')
    hold on
    x = linspace(parameters_pl_pos(j,2), 10, 1000);
    pdf_values = (parameters_pl_pos(j,1) - 1) * (x / parameters_pl_pos(j,2)).^(-parameters_pl_pos(j,1));
    pdf_values = pdf_values / trapz(x, pdf_values);
    plot(x, pdf_values, 'k-', 'LineWidth', 2);
    axis([0 0.3 0 30])
    hold off
end
pause
close

%% Additional and Optional: plot della ccdf power law teorica con i rendimenti positivi
%right tail
%for j=1:n
%    h=plplot(Rj(Rj(:, j) > 0, j), parameters_pl_pos(j,2), parameters_pl_pos(j,1));
%    pause
%    close
%end
%left tail
%for j=1:n
%    plplot(Rj2(Rj2(:, j) > 0, j), parameters_pl_neg(j,2), parameters_pl_neg(j,1));
%    pause
%    close
%end

%% estimating parameters of exponential distribution for each j for positive returns

parameters_pl_pos_exp=zeros(n,2);
parameters_pl_neg_exp=zeros(n,2);

%right tail
for j=1:n
    positive_returns_1 = Rj(Rj(:, j) > 0, j);
    pd = fitdist(positive_returns_1,'Exponential');
    parameters_pl_pos_exp(j,:) = [pd.mu,-pd.NLogL];
end

%left tail
for j=1:n
    negative_returns_1 = Rj2(Rj2(:, j) > 0, j);
    pd = fitdist(negative_returns_1,'Exponential');
    parameters_pl_neg_exp(j,:) = [pd.mu,-pd.NLogL];
end

%%plotting the exponential distribution
% plotting the exponential distribution for positive returns against the empirical data

for j = 1:n
    lambda = parameters_pl_pos_exp(j, 1);
    subplot((n/5),5,j)
    histogram(Rj(Rj(:, j)>0,j), 100, 'normalization', 'pdf', 'FaceColor', 'c', 'EdgeColor', 'black');
    hold on
    x = linspace(min(Rj(:, j)), 10, 10000);
    y = exppdf(x, lambda);
    plot(x, y, 'k-', 'LineWidth', 2);
    axis([0 0.3 0 30])
    hold off
end
pause
close


% plotting the exponential distribution for negative returns against the empirical data
for j = 1:n
    lambda = parameters_pl_neg_exp(j, 1);
    subplot((n/5),5,j)
    histogram(Rj2(Rj2(:, j)>0,j), 100, 'normalization', 'pdf', 'FaceColor', 'c', 'EdgeColor', 'black');
    hold on
    x = linspace(min(Rj2(:, j)), 10, 1000);
    y = exppdf(x, lambda);
    plot(x, y, 'k-', 'LineWidth', 2);
    axis([0 0.3 0 30])
    hold off
end
pause
close

%% Information Criterion
% AIC and BIC
info_par_exp_pos= zeros(n,4);
info_par_exp_neg= zeros(n,4);

for j=1:n
    xmin = parameters_pl_pos(j,2);
    r_exp = Rj(Rj(:, j)>0,j);
    r_par = Rj(Rj(:, j)>xmin,j);
    [aic_pl,bic_pl] = aicbic(parameters_pl_pos(j,3),2,length(r_par));
    info_par_exp_pos(j,[1 2]) =  [aic_pl,bic_pl];
    [aic_pl,bic_pl] = aicbic(parameters_pl_pos_exp(j,2),1,length(r_exp));
    info_par_exp_pos(j,[3 4]) =  [aic_pl,bic_pl];
end

for j=1:n
    xmin = parameters_pl_neg(j,2);
    r_exp = Rj2(Rj2(:, j)>0,j);
    r_par = Rj2(Rj2(:, j)>xmin,j);
    [aic_pl,bic_pl] = aicbic(parameters_pl_neg(j,3),2,length(r_par));
    info_par_exp_neg(j,[1 2]) =  [aic_pl,bic_pl];
    [aic_pl,bic_pl] = aicbic(parameters_pl_neg_exp(j,2),1,length(r_exp));
    info_par_exp_neg(j,[3 4]) =  [aic_pl,bic_pl];
end

% table AIC and BIC
info_pos_table = array2table(info_par_exp_pos,"VariableNames",["AIC PL","BIC PL", "AIC EXP","BIC EXP"]);
info_neg_table = array2table(info_par_exp_neg,"VariableNames",["AIC PL","BIC PL", "AIC EXP","BIC EXP"]);
info_neg_table = horzcat(table(linspace(1,n,n)', 'VariableNames',{'j'}), info_neg_table);
info_table = table ( info_neg_table , info_pos_table, 'VariableNames', {'Left', 'Right'});
disp(info_table);


%% KS test
% unique loop for k-stat right tail (test  + plot) 
kstat_pos = zeros(n,6);  

for j= 1:n
    
    alpha = parameters_pl_pos(j,1);
    xmin = parameters_pl_pos(j,2);
    x2 = Rj(Rj(:, j)>xmin,j);
    alpha_exp = parameters_pl_pos_exp(j,1); 
    u = rand(1000,1);
    y = plrnd(u,xmin,alpha); %generate rn from a desired continuos distribution by mean of inverse CDF apporch
    z = exprnd(alpha_exp,1000,1);
    
    [h, p, ksstat] = kstest2(y,x2); %accept h0
    kstat_pos(j,[1 2 3]) = [h, p, ksstat];
    [h, p, ksstat] = kstest2(z,x2); %reject h0
    kstat_pos(j,[4 5 6]) = [h, p, ksstat];
    
    % plot cdf pareto vs empirical cdf
    subplot(1, 2, 1);
    cdfplot(x2);
    hold on
    cdfplot(y)
    xlim([min(x2), max(x2)])
    hold off
    title(['CDF Pareto (Right tail)'  num2str(j)   "offset"], 'HorizontalAlignment', 'center')
    
    % plot cdf exponential vs empirical cdf
    subplot(1, 2, 2);
    cdfplot(x2)
    hold on
    cdfplot(z)
    xlim([min(x2), max(x2)])
    hold off
    title(['CDF Exp (Right tail)'  num2str(j)   "offset"], 'HorizontalAlignment', 'center')

    pause

end

% unique loop for k-stat left tail (test  + plot) 

kstat_neg = zeros(n,6); 
    
for j= 1:n
    offset= j ;
    alpha = parameters_pl_neg(j,1);
    xmin = parameters_pl_neg(j,2);
    x3 = Rj2(Rj2(:, j)>xmin,j);
    alpha_exp = parameters_pl_neg_exp(j,1); 
    u = rand(100000,1);
    y = plrnd(u,xmin,alpha);
    z = exprnd(alpha_exp,10000,1);
    
    [h, p, ksstat] = kstest2(y,x3); %accept h0
    kstat_neg(j,[1 2 3]) = [h, p, ksstat];
    [h, p, ksstat] = kstest2(z,x2); %reject h0
    kstat_neg(j,[4 5 6]) = [h, p, ksstat];

    % plot cdf pareto vs empirical cdf
    subplot(1, 2,1);
    cdfplot(x3);
    hold on
    cdfplot(y)
    xlim([min(x2), max(x2)])
    hold off
    title(['CDF Pareto (Left tail)' num2str(j)  "offset"], 'HorizontalAlignment', 'center');
    
    % plot cdf exponential vs empirical cdf
    subplot(1, 2, 2);
    cdfplot(x3)
    hold on
    cdfplot(z)
    xlim([min(x2), max(x2)])
    hold off
    title(['CDF Exp (Left tail)'  num2str(j)   "offset"], 'HorizontalAlignment', 'center')
    pause

end

% table
left_table = table(linspace(1,n,n)', kstat_neg(:,1), kstat_neg(:,4), kstat_neg(:,2),kstat_neg(:,5), 'VariableNames', {'j', 'h_pareto', 'p_pareto', 'h_exp', 'p_exp'});
right_table = table( kstat_pos(:,1), kstat_pos(:,4), kstat_pos(:,2),kstat_pos(:,5), 'VariableNames', { 'h_pareto', 'p_pareto', 'h_exp', 'p_exp'});
KS_table = table (left_table, right_table, 'VariableNames', {'Left', 'Right'});
disp(KS_table);


%% AVG
ks_a_pos = length(kstat_pos(:,1))
uffa = ks_a_pos - sum(kstat_pos(:,1))
uffi = ks_a_pos - sum(kstat_pos(:,4))
table(uffa, uffi,'VariableNames', {'PL better', 'Exp better'})


ks_a_pos = length(kstat_neg(:,1))
uff = ks_a_pos - sum(kstat_neg(:,1))
uf = ks_a_pos - sum(kstat_neg(:,4))
table(uff, uf,'VariableNames', {'PL better', 'Exp better'})




