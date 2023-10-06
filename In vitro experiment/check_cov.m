%% Define the DATA and its dosage level

i = 1;
DATA_temp = squeeze(DATA(:,i,:));
conc = Conc(i);
num_sub = 2;


%% Check the DATA correlation

load('Mixture_Data.mat')

DATA = squeeze(BF_41(:,1,:));
c = corrcoef(DATA);

load('MONOCLONAL_DATA.mat')
DATA = squeeze(SENSITIVE_500_BF(:,1,:));

d = corrcoef(DATA);
%% Sto
theta = opt_xx_sto;
p_vec    = theta(1:6:6*num_sub-6);
p_vec    = [p_vec',1-sum(p_vec)];
theta_j = p_vec;
for s = 1:num_sub % Computing the birth rate and death rate under drug d.
   if s ~=num_sub
       r_s  = theta(6*s-4);
       nu_s = theta(6*s-3);  
       b_s  = theta(6*s-2);
       E_s  = theta(6*s-1);
       n_s  = theta(6*s);
       nu_s = nu_s - log(b_s + (1 - b_s)/(1 + (conc/ E_s)^n_s));
   else
       r_s  = theta(6*s-5);
       nu_s = theta(6*s-4);  
       b_s  = theta(6*s-3);
       E_s  = theta(6*s-2);
       n_s  = theta(6*s-1);
       nu_s = nu_s - log(b_s + (1 - b_s)/(1 + (conc/ E_s)^n_s));
   end
   theta_j = [theta_j,[r_s,nu_s]];
end
lam    = [];
r_nu   = [];
for i = 1:num_sub
    r_nu = [r_nu,theta_j(num_sub + 2*i - 1) + theta_j(num_sub + 2*i)];
    lam  = [lam,theta_j(num_sub + 2*i - 1) - theta_j(num_sub + 2*i)];
end

p      = theta_j(1:num_sub);
X_init = round(mean(DATA_temp(:,1)) * p);

theo_Var_sto  = get_Var(X_init, Time, lam, r_nu,0);

emp_Var = cov(DATA_temp);

%% Dyn

theta = opt_xx_dyn;
p_vec    = theta(1:6:6*num_sub-6);
p_vec    = [p_vec',1-sum(p_vec)];
theta_j = p_vec;
for s = 1:num_sub % Computing the birth rate and death rate under drug d.
   if s ~=num_sub
       r_s  = theta(6*s-4);
       nu_s = theta(6*s-3);  
       b_s  = theta(6*s-2);
       E_s  = theta(6*s-1);
       n_s  = theta(6*s);
       nu_s = nu_s - log(b_s + (1 - b_s)/(1 + (conc/ E_s)^n_s));
   else
       r_s  = theta(6*s-5);
       nu_s = theta(6*s-4);  
       b_s  = theta(6*s-3);
       E_s  = theta(6*s-2);
       n_s  = theta(6*s-1);
       nu_s = nu_s - log(b_s + (1 - b_s)/(1 + (conc/ E_s)^n_s));
   end
   theta_j = [theta_j,[r_s,nu_s]];
end
lam    = [];
r_nu   = [];
for i = 1:num_sub
    r_nu = [r_nu,theta_j(num_sub + 2*i - 1) + theta_j(num_sub + 2*i)];
    lam  = [lam,theta_j(num_sub + 2*i - 1) - theta_j(num_sub + 2*i)];
end

p      = theta_j(1:num_sub);
X_init = round(mean(DATA_temp(:,1)) * p);

theo_Var_dyn  = diag(diag(get_Var(X_init, Time, lam, r_nu,0)));

emp_Var = cov(DATA_temp);

