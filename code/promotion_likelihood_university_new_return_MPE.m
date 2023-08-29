function [diff,pseudo_R2]= promotion_likelihood_university_new_return_MPE(x,start_year,end_year)
% new promotion rate and attrition rate estimation, with attrition rate
% estimate for each level and each gender. 
% the order of x:[promotion rate male:12,promotion rate female:34,attrition rate male:567,
% attrition rate female:8910,ratio of new hired tenured as associate:11]

global assi_professor_year asso_professor_year professor_year nh_tenured nh_on_tenure_track 

nt_f_0 = [assi_professor_year(:,2) asso_professor_year(:,2) professor_year(:,2)];
nt_m_0 = [assi_professor_year(:,1) asso_professor_year(:,1) professor_year(:,1)];


k = 3; % number of levels in the organization
n0_f = [nt_f_0(1,1:3) 0]; %, # female in each level, the last element is attrition
%n0_all = [,0];
%n0_m = n0_all-n0_f; % # male
n0_m = [nt_m_0(1,1:3),0];

%simulation time length
T=19;
nt_f = zeros(T,k+1);
nt_m = zeros(T,k+1);

nt_f(1,:) = n0_f;
nt_m(1,:) = n0_m;


% promotion rate and attrition rate at different levels

attrition_rate_m = [x(5),x(5),x(5)]; % attrition rate for male
attrition_rate_f = [x(8),x(8),x(8)];
promotion_rate_m = [x(1),x(2)];
lambda = 1;%x(11);


% initialize state dynamics
p_m = zeros(k+1,k+1);

for i =1:k-1
    p_m(i,i+1)=promotion_rate_m(i);%-0.001*i;
    p_m(i,k+1)=attrition_rate_m(i);
    p_m(i,i) = 1-attrition_rate_m(i)-p_m(i,i+1);
end
p_m(k,k+1) = attrition_rate_m(k);
p_m(k,k) = 1-attrition_rate_m(k);
p_m(k+1,k+1) = 1;

%alpha = x(3);
promotion_rate_f = [x(3) x(4)];%promotion_rate*alpha;
p_f = zeros(k+1,k+1);

for i =1:k-1
    p_f(i,i+1)=promotion_rate_f(i);%-0.001*i;
    p_f(i,k+1)=attrition_rate_f(i);
    p_f(i,i) = 1-attrition_rate_f(i)-p_f(i,i+1);
end
p_f(k,k+1) = attrition_rate_f(k);
p_f(k,k) = 1-attrition_rate_f(k);
p_f(k+1,k+1) = 1;





%% run simulation
% assume we only hire female workers from the bottom, promotion and
% attrition dynamics same for male and female;
% time_equate = T+1;
size_organization = zeros(T+1,1);
total_male = zeros(T+1,1);
total_female = zeros(T+1,1);
% number_leave = zeros(T+1,1);

size_organization(1) = sum(n0_f)+sum(n0_m);
total_male(1) = sum(n0_m);
total_female(1) = sum(n0_f);
% number_leave(1) = 0;

for t = start_year:end_year-1

%     number_leave(t+1) = sum((nt_f(t,1:k).*attrition_rate)+(nt_m(t,1:k).*attrition_rate));
    
    %n_f = nt_f(t,:)*p_f+v_f*(number_leave(t+1)+(size_target(t+1)-size_target(t)));
    %n_m = nt_m(t,:)*p_m+v_m*(number_leave(t+1)+(size_target(t+1)-size_target(t)));
    
    
    n_f = [nt_f_0(t,:),0]*p_f+[nh_on_tenure_track(t,2),lambda*nh_tenured(t,2),(1-lambda)*nh_tenured(t,2),0];
    n_m = [nt_m_0(t,:),0]*p_m+[nh_on_tenure_track(t,1),lambda*nh_tenured(t,1),(1-lambda)*nh_tenured(t,1),0];
    
%     n_f = [nt_f_0(t,:),0]*p_f+[nh_on_tenure_track(t+1,2),nh_on_tenure_track(t+1,2),0,0];
%     n_m = [nt_f_0(t,:),0]*p_m+[nh_tenured(t+1,1),nh_tenured(t+1,1),0,0];
     
    nt_f(t+1,:) = n_f;
    nt_m(t+1,:) = n_m;

    size_organization(t+1) = sum(nt_f(t+1,1:k))+sum(nt_m(t+1,1:k));
    total_male(t+1) = sum(nt_m(t+1,1:k));
    total_female(t+1) = sum(nt_f(t+1,1:k));
end



diff_m = (nt_m_0(start_year:end_year,:)-nt_m(start_year:end_year,1:k))./nt_m_0(start_year:end_year,:);
diff_f = (nt_f_0(start_year:end_year,:)-nt_f(start_year:end_year,1:k))./nt_f_0(start_year:end_year,:);

diff = [diff_m diff_f];

diff = mean(mean(abs(diff)));

data_m = nt_m(start_year:end_year,1:k);
data_f = nt_f(start_year:end_year,1:k);


error_model_m = sum(sum((nt_m_0(start_year:end_year,:)-nt_m(start_year:end_year,1:k)).^2));
error_model_f = sum(sum((nt_f_0(start_year:end_year,:)-nt_f(start_year:end_year,1:k)).^2));
error_average_m = sum(sum((data_m-mean(data_m)).^2));
error_average_f = sum(sum((data_f-mean(data_f)).^2));

pseudo_R2 = [1-error_model_m/error_average_m 1-error_model_f/error_average_f];
% plot(nt_m_0,'--')
% hold on
% plot(nt_m(1:T,1:k))
% legend('1','2','3','11','22','33')
end