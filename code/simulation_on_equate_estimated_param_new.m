function [time_equate] = simulation_on_equate_estimated_param_new(professor_year,asso_professor_year,assi_professor_year,female_ratio,attrition_rate_m,attrition_rate_f,promotion_rate_f,promotion_rate_m,all_level,same_promotion)
k = 3; % number of levels in the organization
n0_f = [assi_professor_year(19,2),asso_professor_year(19,2),professor_year(19,2),0]; %, # female in each level, the last element is attrition
n0_m = [assi_professor_year(19,1),asso_professor_year(19,1),professor_year(19,1),0]; % # male
% max(n0_m./n0_f);
% size = sum(n0_m+n0_f)
% bottom = n0_m(1)/n0_f(1)
% next_bottom = n0_m(2)/n0_f(2)
% top = n0_m(k)/n0_f(k)
% promotion rate and attrition rate at different levels


% initialize state dynamics
p_m = zeros(k+1,k+1);
p_f = zeros(k+1,k+1);

if same_promotion==0
    for i =1:k-1
        p_m(i,i+1)=promotion_rate_m(i);
        p_m(i,k+1)=attrition_rate_m(i);
        p_m(i,i) = 1-attrition_rate_m(i)-promotion_rate_m(i);
        p_f(i,i+1)=promotion_rate_f(i);
        p_f(i,k+1)=attrition_rate_f(i);
        p_f(i,i) = 1-attrition_rate_f(i)-promotion_rate_f(i);
    end
    p_m(k,k+1) = attrition_rate_m(k);
    p_m(k,k) = 1-attrition_rate_m(k);
    p_m(k+1,k+1) = 1;

    p_f(k,k+1) = attrition_rate_f(k);
    p_f(k,k) = 1-attrition_rate_f(k);
    p_f(k+1,k+1) = 1;
else
    promotion_rate = (promotion_rate_m+promotion_rate_f)/2;
    for i =1:k-1
        p_m(i,i+1)=promotion_rate(i);
        p_m(i,k+1)=attrition_rate_m(i);
        p_m(i,i) = 1-attrition_rate_m(i)-promotion_rate(i);
        p_f(i,i+1)=promotion_rate(i);
        p_f(i,k+1)=attrition_rate_f(i);
        p_f(i,i) = 1-attrition_rate_f(i)-promotion_rate(i);
    end
    p_m(k,k+1) = attrition_rate_m(k);
    p_m(k,k) = 1-attrition_rate_m(k);
    p_m(k+1,k+1) = 1;

    p_f(k,k+1) = attrition_rate_f(k);
    p_f(k,k) = 1-attrition_rate_f(k);
    p_f(k+1,k+1) = 1;
end
% initialize new hiring proportion in level and sex, only hire from the bottom

v_f = zeros(1,k+1); % female hire rate
v_m = zeros(1,k+1);


if all_level==0
    
    v_f(1) = female_ratio;
    v_m(1) = 1-female_ratio;
else
    for i =1:k
        v_f(i) = female_ratio*(n0_f(i)+n0_m(i))/(sum(n0_f)+sum(n0_m));
        v_m(i) = (1-female_ratio)*(n0_f(i)+n0_m(i))/(sum(n0_f)+sum(n0_m));
    end
end

% should fix sum(v_f)+sum(v_m)=1

%simulation time length
T=1000;
nt_f = zeros(T,k+1);
nt_m = zeros(T,k+1);

nt_f(1,:) = n0_f;
nt_m(1,:) = n0_m;

size_target = (sum(n0_f)+sum(n0_m))*ones(T+1,1);

time_equate = T+1;
size_organization = zeros(T+1,1);
total_male = zeros(T+1,1);
total_female = zeros(T+1,1);
number_leave_m = zeros(T+1,1);
number_leave_f = zeros(T+1,1);

size_organization(1) = sum(n0_f)+sum(n0_m);
total_male(1) = sum(n0_m);
total_female(1) = sum(n0_f);
number_leave_m(1) = 0;
number_leave_f(1) = 0;

for t = 1:T
    
    number_leave_m(t+1) = sum((nt_m(t,1:k).*attrition_rate_m));
    number_leave_f(t+1) = sum((nt_f(t,1:k).*attrition_rate_f));
    
    n_f = nt_f(t,:)*p_f+v_f*(number_leave_f(t+1)+(size_target(t+1)-size_target(t)));
    n_m = nt_m(t,:)*p_m+v_m*(number_leave_m(t+1)+(size_target(t+1)-size_target(t)));
    nt_f(t+1,:) = n_f;
    nt_m(t+1,:) = n_m;
    if nt_f(t+1,1)+1>=nt_m(t+1,1) && nt_f(t+1,2)+1>=nt_m(t+1,2) && nt_f(t+1,3)+1>=nt_m(t+1,3) && t<=time_equate
        time_equate = t;
    end
    size_organization(t+1) = sum(nt_f(t+1,1:k))+sum(nt_m(t+1,1:k));
    total_male(t+1) = sum(nt_m(t+1,1:k));
    total_female(t+1) = sum(nt_f(t+1,1:k));
    
end
% figure();
% plot(nt_f(1:time_equate,1),'--')
% hold on
% plot(nt_f(1:time_equate,2),'--')
% plot(nt_f(1:time_equate,3),'--')
% 
% plot(nt_m(1:time_equate,1))
% plot(nt_m(1:time_equate,2))
% plot(nt_m(1:time_equate,3))
% 
% legend('f1','f2','f3','m1','m2','m3')
% 
% disp(time_equate)
% disp(v_f)
% disp(v_m)
% disp(p_f)
% disp(p_m)

end

