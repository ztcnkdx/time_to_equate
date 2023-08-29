clear;

% arank=1: professor
% arank=2: associate professor
% ARANK=3: assistant professor
% HRTOTLT: grand total
% HRTOTLM: grand total men
% HRTOTLW: grand total women

university_list = readtable('university_list_r1_r2.csv');

% harvard deleted due to irregular data/tenure system

%
table01 = readtable('s2001_f.csv');
table02 = readtable('s2002_f.csv');
table03 = readtable('s2003_f.csv');
table04 = readtable('s2004_f.csv');
table05 = readtable('s2005_f.csv');
table06 = readtable('s2006_f.csv');
table07 = readtable('s2007_f.csv');
table08 = readtable('s2008_f.csv');
table09 = readtable('s2009_f.csv');
table10 = readtable('s2010_f.csv');
table11 = readtable('s2011_f.csv');

tables01 = {table01,table02,table03,table04,table05,table06,table07,table08,table09,table10,table11};

table12 = readtable('s2012_is.csv');
table13 = readtable('s2013_is.csv');
table14 = readtable('s2014_is.csv');
table15 = readtable('s2015_is.csv');
table16 = readtable('s2016_is.csv');
table17 = readtable('s2017_is.csv');
table18 = readtable('s2018_is.csv');
table19 = readtable('s2019_is.csv');

tables = {table01,table02,table03,table04,table05,table06,table07,table08,table09,table10,table11,table12,table13,table14,table15,table16,table17,table18,table19};

nh01 = readtable('s2001_g.csv');
nh02 = readtable('s2002_g.csv');
nh03 = readtable('s2003_g.csv');
nh04 = readtable('s2004_g.csv');
nh05 = readtable('s2005_g.csv');
nh06 = readtable('s2006_g.csv');
nh07 = readtable('s2007_g.csv');
nh08 = readtable('s2008_g.csv');
nh09 = readtable('s2009_g.csv');
nh10 = readtable('s2010_g.csv');
nh11 = readtable('s2011_g.csv');

nhtables01 = {nh01,nh02,nh03,nh04,nh05,nh06,nh07,nh08,nh09,nh10,nh11};

nh12 = readtable('s2012_nh.csv');
nh13 = readtable('s2013_nh.csv');
nh14 = readtable('s2014_nh.csv');
nh15 = readtable('s2015_nh.csv');
nh16 = readtable('s2016_nh.csv');
nh17 = readtable('s2017_nh.csv');
nh18 = readtable('s2018_nh.csv');
nh19 = readtable('s2019_nh.csv');

nhtables = {nh01,nh02,nh03,nh04,nh05,nh06,nh07,nh08,nh09,nh10,nh11,nh12,nh13,nh14,nh15,nh16,nh17,nh18,nh19};



T01 = 11;
T12 = 8;
T = T01+T12;

professor_year_pool = zeros(T,2);
asso_professor_year_pool = zeros(T,2);
assi_professor_year_pool = zeros(T,2);
nh_on_tenure_track_pool = zeros(T,2);

expected_time_to_tenure = zeros(T,2);


global assi_professor_year asso_professor_year professor_year nh_on_tenure_track nh_tenured
start_year = 1;
end_year = 19;
time_to_equate = zeros(height(university_list),4);
best_x_list = zeros(height(university_list),11);
min_fval_list = zeros(height(university_list),1);
nh_tenure_m = zeros(height(university_list),1);
nh_tenure_f = zeros(height(university_list),1);
nh_tenure_track_m = zeros(height(university_list),1);
nh_tenure_track_f = zeros(height(university_list),1);
available_unit = zeros(height(university_list),1);
no_data_unit = zeros(height(university_list),1);
pseudo_R2_all = zeros(height(university_list),2);
diff_all = zeros(height(university_list),1);

attrition_rate_m_list = zeros(height(university_list),3);
attrition_rate_f_list = zeros(height(university_list),3);
promotion_rate_m_list = zeros(height(university_list),2);
promotion_rate_f_list = zeros(height(university_list),2);

sim_initial_state_m_list = zeros(height(university_list),3);
sim_initial_state_f_list = zeros(height(university_list),3);

for i=1:height(university_list)
    disp(i)
    try
        professor_year = zeros(T,2); % male in column 1, female in column 2
        asso_professor_year = zeros(T,2);
        assi_professor_year = zeros(T,2);
        
        nh_tenured = zeros(T,2);
        nh_on_tenure_track = zeros(T,2);
        
        unitid = table2array(university_list(i,1));
        
        
        
        for j=1:T
            
            if j>11 % for year 2012-2019
                table_temp = tables{j};
                unit = table_temp(table_temp.UNITID==unitid,:);
                
                nh_temp = nhtables{j};
                unit_nh = nh_temp(nh_temp.UNITID==unitid,:);
                
                if height(unit)>10 && height(unit_nh)>3
                    professor_year(j,1) = unit(logical((unit.FACSTAT==20).*(unit.ARANK==1)),:).HRTOTLM;
                    professor_year(j,2) = unit(logical((unit.FACSTAT==20).*(unit.ARANK==1)),:).HRTOTLW;
                    try
                        asso_professor_year(j,1) = unit(logical((unit.FACSTAT==20).*(unit.ARANK==2)),:).HRTOTLM;
                        asso_professor_year(j,2) = unit(logical((unit.FACSTAT==20).*(unit.ARANK==2)),:).HRTOTLW;
                    catch
                        asso_professor_year(j,1) = unit(logical((unit.FACSTAT==30).*(unit.ARANK==2)),:).HRTOTLM;
                        asso_professor_year(j,2) = unit(logical((unit.FACSTAT==30).*(unit.ARANK==2)),:).HRTOTLW;
                    end
                    try
                        assi_professor_year(j,1) = unit(logical((unit.FACSTAT==30).*(unit.ARANK==3)),:).HRTOTLM;
                        assi_professor_year(j,2) = unit(logical((unit.FACSTAT==30).*(unit.ARANK==3)),:).HRTOTLW;
                    catch
                        assi_professor_year(j,1) = unit(logical((unit.FACSTAT==20).*(unit.ARANK==3)),:).HRTOTLM;
                        assi_professor_year(j,2) = unit(logical((unit.FACSTAT==20).*(unit.ARANK==3)),:).HRTOTLW;
                    end
                    
                    try
                        
                        nh_tenured(j,1) = unit_nh(logical((unit_nh.FACSTAT==20)),:).HRTOTLM;
                        nh_tenured(j,2) = unit_nh(logical((unit_nh.FACSTAT==20)),:).HRTOTLW;
                        
                    catch
                        nh_tenured(j,1)=0;
                        nh_tenured(j,2)=0;
                    end
                    
                    try
                        nh_on_tenure_track(j,1) = unit_nh(logical((unit_nh.FACSTAT==30)),:).HRTOTLM;
                        nh_on_tenure_track(j,2) = unit_nh(logical((unit_nh.FACSTAT==30)),:).HRTOTLW;
                    catch
                        nh_on_tenure_track(j,1)=0;
                        nh_on_tenure_track(j,2)=0;
                    end
                end
            elseif j>7
                table_temp = tables{j};
                unit = table_temp(table_temp.UNITID==unitid,:);
                
                nh_temp = nhtables{j};
                unit_nh = nh_temp(nh_temp.UNITID==unitid,:);
                
                if height(unit)>10 && height(unit_nh)>3
                    professor_year(j,1) = unit(logical((unit.ARANK==1)),:).HRTOTLM;
                    professor_year(j,2) = unit(logical((unit.ARANK==1)),:).HRTOTLW;
                    try
                        asso_professor_year(j,1) = unit(logical((unit.ARANK==2)),:).HRTOTLM;
                        asso_professor_year(j,2) = unit(logical((unit.ARANK==2)),:).HRTOTLW;
                    catch
                        asso_professor_year(j,1) = unit(logical((unit.ARANK==9)),:).HRTOTLM;
                        asso_professor_year(j,2) = unit(logical((unit.ARANK==9)),:).HRTOTLW;
                    end
                    try %for assistant profs, we first try get data for assistant on tenure track, then get for tenured
                        assi_professor_year(j,1) = unit(logical((unit.ARANK==10)),:).HRTOTLM;
                        assi_professor_year(j,2) = unit(logical((unit.ARANK==10)),:).HRTOTLW;
                    catch
                        assi_professor_year(j,1) = unit(logical((unit.ARANK==3)),:).HRTOTLM;
                        assi_professor_year(j,2) = unit(logical((unit.ARANK==3)),:).HRTOTLW;
                    end
                    
                    try
                        
                        nh_tenured(j,1) = unit_nh(logical((unit_nh.FUNCTCD==1)),:).HRTOTLM;
                        nh_tenured(j,2) = unit_nh(logical((unit_nh.FUNCTCD==1)),:).HRTOTLW;
                        
                    catch
                        nh_tenured(j,1)=0;
                        nh_tenured(j,2)=0;
                    end
                    
                    try
                        nh_on_tenure_track(j,1) = unit_nh(logical((unit_nh.FUNCTCD==2)),:).HRTOTLM;
                        nh_on_tenure_track(j,2) = unit_nh(logical((unit_nh.FUNCTCD==2)),:).HRTOTLW;
                    catch
                        nh_on_tenure_track(j,1)=0;
                        nh_on_tenure_track(j,2)=0;
                    end
                end
            elseif j>3 % for year 2004-2007
                table_temp = tables{j};
                unit = table_temp(table_temp.UNITID==unitid,:);
                
                nh_temp = nhtables{j};
                unit_nh = nh_temp(nh_temp.UNITID==unitid,:);
                
                if height(unit)>10 && height(unit_nh)>3
                    professor_year(j,1) = unit(logical((unit.ARANK==1)),:).STAFF15;
                    professor_year(j,2) = unit(logical((unit.ARANK==1)),:).STAFF16;
                    try %for associate profs, we first try get data for tenured associate, then get for on tenure track associate
                        asso_professor_year(j,1) = unit(logical((unit.ARANK==2)),:).STAFF15;
                        asso_professor_year(j,2) = unit(logical((unit.ARANK==2)),:).STAFF16;
                    catch
                        asso_professor_year(j,1) = unit(logical((unit.ARANK==9)),:).STAFF15;
                        asso_professor_year(j,2) = unit(logical((unit.ARANK==9)),:).STAFF16;
                    end
                    try %for assistant profs, we first try get data for assistant on tenure track, then get for tenured
                        assi_professor_year(j,1) = unit(logical((unit.ARANK==10)),:).STAFF15;
                        assi_professor_year(j,2) = unit(logical((unit.ARANK==10)),:).STAFF16;
                    catch
                        assi_professor_year(j,1) = unit(logical((unit.ARANK==3)),:).STAFF15;
                        assi_professor_year(j,2) = unit(logical((unit.ARANK==3)),:).STAFF16;
                    end
                    
                    try
                        
                        nh_tenured(j,1) = unit_nh(logical((unit_nh.FUNCTCD==1)),:).STAFF15;
                        nh_tenured(j,2) = unit_nh(logical((unit_nh.FUNCTCD==1)),:).STAFF16;
                        
                    catch
                        nh_tenured(j,1)=0;
                        nh_tenured(j,2)=0;
                    end
                    
                    try
                        nh_on_tenure_track(j,1) = unit_nh(logical((unit_nh.FUNCTCD==2)),:).STAFF15;
                        nh_on_tenure_track(j,2) = unit_nh(logical((unit_nh.FUNCTCD==2)),:).STAFF16;
                    catch
                        nh_on_tenure_track(j,1)=0;
                        nh_on_tenure_track(j,2)=0;
                    end
                end
            else % year 2001-2003
                table_temp = tables{j};
                unit = table_temp(table_temp.unitid==unitid,:);
                
                nh_temp = nhtables{j};
                unit_nh = nh_temp(nh_temp.unitid==unitid,:);
                
                if height(unit)>10 && height(unit_nh)>3
                    professor_year(j,1) = unit(logical((unit.arank==1)),:).staff15;
                    professor_year(j,2) = unit(logical((unit.arank==1)),:).staff16;
                    try
                        asso_professor_year(j,1) = unit(logical((unit.arank==2)),:).staff15;
                        asso_professor_year(j,2) = unit(logical((unit.arank==2)),:).staff16;
                    catch
                        asso_professor_year(j,1) = unit(logical((unit.arank==9)),:).staff15;
                        asso_professor_year(j,2) = unit(logical((unit.arank==9)),:).staff16;
                    end
                    try %for assistant profs, we first try get data for assistant on tenure track, then get for tenured
                        assi_professor_year(j,1) = unit(logical((unit.arank==10)),:).staff15;
                        assi_professor_year(j,2) = unit(logical((unit.arank==10)),:).staff16;
                    catch
                        assi_professor_year(j,1) = unit(logical((unit.arank==3)),:).staff15;
                        assi_professor_year(j,2) = unit(logical((unit.arank==3)),:).staff16;
                    end
                    
                    try
                        
                        nh_tenured(j,1) = unit_nh(logical((unit_nh.functcd==1)),:).staff15;
                        nh_tenured(j,2) = unit_nh(logical((unit_nh.functcd==1)),:).staff16;
                        
                    catch
                        nh_tenured(j,1)=0;
                        nh_tenured(j,2)=0;
                    end
                    
                    try
                        nh_on_tenure_track(j,1) = unit_nh(logical((unit_nh.functcd==2)),:).staff15;
                        nh_on_tenure_track(j,2) = unit_nh(logical((unit_nh.functcd==2)),:).staff16;
                    catch
                        nh_on_tenure_track(j,1)=0;
                        nh_on_tenure_track(j,2)=0;
                    end
                end
                
            end
            
        end
        
        nh_tenure_m(i) = sum(nh_tenured(:,1));
        nh_tenure_f(i) = sum(nh_tenured(:,2));
        
        nh_tenure_track_m(i) = sum(nh_on_tenure_track(:,1));
        nh_tenure_track_f(i) = sum(nh_on_tenure_track(:,2));
        
        flag = sum(sum((professor_year>0)==0));%+sum(sum((nh_on_tenure_track>0)==0));
        
        
        if flag<=8    % set a threshold, if there are more than 4 years (both gender, so 8 entries) of missing data, then we don't use the data for this university
            
            
            if flag>0  %linearly interpolate missing values
                nh_tenured(nh_tenured == 0) = NaN;
                nh_tenured = fillmissing(nh_tenured,'linear');
                
                nh_on_tenure_track(nh_on_tenure_track == 0) = NaN;
                nh_on_tenure_track = fillmissing(nh_on_tenure_track,'linear');
                
                professor_year(professor_year == 0) = NaN;
                professor_year = fillmissing(professor_year,'linear');
                
                assi_professor_year(assi_professor_year == 0) = NaN;
                assi_professor_year = fillmissing(assi_professor_year,'linear');
                
                asso_professor_year(asso_professor_year == 0) = NaN;
                asso_professor_year = fillmissing(asso_professor_year,'linear');
                
            end
            
            professor_year_pool = professor_year_pool+professor_year;
            asso_professor_year_pool = asso_professor_year_pool+asso_professor_year;
            assi_professor_year_pool = assi_professor_year_pool+assi_professor_year;
            %nh_tenured_pool = nh_tenured_pool+nh_tenured;
            nh_on_tenure_track_pool = nh_on_tenure_track_pool+nh_on_tenure_track;
            
            col_count = 1;
            col_count_equate = 1;
            col_count_a = 1;
            for seperate = 1:1%4:15
                for seperate_2 = 1:1%1:2
                    % if seperate_2==1
                    %     start_year = 1;
                    %     end_year = 19;%seperate;
                    % else
                    %     start_year = 1;%seperate;
                    %     end_year = 19;
                    % end
                    
                    min_fval = 1000000000000;
                    
                    for ii = 1:10
                        
                        fun = @(x)promotion_likelihood_university_new(x,start_year,end_year);
                        
                        x0 = rand(1,11);
                        
                        options = optimoptions('fmincon','MaxFunctionEvaluations',10000,'StepTolerance', 1e-12,'Display','off');
                        
                        [x,fval] = fmincon(fun,x0,[],[],[],[],zeros(1,11),ones(1,11),[],options);
                        x;
                        if fval<min_fval
                            best_x0 = x0;
                            min_fval = fval;
                            best_x = x;
                        end
                        
                    end
                    
                    
                    min_fval_list(i,:) = min_fval;
                    best_x_list(i,:) = best_x;
                    
                    attrition_rate_m = [best_x(5) best_x(5) best_x(5)];
                    attrition_rate_f = [best_x(8) best_x(8) best_x(8)];
                    promotion_rate_m = [best_x(1) best_x(2)];
                    promotion_rate_f = [best_x(3) best_x(4)];
                    
                    expected_time_to_tenure(i,1) = first_passage_time(1000,promotion_rate_m(1),attrition_rate_m(1));
                    expected_time_to_tenure(i,2) = first_passage_time(1000,promotion_rate_f(1),attrition_rate_f(1));
                    
                    [diff,pseudo_R2] = promotion_likelihood_university_new_return_MPE(best_x,start_year,end_year);
                    diff_all(i,1) = diff;
                    if pseudo_R2(1)<0 || pseudo_R2(2)<0
                        disp(i);      % 26 92 102 126 138 159 174 199 230
                    end
                    pseudo_R2_all(i,:) = pseudo_R2;
                    available_unit(i,1) = unitid;
                    
                    attrition_rate_m_list(i,col_count_a:col_count_a+2) = attrition_rate_m;
                    attrition_rate_f_list(i,col_count_a:col_count_a+2) = attrition_rate_f;
                    promotion_rate_m_list(i,col_count:col_count+1) = promotion_rate_m;
                    promotion_rate_f_list(i,col_count:col_count+1) = promotion_rate_f;
                    col_count = col_count+2;
                    col_count_a = col_count_a+3;
                    
                    
                    all_level=0;
                    female_ratio = 1;
                    same_promotion = 0;
                    sim_initial_state_m_list(i,:) = [assi_professor_year(19,1),asso_professor_year(19,1),professor_year(19,1)];
                    sim_initial_state_f_list(i,:) = [assi_professor_year(19,2),asso_professor_year(19,2),professor_year(19,2)];


                    sim_result = simulation_on_equate_estimated_param_new(professor_year,asso_professor_year,assi_professor_year,female_ratio,attrition_rate_m,attrition_rate_f,promotion_rate_f,promotion_rate_m,all_level,same_promotion);
                    time_to_equate(i,col_count_equate) = sim_result;
                    
                    
                    
                    female_ratio = 0.75;
                    
                    sim_result = simulation_on_equate_estimated_param_new(professor_year,asso_professor_year,assi_professor_year,female_ratio,attrition_rate_m,attrition_rate_f,promotion_rate_f,promotion_rate_m,all_level,same_promotion);
                    time_to_equate(i,col_count_equate+1) = sim_result;
                    
                    all_level=1;
                    female_ratio = 1;
                    sim_result = simulation_on_equate_estimated_param_new(professor_year,asso_professor_year,assi_professor_year,female_ratio,attrition_rate_m,attrition_rate_f,promotion_rate_f,promotion_rate_m,all_level,same_promotion);
                    time_to_equate(i,col_count_equate+2) = sim_result;
                    
                    
                    
                    female_ratio = 0.75;
                    sim_result = simulation_on_equate_estimated_param_new(professor_year,asso_professor_year,assi_professor_year,female_ratio,attrition_rate_m,attrition_rate_f,promotion_rate_f,promotion_rate_m,all_level,same_promotion);
                    time_to_equate(i,col_count_equate+3) = sim_result;

                    same_promotion = 1;
                    all_level=0;
                    female_ratio = 1;
                    sim_result = simulation_on_equate_estimated_param_new(professor_year,asso_professor_year,assi_professor_year,female_ratio,attrition_rate_m,attrition_rate_f,promotion_rate_f,promotion_rate_m,all_level,same_promotion);
                    time_to_equate(i,col_count_equate+4) = sim_result;
                    
                    
                    
                    female_ratio = 0.75;
                    
                    sim_result = simulation_on_equate_estimated_param_new(professor_year,asso_professor_year,assi_professor_year,female_ratio,attrition_rate_m,attrition_rate_f,promotion_rate_f,promotion_rate_m,all_level,same_promotion);
                    time_to_equate(i,col_count_equate+5) = sim_result;
                    
                    all_level=1;
                    female_ratio = 1;
                    sim_result = simulation_on_equate_estimated_param_new(professor_year,asso_professor_year,assi_professor_year,female_ratio,attrition_rate_m,attrition_rate_f,promotion_rate_f,promotion_rate_m,all_level,same_promotion);
                    time_to_equate(i,col_count_equate+6) = sim_result;
                    
                    
                    
                    female_ratio = 0.75;
                    sim_result = simulation_on_equate_estimated_param_new(professor_year,asso_professor_year,assi_professor_year,female_ratio,attrition_rate_m,attrition_rate_f,promotion_rate_f,promotion_rate_m,all_level,same_promotion);
                    time_to_equate(i,col_count_equate+7) = sim_result;

                    all_level=1;
                    female_ratio = 0.5;

                    sim_result = simulation_on_equate_estimated_param_new(professor_year,asso_professor_year,assi_professor_year,female_ratio,attrition_rate_m,attrition_rate_f,promotion_rate_f,promotion_rate_m,all_level,same_promotion);
                    time_to_equate(i,col_count_equate+8) = sim_result;

                    all_level=0;
                    female_ratio = 0.5;

                    sim_result = simulation_on_equate_estimated_param_new(professor_year,asso_professor_year,assi_professor_year,female_ratio,attrition_rate_m,attrition_rate_f,promotion_rate_f,promotion_rate_m,all_level,same_promotion);
                    time_to_equate(i,col_count_equate+9) = sim_result;


                    col_count_equate = col_count_equate+10;
                    
                    
                end
            end
            
        end
    catch
        no_data_unit(i,1) = unitid;
    end
    
    
    
end
professor_year = professor_year_pool;
asso_professor_year = asso_professor_year_pool;
assi_professor_year = assi_professor_year_pool;
nh_on_tenure_track = nh_on_tenure_track_pool;
nh_tenured = zeros(T,2);
min_fval_pool = 1000000000000;

for ii = 1:10
    
    fun = @(x)promotion_likelihood_university_new(x,start_year,end_year);
    
    x0 = rand(1,11);
    
    options = optimoptions('fmincon','MaxFunctionEvaluations',10000,'StepTolerance', 1e-12,'Display','off');
    
    [x,fval] = fmincon(fun,x0,[],[],[],[],zeros(1,11),ones(1,11),[],options);
    x;
    if fval<min_fval_pool
        best_x0 = x0;
        min_fval_pool = fval;
        best_x = x;
    end
    
end



attrition_rate_m_pool = [best_x(5) best_x(5) best_x(5)];
attrition_rate_f_pool = [best_x(8) best_x(8) best_x(8)];
promotion_rate_m_pool = [best_x(1) best_x(2)];
promotion_rate_f_pool = [best_x(3) best_x(4)];

[diff_pool,pseudo_R2_pool] = promotion_likelihood_university_new_return_MPE(best_x,start_year,end_year);


time_to_equate_pool = zeros(1,10);

all_level=0;
female_ratio = 1;
same_promotion = 0;
sim_result = simulation_on_equate_estimated_param_new(professor_year,asso_professor_year,assi_professor_year,female_ratio,attrition_rate_m_pool,attrition_rate_f_pool,promotion_rate_f_pool,promotion_rate_m_pool,all_level,same_promotion);
time_to_equate_pool(1) = sim_result;



female_ratio = 0.75;
sim_result = simulation_on_equate_estimated_param_new(professor_year,asso_professor_year,assi_professor_year,female_ratio,attrition_rate_m_pool,attrition_rate_f_pool,promotion_rate_f_pool,promotion_rate_m_pool,all_level,same_promotion);
time_to_equate_pool(2) = sim_result;

all_level=1;
female_ratio = 1;
sim_result = simulation_on_equate_estimated_param_new(professor_year,asso_professor_year,assi_professor_year,female_ratio,attrition_rate_m_pool,attrition_rate_f_pool,promotion_rate_f_pool,promotion_rate_m_pool,all_level,same_promotion);
time_to_equate_pool(3) = sim_result;



female_ratio = 0.75;
sim_result = simulation_on_equate_estimated_param_new(professor_year,asso_professor_year,assi_professor_year,female_ratio,attrition_rate_m_pool,attrition_rate_f_pool,promotion_rate_f_pool,promotion_rate_m_pool,all_level,same_promotion);
time_to_equate_pool(4) = sim_result;

same_promotion = 1;

all_level=0;
female_ratio = 1;
sim_result = simulation_on_equate_estimated_param_new(professor_year,asso_professor_year,assi_professor_year,female_ratio,attrition_rate_m_pool,attrition_rate_f_pool,promotion_rate_f_pool,promotion_rate_m_pool,all_level,same_promotion);
time_to_equate_pool(5) = sim_result;



female_ratio = 0.75;
sim_result = simulation_on_equate_estimated_param_new(professor_year,asso_professor_year,assi_professor_year,female_ratio,attrition_rate_m_pool,attrition_rate_f_pool,promotion_rate_f_pool,promotion_rate_m_pool,all_level,same_promotion);
time_to_equate_pool(6) = sim_result;

all_level=1;
female_ratio = 1;
sim_result = simulation_on_equate_estimated_param_new(professor_year,asso_professor_year,assi_professor_year,female_ratio,attrition_rate_m_pool,attrition_rate_f_pool,promotion_rate_f_pool,promotion_rate_m_pool,all_level,same_promotion);
time_to_equate_pool(7) = sim_result;



female_ratio = 0.75;
sim_result = simulation_on_equate_estimated_param_new(professor_year,asso_professor_year,assi_professor_year,female_ratio,attrition_rate_m_pool,attrition_rate_f_pool,promotion_rate_f_pool,promotion_rate_m_pool,all_level,same_promotion);
time_to_equate_pool(8) = sim_result;


all_level=1;
female_ratio = 0.5;

sim_result = simulation_on_equate_estimated_param_new(professor_year,asso_professor_year,assi_professor_year,female_ratio,attrition_rate_m_pool,attrition_rate_f_pool,promotion_rate_f_pool,promotion_rate_m_pool,all_level,same_promotion);
time_to_equate_pool(9) = sim_result;

all_level=0;
female_ratio = 0.5;

sim_result = simulation_on_equate_estimated_param_new(professor_year,asso_professor_year,assi_professor_year,female_ratio,attrition_rate_m_pool,attrition_rate_f_pool,promotion_rate_f_pool,promotion_rate_m_pool,all_level,same_promotion);
time_to_equate_pool(10) = sim_result;

