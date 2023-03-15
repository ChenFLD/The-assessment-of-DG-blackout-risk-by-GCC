% function Final_code_workingV2(temperature)


clc
clear all
close all
tic
define_constants;
%% CPF option setting
mpopt = mpoption('out.all', 0, 'verbose', 2);
mpopt = mpoption(mpopt, 'cpf.stop_at', 'NOSE');
mpopt = mpoption(mpopt, 'cpf.parameterization',3);
mpopt = mpoption(mpopt, 'cpf.adapt_step',1,'cpf.step', 0.1,...
    'cpf.step_min', 0.001,'cpf.step_max', 0.5);
mpopt = mpoption(mpopt, 'cpf.enforce_q_lims', 0);
%% Create a m case file for the vader topology.

mpcb = loadcase('case123.m');
mpcb.baseMVA = 100;

% Load the distribution feeder data from csv file and modify "mpcb".
alloy_data = csvread("converted_alloy_glm_data.csv");
total_buses_in_alloy = max(max(alloy_data(:,1)),max(alloy_data(:,2)));
% modify the mpcb.bus first.
bus_data_format = mpcb.bus;
alloy_bus_data = [];
% update bus number and bus type data.
alloy_bus_data = [alloy_bus_data transpose(1:total_buses_in_alloy) [3;ones(total_buses_in_alloy-1,1)]];


bus_data_format(:,3) = bus_data_format(:,3)/10;
bus_data_format(:,4) = bus_data_format(:,4)/100;
% update the PD, QD.
pd_qd = [];
vm_va = [];
for i = 1:total_buses_in_alloy
    if i <=123
        pd_qd(i,1) = bus_data_format(i,3);
        pd_qd(i,2) = bus_data_format(i,4);
        
        vm_va(i,1) = bus_data_format(i,8);
        vm_va(i,2) = bus_data_format(i,9);
    else
        try
            % "122" because we exclude slack bus.
            pd_qd(i,1) = bus_data_format(i-122,3);
            pd_qd(i,2) = bus_data_format(i-122,4);
            
            vm_va(i,1) = bus_data_format(i-122,8);
            vm_va(i,2) = bus_data_format(i-122,9);
        catch ME % when i >= 246
            pd_qd(i,1) = bus_data_format(i-244,3);
            pd_qd(i,2) = bus_data_format(i-244,4);
            
            vm_va(i,1) = bus_data_format(i-244,8);
            vm_va(i,2) = bus_data_format(i-244,9);
        end
    end
end

% Initial guess of the power flow. Made the intial guess 1+j0.
vm_va(:,1) = ones(total_buses_in_alloy,1);
vm_va(:,2) = zeros(total_buses_in_alloy,1);

alloy_bus_data = [alloy_bus_data pd_qd];
% add all other extra rows unrelated to typical power flow program.
alloy_bus_data = [alloy_bus_data zeros(total_buses_in_alloy, 1)... 
    zeros(total_buses_in_alloy,1) ones(total_buses_in_alloy,1) vm_va... 
    1.1*ones(total_buses_in_alloy,1) ones(total_buses_in_alloy,1)... 
    1.1*ones(total_buses_in_alloy,1) 0.9*ones(total_buses_in_alloy,1)];

% Let us make the gen data now. (bus data above is complete).
gen_bus_format = mpcb.gen;
alloy_gen_data = gen_bus_format;
alloy_gen_data(1,2) = gen_bus_format(1,2)*2;%*100;
alloy_gen_data(1,3) = gen_bus_format(1,3)*2;%*100;
alloy_gen_data(1,4) = gen_bus_format(1,4)*2;%*100; % qmax
alloy_gen_data(1,5) = gen_bus_format(1,5)*2;%*100; % qmin


% Let us make the branch data now. (gen data above is complete).
branch_format = mpcb.branch;
alloy_branch_data = [];
total_branch_elements = size(alloy_data, 1);
alloy_branch_data = [alloy_branch_data alloy_data(:,1)]; % from bus.
alloy_branch_data = [alloy_branch_data alloy_data(:,2)]; % to bus.
alloy_branch_data = [alloy_branch_data alloy_data(:,3)]; 
alloy_branch_data = [alloy_branch_data alloy_data(:,4)];
alloy_branch_data = [alloy_branch_data alloy_data(:,5)];
alloy_branch_data = [alloy_branch_data 250*ones(total_branch_elements,1)... 
    250*ones(total_branch_elements,1) 250*ones(total_branch_elements,1)... 
    zeros(total_branch_elements,1)...
     zeros(total_branch_elements,1) ones(total_branch_elements,1)... 
     -360*ones(total_branch_elements,1) 360*ones(total_branch_elements,1)];
 
% Update the mpcb with the alloy data now!
mpcb.bus = alloy_bus_data;
mpcb.gen = alloy_gen_data;
mpcb.branch = alloy_branch_data;

% Assign some to zero injection buses as well.
%--------------------------------------------------------------------------
a = 1;
b = size(mpcb.bus,1);

%rng(1);
for i = 1:b
    if mpcb.bus(i,PD) == 0
        rand_bus_no = round((b-a).*rand(1,1) + a);
        while mpcb.bus(rand_bus_no,PD) == 0
            rand_bus_no = round((b-a).*rand(1,1) + a);
        end
        mpcb.bus(i,PD) = mpcb.bus(rand_bus_no,PD);
    end
    if mpcb.bus(i,QD) == 0
        rand_bus_no = round((b-a).*rand(1,1) + a);
        while mpcb.bus(rand_bus_no,PD) == 0
            rand_bus_no = round((b-a).*rand(1,1) + a);
        end
        mpcb.bus(i,QD) = mpcb.bus(rand_bus_no,QD);
    end
end

% Make Ybus. Note: The below code ignores shunts since we removed the
% diagonal elements in it. This discrepency is compensated by above code.
mpcb=ext2int(mpcb);
mpcb=rmfield(mpcb,'order');
% mpcb.baseMVA = 1;

[Ybus, ~, ~] = makeYbus(mpcb.baseMVA, mpcb.bus, mpcb.branch);
Ybus = full(Ybus);
output_req_Ybus = Ybus - diag(diag(Ybus));
bus_matrix_G = real(output_req_Ybus);
bus_matrix_B = imag(output_req_Ybus);

num_bus = size(mpcb.bus,1);

% STEP1: Find the lambda of the basecase file loaded using CPF.

% Find the maximum load multiplier a specific ieee system can take to have
% a valid power flow solution.
mpcb_noload=mpcb;
mpcb_baseload=mpcb;
% no generation
mpcb_noload.gen(:, [PG ]) = mpcb_noload.gen(:, [PG ]) * 0;
% and no load
mpcb_noload.bus(:, [PD QD]) = mpcb_noload.bus(:, [PD QD]) * 0;
%%

G = graph(mpcb_baseload.branch(:,1),mpcb_baseload.branch(:,2));
plot(G)
deg_list = degree(G);
%% Code below intially done to test if test system converges.

% target load
mpcb_baseload.bus(:, [PD QD]) = mpcb_baseload.bus(:, [PD QD])*1;

% Run CPF.
results = runcpf(mpcb_noload, mpcb_baseload, mpopt);
feasible_max_load_multiplier_for_loaded_case = results.cpf.lam(1,end);
var1 = feasible_max_load_multiplier_for_loaded_case;
%% Code below to characterize the boundary.

pf_flag = 1;
multiplier = 0;
theta = acos(0.98);
mpcb_original = mpcb;
% volt_save = [];
while pf_flag == 1
    mpcb.bus(:,3) = mpcb_original.bus(:,3)*multiplier;
    mpcb.bus(:,4) = tan(theta)*mpcb.bus(:,3);
    resultss = runpf(mpcb, mpoption('out.all', 0, 'verbose', 0));
    pf_flag = resultss.success;
    multiplier = multiplier +1;
end
boundary_multiplier_consistent_approach = multiplier - 1;

% STEP2: Use optimization problem to characterize the power bundles.

mat_P_store = [];
store = [];

previous_iter_cvx_value = 72638; % initalization dummy number.

all_multiplier_vector = linspace(0,boundary_multiplier_consistent_approach,360);
for kool = 1:length(results.cpf.V_hat(1,:))%1:length(all_multiplier_vector)
    
    % Solve power flow and obtain the voltage data.
    mpcb.bus(:,3) = mpcb_original.bus(:,3)*all_multiplier_vector(kool);
    rect_volt = results.cpf.V_hat(:,kool);
    v_kr = real(rect_volt);
    v_ki = imag(rect_volt);
    
    % calculate the jacobian.
    [Jacob]=CalculateJacobian(num_bus,bus_matrix_G,bus_matrix_B,...
        v_kr,v_ki);
    
    % optimization problem below.
    cvx_begin
    variable y
    maximize(sum(sum(Jacob,1)))
    subject to
    Jacob >= 0;
    norm(y,2) <= 1;
    cvx_end
    
    store(kool,1) = cvx_optval;
    previous_iter_cvx_value = cvx_optval;
    mat_P_store(kool,:) =  transpose(results.bus(:,3)*results.cpf.lam_hat(kool));
end

%% Calculate means for all buses for 10 regions.
% % Normalize the margin matrix.
ieee_min_store_value = min(store(:));
NormStore = store - min(store(:));
ieee_max_partial_norm_store_value = max(NormStore(:));
NormStore = NormStore./max(NormStore(:));

t_iter = 0;
mean_mat = [];
for margin = 0.9:-0.1:0
    t_iter = t_iter +1;
    % finds the data relevant to the required margin specified.
    Indices_req_margin = find((NormStore>=margin)&(NormStore<...
        margin+0.1));
    for busno = 1:num_bus
        % It is first column because I assumed that powers at every bus
        % are same.
        mean_value = mean(mat_P_store(Indices_req_margin,busno));
        mean_mat(busno,t_iter) = mean_value;
    end
end

%% Working with SCE customer data taken after using py file.

load('pnas_work_savepoint2')

% Making certain cutomers zero... because I couldn't find their latitudes
% and longitudes so for consistency sake, I am removing the data
% corresponding those two loads present in the cobalt feeder.
cobalt_customersData(46) = [];
cobalt_customersData(262) = [];

alloy_customersData = alloy_customersData(~cellfun('isempty',alloy_customersData));
uranium_customersData = uranium_customersData(~cellfun('isempty',uranium_customersData));
cobalt_customersData = cobalt_customersData(~cellfun('isempty',cobalt_customersData));
%% Working with temperature data and a lot of data cleaning!

load("xls_option.mat")
% remove unnecessary columns.
data(:,[1 4:43 45:end]) = [];
% remove the Nan and also preserve the data quality and resolution.
data = data_cleaner(data);
% remove the headings.
data(1,:) = [];
% remove the time stamps and only keep day-month-year.
data(:,1) = strtok(data(:,1), 'T');
% convert the temperature data into a matrix so we can use it!
[mat_day_hoursOfaDay,temperature_data_mat,unique_days] = Convert2Matrix(data);

%% Combine the temperature data and load data together.

% remove day 343 from everwhere since the temperature data is missing
% during this day for 5 hours!!!

temperature_data_mat(343,:) = [];

% Remove NaNs in the temperature dataset and replace them with neighboring
% temperature value.
for row = 1:size(temperature_data_mat,1)
    for col = 1:size(temperature_data_mat,2)
        if isnan(temperature_data_mat(row,col))
            temperature_data_mat(row,col) = temperature_data_mat(row-1,col);
        end
    end
end
unique_days(343) = [];
for i = 1:3
    if i == 1
        for j = 1:size(alloy_customersData,1)
            alloy_customersData{j,1}(343,:) = [];
        end
    elseif i == 2
        for j = 1:size(uranium_customersData,1)
            uranium_customersData{j,1}(343,:) = [];
        end
    elseif i == 3
        for j = 1:size(cobalt_customersData,1)
            cobalt_customersData{j,1}(343,:) = [];
        end
    end
end

total_customer_list = cell(300,1);
cust_index = 0;
for i = 1:3 % for alloy, uranium and cobalt.
    if i ==1
        not_segregated_customer_set = alloy_customersData;
    elseif i ==2
        not_segregated_customer_set = uranium_customersData;
    elseif i ==3
        not_segregated_customer_set = cobalt_customersData;
    end
    
    total_customers_no = length(not_segregated_customer_set);
    
    total_valid_customers_5_in_set = fix(total_customers_no/5);
    
    count = 0;
    
    for j = 1:total_valid_customers_5_in_set
        count = count +1;
        customer_data1 = not_segregated_customer_set{count,1};
        count = count +1;
        customer_data2 = not_segregated_customer_set{count,1};
        count = count +1;
        customer_data3 = not_segregated_customer_set{count,1};
        count = count +1;
        customer_data4 = not_segregated_customer_set{count,1};
        count = count +1;
        customer_data5 = not_segregated_customer_set{count,1};
        mean_customer_data = (customer_data1 + customer_data2 + ...
            customer_data3 + customer_data4 + customer_data5)/5;
        cust_index = cust_index + 1;
        total_customer_list{cust_index,1} = mean_customer_data;
    end
end
total_customer_list = total_customer_list(~cellfun('isempty',total_customer_list));
%% OUTLIER
% Concept references.
% ref: https://statisticsbyjim.com/basics/outliers/
% ref: https://machinelearningmastery.com/how-to-use-statistics-to-identify-outliers-in-data/


% Standard deviation method: (Z-scores).
% Let us check percent of the given one year data for each customer that
% has outliers that are greater than +/- 3 standard deviations.
z_scores_of_observations = cell(size(total_customer_list,1),1);
percentage_outliers_customers = cell(size(total_customer_list,1),1);
for i = 1:size(total_customer_list,1) % for each customer
    % calculate mean and standard deviation for each hour of the year.
    mean_values = (mean(total_customer_list{i,1},1));
    std_value = (std(total_customer_list{i,1}));
    z_scores_of_observations{i,1} = [];
    percentage_outliers_customers{i,1} = [];
    for j = 1:24 % each hour (column i.e., random variable)
        count_outliers = 0;
        for k = 1:length(total_customer_list{i,1}(:,j)) % each observation
            obs_value = total_customer_list{i,1}(k,j);
            mu_value = mean_values(1,j);
            sigma_value = std_value(1,j);
            z_value = (obs_value-mu_value)/sigma_value;
            z_scores_of_observations{i,1}(k,j) = z_value;
            % Calculate the total outliers as counter for each hour of each
            % customer.
            if z_value > 3 || z_value < -3
                count_outliers = count_outliers + 1;
            end
        end
        % percentage of outliers in each hour.
        percent_value = count_outliers/length(total_customer_list{i,1}(:,j));
        percentage_outliers_customers{i,1} = [percentage_outliers_customers{i,1}; percent_value*100];
    end
end

% Interquantile method:
iqr_percentage_outliers_customers = cell(size(total_customer_list,1),1);
% get iqr values for each hour of each customer
for i = 1:size(total_customer_list,1) % for each customer
    % iqr values for each hour of a specific customer.
    iqr_values = iqr(total_customer_list{i,1})*3;
    for j = 1:24 % each hour (column i.e., random variable)
        if i == 9
            k = 1;
        end
        count_outliers = 0;
        q25 = prctile(total_customer_list{i,1}(:,j),25);
        q75 = prctile(total_customer_list{i,1}(:,j),75);
        lower_value =q25-iqr_values(1,j);
        upper_value =q75+iqr_values(1,j);
        for k = 1:length(total_customer_list{i,1}(:,j)) % each observation
            obs_value = total_customer_list{i,1}(k,j);
            if obs_value < lower_value || obs_value > upper_value
                count_outliers = count_outliers + 1;
            end
        end
        % percentage of outliers in each hour.
        percent_value = count_outliers/length(total_customer_list{i,1}(:,j));
        iqr_percentage_outliers_customers{i,1} = [iqr_percentage_outliers_customers{i,1}; percent_value*100];
    end
end

Q = cell(size(total_customer_list,1),1);
threshold_z_score = 2;
customer_hour_to_be_removed = [];
for i = 1:size(total_customer_list,1) % for each customer
    for j = 1:24 % each hour (column i.e., random variable)
        count_outliers = 0;
        if i == 9
            kdwsasa = 1;
        end
        q10 = prctile(total_customer_list{i,1}(:,j),10);
        q20 = prctile(total_customer_list{i,1}(:,j),20);
        q30 = prctile(total_customer_list{i,1}(:,j),30);
        q40 = prctile(total_customer_list{i,1}(:,j),40);
        q50 = prctile(total_customer_list{i,1}(:,j),50);
        q60 = prctile(total_customer_list{i,1}(:,j),60);
        q70 = prctile(total_customer_list{i,1}(:,j),70);
        q80 = prctile(total_customer_list{i,1}(:,j),80);
        q90 = prctile(total_customer_list{i,1}(:,j),90);
        q100 = prctile(total_customer_list{i,1}(:,j),100);
        q_vec = [q10 q20 q30 q40 q50 q60 q70];% q80 q90 q100];
        
        percent_q_change_value = [];
        for k = 1:length(q_vec)-1
            percent_q_change_value = [percent_q_change_value (q_vec(1,k+1)-q_vec(1,k))/q_vec(1,k)];
        end
        
        mu_data = mean(percent_q_change_value);
        std_data = std(percent_q_change_value);
        z_data = (percent_q_change_value-mu_data)/std_data;
        truth_array = z_data>= threshold_z_score;
        decision_var = any(truth_array);
        t1= percent_q_change_value>=1;
        d2 = any(t1);
        if d2 == 1
            if decision_var == 1
                customer_hour_to_be_removed = [customer_hour_to_be_removed; i j];
            end
        end
    end
end

list_of_customer_to_be_removed = unique(customer_hour_to_be_removed(:,1));

% remove the above customer data.
for i = 1:length(list_of_customer_to_be_removed)
    total_customer_list{list_of_customer_to_be_removed(i,1), 1} = [];
end
total_customer_list = total_customer_list(~cellfun('isempty',total_customer_list));
%% Assign loads from real-data to buses in the grid
%rng(1);
a = 1;
b = length(total_customer_list);
cust_indx = fix(a + (b-a).*rand(num_bus,1));

% Rearrange the customer data as ieee data with respective hour.
p_bundle = cell(24,1);
for i = 1:24
    p_bundle_mat = [];
    for j = 1:364
        temp = [];
        for k = 1:num_bus
            L = cust_indx(k);
            temp = [temp total_customer_list{L,1}(j,i)]; % 1x123
        end
        p_bundle_mat = [p_bundle_mat;temp]; % 364x123
    end
    p_bundle{i,1} = p_bundle_mat;
end
%% Create a different climate change scenario for risk comparison results.
load('test_data.mat');  
load("losangeles_load_45.mat")
load("losangeles_load_85.mat")


%temperature_years  = cell(12,1);%cell(4+4+4,1); 
temperature_years  = cell(18,1);%3*6
% First 6 - max
% Second 6 - min
% Third 6 - mean
% The latter 18 is RCP 8.5


%% Let us predict the consumer load behavior for future climate change scenarios.

% % PARALLEL PROCESSING CODE
% % 
% % % get p_mat and S_mat for ML fit 
% % p_mat = [];
% % S_mat = [];
% % for i = 1:size(ML_fit_data,1)
% %     for j = 1:busno
% %         
% %     end
% % end

p_bundle_list = cell(size(temperature_years,1),1);
for k = 1:size(temperature_years,1) % total climate change scenarios.
    %temperature_data_mat_22 = temperature_years{i,1};
    temperature_data_mat_22 = zeros(3,1);
    p_bundle_list{k,1} = cell(size(temperature_data_mat_22,2),1);
    for i=1:3%max,min,mean
       for j=1:3%year
          for each_sce=1:6
              p_bundle_list{6*(i-1)+each_sce,1}{j,1} = losangeles_load_45{each_sce,2}{j,i};
%               p_bundle_list{18+6*(i-1)+each_sce,1}{j,1} = losangeles_load_85{each_sce,2}{j,i};
          end
       end
    end
   
%        i
%        each_year
%         time_taken_for_one_year = toc
%    end
end
%% Normalize the real-time power data to ieee range.

convert_to_old_mean_bundle_list = cell(size(temperature_years,1),1);
%rng(1);
convert_to_old_mean_bundle = cell(size(temperature_years,1));
for i = 1:size(temperature_years,1)
    
    convert_to_old_mean_bundle{i,1} = cell(size(p_bundle_list{i,1},1),1); % 150
    
    % calculate the mean and std of the basecase system.
    hour_i = 1;
    mu_target_all_buses = mean(p_bundle{hour_i,1}); %1x123 vector.
    std_target_all_buses = std(p_bundle{hour_i,1});
    
    mu_ref_all_buses = transpose(mean(mean_mat,2)); %1x123 vector.
    std_ref_all_buses = std(transpose(mean_mat));
    
    mat1 = p_bundle_list{i,1};
    for each_year = 1:size(p_bundle_list{i,1},1)
        mat2 = mat1{each_year,1};
        logic_vec = mean(p_bundle{hour_i,1}) <= mean(mat2);
        
        z_score_of_target = zeros(size(mat2));
        % calculate the z score based on base case.
        for each_bus = 1:size(mat2,2)
            actual_value = mat2(:,each_bus);
            zscore_value = (actual_value - mu_target_all_buses(1,each_bus))/std_target_all_buses(1,each_bus);
            if logic_vec(each_bus) == 0 % case mentioned in above description.
                zscore_value = abs(zscore_value);
            end
            
            z_score_of_target(:,each_bus) = zscore_value;
        end
        % normalize the data
        convert_to_old_mean_mat{i,1}{each_year,1} = [];
        convert_to_old_mean_mat{i,1}{each_year,1} = (mu_ref_all_buses) + z_score_of_target.*((std_ref_all_buses));
    end
%     convert_to_old_mean_bundle{i,1}{}
end

%% Reduce the load consumption of consumers for the years between X to Y. Alternatively scale down all consumption values of ALL years.

percent = 0.1;
reduce_percent = 1-percent;
for i = 1:size(convert_to_old_mean_mat,1)
    for j = 1:size(convert_to_old_mean_mat{i,1},1)
        convert_to_old_mean_mat{i,1}{j,1} = convert_to_old_mean_mat{i,1}{j,1}*reduce_percent;
    end
end

%%

normalized_store_conditional_cvx_values_list = cell(size(temperature_years,1),1);
store_conditional_cvx_values = cell(size(temperature_years,1),1);

%rng(1)
define_constants;
power_factor = 0.98;
theta = acos(power_factor);
mpopt = mpoption('out.all', 0, 'verbose', 0);
mpopt = mpoption(mpopt, 'cpf.stop_at', 1);
num_bus = size(mpcb_original.bus(:,1),1);

for i = 1:size(temperature_years,1)
    temperature_data_mat = temperature_years{i,1};
    convert_to_old_current_season = convert_to_old_mean_mat{i,1};
    
    store_conditional_cvx_values{i,1} = cell(size(convert_to_old_current_season,1),1);
    normalized_store_conditional_cvx_values_list{i,1} = cell(size(convert_to_old_current_season,1),1);
    
    for each_year = 1:size(convert_to_old_current_season,1) % 150 years
        current_year_data = convert_to_old_current_season{each_year,1};
        
        
        store_conditional_cvx_values{i,1}{each_year,1} = ones(size(current_year_data,1),1)*99999;
        normalized_store_conditional_cvx_values_list{i,1}{each_year,1} = ones(size(current_year_data,1),1)*99999;

        store_vec = ones(size(current_year_data,1),1)*99999;
        normalized_store_vec = ones(size(current_year_data,1),1)*99999;
        
        for each_sample = 1:size(current_year_data,1)

            p_values = current_year_data(each_sample,:);
            
            mpcb_work = mpcb_original;
            mpcb_work.bus(:,PD) = transpose(p_values)/mpcb_work.baseMVA;
%             mpcb_work.bus(:,QD) = transpose(q_values);
            
            % solve power flow for voltage measurements.
            resulttt = runcpf(mpcb_noload, mpcb_work, mpopt);
            fprintf('Solving power flow solution corresponding to season %d at year %d of day %d \n\n\n', i, 2020+each_year*10,each_sample)
            
            if resulttt.success == 1
            v_kr=real(resulttt.cpf.V(:,end));
            v_ki=imag(resulttt.cpf.V(:,end));
            
            % calculate the jacobian.
            [Jacob]=CalculateJacobian(num_bus,bus_matrix_G,bus_matrix_B,v_kr,v_ki);
            
            cvx_optval = implement_cvx(Jacob);
            
            store_vec(each_sample,1) = cvx_optval;
            normalized_store_vec(each_sample,1) = (cvx_optval - ieee_min_store_value)/...
                ieee_max_partial_norm_store_value;
            if (cvx_optval - ieee_min_store_value)/ieee_max_partial_norm_store_value > 1
                normalized_store_vec(each_sample,1) = 1;
            elseif (cvx_optval - ieee_min_store_value)/ieee_max_partial_norm_store_value < 0
                normalized_store_vec(each_sample,1) = 0;
            end
            else
                store_vec(each_sample,1) = 9999; % This indicates that the location where there are "9999" in "store_vec(each_sample,1)" are the location when there is no power flow solution which means no voltages to calculate the cvx value. Hence I will just make the risk = 100%.
                normalized_store_vec(each_sample,1) = 0;
                fprintf("Uh-oh! no power flow solution, predicted data is way beyond boundary of the system.")
            end
        end
        store_conditional_cvx_values{i,1}{each_year,1} = store_vec;
        normalized_store_conditional_cvx_values_list{i,1}{each_year,1} = normalized_store_vec;
    end
end

%% Start counting and extracting the probabilities for all margin layer given hour.
%clc
%clear all

data_mat = [];

for each_temperature_scenario = 1:length(normalized_store_conditional_cvx_values_list)
    data_to_plot = [];
    for each_year_in_scenario = 1:length(normalized_store_conditional_cvx_values_list{each_temperature_scenario,1})
        data_to_plot = [data_to_plot, mean(normalized_store_conditional_cvx_values_list{each_temperature_scenario,1}{each_year_in_scenario,1})];
    end
    data_mat(:,each_temperature_scenario) = data_to_plot;
end
data_season = [];
d=[0,90,181,273,365];
for each_temperature_scenario = 1:length(normalized_store_conditional_cvx_values_list)
    data_plot = [];
    for each_year_in_scenario = 1:length(normalized_store_conditional_cvx_values_list{each_temperature_scenario,1})
        for i=1:4
        data_plot = [data_plot, mean(normalized_store_conditional_cvx_values_list{each_temperature_scenario,1}{each_year_in_scenario,1}(d(i)+1:d(i+1),1))];
        end
    end
    data_season(:,each_temperature_scenario) = data_plot;
end
data_mat_45=1-data_mat(:,1:18);

data_season_45=1- data_season(:,1:18);

figure
% means
% C = linspecer(length(normalized_store_conditional_cvx_values_list));
% colororder(C);
plot(2030:10:2050,data_mat_45(:,1:6),'LineWidth',4)
xlabel("Year")
ylabel("Risk to system collapse in RCP4.5(max)")
set(gcf,'Position',  [100, 100, 755, 550])
set(gca,'fontsize',20)

figure
% maximum upper bound
% C = linspecer(length(normalized_store_conditional_cvx_values_list));
% colororder(C);
plot(2030:10:2050,data_mat_45(:,7:12),'LineWidth',4)
xlabel("Year")
ylabel("Risk to system collapse in RCP4.5(min)")
set(gcf,'Position',  [100, 100, 755, 550])
set(gca,'fontsize',20)

figure
% minimum lower bound
% C = linspecer(length(normalized_store_conditional_cvx_values_list));
% colororder(C);
plot(2030:10:2050,data_mat_45(:,13:18),'LineWidth',4)
xlabel("Year")
ylabel("Risk to system collapse in RCP4.5(mean)")
set(gcf,'Position',  [100, 100, 755, 550])
set(gca,'fontsize',20)

%% below function
function [cvx_optval] = implement_cvx(Jacob)
% optimization problem below.
    cvx_begin
    variable y
    maximize(sum(sum(Jacob,1)))
    subject to
    Jacob >= 0;
    norm(y,2) <= 1;
    cvx_end
end

function convert_to_old_mean_bundle = norm_data_to_ieee_range(convert_to_old_mean_bundle,p_bundle,mean_mat,p_bundle_22_n,flag)
for j = 1:24
    convert_to_old_mean_mat = [];
    mu_target_all_buses = mean(p_bundle{j,1}); %1x123 vector.
    std_target_all_buses = std(p_bundle{j,1});
    if flag == 1
        logic_vec = mean(p_bundle{j,1}) <= mean(p_bundle_22_n{j,1});
    elseif flag == 0
        logic_vec = ones(1,size(p_bundle{1,1},2));
    end
    
    mu_ref_all_buses = transpose(mean(mean_mat,2)); %1x123 vector.
    std_ref_all_buses = std(transpose(mean_mat));
    for each_dayy = 1:size(p_bundle{j,1},1)
        for bus_nummmber = 1:size(mean_mat,1) % bus number
            z_score_of_target = (p_bundle_22_n{j,1}(each_dayy,bus_nummmber) - mu_target_all_buses(1,bus_nummmber))/std_target_all_buses(1,bus_nummmber);
            if logic_vec(bus_nummmber) == 0 % case mentioned in above description.
                z_score_of_target = abs(z_score_of_target);
            end
            convert_to_old_mean_mat(each_dayy,bus_nummmber) = (mu_ref_all_buses(1,bus_nummmber) + std_ref_all_buses(1,bus_nummmber)*z_score_of_target);
            %             z_score_of_target
        end
        convert_to_old_mean_mat(each_dayy,bus_nummmber+1) = p_bundle_22_n{j,1}(each_dayy,bus_nummmber+1);
    end
    convert_to_old_mean_bundle{j,1} = convert_to_old_mean_mat;
end
end

function p_bundle_22 = predict_behavior(p_bundle, temperature_data_mat_2, ML_fit_data,total_system_load_bad_point)
p_bundle_22 = cell(24,1);
for hour_i = 1:24 % 24 hours.
    p_bundle_22{hour_i,1} = [];
    for bus_j = 1:size(p_bundle{hour_i,1},2) % 123 buses
        % get mean and std of the original load data to find the z score
        % corresponding to the day of the year.
        muu = mean(p_bundle{hour_i,1}(:,bus_j));
        stdd = std(p_bundle{hour_i,1}(:,bus_j));
        for each_day_k = 1:364
            x = temperature_data_mat_2(each_day_k,hour_i);
            z_score_value = (p_bundle{hour_i,1}(each_day_k,bus_j) - muu)/stdd;
            if z_score_value > 2 % THIS PART OF CODE CAN CHANGE THE FINAL RISK..LOOK INTO THIS BEFORE GETTING FINAL RESULTS FOR THE PAPER.
                z_score_value = 1;
            elseif z_score_value < -1
                z_score_value = 0;
            end
            p = ML_fit_data{hour_i,1}{bus_j,1}{1,1};
            S = ML_fit_data{hour_i,1}{bus_j,1}{2,1};
            [y_fit,delta] = polyval(p,x,S);
            y_hat = y_fit + (z_score_value)*delta;
            p_bundle_22{hour_i,1}(each_day_k,bus_j) = y_hat;
        end     
    end
    % THE BELOW DOES NOT MATTER BECAUSE WE ARE NORMALIZING THE P_BUNDLE TO
    % IEEE BOUNDARY RANGE ANYWAY IN THE LATER SECTIONS OF THE CODE.
%     % The below code to make sure that I do not enter
%     % to the non-linearity region of the linear programing optimization
%     % model.
%     
%     vals = sum(p_bundle_22{hour_i,1},2);
%     for i = 1:length(vals)
%         if vals(i) >= total_system_load_bad_point
%             kkk = 1;
%             print("omg check here!!!")
% %             break;
%         end
%     end
    
     
end
end

function new_p_bundle = generate_missing_data(temperature_data_mat,p_bundle,no_samples)
new_p_bundle = cell(24,1);
for i = 1:24
    new_p_bundle{i,1} = [];
    new_p_bundle_mat_all_temps_given_hour = [];
    % find the unique temperatures in hour "i".
    temp_valuesss = unique(temperature_data_mat(:,i));
    % collect the day indices for each unique temperature and fit the data
    % and sample from the distribution.
    for j = 1:length(temp_valuesss)
        generated_data = [];
        new_p_bundle_mat = [];
        % select a unique temperature and its data corresponding days.
        daysss_indexx = temperature_data_mat(:,i) == temp_valuesss(j);
        data_to_fit = p_bundle{i,1}(daysss_indexx,:);
        [std_array,data_to_fit] = find_std_loop(p_bundle{i,1},data_to_fit,temperature_data_mat(:,i),temp_valuesss(j));
        % below randomly generates the z-score. but it only generate the
        % +ve z scores and hence the minimum value in the generated data
        % will be mean. which is the predicted value of ML_fit_data present
        % in the "p_bundle" variable. This mean is atleast greater than the
        % previous temperature scenario (I know its bit hard to catch up
        % but it is what it is.).
        a = 1;
        b = 1.05; % This means that the data mostly lies around 1 std from the mean. we do not want it to fluctuatet too much as the results will be not explainable then in risk versus time/temperature!!!.
        r = a + (b-a).*rand(no_samples,1); % THIS CAN ALSO CHANGE RISK. SO BEFORE FINAL RESULTS MAKE SURE TO CHANGE THIS.
        
        mean_array = mean(data_to_fit);
        for each_sample = 1:length(r)
            data_sample = mean_array + r(each_sample)*std_array;
            generated_data = [generated_data;data_sample];
        end
        new_p_bundle_mat = [new_p_bundle_mat;[generated_data ones(size(generated_data,1),1)*j]];
        new_p_bundle_mat_all_temps_given_hour = [new_p_bundle_mat_all_temps_given_hour;new_p_bundle_mat];
    end
    
    new_p_bundle{i,1} = new_p_bundle_mat_all_temps_given_hour;
end
end

function [std_out,data_to_fit] = find_std_loop(p_bundle_source_mat,data_to_fit,temp_array,temp_value)

% validate if a different std from next day is needed due to lack of sample
% size.

% calculate the std_out.
found_flag = 0;
if size(data_to_fit,1) >= 5
    % sample size is enough
    std_out = std(data_to_fit);
else
    % start searching for a temperature values in "temp_value" neighborhood
    % until a new temperature value with good ample size is found to
    % calculate the standard deviation.
    delta = 1; % +/- 2 degrees is acceptable neighborhood to check for sample size.
    while found_flag == 0
    delta = delta +1;
    for i = 1:length(temp_array)
        if temp_array(i) <= temp_value+delta & temp_array(i) >= temp_value-delta
            % check if samples are enough
            new_days = temp_array == temp_array(i);
            data_to_fit = p_bundle_source_mat(new_days,:);
            if size(data_to_fit,1) >= 5
                std_out = std(data_to_fit);
                found_flag = 1;
                break
            end
        end
    end
    end
end
end

function [output_data] = calculate_avg_risk(risk_data)
output_data = cell(24,1);
for i = 1:24
    output_data{i,1} = [];
    for j = 1:10
        avg_risk_value = mean(risk_data{j,1}{i,1}(:,2));
        output_data{i,1}(j,1) = avg_risk_value;
    end
end
end

function [] = plot_data(temperature_probability_toggle,selected_hour,final_distribution_data)
figure
for margins = 1:10
    subplot(2,5,margins)
    data_lol = final_distribution_data{margins,1}{selected_hour,1}(:,2);
    h1 = histogram(data_lol,length(data_lol));
    h1.BinWidth = 0.1;
    if temperature_probability_toggle== 2
        str_1 = sprintf("Margin layer = %d",margins);
    else
        str_1 = sprintf("Temperature");
    end
    ylim([0 400]);
    xlim([0 1.1]);
    xlabel(str_1);
    ylabel('Frequency');
end
if temperature_probability_toggle == 2
    str = sprintf("Probability of crossing various margin layers at hour = %d",selected_hour);
    sgtitle(str)
else
    str = sprintf("Temperature of various margin layers at hour = %d",selected_hour);
    sgtitle(str)
end
end
function [final_distribution_data] = count_extract_probabilities_to_standard_representation(temperature_data_mat,normalized_store_conditional_cvx_values)
final_distribution_data = cell(10,1);
for i = 1:10
    final_distribution_data{i,1} = cell(24,1);
    for j = 1:24
        final_distribution_data{i,1}{j,1} = [];
    end
end
for i = 1:24
    unique_temp_list = unique(temperature_data_mat(:,i));
    for j = 1:364
        current_temperatureee = temperature_data_mat(j,i);
        % find the index of current_temperature from unique list.
        index_req = find(unique_temp_list == current_temperatureee);
        margin_set = normalized_store_conditional_cvx_values{i,1}{index_req,1};
        % loop for each margin and identify the probability of that margin
        % layer.
        margin_counter_index = 10;
        temp_count_store = zeros(10,1);
        for k = 0.9:-0.1:0 % 0:0.1:0.9
            % finds the data relevant to the required margin specified.
            if k == 0.9
                Indices_req_margin = find((margin_set>=k)&(margin_set<=k+0.1));
            else
                Indices_req_margin = find((margin_set>=k)&(margin_set<k+0.1));
            end
            % If the margin layer value is > 0.1 then it crossed layer 2.
            % similarly if crossed 0 then crossed layer 1, if crossed 0.9
            % then crossed layer 10 i.e., black out.
            current_crossing_layer = margin_counter_index;
            all_layers_beneath_current_layer = 1:current_crossing_layer;
            count_of_points_crossed_current_layer = length(Indices_req_margin);
            %update the count vector for 10 margins
            temp_count_store(all_layers_beneath_current_layer,1) = temp_count_store(all_layers_beneath_current_layer,1) + count_of_points_crossed_current_layer;
            margin_counter_index = margin_counter_index - 1;
        end
        for each_mar = 1:10
            final_distribution_data{each_mar,1}{i,1}(j,2) = temp_count_store(each_mar,1)/length(margin_set);
            final_distribution_data{each_mar,1}{i,1}(j,1) = current_temperatureee;
        end
        
    end
end

end

function [store_conditional_cvx_values,normalized_store_conditional_cvx_values,mat_conditional_P_store] = calculate_and_count_risk(temperature_data_mat,convert_to_old_mean_bundle,mpcb_original,bus_matrix_G,bus_matrix_B,ieee_max_partial_norm_store_value,ieee_min_store_value)
rng(1)
define_constants;
power_factor = 0.98;
theta = acos(power_factor);
store_conditional_cvx_values = cell(24,1);
normalized_store_conditional_cvx_values = cell(24,1);
mat_conditional_P_store = cell(24,1);
mpopt = mpoption('out.all', 0, 'verbose', 0);
num_bus = size(mpcb_original.bus(:,1),1);
% for every hour
for i = 1:24
    temp_valuesss = unique(temperature_data_mat(:,i));
    store_conditional_cvx_values{i,1} = cell(length(temp_valuesss),1);
    normalized_store_conditional_cvx_values{i,1} = cell(length(temp_valuesss),1);
    mat_conditional_P_store{i,1} = cell(length(temp_valuesss),1);
    % for every day of the year at specific hour "i".
    for j = 1:length(temp_valuesss) % length of unique temperature given current hour "i".
        store_conditional_cvx_values{i,1}{j,1} = [];
        normalized_store_conditional_cvx_values{i,1}{j,1} = [];
        mat_conditional_P_store{i,1}{j,1} = [];
        daysss_indexx = convert_to_old_mean_bundle{i,1}(:,end) == j;
        % given temperature value, find the samples and calculate cvx
        % values and save the information.
        p_set = convert_to_old_mean_bundle{i,1}(daysss_indexx,1:end-1);
        for iter_loop = 1:size(p_set,1)
            p_values = p_set(iter_loop,:);
%             q_values = tan(theta)*p_values;
            
            mpcb_work = mpcb_original;
            mpcb_work.bus(:,PD) = transpose(p_values);
%             mpcb_work.bus(:,QD) = transpose(q_values);
            
            resulttt = runpf(mpcb_work,mpopt);
            rect_volt = resulttt.bus(:,8).*(cos(resulttt.bus(:,9)+1j*...
                sin(resulttt.bus(:,9))));
            
            v_kr = real(rect_volt);
            v_ki = imag(rect_volt);
            % calculate the jacobian.
            [Jacob]=CalculateJacobian(num_bus,bus_matrix_G,bus_matrix_B,v_kr,v_ki);
            
            % optimization problem below.
            cvx_begin
            variable y
            maximize(sum(sum(Jacob,1)))
            subject to
            Jacob >= 0;
            norm(y,2) <= 1;
            cvx_end
            
            store_conditional_cvx_values{i,1}{j,1}(iter_loop,1) = ...
                cvx_optval;
            mat_conditional_P_store{i,1}{j,1}(iter_loop,:) = ...
                p_values;
            normalized_store_conditional_cvx_values{i,1}{j,1}(...
                iter_loop,1) = (cvx_optval - ieee_min_store_value)/...
                ieee_max_partial_norm_store_value;
            if (cvx_optval - ieee_min_store_value)/ieee_max_partial_norm_store_value > 1
                normalized_store_conditional_cvx_values{i,1}{j,1}(...
                    iter_loop,1) = 1;
            elseif (cvx_optval - ieee_min_store_value)/ieee_max_partial_norm_store_value < 0
                normalized_store_conditional_cvx_values{i,1}{j,1}(...
                    iter_loop,1) = 0;
            end
        end
    end
end
end

function [data_to_fit,temp_current_index] = calculate_corr_matrix(temp_valuesss,temperature_data_mat,temp_current_index,current_hour_index,p_bundle)
error_flag = 1;
while error_flag == 1
    % check the position of index of temperature. always move towards mean.
    if temp_current_index/length(temp_valuesss) > 0.5
        daysss_indexx = temperature_data_mat(:,current_hour_index) == temp_valuesss(temp_current_index-1);
        data_to_fit = p_bundle{current_hour_index,1}(daysss_indexx,:);
        temp_current_index = temp_current_index -1;
    else
        daysss_indexx = temperature_data_mat(:,current_hour_index) == temp_valuesss(temp_current_index+1);
        data_to_fit = p_bundle{current_hour_index,1}(daysss_indexx,:);
        temp_current_index = temp_current_index +1;
    end
    % check for errors and if there is error repeat. otherwise exit witht
    % he data.
    R = corrcoef(data_to_fit);
    if R == 1
        error_flag = 1;
    elseif sum(isnan(R(:))) ~= 0
        error_flag = 1;
    else
        error_flag = 0;
    end
end
end

% Read the SCE data.
function [alloy_loc_mat,uranium_loc_mat,cobalt_loc_mat,alloy_customersData,uranium_customersData,cobalt_customersData,alo_info_mat,ura_info_mat,cob_info_mat] = read_SCE_data()

% We have three feeders and each feeder has several customers.
% First we complete the processing of alloy then uranium and finally cobalt

% Alloy
[alloy_loc_mat,alloy_customersData,alo_info_mat] = process_feeder_all_customers(1);
%Uranium
[uranium_loc_mat,uranium_customersData,ura_info_mat] = process_feeder_all_customers(2);
%Cobalt
[cobalt_loc_mat,cobalt_customersData,cob_info_mat] = process_feeder_all_customers(3);
end


function [loc_mat,feeder_all_customer_cells,info_mat] = process_feeder_all_customers(feeder_flag)
%feeder_flag = 1 implies alloy, 2 implies uranium and 3 implies cobalt.
info_mat = [];
% Select specific feeder.
if feeder_flag ==1
    % Alloy has 43 customers.
    x = 43;
    f1 = 'C:\Users\kisha\Desktop\Temp_load_paperV4\Temp_loa_paperV4\Distribution test case method\Distribution test case data generation for margin curves\Working\pnas data\alloy\';
    feeder_all_customer_cells = cell(300,1);
    k = 'C:\Users\kisha\Desktop\Temp_load_paperV4\Temp_loa_paperV4\Distribution test case method\Distribution test case data generation for margin curves\Working\pnas data\';
    [~,~,loc_data]  = xlsread([k 'alloy_loc_info.xlsx']);
    loc_mat = [];
elseif feeder_flag ==2
    % Uranium has 146 customers.
    x = 146;
    f1 = 'C:\Users\kisha\Desktop\Temp_load_paperV4\Temp_loa_paperV4\Distribution test case method\Distribution test case data generation for margin curves\Working\pnas data\uranium\';
    feeder_all_customer_cells = cell(300,1);
    k = 'C:\Users\kisha\Desktop\Temp_load_paperV4\Temp_loa_paperV4\Distribution test case method\Distribution test case data generation for margin curves\Working\pnas data\';
    [~,~,loc_data]  = xlsread([k 'uranium_loc_info.xlsx']);
    loc_mat = [];
elseif feeder_flag ==3
    % Cobalt has 292 customers.
    x = 292;
    f1 = 'C:\Users\kisha\Desktop\Temp_load_paperV4\Temp_loa_paperV4\Distribution test case method\Distribution test case data generation for margin curves\Working\pnas data\cobalt\';
    feeder_all_customer_cells = cell(300,1);
    k = 'C:\Users\kisha\Desktop\Temp_load_paperV4\Temp_loa_paperV4\Distribution test case method\Distribution test case data generation for margin curves\Working\pnas data\';
    [~,~,loc_data]  = xlsread([k 'cobalt_loc_info.xlsx']);
    loc_mat = [];
end

loc_data = loc_data(2:end,[5 21 22]);
loc_data(end,:) = [];
koo = 0;
for ix = 1:x
    customer_Power = ones(366,24)*999999;
    % Load a specific customer's year long load information.
    load([f1 'customer_' num2str(ix) '.mat']);
    koo = koo + 1;
    % To skip bad data customer consumption files.
    if size(customerData,1) == 8784
        
        % row number = 4 because all the data in "customerData" belongs to
        % a unique customer.. That is how I processed the data in python
        % code which is given as input to this function.
        customer_ID = string(cell2mat(customerData(4,2)));
        lat_long = [];
        for i = 1:size(loc_data,1)
            if loc_data{i,1}(1,1) == "P"
                loc_data{i,1} = loc_data{i,1}(1,2:end);
            end
            if loc_data{i,1}(1,1) == "B"
                loc_data{i,1} = loc_data{i,1}(1,2:end);
            end
            if loc_data{i,1}(1,1) == "E"
                loc_data{i,1} = loc_data{i,1}(1,2:end);
            end
            if loc_data{i,1}(1,1) == "V"
                loc_data{i,1} = loc_data{i,1}(1,2:end);
            end
            if loc_data{i,1}(1,1) == "S"
                loc_data{i,1} = loc_data{i,1}(1,2:end);
            end
            if customer_ID == string(cell2mat(loc_data(i,1)))
                lat_long = cell2mat(loc_data(i,[2 3]));
                break;
            elseif (customer_ID == "5407570") && (feeder_flag == 3)
                lat_long = cell2mat(loc_data(10,[2 3]));
            elseif (customer_ID == "1087853E") && (feeder_flag == 3)
                lat_long = cell2mat(loc_data(4,[2 3]));
            elseif (customer_ID == "1335106E") && (feeder_flag == 3)
                lat_long = cell2mat(loc_data(268,[2 3]));
            end
        end
        loc_mat = [loc_mat;lat_long];
        
        if isempty(lat_long)
            fprintf("OMG, I did not find anything!!!")
            customer_ID
            customer_no = ix
            feeder_no = feeder_flag
        end
        % Let us sort the data first based on the time stamps.
        [sorted_dates, sorted_dates_index] = sort(datenum(customerData(:,3)));
        sorted_customerData = customerData(sorted_dates_index,:);
        
        % unique days in an year.
        customerData_day_count = unique(sorted_dates);
        
        for i = 1:length(customerData_day_count) % every day in an year.
            specificDay_indices_inData = (sorted_dates ==...
                customerData_day_count(i));
            Specific_day_24hour_data = sorted_customerData(specificDay_indices_inData,:);
            for j = 1:24 % one to 24 hours
                % find hour indices.
                hour_index = (cell2mat(Specific_day_24hour_data(:,4))==j);
                Power_hourWise = cell2mat(Specific_day_24hour_data...
                    (hour_index,6));
                
                % juts for verification.
                %          ALLOY_hours(i,j) = j;
                
                % actual final data with time series of 24 hours and several days.
                %                 disp(i)
                if isempty(Power_hourWise)
                    %do nothing and just exit this customer.
                    
                    % THIS WILL MAKE 999999 STAY IN THE ROW IF THE LOOP IS
                    % SATISFIED AND THEN LATER I CAN REMOVE IT USING "FIND"
                    % [which feeder customer number day hour]
                    info_mat = [info_mat; feeder_flag ix i j];
                    
                    % I am going to create a fake hour data here!!!
                    hour_index = (cell2mat(Specific_day_24hour_data(:,4))==j-1);
                    Power_hourWise = cell2mat(Specific_day_24hour_data...
                        (hour_index,6));
                    customer_Power(i,j) = Power_hourWise;
                else
                    customer_Power(i,j) = Power_hourWise;
                end
                %                 disp(j)
                %                 disp(Power_hourWise)
                %                 disp(customer_Power(i,j))
                %                 customer_Power(i,j) = Power_hourWise;
            end
        end
        customer_Power(end,:) = [];
        feeder_all_customer_cells{ix,1} = customer_Power;
    else
        % do nothing.. Skip the customer...Bad data...does not have full
        % information.
        a = 1;
    end
end
end

% Function to remove the unnecessary data points from the temperature data.

function [data] = data_cleaner(data)
count = 1;
for i = 2:size(data,1)
    count = count + 1;
    var = string(cell2mat(data(count,2)));
    
    if var ~= "FM-15"
        data(count,:) = [];
        count = count - 1;
    end
end
end

% Function to convert the temperature data into a mtrix form so it is
% readable.
function [mat_day_hoursOfaDay,data_mat,unique_days] = Convert2Matrix(data)

% all days in the given data set. each day should ahve 24 hours data.
all_days = datenum(data(:,1));
unique_days = unique(all_days);
no_unique_days = length(unique_days);

mat_day_hoursOfaDay = [];
data_mat = [];
% We need to find the missing specific hours in a day from the data.
for unique_day_no = 1:no_unique_days
    [day_row_index, ~, ~] = find(all_days == unique_days(unique_day_no));
    total_hours_in_a_day = length(day_row_index);
    mat_day_hoursOfaDay = [mat_day_hoursOfaDay; unique_days(unique_day_no) total_hours_in_a_day];
    % BY OBSERVING "mat_day_hoursOfaDay", WE ARE MISSING 5 HOURS OF TEMPERATURE
    % DATA FROM 9TH DECEMBER 2015 and also couple of other days.
    
    % bad data fixing..missing data point.
    if unique_day_no == 70
        data(1666,3) = {[68]};
    end
    if unique_day_no == 117
        data(2794,3) = {[87]};
    end
    if unique_day_no == 252
        data(6039,3) = {[83]};
    end
    if unique_day_no == 282
        data(6756,3) = {[103]};
    end
    if unique_day_no == 336
        data(8052,3) = {[76]};
    end
    
    % Code to make the data_mat.
    temprature_data_of_a_day = cell2mat(data(day_row_index,3));
    % bad data fixing..missing data point.
    if unique_day_no == 92
        temprature_data_of_a_day(24,1) = 60;
    end
    % NOTE: KISHAN REMOVE DAY 343 FROM EVERYWHERE!!!!!!!!!!!
    if unique_day_no == 343
        temprature_data_of_a_day(21:24,1) = 999999;
    end
    data_mat = [data_mat; transpose(temprature_data_of_a_day)];
end


end

function [new_mean_view_instance] = handle_customer_average_model(instance)
new_mean_view_instance = [];
% find the mean of Watts for a unique temperature and create a unique data
% point for a unique temperature. The unique data point is the mean value
% of all the data points(watts) at that temperature value.
unique_temperature_values = unique(instance(:,1));
for i = 1:length(unique_temperature_values)
    index = (instance(:,1) == unique_temperature_values(i));
    g = instance(:,2);
    mean_watts_for_a_given_temp = mean(g(index,1));
    new_mean_view_instance = [new_mean_view_instance; unique_temperature_values(i) mean_watts_for_a_given_temp];
end
end

function all_coefficients = create_Ucurves(no_of_Ushapes)

    t_max = 100;
    t_min = 50;

    p_max = 1e4;
    p_min = 1e3;

    T = linspace(t_min,t_max,1000);

    all_square_terms = zeros(no_of_Ushapes,1);
    all_linear_terms = zeros(no_of_Ushapes,1);
    all_constant_terms = zeros(no_of_Ushapes,1);
    
    for i=1:no_of_Ushapes
        t_sweet = 70 + (t_max - t_min)/5*(rand-0.5); %25C is the mean, max 30, min 20
        p_sweet = p_min + (p_max - p_min)/2 * rand; %3500 Watts is the mean power
        P = (T-t_sweet).^2 + p_sweet; % ,max 7000, min 1000
        square_term = 1;
        linear_term = -2*t_sweet;
        constant_term = t_sweet^2 + p_sweet;

        %Window y Scaling
        coefficient_scaling = 1;
        c = 1.01;
        while (all(P<p_max))
            P = c*P;
            coefficient_scaling = c*coefficient_scaling;
        end

        c=0.99;
        while (any(P>p_max))
            P = c*P;
            coefficient_scaling = c*coefficient_scaling;
        end

        %plot(T,P) %plot a single U shape

        all_square_terms(i) = coefficient_scaling * square_term;
        all_linear_terms(i) = coefficient_scaling * linear_term;
        all_constant_terms(i) = coefficient_scaling * constant_term;
    end

    all_coefficients = [all_square_terms all_linear_terms all_constant_terms];
end

function predicted_power = getPower(c_index, temp, all_coefficients)
    square_term = all_coefficients(c_index,1); 
    linear_term = all_coefficients(c_index,2); 
    constant_term = all_coefficients(c_index,3);
    predicted_power = square_term * temp.^2 + linear_term * temp + constant_term;
end

