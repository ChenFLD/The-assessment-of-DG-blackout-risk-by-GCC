% 清空环境变量
clear
clc
tic

%第一步 读取数据
load('test_data.mat');  %载入训练数据
load("los_angeles_45.mat") %载入预测数据
load("los_angeles.mat")
k=1;c=[7,8,9];
%temp_pretict_45=cell(6,1)(9,365);

for i=1:3 %max,min,mean
    for j=1:3 %year
        for n=1:6
        temp_pretict_45_m{n,1}{j,1}(i,:)=los_angeles_45{1+j,n}(:,c(i))';%
        end
        temp_pretict_85_m{1,1}{j,1}(i,:)=los_angeles_NorESM1{1,j}(:,c(i))';
        temp_pretict_85_m{2,1}{j,1}(i,:)=los_angeles_MPI_ESM_MR{1,j}(:,c(i))';
        temp_pretict_85_m{3,1}{j,1}(i,:)=los_angeles_INMCM4{1,j}(:,c(i))';
        temp_pretict_85_m{4,1}{j,1}(i,:)=los_angeles_GFDL_CM3{1,j}(:,c(i))';
        temp_pretict_85_m{5,1}{j,1}(i,:)=los_angeles_CSIRO{1,j}(:,c(i))';
        temp_pretict_85_m{6,1}{j,1}(i,:)=los_angeles_bcc_csm{1,j}(:,c(i))';
    end
end

n = 24; 
[x_sort s_index1]=sort(temperature_data_mat(10,:));
[x_s2 s_index]=sort(s_index1);
for i=1:6
    for year=1:3
        for day=1:365
            xmean = temp_pretict_85_m{i,1}{year,1}(3,day); 
            xmin = temp_pretict_85_m{i,1}{year,1}(2,day); 
            xmax = temp_pretict_85_m{i,1}{year,1}(1,day); 
            xeps = 0.01; 
            x = unifrnd(xmin,xmax,[n,1]); 
            while abs(xmean - mean(x)) >= xeps 
                if xmean > mean(x) 
                 x(find(x < xmean,1)) = unifrnd(xmean,xmax); 
                elseif xmean < mean(x) 
                 x(find(x > xmean,1)) = unifrnd(xmin,xmean); 
                end 
            end 
            temp_pretict_85_n{i,1}{year,1}(:,day)=x;
        end
            temp_pretict_85_n{i,1}{year,1}=sort(temp_pretict_85_n{i,1}{year,1},2);
            temp_pretict_85_n{i,1}{year,1}=[temp_pretict_85_n{i,1}{year,1},s_index'];
            temp_pretict_85{i,1}{year,1}=sortrows(temp_pretict_85_n{i,1}{year,1},366);
    end
end
for i=1:6
    for year=1:3
        for day=1:365
            xmean = temp_pretict_45_m{i,1}{year,1}(3,day); 
            xmin = temp_pretict_45_m{i,1}{year,1}(2,day); 
            xmax = temp_pretict_45_m{i,1}{year,1}(1,day); 
            xeps = 0.01; 
             x = unifrnd(xmin,xmax,[n,1]);  
            while abs(xmean - mean(x)) >= xeps 
                if xmean > mean(x) 
                 x(find(x < xmean,1)) = unifrnd(xmean,xmax);  
                elseif xmean < mean(x) 
                 x(find(x > xmean,1)) = unifrnd(xmin,xmean); 
                end 
            end 
            temp_pretict_45_n{i,1}{year,1}(:,day)=x;
        end
        temp_pretict_45_n{i,1}{year,1}=sort(temp_pretict_45_n{i,1}{year,1},2);
        temp_pretict_45_n{i,1}{year,1}=[temp_pretict_45_n{i,1}{year,1},s_index'];
        temp_pretict_45{i,1}{year,1}=sortrows(temp_pretict_45_n{i,1}{year,1},366);
    end
end



%% 第二步 设置训练数据和预测数据
numcom=293;
for i=1:364
    for h=1:24
    temp_data(i,h)=convtemp(temperature_data_mat(i,h),'F','C');
%     temp_data(i,2)=convtemp(min(temperature_data_mat(i,:)),'F','C');
%     temp_data(i,3)=convtemp(mean(temperature_data_mat(i,:)),'F','C');
    end
    for j=1:numcom
        for l=1:24
            load_original{i,1}(j,l)=p_bundle{l,1}(i,j);
            load_data{l,1}(i,j)=load_original{i,1}(j,l);
        end
%             load_data{2,1}(i,j)=min(load_original{i,1}(j,:));
%             load_data{3,1}(i,j)=mean(load_original{i,1}(j,:));
    end
end
% for j=1:numcom
%     for h=1:2
%         x=temperature_data_mat(1:364,a(h))';
%     input_train(1,1:364)=convtemp(x,'F','C');
%     output_train(1,1:364)=p_bundle{a(h)}(1:364,j)';
%% Let us create poly model that can extrapolate the user behaviour to different climate change values.
    
    ML_fit_data = cell(24,1);
%    figure 
    for bus_j = 1:1 % 123 buses
        for h = 1:length(load_data) % max,min,mean
%        ML_fit_data{h,1} = {};
            x = temp_data(:,h);
            y = load_data{h,1}(:,bus_j);
            [p,S] = polyfit(x,y,2);
            
%                     % predict output and see a plot for my verification.
%                     % this y_fit is like a mean and delta is like std for every unique
%                     % temperature. y_fit gives the load consumption values.
%                     [y_fit,delta] = polyval(p,x,S);
%                     subplot(4,6,h)
%                     plot(x,y,'bo')
%                     hold on
%                     plot(x,y_fit,'r-')
%                     plot(x,y_fit+2*delta,'m--',x,y_fit-2*delta,'m--')
%                     xlabel('Temperature(^。C)')
%                     ylabel('Load(kWh)')
%                     title('Linear Fit of Data with 95% Prediction Interval')
%                     legend('Data','Parabola "U" fit','95% Prediction Interval')
            
            % Save the ML fit variables for each hour and consumer.
            ML_fit_data{h,1}{bus_j,1} = cell(2,1);
            ML_fit_data{h,1}{bus_j,1}{1,1} = p;
            ML_fit_data{h,1}{bus_j,1}{2,1} = S;
        end
    end

    %% Let us predict the consumer load behavior for future climate change scenarios.
    
    % PARALLEL PROCESSING CODE
    % 
    % % get p_mat and S_mat for ML fit 
    % p_mat = [];
    % S_mat = [];
    % for i = 1:size(ML_fit_data,1)
    %     for j = 1:busno
    %         
    %     end
    % end
%    losangeles_load=cell(6,1)(3,3);
for i = 1:size(los_angeles_45,2) % total climate change scenarios.
   for each_type = 1:3%max,min,mean
                % Calculate for each bus, p values.
       for each_year = 1:3
           for each_bus = 1:293
               p = ML_fit_data{1,1}{each_bus,1}{1,1};
               S = ML_fit_data{1,1}{each_bus,1}{2,1};
                for each_day = 1:365
                    temp_value_45=temp_pretict_45{i,1}{each_year,1}(each_type,each_day);
                    temp_value_85=temp_pretict_85{i,1}{each_year,1}(each_type,each_day);
                    [y_fit_45,delta_45] = polyval(p,temp_value_45,S);
                    [y_fit_85,delta_85] = polyval(p,temp_value_85,S);
%                     temp_mat_45(each_day,each_bus) = y_fit_45;
%                     temp_mat_85(each_day,each_bus) = y_fit_85;
                    losangeles_load_45{i,1}{each_year,each_type}(each_day,each_bus) = y_fit_45;
                    losangeles_load_85{i,1}{each_year,each_type}(each_day,each_bus) = y_fit_85;
                end
            end
        end
    end
end
scero={'NorESM1-M';'MPI-ESM-MR';'inmcm4';'GFDL-CM3';'CSIRO';'bcc-csm1-1'};
losangeles_load_45=[scero,losangeles_load_45];
losangeles_load_85=[scero,losangeles_load_85];

save('losangeles_load_45.mat','losangeles_load_45');
save('losangeles_load_85.mat','losangeles_load_85');
toc