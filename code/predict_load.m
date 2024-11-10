%% To run this program successfully you need to comment out steps 2) and 3).

%% 1) input historical temperature and electricity load
% The local DG system data is sensitive information regarding the infrastructure security of 
% the urban area and is subject to the non-disclosure agreement (NDA).
% We can not disclose all utility's load curve. Here, 'p_bundle' can not be loaded. 
% But we will disclose the fitted curves for electricity load and temperature for 10 utilities.
numcom=33;
for i=1:364
    for h=1:24
    temp_data(i,h)=temperature_data_mat(i,h);
%     temp_data(i,2)=convtemp(min(temperature_data_mat(i,:)),'F','C');
%     temp_data(i,3)=convtemp(mean(temperature_data_mat(i,:)),'F','C');
    end
    for j=1:numcom
        for h=1:24
            load_original{i,1}(j,h)=p_bundle{h,1}(i,j);
            load_data{h,1}(i,j)=load_original{i,1}(j,h);
        end
%             load_data{2,1}(i,j)=min(load_original{i,1}(j,:));
%             load_data{3,1}(i,j)=mean(load_original{i,1}(j,:));
    end
end

%% 2) Let us create poly model that can extrapolate the user behaviour to different climate change values.
    
%     ML_fit_data = cell(24,1);
%    figure 
    for bus_j = 1:33 % 123 buses
        for h = 1:length(load_data) % max,min,mean
%        ML_fit_data{h,1} = {};
            x = temp_data(:,h);
            y = load_data{h,1}(:,bus_j);
            [p,S] = polyfit(x,y,2);
            
            % Save the ML fit variables for each hour and consumer.
            ML_fit_data{h,1}{bus_j,1} = cell(2,1);
            ML_fit_data{h,1}{bus_j,1}{1,1} = p;
            ML_fit_data{h,1}{bus_j,1}{2,1} = S;
        end
    end
    

%% 3) Create a different climate change scenario for risk comparison results.
clear
tic


load('historical_temperature.mat');  
temp_list=csvread('US06037_MeanModel_english.csv'); % Climate prediction data from USGS NCCV database

n = 24; 
[x_sort s_index1]=sort(temperature_data_mat(14,:));
[x_s2 s_index]=sort(s_index1);
for year=1:150
        for month=1:12
            xmean = temp_list(12*(year-1)+month,1);
            xmax = temp_list(12*(year-1)+month,2); 
            xmin = temp_list(12*(year-1)+month,3); 
            xeps = 1; 
            x = unifrnd(xmin,xmax,[1,n]); 
            while abs(xmean - mean(x)) >= xeps 
                if xmean > mean(x) 
                 x(find(x < xmean,1)) = unifrnd(xmean,xmax); 
                elseif xmean < mean(x) 
                 x(find(x > xmean,1)) = unifrnd(xmin,xmean); 
                end 
            end 
            temp_pretict_45_n{year,1}(month,:)=sort(x);
           
        for i=1:24
            temp_pretict_45{year,1}{month,i}=temp_pretict_45_n{year,1}(month,s_index(i));
        end
        end
end
for year=1:150
        for month=1:12
            xmean = temp_list(12*(year-1)+month,4);
            xmax = temp_list(12*(year-1)+month,5); 
            xmin = temp_list(12*(year-1)+month,6); 
            xeps = 1; 
            x = unifrnd(xmin,xmax,[1,n]); 
            while abs(xmean - mean(x)) >= xeps 
                if xmean > mean(x) 
                 x(find(x < xmean,1)) = unifrnd(xmean,xmax); 
                elseif xmean < mean(x) 
                 x(find(x > xmean,1)) = unifrnd(xmin,xmean); 
                end 
            end 
            temp_pretict_85_n{year,1}(month,:)=sort(x);
           
        for i=1:24
            temp_pretict_85{year,1}{month,i}=temp_pretict_85_n{year,1}(month,s_index(i));
        end
        end
end

    %% 4) Let us predict the consumer load behavior for future climate change scenarios.
    
    % PARALLEL PROCESSING CODE
    % Therefore, we can only disclose fitted curves for electricity load
    % and temperature for 10 utilities. They are then randomly distributed over the 123 nodes.
    load ('fitted_curve.mat')
   for each_hour = 1:24%max,min,mean
                % Calculate for each bus, p values.
       for each_year = 1:150
           for each_bus=1:33
               selet_user=randperm(10,1);
               p = ML_fit_data{each_hour,1}{selet_user,1}{1,1};
               S = ML_fit_data{each_hour,1}{selet_user,1}{2,1};
                for each_month = 1:12
                    temp_value_45=temp_pretict_45{each_year,1}{each_month,each_hour};
                    [y_fit_45,delta_45] = polyval(p,temp_value_45,S);
                    load_1950_2099{1,1}{each_year,each_month}(each_hour,each_bus) = y_fit_45;
                     temp_value_85=temp_pretict_85{each_year,1}{each_month,each_hour};
                    [y_fit_85,delta_85] = polyval(p,temp_value_85,S);
                    load_1950_2099{2,1}{each_year,each_month}(each_hour,each_bus) = y_fit_85;
                end
           end
        end
   end
scero={'RCP4.5';'RCP8.5'};
load_1950_2099=[scero,load_1950_2099];
save('load_1950_2099.mat','load_1950_2099');
for type=1:2
    for year=1:150
        for season=1:3
            for m=1:3
                load_1950_2099_org{type,season}{year,1}(24*(m-1)+1:24*m,:)=load_1950_2099{type,2}{year,3*(season-1)+m+2};
            end
        end
        for m=1:2
        load_1950_2099_org{type,4}{year,1}(24*(m-1)+1:24*m,:)=load_1950_2099{type,2}{year,m};
        end
        load_1950_2099_org{type,4}{year,1}(49:72,:)=load_1950_2099{type,2}{year,12};
    end
end
season={'','spring','summer','fall','winter'};
load_1950_2099_org=[scero,load_1950_2099_org];
load_1950_2099_org=[season;load_1950_2099_org];
save('load_1950_2099_org.mat','load_1950_2099_org')

for type=1:2
for year=1:150
    for month=1:12
        for bus=1:33
            load_1950_2099_maxminmean{type,1}{year,1}(month,bus)=max(load_1950_2099{type,2}{year,month}(:,bus));
            load_1950_2099_maxminmean{type,2}{year,1}(month,bus)=min(load_1950_2099{type,2}{year,month}(:,bus));
            load_1950_2099_maxminmean{type,3}{year,1}(month,bus)=mean(load_1950_2099{type,2}{year,month}(:,bus));
        end
    end
end
end
load_1950_2099_maxminmean=[scero,load_1950_2099_maxminmean];
save('load_1950_2099_maxminmean.mat','load_1950_2099_maxminmean')
toc
