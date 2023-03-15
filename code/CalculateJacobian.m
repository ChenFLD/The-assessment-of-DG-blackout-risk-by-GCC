
function [Jacob] = CalculateJacobian(no_bus,bus_matrix_G,bus_matrix_B,v_kr,v_ki)
%Forming dP_dVr nXn Matrix
Apr = zeros(no_bus,no_bus);

for i = 1:1:no_bus
    for j = 1:1:no_bus
        if i == j
            t_d1 = sum(-bus_matrix_G(i,:));
            V_real_G_dk = transpose(v_kr)*transpose(bus_matrix_G(i,:));
            V_imag_B_dk = transpose(v_ki)*transpose(bus_matrix_B(i,:));
            t_d2 = (V_real_G_dk - V_imag_B_dk);
            Apr(i,j) = 2*t_d1*v_kr(i)+ t_d2; % d=k condition
        else
            V_real_G_dk = v_kr(i)*bus_matrix_G(i,j);
            V_imag_B_dk = v_ki(i)*bus_matrix_B(i,j);
            Apr(i,j) = V_real_G_dk + V_imag_B_dk; %d~= k condition
        end
    end
end

%Forming dP_dVi nXn Matrix
Api = zeros(no_bus,no_bus);

for i = 1:1:no_bus
    for j = 1:1:no_bus
        if i == j
            t_d1 = sum(-bus_matrix_G(i,:));
            V_real_B_dk = transpose(v_kr)*transpose(bus_matrix_B(i,:));
            V_imag_G_dk = transpose(v_ki)*transpose(bus_matrix_G(i,:));
            t_d3 = (V_real_B_dk + V_imag_G_dk);
            Api(i,j) = 2*t_d1*v_ki(i)+ t_d3; % d=k condition
        else
            V_real_B_dk = v_kr(i)*bus_matrix_B(i,j);
            V_imag_G_dk = v_ki(i)*bus_matrix_G(i,j);
            Api(i,j) = V_imag_G_dk - V_real_B_dk; %d~= k condition
        end
    end
end

%Forming dQ_dVr nXn Matrix
Aqr = zeros(no_bus,no_bus);

for i = 1:1:no_bus
    for j = 1:1:no_bus
        if i == j
            t_d4 = sum(bus_matrix_B(i,:));
            V_real_B_dk = transpose(v_kr)*transpose(bus_matrix_B(i,:));
            V_imag_G_dk = transpose(v_ki)*transpose(bus_matrix_G(i,:));
            t_d3 = (V_real_B_dk + V_imag_G_dk);
            Aqr(i,j) = 2*t_d4*v_kr(i) - t_d3; % d=k condition
        else
            V_real_B_dk = v_kr(i)*bus_matrix_B(i,j);
            V_imag_G_dk = v_ki(i)*bus_matrix_G(i,j);
            Aqr(i,j) = V_imag_G_dk - V_real_B_dk; %d~= k condition
        end
    end
end

%Forming dQ_dVi nXn Matrix
Aqi = zeros(no_bus,no_bus);

for i = 1:1:no_bus
    for j = 1:1:no_bus
        if i == j
            t_d4 = sum(bus_matrix_B(i,:));
            V_real_G_dk = transpose(v_kr)*transpose(bus_matrix_G(i,:));
            V_imag_B_dk = transpose(v_ki)*transpose(bus_matrix_B(i,:));
            t_d2 = (V_real_G_dk - V_imag_B_dk);
            Aqi(i,j) = 2*t_d4*v_ki(i) + t_d2; % d=k condition
        else
            V_real_G_dk = v_kr(i)*bus_matrix_G(i,j);
            V_imag_B_dk = v_ki(i)*bus_matrix_B(i,j);
            Aqi(i,j) = -V_real_G_dk - V_imag_B_dk; %d~= k condition
        end
    end
end

%Forming Jacobian matrix
Jacob = zeros(2*no_bus,2*no_bus);
for i = 1:1:2*no_bus
    for j = 1:1:2*no_bus
        if (i <= no_bus)&&(j <= no_bus)
            Jacob(i,j) = Apr(i,j);
        elseif (i <= no_bus)&&(j > no_bus)
            Jacob(i,j) = Api(i,j - no_bus);
        elseif (i > no_bus) && (j <= no_bus)
            Jacob(i,j) = Aqr(i - no_bus,j);
        elseif (i > no_bus)&&(j > no_bus)
            Jacob(i,j) = Aqi(i - no_bus,j - no_bus);
        end
    end
end


% Reshaping the Jacobian matrix to a shape of 2*(n-1) X 2*(n-1).
% Removing the slack bus from the Jacobian.
Jacob(1,:) = [];
Jacob(no_bus,:) = [];
Jacob(:,1) = [];
Jacob(:,no_bus) = [];
end