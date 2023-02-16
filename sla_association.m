%% Initializations

sum_ov_num1 = 0;
sum_ov_num2 = 0;
sum_ov_num3 = 0;
sum_ov_num4 = 0;

mean_dist1 = 0;
mean_dist2 = 0;
mean_dist3 = 0;
mean_dist4 = 0;

% number of iterations over which results are averaged
repeats = 1;

for m = 1:repeats

    load('topology4_uni.mat')    

    b = 0.5;                                            % SLA learning rate
    sla_ite = 0;                                        % SLA iterations counter
    convergence = 0;                                    % SLA convergence flag
    action_prob(1:users,1:femto) = 1/femto;             % SLA initial action probabilities
    user_feedback = zeros(1,users);                     % user's feedback                   

    chosen_femto = zeros(1,users);                      % users' chosen femto cell
    sum_user_contract_rewards = zeros(users,femto);     % users' sum of contract rewards for each femto
    norm_user_contract_rewards = zeros(users,femto);    % normalized users' sum of contract rewards for each femto

    % useful variables for results evaluation
    plot_ov_number = [];
    plot_ov_distances = [];
    plot_user_feedback = [];
    plot_mean_feedback = [];

    sum_band2 = 0;

    %% SLA Loop - Users choose the femto in order to associate with

    while convergence == 0

        % Users choose their femto cell according to their action probabilities
        for i = 1:users
            chosen_femto(i) = randsample(1:femto,1,true,action_prob(i,:));
        end

        % System's scatter plot with initial associations
        if sla_ite == 0
            scatter_plot(femto,x_femto,y_femto,z_femto,users,x_users,y_users,z_users,chosen_femto,sla_ite);
        end

        % -------- Contract theoretic power allocation in each cell ----------

        overall_efforts = zeros(1,users);                       % users' effort
        overall_powers = zeros(1,users);                        % users' power
        overall_rewards = zeros(1,users);                       % users' reward 

        overall_types = zeros(1,users);                         % user's contract type
        overall_gains = zeros(1,users);                         % user's channel gain
        overall_distances = zeros(1,users);                     % user's distance from the BS

        overall_cell_utilities = zeros(1,users);                % cell utilities
        overall_user_utilities = zeros(1,users);                % user utilities

        for i = 1:femto 
            users_femto = [];                           
            channels_femto = [];
            for j = 1:users
                if chosen_femto(j) == i
                    users_femto = [users_femto j];
                    channels_femto = [channels_femto channel_gains(j,i)];
                end
            end

            % execute the contract only if there are associated users
            if length(users_femto) > 0
                % sort users according to their channel gain 
                [channels_femto,indeces] = sort(channels_femto);      
                users_femto = users_femto(indeces);

                % find the optimal contract items for the users 
                [p,q,r,c_util,u_util,types] = continuum_contract(length(users_femto),channels_femto,noise(i),bandwidth(i));

                % append the contract-related variables for all cells in concentrated vectors
                for j = 1:length(users_femto)
                    overall_powers(users_femto(j)) = p(j);
                    overall_efforts(users_femto(j)) = q(j);
                    overall_rewards(users_femto(j)) = r(j);
                    overall_cell_utilities(users_femto(j)) = c_util(j);
                    overall_user_utilities(users_femto(j)) = u_util(j);
                    overall_types(users_femto(j)) = types(j);
                    overall_gains(users_femto(j)) = channels_femto(j);
                    overall_distances(users_femto(j)) = distances(users_femto(j),i);

                    sum_user_contract_rewards(users_femto(j),i) = sum_user_contract_rewards(users_femto(j),i) + r(j);
                end  
            end
            norm_user_contract_rewards(:,i) = (sum_user_contract_rewards(:,i) - min(sum_user_contract_rewards(:,i)))*1 / (max(sum_user_contract_rewards(:,i)) - min(sum_user_contract_rewards(:,i))) + 1;
        end

        % ----------- End of contract theoretic power allocation -------------

        % Useful variables
        user_distance = zeros(1,users);                 % users' distance from the associated cell
        ov_number = zeros(1,femto);                     % number of users per femto cell
        ov_distances = zeros(1,femto);                  % sum of users' distances per femto cell
        ov2_distances = zeros(1,femto);                 % sum of square users' distances per femto cell
        max_distances = zeros(1,femto);                 % max of users' distances per femto cell
        min_distances = 1000*zeros(1,femto);            % min of users' distances per femto cell
        for i = 1:femto
            for j = 1:users
                if chosen_femto(j) == i
                    user_distance(j) = distances(j,i);
                    ov_number(i) = ov_number(i) + 1;
                    ov_distances(i) = ov_distances(i) + distances(j,i);
                    ov2_distances(i) = ov2_distances(i) + (1/distances(j,i));
                    if distances(j,i) > max_distances(i)
                       max_distances(i) = distances(j,i);
                    end
                   if distances(j,i) < min_distances(i)
                       min_distances(i) = distances(j,i);
                   end
                end
            end
        end   

        for j = 1:femto
            sum_band2 = sum_band2 + bandwidth(j)^2;
        end

        % Update users' SLA feedback
        sum_user_feedback = zeros(1,femto);
        max_user_feedback = zeros(1,femto);
        for i = 1:users
            for j = 1:femto
                if chosen_femto(i) == j
                   user_feedback(i) = (norm_user_contract_rewards(i,j)) * sqrt((1/user_distance(i))/ov2_distances(j)) * sqrt(bandwidth(j)/sum(bandwidth));
                   sum_user_feedback(j) = sum_user_feedback(j) + user_feedback(i);
                   if user_feedback(i) > max_user_feedback(j)
                       max_user_feedback(j) = user_feedback(i);
                   end
                end
            end
        end
        
        % Update the SLA action probabilities
        for i = 1:users
            for j = 1:femto
               if chosen_femto(i) == j
                   action_prob(i,j) = action_prob(i,j) + b*user_feedback(i)*(1-action_prob(i,j));
               else
                   action_prob(i,j) = action_prob(i,j) - b*user_feedback(i)*action_prob(i,j);
               end
            end    
        end

        % Check for SLA convergence 
        count = 0;
        for i = 1:users
            if max(action_prob(i,:)) >= 0.95
               count = count + 1;
            end
        end
  
        if count == users
            convergence = 1;
            scatter_plot(femto,x_femto,y_femto,z_femto,users,x_users,y_users,z_users,chosen_femto,sla_ite);
        else
            sla_ite = sla_ite + 1;
        end

        % Save data for plots
        plot_ov_number = [plot_ov_number ; ov_number];
        plot_ov_distances = [plot_ov_distances ; ov_distances];
        plot_user_feedback = [plot_user_feedback ; user_feedback];
        plot_mean_feedback = [plot_mean_feedback ; mean(user_feedback)];

    end
    
    sum_ov_num1 = sum_ov_num1 + ov_number(1);
    sum_ov_num2 = sum_ov_num2 + ov_number(2);
    sum_ov_num3 = sum_ov_num3 + ov_number(3);
    sum_ov_num4 = sum_ov_num4 + ov_number(4);
    
    mean_dist1 = mean_dist1 + ov_distances(1)/ov_number(1);
    mean_dist2 = mean_dist2 + ov_distances(2)/ov_number(2);
    mean_dist3 = mean_dist3 + ov_distances(3)/ov_number(3);
    mean_dist4 = mean_dist4 + ov_distances(4)/ov_number(4);

end

%% SLA evaluation results

% Plot the change in the number of associated users with each cell at each
% iteration
figure();
plot(1:sla_ite+1,plot_ov_number(1:sla_ite+1,:));
xlabel('SLA iteration')
ylabel('Number of associated users')
legend('UAV 1','UAV 2','UAV 3','UAV 4')
xlim([1 inf])

% Plot the change in the sum of users' distances from the center of each cell
% at each iteration
figure();
plot(1:sla_ite+1,plot_ov_distances(1:sla_ite+1,:));
xlabel('SLA iteration')
ylabel('Sum of user distances')
legend('UAV 1','UAV 2','UAV 3','UAV 4')
xlim([1 inf]) 

% Plot the change in each user's feedback at each iteration
figure();
plot(1:sla_ite+1,plot_user_feedback(1:sla_ite+1,:));
xlabel('SLA iteration')
ylabel('User''s feedback')
xlim([1 inf]) 

% Plot the change in the overall users mean feedback at each iteration 
% (for better visualization purposes)
figure();
plot(1:sla_ite+1,plot_mean_feedback(1:sla_ite+1));
xlabel('SLA iteration')
ylabel('Overall users mean feedback')
xlim([1 inf]) 

% Plot the sum of users' distances from the center of each cell
figure();
bar(ov_distances)
xlabel('Femto cell')
ylabel('Sum of user distances')

% Plot the number of associated users with each cell
figure();
bar(ov_number);
xlabel('Femto cell')
ylabel('Number of associated users')
set(gca,'FontSize',20);


% Plot the number of associated users with each cell
color = [0 0.4470 0.7410 ; 0.8500 0.3250 0.0980 ; 0.9290 0.6940 0.1250 ; 0.4940 0.1840 0.5560];

figure();
x = categorical({'UAV1','UAV2','UAV3','UAV4'});
x = reordercats(x,{'UAV1','UAV2','UAV3','UAV4'});
b = bar(x,ov_number,'FaceColor','flat');
b.CData(1,:) = color(1,:);
b.CData(2,:) = color(2,:);
b.CData(3,:) = color(3,:);
b.CData(4,:) = color(4,:);
ylabel('Associated users'' number','FontSize',20);
set(gca,'FontSize',20);


%% Contract evaluation results

% Create a struct to hold the results for all cells
% Then, execute contract_plots.m to plot the results for each cell

cell = struct;
cell(femto).distances = [];
cell(femto).types = [];
cell(femto).gains = [];
cell(femto).powers = [];
cell(femto).efforts = [];
cell(femto).rewards = [];
cell(femto).cell_utilities = [];
cell(femto).user_utilities = [];

for i = 1:femto
    for j = 1:users
        if chosen_femto(j) == i
            cell(i).distances = [cell(i).distances overall_distances(j)];
            cell(i).types = [cell(i).types overall_types(j)];
            cell(i).gains = [cell(i).gains overall_gains(j)];        
            cell(i).powers = [cell(i).powers overall_powers(j)];
            cell(i).efforts = [cell(i).efforts overall_efforts(j)];
            cell(i).rewards = [cell(i).rewards overall_rewards(j)];           
            cell(i).cell_utilities = [cell(i).cell_utilities overall_cell_utilities(j)];
            cell(i).user_utilities = [cell(i).user_utilities overall_user_utilities(j)];     
        end
    end
    
    % Sort users within the cell according to their types
    [cell(i).types,indeces] = sort(cell(i).types); 
    cell(i).distances = cell(i).distances(indeces);
    cell(i).gains = cell(i).gains(indeces);
    cell(i).powers = cell(i).powers(indeces); 
    cell(i).efforts = cell(i).efforts(indeces);    
    cell(i).rewards = cell(i).rewards(indeces); 
    cell(i).cell_utilities = cell(i).cell_utilities(indeces);
    cell(i).user_utilities = cell(i).user_utilities(indeces);    
end


%% Monte carlo results

color = [0 0.4470 0.7410 ; 0.8500 0.3250 0.0980 ; 0.9290 0.6940 0.1250 ; 0.4940 0.1840 0.5560];

b = [mean_dist1/repeats mean_dist2/repeats mean_dist3/repeats mean_dist4/repeats];

figure();
x = categorical({'UAV1','UAV2','UAV3','UAV4'});
x = reordercats(x,{'UAV1','UAV2','UAV3','UAV4'});
b = bar(x,b,'FaceColor','flat');
b.CData(1,:) = color(1,:);
b.CData(2,:) = color(2,:);
b.CData(3,:) = color(3,:);
b.CData(4,:) = color(4,:);
ylabel('Mean users'' distance [m]','FontSize',20);
set(gca,'FontSize',20);

color = [0 0.4470 0.7410 ; 0.8500 0.3250 0.0980 ; 0.9290 0.6940 0.1250 ; 0.4940 0.1840 0.5560];

b = [round(sum_ov_num1/repeats) round(sum_ov_num2/repeats) round(sum_ov_num3/repeats) round(sum_ov_num4/repeats)];

figure();
x = categorical({'UAV1','UAV2','UAV3','UAV4'});
x = reordercats(x,{'UAV1','UAV2','UAV3','UAV4'});
b = bar(x,b,'FaceColor','flat');
b.CData(1,:) = color(1,:);
b.CData(2,:) = color(2,:);
b.CData(3,:) = color(3,:);
b.CData(4,:) = color(4,:);
ylabel('Associated users'' number','FontSize',20);
set(gca,'FontSize',20);



