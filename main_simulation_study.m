%% Generate weakly stationary paths
% define simulation parameters
n_clusters = 5;
obs_num_per_cluster = 10;
alpha = 0.31:0.02:0.39;
total_time_steps = 100;
obs_num_per_step = 5;
total_num_observations = total_time_steps * obs_num_per_step;
total_num_paths = 100;

% simulation weakly stationary stochastic processes
[obs_chain, cluster_ind] = sim_wssp_paths(n_clusters, obs_num_per_cluster, alpha, ...
                         total_num_observations, total_num_paths);
                                         
%% Offline dataset experiments
test_time_steps = 50;  
test_num_sims = 100; 
miscls_rate_offline_algo1 = zeros(test_time_steps, test_num_sims);
miscls_rate_offline_algo2 = zeros(test_time_steps, test_num_sims);
avg_miscls_rate_offline_algo1 = zeros(test_time_steps,1);
avg_miscls_rate_offline_algo2 = zeros(test_time_steps,1);

for t = 1:test_time_steps
    for sim = 1:test_num_sims;
    % parfor sim = 1:test_num_sims;  % parallel computing if necessary
        
        % scale the obsersed times series to be mean 0
        obs = obs_chain(1:(t * obs_num_per_step), :, sim)';
        obs = scale_mean(obs, 0);
        
        % full matrix is observed under offline dataset
        obs_idx = ones(size(obs));
        
        % clustering the observed time series
        [I_chain_algo1, dm] = unsup_wssp_offline_algo(obs, obs_idx, n_clusters);
        [I_chain_algo2, ~] = unsup_wssp_online_algo(obs, obs_idx, n_clusters, dm);
       
        % calculate misclassification rate
        miscls_rate_offline_algo1(t,sim) = misclassify_rate(I_chain_algo1, cluster_ind);
        miscls_rate_offline_algo2(t,sim) = misclassify_rate(I_chain_algo2, cluster_ind);

        fprintf('simluation iter %i for time step %i. \n', sim, t)
    end
    avg_miscls_rate_offline_algo1(t) = mean(miscls_rate_offline_algo1(t,:));
    avg_miscls_rate_offline_algo2(t) = mean(miscls_rate_offline_algo2(t,:));
end

% plot of clsutering results
x = 1:test_time_steps;
figure
plot(x, avg_miscls_rate_offline_algo1, 'b', 'LineWidth', 2)
hold on
plot(x, avg_miscls_rate_offline_algo2, '-.r', 'LineWidth', 2)
hold off
title('Offline Dataset with Covariance Distance Clustering')
xlabel('time step')
ylabel('misclassification rate')
legend('Algorithm 1', 'Algorithm 2')        


%% Online dataset experiments
test_time_steps = 20;  
test_num_sims = 5; 
miscls_rate_online_algo1 = zeros(test_time_steps, test_num_sims);
miscls_rate_online_algo2 = zeros(test_time_steps, test_num_sims);
avg_miscls_rate_online_algo1 = zeros(test_time_steps,1);
avg_miscls_rate_online_algo2 = zeros(test_time_steps,1);

for t = 1:test_time_steps
    parfor sim = 1:test_num_sims;
    % parfor sim = 1:test_num_sims;  % parallel computing if necessary
        
        % scale the obsersed times series to be mean 0
        obs = obs_chain(1:(t * obs_num_per_step), :, sim)';
        
        % full matrix is observed under offline dataset
        obs_idx = zeros(size(obs));
        for i = 1:obs_num_per_cluster:(n_clusters * obs_num_per_cluster)
            obs_idx(i:(i+4), :) = 1;
            for j = 1:(obs_num_per_cluster - 5)
                if t > j * 10
                    obs_idx(5+j, (j * 10 * obs_num_per_step + 1):end) = 1;
                end
            end
        end
        
        keep_idx = find(sum(obs_idx,2) ~= 0);
        obs = obs(keep_idx, :);
        obs_idx = obs_idx(keep_idx, :);
        cluster_ind_online = cluster_ind(keep_idx);
        
        % clustering the observed time series
        [I_chain_algo1, dm] = unsup_wssp_offline_algo(obs, obs_idx, n_clusters);
        [I_chain_algo2, ~] = unsup_wssp_online_algo(obs, obs_idx, n_clusters, dm);
       
        % calculate misclassification rate
        miscls_rate_offline_algo1(t,sim) = misclassify_rate(I_chain_algo1, cluster_ind_online);
        miscls_rate_offline_algo2(t,sim) = misclassify_rate(I_chain_algo2, cluster_ind_online);

        fprintf('simluation iter %i for time step %i. \n', sim, t)
    end
    avg_miscls_rate_offline_algo1(t) = mean(miscls_rate_offline_algo1(t,:));
    avg_miscls_rate_offline_algo2(t) = mean(miscls_rate_offline_algo2(t,:));
end

% plot of clsutering results
x = 1:test_time_steps;
figure
plot(x, avg_miscls_rate_offline_algo1, 'b', 'LineWidth', 2)
hold on
plot(x, avg_miscls_rate_offline_algo2, '-.r', 'LineWidth', 2)
hold off
title('Online Dataset with Covariance Distance Clustering')
xlabel('time step')
ylabel('misclassification rate')
legend('Algorithm 1', 'Algorithm 2')    