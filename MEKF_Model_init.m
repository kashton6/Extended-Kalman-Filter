%%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Set Initial Conditions
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 
t = [0:0.1:300];
w = length(t);
 %%% MEKF initialize values %%%
 sig_v = sqrt(10)*1e-10;% gyro sensor noise
sig_u = sqrt(10)*1e-3;
sun_sensor_var = 0.04; %sun sensor noise variance

    del_t = 0.1;
    Beta_init = 0.01*[6;5;4]; %Initialize gyro Bias

    Q_k = [(sig_v^2*del_t + 1/3*sig_u^2*del_t^3)*eye(3), (0.5*sig_u^2*del_t^2)*eye(3);
        (0.5*sig_u^2*del_t^2)*eye(3), (sig_u^2*del_t)*eye(3)];
    gamma = [-eye(3),zeros(3,3);zeros(3,3),eye(3)];
    Q_sq = chol((gamma*Q_k*gamma'),'lower');
    
    %initilize the quaternion value
    qe_init = [0.25;0.25;0.25;0.98];
   
    %initialize covariance matrix and its sqrt form
    P = 0.001*eye(6);
    P_sq = chol(P,'lower');
    P_sq_init = P_sq;
    
    %initialize Sensor Noise matrix and its sqrt form
    R = sun_sensor_var^2*eye(3);
    R_sq = chol(R,'lower');
    
    [np,mp] = size(P);
    [nr,mr] = size(R_sq);

    
tfinal = t(w);

% Simulation with digital controller
simout1=sim('MEKF_Model','StopTime','tfinal', ...
    'SaveTime','on','TimeSaveName','timeoutNew',...
    'SaveOutput','on','OutputSaveName','youtNew');
time2=simout1.get('timeoutNew');
y=simout1.get('youtNew');
qest_simu=y{1}.Values.Data';

figure
hold on
subplot(2,2,1)
plot(t,qest_simu(1,:),'g')
subplot(2,2,2)
plot(t,qest_simu(2,:),'r')
subplot(2,2,3)
plot(t,qest_simu(3,:),'b')
subplot(2,2,4)
plot(t,qest_simu(4,:),'m')
title('q1 Evolution/Estimates')
xlabel('Time (min)')
ylabel('q1 (rad)')
grid on
legend('Estimated q1 (simulink)','Estimated q2 (simulink)','Estimated q3 (simulink)','Estimated q4 (simulink)''Location','southwest')
    

% Save Workspace.
filename = 'simData';
save('simData')
