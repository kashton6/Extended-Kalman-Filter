clear all
close all
% clc
%%%%Create the q and omega vectors simulated data

tspan = [0:0.1:300];
% for qqq = 1:1
J = [3,-2,-1;-2,3,-1;-1,-1,4];
p0 = [0;0;0;1;0.05;0;0]; % initial conditions q0 concatinated with w0
K = @(t,x) [(0.5.*([x(4), -x(3), x(2); x(3), x(4), -x(1); -x(2), x(1), x(4); -x(1), -x(2), -x(3)]*[x(5); x(6); x(7)]));
    (inv(J)*[0.025.*sin(0.2.*t);0;0] - inv(J)*[0, -x(7), x(6); x(7), 0 , -x(5); -x(6), x(5), 0]*J*[x(5); x(6); x(7)])];
% tic
[t,x1] = ode45(K, tspan, p0);
% toc
Mag = [];
Mag = sqrt((x1(:,1)).^2 + (x1(:,2)).^2 + (x1(:,3)).^2 + (x1(:,4)).^2);
w = length(t);
Rb = zeros(3,length(t));
sig_v = sqrt(10)*1e-10;% gyro sensor noise
sig_u = sqrt(10)*1e-9;
sun_sensor_var = 0.04;
w = length(t);
M = 1; %number of monte carlo iterations
error = zeros(M,w);
% tic
% for j = 1:M
    % for loop below creates matrix CIB for each q vector at each time interval
    %generate quaternions and data, both noisy and true
    tic
    for ii = 1:length(t)
        qq1 = x1(ii,1);
        qq2 = x1(ii,2);
        qq3 = x1(ii,3);
        qq4 = x1(ii,4);
        CIB = [(1-2.*(qq2.^2 + qq3.^2)) ,(2.*qq1.*qq2 + 2.*qq3.*qq4), (2.*qq1.*qq3 - 2.*qq2.*qq4);
            (2.*qq1.*qq2 - 2.*qq3.*qq4) ,(1 - 2.*(qq1.^2 + qq3.^2)), 2.*qq2.*qq3 + 2.*qq1.*qq4;
            (2.*qq1.*qq3 + 2.*qq2.*qq4) ,(2.*qq2.*qq3 -2.*qq1.*qq4), (1-2.*(qq1.^2 + qq2.^2))];
        
        
        eps = (sun_sensor_var)*randn(3,w);%Gausian noise
        Ri = [0;1;0];
        Rb(:,ii) = (CIB*Ri) + eps(:,1); %noisy values
        Rb_hat(:,ii) = (CIB*Ri); %true values
        noise = (sig_v)*randn(3,w);
        vk(:,ii) = noise(:,1);
         
%         ii = ii+1;
%     Qx = [0 -qe(3) qe(2);qe(3) 0 -qe(1);-qe(2) qe(1) 0]
    end
    Ang = acosd(Rb);
    omega_tilda = x1(:,5:7)+ [0.0625;0.0625;0.0625]' + vk'; %measured omega
    
%  for j = 1:M        
    %%%initialize stuff
for j = 1:1
    %%% EKF initialize values %%%
    eyez = sparse(eye(3)); %sparse 3x3 identity matrix
    zero = zeros(3,3);
    del_t = t(2)-t(1); %t(2)-t(1);
    R = sun_sensor_var^2*eye(3);
    Beta = 0.01*[6;5;4]; %Initialize gyro Bias
    Beta_init = Beta;
%     del_x = 0; %initializing delta x_k minus
%     r_ix = [0 -Ri(3) Ri(2);Ri(3) 0 -Ri(1);-Ri(2) Ri(1) 0]; %the direction to
    %the Sun is given by the following unit vector in the inertial-frame:
    Q_k = [(sig_v^2*del_t + 1/3*sig_u^2*del_t^3)*eye(3), (0.5*sig_u^2*del_t^2)*eye(3);
        (0.5*sig_u^2*del_t^2)*eye(3), (sig_u^2*del_t)*eye(3)];
    gamma = [-eye(3),zeros(3,3);zeros(3,3),eye(3)];
    q1 = 0.01; %initilize the quaternion value
    q2 = 0.01;
    q3 = -0.02;
    q4 = 1;
    qe = [q1;q2;q3;q4];
    qe_init = qe;
    P = 0.0001*eye(6);%initialize covariance matrix
%     Pinv = inv(P);
    P_sq = chol(P,'lower');
    P_sq_init = P_sq;
    R_sq = chol(R,'lower');
    [np,mp] = size(P);
    [nr,mr] = size(R_sq);
     Q_sq = chol((gamma*Q_k*gamma'),'lower');
%     Rinv = inv(R);
%     Vc = chol(Rinv,'lower');
%     
%     [Z,E] = eig(full(Q_k));
end
% Z*E*Z'-Q_k;

%
%%% Squaqre root Form P matrix
    %%%EKF routine
   
%     Q_sq = chol(Q_k);
%     tic
    
    for i = 1:w
        Q = [qe(1:3,1)]; %change these q's to be q(1,1)
        Qx = [0 -qe(3) qe(2);qe(3) 0 -qe(1);-qe(2) qe(1) 0]; %cross product matrix of q
        Xi_q = [qe(4)*eye(3)+Qx;-Q'];
        Phi_q = [qe(4)*eye(3)-Qx;-Q'];
        A_q = Xi_q'*Phi_q;   %A(q_k minus) matrix
        
        mat = A_q*Ri;
        mat_x = [0 -mat(3) mat(2);mat(3) 0 -mat(1);-mat(2) mat(1) 0]; %cross product matrix of q
        %%%Gain
        H = [mat_x, zero]; %sensitivity matrix
        
%         P_sq = chol(P,'lower');
        A = [R_sq',zeros(nr,np);(P_sq')*H', P_sq'];
        [TA,UA] = qr(A);
        P_sq = UA(4:9,4:9)'; %Updated P_sq
        W = UA(1:3,4:9)';
%        
        sq_Bk = UA(1:3,1:3)';
        K = W*inv(sq_Bk);
%         K = (P*H')*inv(H*P*H' + R);
        
        
        %%%Update
        full_P  = P_sq*P_sq';
        vestt(:,i) = diag(full_P).^(0.5);
%        

%         rb = Rb(:,i);
        del_x = K*(Rb(:,i)-A_q*Ri);
        ytil(:,i) = A_q*Ri;
        
        % q = [q1;q2;q3;q4];
        del_alpha = del_x(1:3,1);
        del_beta = del_x(4:6,1);
        
        qe = qe + 0.5*Xi_q*del_alpha;
        
        qe = 1/(norm(qe))*qe; %renormalize
        qest(:,i) = qe;
         %store the quaternion
%         constraint = qe'*qe; %should equal 1 or very close, normalize
        Beta = Beta + del_beta;
        
        
        %%%Propagation%%%
        omega_meas = omega_tilda(i,:);
        omega = omega_tilda(i,:)' - Beta;
        omega_est(:,i) = omega;
        om_mag = norm(omega); %norm of omega estimate
        sin_omega = sin(0.5*om_mag*del_t);
        
        Psi = (sin_omega*omega)/om_mag;
        Psi_x = [0 -Psi(3) Psi(2);Psi(3) 0 -Psi(1);-Psi(2) Psi(1) 0];
        
        Omega_omega = [(cos(0.5*om_mag*del_t)*eyez - Psi_x), Psi; -Psi', cos(0.5*om_mag*del_t)];
        cos_omega = (1-cos(om_mag*del_t))/om_mag^2;
        
        omega_x = [0 -omega(3) omega(2);omega(3) 0 -omega(1);-omega(2) omega(1) 0];
        Phi_discrete = [(eyez-omega_x*sin(om_mag*del_t)/om_mag + (omega_x^2)*cos_omega),omega_x*cos_omega-eyez*del_t-omega_x^2*(om_mag*del_t-sin(om_mag*del_t))/om_mag^3;
            zero,eyez;];
        
          qe = Omega_omega*qe; %quaternion propagation
          
          
%         P_sq_plus = chol(P,'lower');
       
          C = [P_sq'*Phi_discrete';Q_sq'];
          [Tc,Uc] = qr(C);
%           Uc = full(Uc);
          sq_P_kp1 = Uc(1:6,1:6)';
%           P = sq_P_kp1'*sq_P_kp1;
          P_sq = sq_P_kp1;
          
%         P = Phi_discrete*P*Phi_discrete' + gamma*Q_k*gamma';
        %covariance Propagation
%         error_norm(:,i) = norm((ytil(:,i)-Rb_hat(:,i)),2); %tracking estimate error
%         error_norm_noise(:,i) = norm((Rb(:,i)-Rb_hat(:,i)),2);
%         
    end
%     toc
%     error(j,:) = error_norm;
%     error_noise(j,:) = error_norm_noise;
    
    
%  end
toc

% toc
% tic
% % [u,d] = udu(P);
% k = 1;
% % toc
% Monte_norm = mean(error,1);
% % erromonte = mean(error_noise,1);
% % semilogx(Monte_norm)
% figure
% hold on
% semilogy(t/60,Monte_norm*180/pi);
% % semilogy(t/60,erromonte*180/pi);
% grid on
% xlabel('Time (min)')
% ylabel('Error Norm')
% set(gca, 'XScale', 'linear')

% qerr=quat_err(qest',x1(:,1:4));erre=asin(qest(1:3,:))*2*180/pi;

figure
subplot(2,1,1)
hold on
% plot(t,x1(:,3),'k*',t,qest(3,:)) %t,x1(:,2),'r',t, x1(:,3),'g',t,x1(:,4),t,Mag)
% plot(t,x1(:,7),'b*',t,omega_est(3,:)) %angular velo
plot(t/60,Rb(1,:)*180/pi,'--m','MarkerSize',3)
plot(t/60,ytil(1,:)*180/pi,'g',t/60,Rb_hat(1,:)*180/pi,'--k')%,t,omega_est(1,:),t,omega_tilda(:,1),t,x1(:,5))
% plot(t,omega_est(1,:),t,omega_tilda(:,1),t,x1(:,5))%attitude z angle in body frame
title('Roll Angle Time Evolution/Estimates')
xlabel('Time (min)')
ylabel('Roll Angle (deg)')
grid on
% title('Components of q and its Magnitude')
legend('Noisy Roll Data', 'Estimated Roll', 'True Roll')


        k1 = ytil(1,:)-Rb_hat(1,:);
        error1 = Rb(1,:)-Rb_hat(1,:);
            subplot(2,1,2)
            hold all;
            plot(t/60,3*180/pi*vestt(1,:),':m','LineWidth',2);
            plot(t/60,k1*180/pi,'k');
        %     plot(t/60,error1*180/pi);
            plot(t/60,-3*180/pi*vestt(1,:),':c','LineWidth',2);
        %     xlim([0,t(w)])
         title('Roll Error 3\sigma Bounds')
         legend('+3\sigma bound','Attitude Error','-3\sigma bound')
            xlabel('Time (min')
            ylabel('3 \sigma bound Error (deg)')
         grid on
 
 figure
subplot(2,1,1)
hold on
% plot(t,x1(:,3),'k*',t,qest(3,:)) %t,x1(:,2),'r',t, x1(:,3),'g',t,x1(:,4),t,Mag)
% plot(t,x1(:,7),'b*',t,omega_est(3,:)) %angular velo
plot(t/60,Rb(2,:)*180/pi,'--m','MarkerSize',3)
plot(t/60,ytil(2,:)*180/pi,'g',t/60,Rb_hat(2,:)*180/pi,'--k')%,t,omega_est(1,:),t,omega_tilda(:,1),t,x1(:,5))
% plot(t,omega_est(1,:),t,omega_tilda(:,1),t,x1(:,5))%attitude z angle in body frame
title('Pitch Angle Time Evolution/Estimates')
xlabel('Time (min)')
ylabel('Pitch Angle (deg)')
grid on
% title('Components of q and its Magnitude')
legend('Noisy Pitch Data', 'Estimated Pitch Angle', 'True Pitch Angle')


        k2 = ytil(2,:)-Rb_hat(2,:);
        error1 = Rb(2,:)-Rb_hat(2,:);
            subplot(2,1,2)
            hold all;
            plot(t/60,3*180/pi*vestt(2,:),':m','LineWidth',2);
            plot(t/60,k2*180/pi,'k');
        %     plot(t/60,error1*180/pi);
            plot(t/60,-3*180/pi*vestt(2,:),':c','LineWidth',2);
        %     xlim([0,t(w)])
         title('Pitch Error 3\sigma Bounds')
         legend('+3\sigma bound','Pitch Error','-3\sigma bound')
            xlabel('Time (min')
            ylabel('3 \sigma bound Error (deg)')
         grid on

 
 figure
subplot(2,1,1)
% plot(t,x1(:,3),'k*',t,qest(3,:)) %t,x1(:,2),'r',t, x1(:,3),'g',t,x1(:,4),t,Mag)
% plot(t,x1(:,7),'b*',t,omega_est(3,:)) %angular velo
plot(t/60,Rb(3,:)*180/pi,'--m',t/60,ytil(3,:)*180/pi,'g',t/60,Rb_hat(3,:)*180/pi,'--k')%,t,omega_est(1,:),t,omega_tilda(:,1),t,x1(:,5))
% plot(t,omega_est(1,:),t,omega_tilda(:,1),t,x1(:,5))%attitude z angle in body frame
title('Yaw Angle Time Evolution/Estimates')
xlabel('Time (min)')
ylabel('Yaw Angle (deg)')
grid on
legend('Noisy Yaw Data', 'Estimated Yaw Angle', 'True Yaw Angle')

        subplot(2,1,2)
        hold all;
        plot(t/60,180/pi*(ytil(3,:)-Rb_hat(3,:)),'k')
        plot(t/60,3*180/pi*(vestt(3,:)),':m','LineWidth',2);
        plot(t/60,-3*180/pi*(vestt(3,:)),':c','LineWidth',2);
        title('Yaw Error 3\sigma Bounds')
    %     xlim([0,t(w)])
        xlabel('Time (min')
        ylabel('3 \sigma bound Error')
        grid on
 
