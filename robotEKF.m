% Code nay co chinh sua va sao chep tu code mau cua co :D 
pos_true_X = [0]; 
pos_true_Y = [0];
pos_true_T = [0];
pos_mes_X = [0];
pos_mes_Y = [0];
pos_mes_T = [0];
posX = [0];
posY = [0];
posT = [0];
pri_est_X = [0];
pri_est_Y = [0];
pri_est_T = [0];
pos_est_X = [0];
pos_est_Y = [0];
pos_est_T = [0];

N = 40;
vL2 = [ones(20,1)*pi/4;ones(20,1)*pi];
vR2 = [ones(20,1)*pi; ones(20,1)*pi];
vL = zeros(N);
vR = zeros(N);
PNoise_L = rand(N,1);
PNoise_R = rand(N,1);
ONoise_x = rand(N,1);
ONoise_y= rand(N,1);
ONoise_t= rand(N,1);
dt = 0.05;
r_wheel = 0.5;
b_wheel = 0.6;
P = zeros(3,3);
i = 2;
while(i<=N)
     %1 Velocites with noise
    vL=vL2(i-1)+PNoise_L(i-1);
    vR=vR2(i-1)+PNoise_R(i-1);
    dSL=dt*r_wheel*vL;
    dSR=dt*r_wheel*vR;
    dS=(dSL+dSR)/2;
    dTheta=(dSR-dSL)/b_wheel;
     %Real path (Kinematic model + velocities with noise)
    pos_true_X(i)=pos_true_X(i-1)+ dS*cos(pos_true_T(i-1)+ dTheta/2);
    pos_true_Y(i)=pos_true_Y(i-1)+ dS*sin(pos_true_T(i-1)+dTheta/2);
    pos_true_T(i)=pos_true_T(i-1)+ dTheta;
     %2. Measuremement
    %z(k)=h(x(k),u(k)) +w(k)), measurement of odometry= real path +measurement noise
    pos_mes_X(i)=pos_true_X(i)+ONoise_x(i);
    pos_mes_Y(i)=pos_true_Y(i)+ONoise_y(i);
    pos_mes_T(i)=pos_true_T(i)+ONoise_t(i);
    pos_mes=[pos_mes_X(i) pos_mes_Y(i) pos_mes_T(i)];%1x3
%3. Estimation x^-(k)=f(x^(k-1),u(k),0)
    % Velocites without noise
    vL=vL2(i-1);
    vR=vR2(i-1);
    dSL=dt*r_wheel*vL;
    dSR=dt*r_wheel*vR;
    dS=(dSL+dSR)/2;
    dTheta=(dSR-dSL)/b_wheel;
    
    % Theory path = kinematic + velocities without noise
    posX(i)=posX(i-1)+ dS*cos(posT(i-1)+ dTheta/2);
    posY(i)=posY(i-1)+ dS*sin(posT(i-1)+dTheta/2);
    posT(i)=posT(i-1)+ dTheta;
    
    pri_est_X(i)=pos_est_X(i-1)+dS*cos(pos_est_T(i-1)+ dTheta/2);
    pri_est_Y(i)=pos_est_Y(i-1)+dS*sin(pos_est_T(i-1)+dTheta/2);
    pri_est_T(i)=pos_est_T(i-1)+dTheta;
    pri_est=[ pri_est_X(i) pri_est_Y(i) pri_est_T(i)]; %1x3
    
%4. P va Kalman gain K
    A=[1 0 -dS*sin(pos_est_T(i-1)+dTheta/2);0 1 dS*cos(pos_est_T(i-1)+dTheta/2); 0 0 1]; % Jacobian of input system
    W=dt*r_wheel/2*[cos(pos_est_T(i-1)+dTheta/2)-(dS/b_wheel)*sin(pos_est_T(i-1)+dTheta/2) ... % Jacobian of input noise
                    cos(pos_est_T(i-1)+dTheta/2)+(dS/b_wheel)*sin(pos_est_T(i-1)+dTheta/2);...
                    sin(pos_est_T(i-1)+dTheta/2)+(dS/b_wheel)*cos(pos_est_T(i-1)+dTheta/2) ...
                    sin(pos_est_T(i-1)+dTheta/2)-(dS/b_wheel)*cos(pos_est_T(i-1)+dTheta/2);...
                    1/b_wheel -1/b_wheel];
                
    sig = sqrt((0.5^2)/12);
    Q = [sig*(vL2(i)^2) 0; 0 sig*(vR2(i)^2)];
    R = [10*sig 0 0; 0 10*sig 0; 0 0 sig];
    P=A*P*A'+W*Q*W'; %3x3
    K=P/(P+R); %3x3
    
    % Correct priori estimate with measurement z
    residual=pos_mes-pri_est;
    pos_est= pri_est' + K*(residual)';   % 1x3+ 3*3(3x1-3x1)=3x1
    
    pos_est_X(i)= pos_est(1); % posteriori estimate
    pos_est_Y(i)= pos_est(2);
    pos_est_T(i)= pos_est(3);
    
    P=P-K*P; % update 
        
    i=i+1;
 end  
plot(posX,posY);
hold on;
plot(pos_est_X,pos_est_Y);
hold on;
plot(pri_est_X,pri_est_Y,'g');
