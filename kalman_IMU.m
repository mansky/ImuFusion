% 这个段程序是用来融合加速度计和陀螺仪的数据计算姿态角的
% 使用卡尔曼滤波对陀螺仪数据进行校准，然后在由角速度计算姿态角
% 首先要创建sensor矩阵，它是一个N行6列的矩阵，分别是gyrox,gyroy,gyroz,accx,accy,accz
% 陀螺仪的单位是rad/s，加速度的单位是9.8m/s?，即当加速度计水平时，Z轴加速度为9.8m/s?，所以accz=1。
% 计算出的陀螺仪的校正值存放在result中的前三列
% 如果手上没有数据，可以自己生成数据，然后加噪声。
% 思路是瞎编陀螺仪数据，使用四元数微分方程生成四元数，在从四元数得出旋转矩阵，通过旋转矩阵直接求加速度就可以了

%生成单位矩阵
E=eye(6);
%协方差矩阵
P=[10    0   0   0   0   0
    0    10   0   0   0   0
    0    0   10   0   0   0
    0    0   0   10   0   0
    0    0   0   0   10   0
    0    0   0   0   0   10];
dt=0.01;%100Hz,dt为两组数据的时间间隔
% X为最优估计解，即输出的数据,从上到下依次是wx,wy,wz,ax,ay,az。
% w代表角速度，单位rad/s，a为单位加速度，单位为 9.8m/s?
X=[0;0;0;0;0;1];


%系统转移过程噪声
Q=[0.01    0    0       0       0       0
   0     0.01   0       0       0       0
   0       0    0.01    0       0       0
   0       0    0       0.01    0       0
   0       0    0       0       0.01    0
   0       0    0       0       0       0.01 ];

%测量噪声
R=[0.1    0    0       0        0       0
   0     0.1   0       0        0       0
   0      0    0.1     0        0       0
   0      0    0       10       0       0
   0      0    0       0        10      0
   0      0    0       0        0       10 ];

%测量矩阵
H=[1      0     0       0        0       0
   0      1     0       0        0       0
   0      0     1       0        0       0
   0      0     0       1        0       0
   0      0     0       0        1       0
   0      0     0       0        0       1 ];

%测量值
Z=[0;0;0;0;0;1];

%需要把数据导入sensor,1-6列分别为gx,gy,gz,ax,ay,az
n=size(sensor,1);%求sensor数组的行数
result=zeros(n,6);
for i =1:n
    % A为状态转移矩阵，实际的系统为非线性
    A=[1            0           0               0               0               0
       0            1           0               0               0               0
       0            0           1               0               0               0
       0           -X(6,1)*dt   X(5,1)*dt       1            X(3,1)*dt       -X(2,1)*dt
      X(6,1)*dt    0            -X(4,1)*dt   -X(3,1)*dt           1            X(1,1)*dt 
      -X(5,1)*dt   X(4,1)*dt    0            X(2,1)*dt        -X(1,1)*dt        1        ];
    %把当前时刻的观察值导入到观察矩阵Z中
    Z=sensor(i,1:6)';
    %计算当前时刻X的估计值
    X_predict=A*X;
    %计算当前时刻协方差矩阵P的估计值
    P_predict=A'*P*A+Q;
    %计算卡尔曼增益，除以某个矩阵等于乘以矩阵的逆矩阵
    K=H'*P_predict/(H*P_predict*H'+R);
    %计算最优估计解
    X=X_predict+K*(Z-H*X_predict); 
    %计算最优估计解对应的协方差矩阵P
    P=(E-K*H)*P_predict;  
    %对加速度进行单位化
    accNorm=X(4,1)^2+X(5,1)^2+X(6,1)^2;
    accNorm=accNorm^0.5;
    X(4,1)=X(4,1)/accNorm;
    X(5,1)=X(5,1)/accNorm;
    X(6,1)=X(6,1)/accNorm;
    %复制结果到result矩阵
    result(i,:)=X';
end

kalmanQuat=zeros(n,4);
kalmanQuat(:,1)=1;
kalmanAngles=zeros(n,3);
%通过四元数微分方程解算姿态
for i =1:n-1
          qDot1 =  kalmanQuat(i,1)+0.5 * dt*(-kalmanQuat(i,2)*result(i,1) - kalmanQuat(i,3)*result(i,2) - kalmanQuat(i,4)*result(i,3));
          qDot2 =  kalmanQuat(i,2)+0.5 * dt*(kalmanQuat(i,1)*result(i,1) + kalmanQuat(i,3)*result(i,3) - kalmanQuat(i,4)*result(i,2));
          qDot3 =  kalmanQuat(i,3)+0.5 * dt*(kalmanQuat(i,1)*result(i,2) - kalmanQuat(i,2)*result(i,3) + kalmanQuat(i,4)*result(i,1));
          qDot4 =  kalmanQuat(i,4)+0.5 * dt*(kalmanQuat(i,1)*result(i,3) + kalmanQuat(i,2)*result(i,2) - kalmanQuat(i,3)*result(i,1));
          
          quatNorm=qDot1^2+qDot2^2+qDot3^2+qDot4^2;
          qDot1=qDot1/quatNorm;
          qDot2=qDot2/quatNorm;
          qDot3=qDot3/quatNorm;
          qDot4=qDot4/quatNorm;
          kalmanQuat(i+1,1)=qDot1;
          kalmanQuat(i+1,2)=qDot2;
          kalmanQuat(i+1,3)=qDot3;
          kalmanQuat(i+1,4)=qDot4;
          
          sinr_cosp = (+2.0 * (kalmanQuat(i+1,1) * kalmanQuat(i+1,2) + kalmanQuat(i+1,3) * kalmanQuat(i+1,4)));
          cosr_cosp = (+1.0 - 2.0 * (kalmanQuat(i+1,2) * kalmanQuat(i+1,2) + kalmanQuat(i+1,3) * kalmanQuat(i+1,3)));
          kalmanAngles(i,2) = atan2(sinr_cosp, cosr_cosp);

          sinp = +2.0 * (kalmanQuat(i+1,1) * kalmanQuat(i+1,3) - kalmanQuat(i+1,4)* kalmanQuat(i+1,2));
        if (abs(sinp) >= 1)
            kalmanAngles(i,1) = copysign(3.1415926 / 2, sinp); % 俯仰为90°
        else
            kalmanAngles(i,1) = asin(sinp);
        end

        % 航向轴
        siny_cosp = +2.0 * (kalmanQuat(i+1,1) * kalmanQuat(i+1,4) + kalmanQuat(i+1,2) * kalmanQuat(i+1,3));
        cosy_cosp = +1.0 - 2.0 * (kalmanQuat(i+1,3) * kalmanQuat(i+1,3) + kalmanQuat(i+1,4) * kalmanQuat(i+1,4)); 
        kalmanAngles(i,1) = kalmanAngles(i,1)*57.3;
        kalmanAngles(i,2) = kalmanAngles(i,2)*57.3;
        kalmanAngles(i,3) = atan2(siny_cosp, cosy_cosp)*57.3;
end
%删除中间变量，保持工作区的清洁
clear A accNorm c cosr_cosp cosy_cosp dt E H i K n P P_predict
clear Q qDot1 qDot2 qDot3  qDot4 quatNorm R sinp sinr_cosp siny_cosp X X_predict Z



