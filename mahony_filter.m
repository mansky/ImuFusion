% 这个段程序是用来融合加速度计和陀螺仪的数据计算姿态角的
% 使用互补滤波对陀螺仪数据进行校准，然后在由角速度计算姿态角
% 首先要创建sensor矩阵，它是一个N行6列的矩阵，分别是gyrox,gyroy,gyroz,accx,accy,accz
% 陀螺仪的单位是rad/s，加速度的单位是9.8m/s?，即当加速度计水平时，Z轴加速度为9.8m/s?，所以accz=1。
% 计算出的陀螺仪的校正值存放在result中的前三列
% 如果手上没有数据，可以自己生成数据，然后加噪声。
% 思路是瞎编陀螺仪数据，使用四元数微分方程生成四元数，在从四元数得出旋转矩阵，通过旋转矩阵直接求加速度就可以了
dt=0.01;  
Kp=1;
Ki=0;
halfT=0.5*dt;
mahonyAngle=zeros(size(sensor,1),3);
mahonyQuat=zeros(size(sensor,1),4);
mahonyQuat(:,1)=1;
exInt=0;
eyInt=0;
ezInt=0;
for i = 2:size(sensor,1)
    q0=mahonyQuat(i-1,1);
    q1=mahonyQuat(i-1,2);
    q2=mahonyQuat(i-1,3);
    q3=mahonyQuat(i-1,4);

    gx=sensor(i,1);
    gy=sensor(i,2);
    gz=sensor(i,3);
    ax=sensor(i,4);
    ay=sensor(i,5);
    az=sensor(i,6);
    norm = sqrt(ax*ax + ay*ay + az*az);      
    ax = ax /norm;
    ay = ay / norm;
    az = az / norm;	

    vx = 2*(q1*q3 - q0*q2);											
    vy = 2*(q0*q1 + q2*q3);
    vz = q0*q0 - q1*q1 - q2*q2 + q3*q3 ;

    ex = (ay*vz - az*vy) ;                           				
    ey = (az*vx - ax*vz) ;
    ez = (ax*vy - ay*vx) ;
    
    exInt = exInt + ex * Ki;							
    eyInt = eyInt + ey * Ki;
    ezInt = ezInt + ez * Ki;

    gx = gx + Kp*ex + exInt;					   						
    gy = gy + Kp*ey + eyInt;
    gz = gz + Kp*ez + ezInt;				   						

    q_out1 = q0 + (-q1*gx - q2*gy - q3*gz)*halfT;
    q_out2 = q1 + (q0*gx + q2*gz - q3*gy)*halfT;
    q_out3 = q2 + (q0*gy - q1*gz + q3*gx)*halfT;
    q_out4 = q3 + (q0*gz + q1*gy - q2*gx)*halfT;

    norm = sqrt(q_out1*q_out1 + q_out2*q_out2 + q_out3*q_out3 + q_out4*q_out4);

    q_out1 = q_out1 / norm;
    q_out2 = q_out2 / norm;
    q_out3 = q_out3 / norm;
    q_out4 = q_out4 / norm;
    mahonyQuat(i,1)=q_out1;
    mahonyQuat(i,2)=q_out2;
    mahonyQuat(i,3)=q_out3;
    mahonyQuat(i,4)=q_out4;
end
%从四元数中求出角度
for i = 1:size(sensor,1)
    q0=mahonyQuat(i,1);
    q1=mahonyQuat(i,2);
    q2=mahonyQuat(i,3);
    q3=mahonyQuat(i,4);
    %roll
    mahonyAngle(i,2)= 57.3*atan2(2 * (q2*q3 + q0*q1), q0*q0 - q1*q1 - q2*q2 + q3*q3);
    %pitch
    mahonyAngle(i,1)=57.3* asin(-2 * (q1*q3 - q0*q2));
    %yaw
    mahonyAngle(i,3) = 57.3*atan2(2 * (q1*q2 + q0*q3), q0*q0 + q1*q1 - q2*q2 - q3*q3);
end
%删除中间变量，保持工作区的清洁
clear  ax ay az dt ex exInt ey eyInt ez ezInt gx gy gz halfT i Ki Kp norm
clear q0 q1 q2 q3  q_out1 q_out2 q_out3 q_out4 qDot1 qDot2 qDot3 qDot4
clear vx vy vz

