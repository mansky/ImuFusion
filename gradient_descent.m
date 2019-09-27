% 这个段程序是用来融合加速度计和陀螺仪的数据计算姿态角的
% 使用梯度下降对陀螺仪数据进行校准，然后在由角速度计算姿态角
% 首先要创建sensor矩阵，它是一个N行6列的矩阵，分别是gyrox,gyroy,gyroz,accx,accy,accz
% 陀螺仪的单位是rad/s，加速度的单位是9.8m/s?，即当加速度计水平时，Z轴加速度为9.8m/s?，所以accz=1。
% 计算出的陀螺仪的校正值存放在result中的前三列
% 如果手上没有数据，可以自己生成数据，然后加噪声。
% 思路是瞎编陀螺仪数据，使用四元数微分方程生成四元数，在从四元数得出旋转矩阵，通过旋转矩阵直接求加速度就可以了
IMU_Dt=0.01;  
gradientAngle=zeros(size(sensor,1),3);
gradientQuat=zeros(size(sensor,1),4);
gradientQuat(:,1)=1;
for i = 2:size(sensor,1)
    q0=gradientQuat(i-1,1);
    q1=gradientQuat(i-1,2);
    q2=gradientQuat(i-1,3);
    q3=gradientQuat(i-1,4);

    gx=sensor(i,1);
    gy=sensor(i,2);
    gz=sensor(i,3);
    ax=sensor(i,4);
    ay=sensor(i,5);
    az=sensor(i,6);


    %角速度模长
    Gyro_Length=sqrt(gx*gx+gy*gy+gz*gz)*57.3;%单位deg/s

    %四元数微分方程计算本次待矫正四元数 
    qDot1 = 0.5 * (-q1*gx - q2*gy - q3*gz);
    qDot2 = 0.5 * (q0*gx + q2*gz - q3*gy);
    qDot3 = 0.5 * (q0*gy - q1*gz + q3*gx);
    qDot4 = 0.5 * (q0*gz + q1*gy - q2*gx);

    %加速度计输出有效时,利用加速度计补偿陀螺仪 
    if(ax*ay*az==0)
     return
    end 

    Anorm=sqrt(ax * ax + ay * ay + az * az);
    ax = ax/Anorm;
    ay = ay/Anorm;
    az = az/Anorm;
    % 避免重复运算 
    m2q0 = 2.0 * q0;
    m2q1 = 2.0 * q1;
    m2q2 = 2.0 * q2;
    m2q3 = 2.0 * q3;
    m4q0 = 4.0 * q0;
    m4q1 = 4.0 * q1;
    m4q2 = 4.0 * q2;
    m8q1 = 8.0 * q1;
    m8q2 = 8.0 * q2;
    q0q0 = q0*q0;
    q1q1 = q1*q1 ;
    q2q2 = q2*q2;
    q3q3 = q3*q3;

    %梯度下降算法,计算误差函数的梯度 
    s0 = m4q0 * q2q2 + m2q2 * ax + m4q0 * q1q1 - m2q1 * ay;
    s1 = m4q1 * q3q3 - m2q3 * ax + 4.0 * q0q0 *q1 - m2q0 * ay - m4q1 + m8q1 * q1q1 + m8q1 * q2q2 + m4q1 * az;
    s2 = 4.0 * q0q0 * q2 + m2q0 * ax + m4q2 * q3q3 - m2q3 * ay - m4q2 + m8q2 * q1q1 + m8q2 * q2q2 + m4q2 * az;
    s3 = 4.0 * q1q1 * q3 - m2q1 * ax + 4.0 * q2q2 *q3 - m2q2 * ay;

    % 梯度归一化
    Snorm=sqrt(s0 * s0 + s1 * s1 + s2 * s2 + s3 * s3);
    s0 = s0 /Snorm;
    s1 = s1 /Snorm;
    s2 = s2 /Snorm;
    s3 = s3 / Snorm;

    BETADEF=IMU_Dt+0.01*Gyro_Length*IMU_Dt;
    qDot1 = qDot1 - BETADEF * s0;
    qDot2 = qDot2 - BETADEF * s1;
    qDot3 = qDot3 - BETADEF * s2;
    qDot4 = qDot4 - BETADEF * s3;

    % 补偿由四元数微分方程引入的姿态误差 */
    %将四元数姿态导数积分,得到当前四元数姿态 */
    %二阶毕卡求解微分方程 */
    delta = (IMU_Dt * gx) * (IMU_Dt * gx) + (IMU_Dt * gy) * (IMU_Dt * gy) + (IMU_Dt * gz) * (IMU_Dt * gz);
    q_out1 = (1.0 - delta / 8.0) * q0 + qDot1 * IMU_Dt;
    q_out2 = (1.0 - delta / 8.0) * q1 + qDot2 * IMU_Dt;
    q_out3 = (1.0 - delta / 8.0) * q2 + qDot3 * IMU_Dt;
    q_out4 = (1.0 - delta / 8.0) * q3 + qDot4 * IMU_Dt;
    %单位化四元数 */
    recipNorm=1/sqrt(q_out1 * q_out1 + q_out2 * q_out2 + q_out3 * q_out3 + q_out4 * q_out4);
    q_out1 = q_out1*recipNorm;
    q_out2 = q_out2*recipNorm;
    q_out3 = q_out3*recipNorm;
    q_out4 = q_out4*recipNorm;
    
    gradientQuat(i,1)=q_out1;
    gradientQuat(i,2)=q_out2;
    gradientQuat(i,3)=q_out3;
    gradientQuat(i,4)=q_out4;
end
%从四元数中求出角度
for i = 1:size(sensor,1)
    q0=gradientQuat(i,1);
    q1=gradientQuat(i,2);
    q2=gradientQuat(i,3);
    q3=gradientQuat(i,4);
    %roll
    gradientAngle(i,2)= 57.3*atan2(2 * (q2*q3 + q0*q1), q0*q0 - q1*q1 - q2*q2 + q3*q3);
    %pitch
    gradientAngle(i,1)=57.3* asin(-2 * (q1*q3 - q0*q2));
    %yaw
    gradientAngle(i,3) = 57.3*atan2(2 * (q1*q2 + q0*q3), q0*q0 + q1*q1 - q2*q2 - q3*q3);
 end
%删除中间变量，保持工作区的清洁
clear Anorm ax ay az BETADEF delta gx gy gz Gyro gz i IMU_Dt
clear Gyro_Length m2q0 m2q1  m2q2  m2q3  m4q0 m4q1 m4q2 m4q3 m8q0 m8q1  m8q2
clear q0 q0q0 q1 q1q1 q2 q2q2 q3 q3q3 q_out1 q_out2 q_out3 q_out4
clear qDot1 qDot2 qDot3 qDot4 recipNorm result s0 s1 s2 s3 Snorm

