clear all;clc;close all;
T =2;%扫描周期
M=5;%蒙特卡洛模拟次数，用于测试，计算滤波误差均值和误差标准差
N = 900/T;%总得采用点数
N1 = 400/T;%第一次转弯处采样起点
N2 = 600/T;%第一次匀速处采用起点
N3 = 610/T;%第二转弯处采样起点
N4 = 660/T;%第二匀速处采用起点
Delta = 100;%测量噪声标准差
Rx = zeros(N,1);
Ry = zeros(N,1);
Zx = zeros(N,M);
Zy = zeros(N,M);
%沿y轴匀速直线运动
t = 2:T:400;
x0 = 2000 + 0*t';
y0 = 10000 - 15*t';
%慢转弯
t= 400+T:T:600;
x1 = x0(N1) + 0.075*((t'-400).^2)/2;
y1 = y0(N1) - 15*(t'-400)+0.075*((t'-400).^2)/2;
%匀速
t = 600+T:T:610;
vx = 0.075*(600-400);
x2 = x1(N2-N1)+vx*(t'-600);
y2 = y1(N2-N1)+0*t';
%快转弯
t = 610+T : T:660;
x3=x2(N3-N2)+(vx*(t'-610)-0.3*((t'-610).^2)/2);
y3=y2(N3-N2)-0.3*((t'-610).^2)/2;
%匀速
t=660+T:T:900;
vy=-0.3*(660-610);
x4=x3(N4-N3)+0*t';
y4=y3(N4-N3)+vy*(t'-660);
%轨迹合成
Rx = [x0;x1;x2;x3;x4];
Ry = [y0;y1;y2;y3;y4];
%每次蒙塔卡洛仿真的滤波估计位置的初始化
mtEstPx=zeros(M,N);
mtEstPy=zeros(M,N);
%产生观测数据，要仿真M次，必须有M次的观测数据
nx = randn(N,M)*Delta;%产生测量噪声
ny = randn(N,M)*Delta;
Zx = Rx*ones(1,M)+nx;%真实的轨迹上叠加噪声
Zy = Ry*ones(1,M)+ny;

for m = 1:M
    %滤波初始化
    mtEstPx(m,1)=Zx(1,m);%初始数据
    mtEstPy(m,1)=Zx(2,m);
    xn(1)=Zx(1,m);%滤波初值
    xn(2)=Zx(2,m);
    yn(1)=Zy(1,m);
    yn(2)=Zy(2,m);
    %非机动模型参数
    phi = [1 T 0 0 ;%运动模型
                0 1 0 0;
                0 0 1 T;
                0 0 0 1];
   h= [1 0 0 0;%测量方程
            0 0 1 0];
   g=[T/2 0;
            1 0 ;
            0 T;
            0 1];
%    q=0.01*[];
   R=[Delta.^2 0 ;%
            0 Delta.^2];
    vx=(Zx(2)-Zx(1,m))/2;
    vy=(Zy(2)-Zy(1,m))/2;
    %初始状态估计
    x_est=[Zx(2,m);vx;Zy(2,m);vy];
    p_est=[Delta^2,Delta^2/T,0,0;
                    Delta^2/T,2*Delta^2/(T^2),0,0;
                    0,0,Delta^2,Delta^2/T;
                    0,0,Delta^2/T,2*Delta^2/(T^2)];
     mtEstPx(m,2)=x_est(1);
     mtEstPy(m,2)=x_est(3);
     %滤波开始
     for r=3:N
        z=[Zx(r,m);Zy(r,m)];
        if r<20%前20次只做非机动模型滤波
            x_pre=phi*x_est;%预测,Warning :该地方未添加过程噪声协方差Q
            p_pre=phi*p_est*phi';%预测误差协方差,Warning :该地方未添加过程噪声协方差Q
            k=p_pre*h'*inv(h*p_pre*h'+R);
            x_est=x_pre+k*(z-h*x_pre);%滤波
            p_est=(eye(4)-k*h)*p_pre;%滤波协方差
            
            xn(r)=x_est(1);%记录采样点滤波数据
            yn(r)=x_est(3);
            mtEstPx(m,r)=x_est(1);%记录第m次仿真滤波估计数据
            mtEstPy(m,r)=x_est(3);
        else
            if r==20
                x_est=[x_est;0;0];%扩维
                p_est=p_est;
                p_est(6,6)=0;%扩维
                for i =1:3
                    xn_est{i,1}=x_est;
                    pn_est{i,1}=p_est;
                end
                u=[0.8,0.1,0.1];%模型概率初始化
            end
            [x_est,p_est,xn_est,pn_est,u]=imm(xn_est,pn_est,T,z,Delta,u);
            xn(r)=x_est(1);
            yn(r)=x_est(3);
            mtEstPx(m,r)=x_est(1);
            mtEstPy(m,r)=x_est(3);
        end
     end%结束一次滤波
% figure(m);
% plot(Rx,Ry,'r',Zx,Zy,'g*',xn,yn,'b*');
% legend('真实轨迹','观察样本','估计轨迹');
end
err_x=zeros(N,1);
err_y=zeros(N,1);
delta_x = zeros(N,1);
delta_y=zeros(N,1);
for r = 1:N
    ex=sum(Rx(r)-mtEstPx(:,r));
    ey=sum(Ry(r)-mtEstPy(:,r));
    err_x(r)=ex/M;
    err_y(r)=ey/M;
    eqx=sum((Rx(r)-mtEstPx(:,r)).^2);
    eqy=sum((Ry(r)-mtEstPy(:,r)).^2);
    
    delta_x(r)=sqrt(abs(eqx/M-(err_x(r)^2)));
    delta_y(r)=sqrt(abs(eqy/M-(err_y(r)^2)));
end
figure(1);
plot(Rx,Ry,'r',Zx,Zy,'g.',xn,yn,'bo');
legend('真实轨迹','观察样本','估计轨迹');
figure(2);
subplot(211);
plot(err_x);
title('average of Error of x dirction ');
subplot(212);
plot(err_y);
title('average of Error of y dirction ');
figure(3);
subplot(211);
plot(delta_x);
title('standard deviation  of Error of x dirction ');
subplot(212);
plot(delta_y);
title('standard deviation  of Error of y dirction ');
