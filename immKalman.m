% 1、imm（交互多模型卡尔曼滤波器）是采用的对个卡尔曼滤波器模型同时做跟踪，
%         然后利用模型概率综合出，对每个模型的可行度，对每个模型输出的预测值和协方差，进行加权求平均
% 2、imm重点的核心在于，多个滤波器的概率的更新的方法，imm采用的是最大似然估计。
% 3、imm滤波器应该考虑的因素：
%     a、选择一定个数的imm滤波器，包括较为精确的模型和较为粗糙的模型，imm滤波算法不仅描述了目标的
%             连续运动状态还描述了目标的机动性
%     b、马尔科夫链状态转移概率的选择，对imm滤波器的性能较大影响。
%     c、imm算法具有模块化的特性，当运动模型较为精确的时候，可以采用比较精确的运动模型。
%            当无法预料目标的运动规律的时候，那就应该选择更一般的模型，使得该模型具有更强的鲁棒性。
function immKalman
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
plot(Rx,Ry,'r.',Zx,Zy,'g+',xn,yn,'-bo');
legend('真实值 ','观测值','估计值');
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%x_est,p_est,记返回第m次仿真，第r个采样点的滤波结果
%xn_est,pn_est,记录每个模型对应的m次仿真，第r个采样点的滤波结果
%T雷达周期
%Z测量值
%Delta过程噪声因子
%u初始模型选择概率矩阵
function [x_est,p_est,xn_est,pn_est,u]=imm(xn_est,pn_est,T,Z,Delta,u)
p=[0.95 0.025 0.025;%控制模型转换的马尔科夫链的转移概率矩阵
        0.025 0.95 0.025;
        0.025 0.025 0.95];
%采用三个不同的模型参数，模型1个非机动模型，模型2、3为机动模型（Q值不同）    
phi{1,1}=[1 T 0 0 ;
                0 1 0 0;
                0 0 1 T;
                0 0 0 1];
phi{1,1}(6,6)=0;%模型1,状态转移方程。

phi{2,1}=[1 T 0 0 T^2/2 0;
                0 1 0 0 T 0;
                0 0 1 T 0 T^2/2;
                0 0 0 1 0 T;
                0 0 0 0 1 0 ;
                0 0 0 0 0 1];%模型2,状态转移方程。
            
phi{3,1}=phi{2,1};%模型3,状态转移方程。

g{1,1}=[T/2,0;
                1,0;
                0 T/2;
                0 1];
g{1,1}(6,2)=0;%模型1

g{2,1}=[T^2/4,0;
                T/2,0;
                0,T^2/4;
                0 T/2;
                1 0 ;
                0 1];%模型2
            
g{3,1}=g{2,1};%模型3

q{1,1}=zeros(2);
q{2,1}=0.001*eye(2);
q{3,1}=0.014*eye(2);
% Q = g{j,1}*q{j,1}*g{j,1}';
H = [1 0 0 0 0 0;
        0 0 1 0 0 0 ];
R = eye(2)*Delta^2;%测量噪声协方差矩阵
mu=zeros(3,3);%混合概率矩阵
c_mean = zeros(1,3);%归一化常数
for i = 1:3
    c_mean=c_mean+p(i,:)*u(i);
end
for i = 1:3
    mu(i,:)=p(i,:)*u(i)./c_mean;
end
%交互输入
for j = 1:3
    x0{j,1}=zeros(6,1);
    p0{j,1}=zeros(6);
    for i = 1:3
        x0{j,1}=x0{j,1}+xn_est{i,1}*mu(i,j);
    end
    for i =1:3
        p0{j,1}=p0{j,1}+mu(i,j)*(pn_est{i,1}...
                    +(xn_est{i,1}-x0{j,1})*(xn_est{i,1}-x0{j,1})');
    end
end
%模型条件滤波
a=zeros(1,3);
for j = 1:3
    x_pre{j,1}=phi{j,1}*x0{j,1};
    p_pre{j,1}=phi{j,1}*p0{j,1}*phi{j,1}'+g{j,1}*q{j,1}*g{j,1}';%  Q = g{j,1}*q{j,1}*g{j,1}';
    k{j,1}=p_pre{j,1}*H'*inv(H*p_pre{j,1}*H'+R);
    xn_est{j,1}=x_pre{j,1}+k{j,1}*(Z-H*x_pre{j,1});
    pn_est{j,1}=(eye(6)-k{j,1}*H)*p_pre{j,1};
end
%模型概率更新
for j = 1:3
    v{j,1}=Z-H*x_pre{j,1};%新息
    s{j,1}=H*p_pre{j,1}*H'+R;%观测协方差矩阵
    n=length(s{j,1})/2;
    a(1,j)=1/((2*pi)^n*sqrt(det(s{j,1})))*exp(-0.5*v{j,1}'...
                *inv(s{j,1})*v{j,1});%观测相对于模型j的似然函数
end
c=sum(a.*c_mean);%归一化常数
u=a.*c_mean./c;%概率模型更新

%交互输出
xn=zeros(6,1);
pn=zeros(6);
for j = 1:3
    xn=xn+xn_est{j,1}.*u(j);
end
for j = 1:3
    pn=pn+u(j).*(pn_est{j,1}+(xn_est{j,1}-xn)*(xn_est{j,1}-xn)');
end
%返回滤波结果
x_est=xn;
p_est=pn;


