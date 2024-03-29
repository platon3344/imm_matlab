function immKalman_test 
clc;clear;close all;
T = 2;%
M = 5;
t = 900;
N = t/T;
delta = 100;
t0 = 0:T:400 - T; n0 = 400/T-1;
t1 = 0 : T : 200- T;n1 = 200/T-1;%400 ~ 600  第一次拐弯起点
t2 = 0  :  T : 10- T;n2 = 10/T-1;%600 ~ 610 第二匀速起点
t3 = 0 : T : 50- T;n3 = 50/T-1;%610 ~ 660 
t4 = 0 : T : 240- T ;n4 = 240/T-1;
Xstart = 2000;Ystart = 10000;
v0x = 0; v0y = -15;
v1x = 0;a1x = 0.075;v1y = -15;a1y = 0.075;
v2x = 200*0.075;v2y = 0;
v3x = v2x ; a3x = -0.3;v3y = 0;a3y = a3x;
v4x = 0;v4y = -0.3*50;
%真是观测值生成
x0 = Xstart + v0x*t0;
y0 = Ystart + v0y*t0;
x1 = x0(n0) + v1x*t1 + 1/2*a1x*t1.^2;
y1 = y0(n0) + v1y*t1 + 1/2*a1y*t1.^2;
x2 = x1(n1) + v2x*t2;
y2 = y1(n1) + v2y*t2;
x3 = x2(n2) + v3x*t3 +1/2*a3x*t3.^2;
y3 = y2(n2) + v3y*t3 + 1/2*a3y*t3.^2;
x4 = x3(n3) + v4x*t4;
y4 = y3(n3) + v4y*t4;
x = [x0,x1,x2,x3,x4]';
y = [y0,y1,y2,y3,y4]';

nx = randn(N,M)*delta;
ny = randn(N,M)*delta;

zx = x*ones(1,M) + nx;
zy = y*ones(1,M) + ny;
plot(x,y);
figure;
hold on ;
plot(x0,y0,'r+');
plot(x1,y1,'b+');
plot(x2,y2,'g+');
plot(x3,y3,'y+');
plot(x4,y4,'k+');


zx_est = zeros(N,M);
zy_est = zeros(N,M);

for k = 1:M
     for i = 3:N
         z = [zx(i,k),zy(i,k)]';
         if i < 20 %采用标准的卡尔曼滤波
             if i == 3 %只在等于i==3赋初值
                 vx = (zx(2,k) - zx(1,k))/T;
                 vy = (zy(2,k) - zy(1,k))/T;
                 x_est=[zx(i,1),vx,zy(i,1),vy]';
                 p_est=[delta^2,           delta^2/T,            0,              0;
                              delta^2/T,        2*delta^2/(T^2),  0,              0;
                              0,                    0,                          delta^2,    delta^2/T;
                              0,                    0,                          delta^2/T,  2*delta^2/(T^2)];
             end
             [x_est,p_est] = kalman_std(x_est,p_est,T,z,delta);
             zx_est(i,k) = x_est(1,1);%保存估计值
             zy_est(i,k) = x_est(3,1);
         else % 采用imm卡尔曼滤波
             if i == 20 %只在i==20赋初值
                 x_est = [x_est;0;0];
                 p_est = p_est;
                 p_est(6,6) = 0;
                 for i = 1:3
                     xn_est{i,1}=x_est;
                     pn_est{i,1}=p_est;
                 end
                 u = [0.8,0.1,0.1];
             end
             [x_est,p_est,xn_est,pn_est,u] = imm_test(xn_est,pn_est,T,z,delta,u);
             zx_est(i,k) = x_est(1,1);%保存估计值
             zy_est(i,k) = x_est(3,1);
         end
     end
end


for i = 1:size(zx_est,1)
    zx_est_means(i) = sum(zx_est(i,:))/M;
    zy_est_means(i) = sum(zy_est(i,:))/M;
end
figure;
hold on;
plot(x,y,'r.')
plot(zx_est_means,zy_est_means,'bo');
legend('真实值','估计值');

function [x_est,p_est] = kalman_std(x_est,p_est,T,z,delta)
%模型 [x,vx,y,vy]
    F = [1,T,0,0;0,T,0,0;0,0,1,T;0,0,0,T];
    H = [1,0,0,0;0,0,1,0];
    R = [delta^2,0;0,delta^2];
    Q = 0.01*eye(2);
    G = [T/2,0;1,0;0,T/2;0,1];              
     x_pre = F*x_est;
     p_pre = F*p_est*F' + G*Q*G';
     kg = p_pre*H'*inv(H*p_pre*H'+R);
     x_est = x_pre + kg*(z-H*x_pre);
     p_est = (eye(4)-kg*H)*p_pre;


function [x_est,p_est,xn_est,pn_est,u] = imm_test(xn_est,pn_est,T,z,delta,u)
%模型 [x,vx,y,vy,ax,ay]
p = [0.95,0.025,0.025;
        0.025,0.95,0.025;
        0.025,0.025,0.95];
%模型1
F{1,1} = [1,T,0,0;0,T,0,0;0,0,1,T;0,0,0,T];
F{1,1}(6,6) = 0;
%模型2
F{2,1} = [1,T,0,0,T^2/2,0;
                0,1,0,0,T,0;
                0,0,1,T,0,T^2;
                0,0,0,1,0,T;
                0,0,0,0,1,0;
                0,0,0,0,0,1];
%模型3
F{3,1} = F{2,1};
%模型1，驱动方程
G{1,1} = [T/2,0;
                1,0;
                 0,T/2;
                 0,1];
G{1,1}(6,2)=0;
%模型2，驱动方程
G{2,1} = [T^2/4,0;T^2,0;0,T^2/4;0,T/2;1,0;0,1];
%模型3，驱动方程
G{3,1}=G{2,1};
H = [1,0,0,0,0,0;
         0,0,1,0,0,0];
R = eye(2)*delta^2;
Q = 0.01*eye(2);
mu = zeros(3,3);
c_means = zeros(1,3);

for i = 1:3
    c_means = c_means + p(i,:)*u(i);
end

for i = 1:3 
    mu(i,:) = p(i,:)*u(i) ./ c_means;
end

for j = 1:3
    x0_est{j,1} = zeros(6,1);
    p0_est{j,1} = zeros(6,6);
    for i = 1:3
        x0_est{j,1} = x0_est{j,1} + mu(i,j)*xn_est{i,1};
    end
    for i = 1:3
        p0_est{j,1} = p0_est{j,1} + mu(i,j)*(pn_est{j,1} + (xn_est{i,1} - x0_est{j,1})*(xn_est{i,1} - x0_est{j,1})');
    end
end

for i = 1:3
    x_pre{i,1} = F{i,1}*x0_est{i,1};
    p_pre{i,1} = F{i,1}*p0_est{i,1}*F{i,1}' + G{i,1}*Q*G{i,1}';
    kg{i,1} = p_pre{i,1}*H'*inv(H*p_pre{i,1}*H' + R);
    xn_est{i,1} = x_pre{i,1} + kg{i,1}*(z - H*x_pre{i,1});
    pn_est{i,1} = (eye(6) - kg{i,1}*H)*p_pre{i,1};
end

%更新概率转移值
for i = 1:3
    v{i,1} = z - H*x_pre{i,1};
    s{i,1} = H*p_pre{i,1}*H' + R;
    n = length(s{i,1})/2;
    a(1,i) = 1./sqrt((2*pi)^n*det(s{i,1}))*exp(-1/2*v{i,1}'*inv(s{i,1})*v{i,1});
end

u = a.*c_means./sum(a.*c_means);
x_est = zeros(6,1);
p_est = zeros(6,6);

for i = 1:3
    x_est = x_est + u(i)*xn_est{i,1};
end
for i = 1:3
    p_est = p_est + u(i)*(pn_est{i,1} + (xn_est{i,1} - x_est)*(xn_est{i,1} - x_est)');
end






