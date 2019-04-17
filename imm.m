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

