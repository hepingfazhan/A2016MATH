clc;clear;
rho_sea=1025;%海水密度kg/m3
rho_steel=7850;%钢的密度
m_hc=1200;%重物球的质量kg
v_hc=m_hc/rho_steel;%球的体积
m_fb=1000;%浮标的质量kg
g=9.80665;%重力加速度
H=18;%水面高度

%%风速12m/s时
%浮标受力分析
v=[12,24];
A=pi*(1^2);%浮标底面积公式
G_fb=m_fb*g;%浮标重量
G_tube=10*g.*ones(1,4);%钢管的重力（向量）
F_tube=rho_sea*g*((0.025^2)*pi*1).*ones(1,4);%钢管的浮力（向量）
G5=100*g;%钢桶的重力
F5=rho_sea*g*((0.15^2)*pi*1);%钢桶的浮力
Gc=m_hc*g;%重物的重力
Fc=rho_sea*g*v_hc;%重物圆的浮力
s=22.05;
G_mao=7*g*s;%锚链的重力
F_mao=rho_sea*g*(7*s/rho_steel);%锚链的浮力，假设锚链为纯钢
G_single=7*22.05*g/210;
F_fusingle=rho_sea*g*(7*22.05/210)/rho_steel;
alpha_1=zeros(1,5);%储存v1速度下alpha所得大小
alpha_2=zeros(1,5);%储存v2速度下alpha所得大小

G_level=zeros(1,9);%%算出对于此系统的重力与浮力的差值
for i = 1:9
    if i==1%浮标
        G_level(i)=G_fb;
    elseif i<=5&&i>=2%钢管
        G_level(i)=G_tube(i-1)-F_tube(i-1);
    elseif i==6%钢桶+重物球
        G_level(i)=(G5-F5)+(Gc-Fc);
    elseif i==7
        G_level(i)=Gc-Fc;
    elseif i==8%锚链总重
        G_level(i)=(G_mao-F_mao);
    elseif i>8
        G_level(i)=G_single-F_fusingle;%单个锚链的准重力
    end
end

%%使用迭代法
%对h初始值进行估计,假设theta0=0°，由悬链线公式可近似得h在1.3~1.32之间
%由公式得，x=1.3时,h=2.4165;x=1.4时,h=-1.8857
x_d0=1.31;
res_d0=5.0084;
x_u0=1.32;
res_u0=-0.16;
eps1=1e-04;
% x_fix=1.31712;
%利用迭代法进行计算

%%若将v作为未知量进行带入运算，可得知：
%当锚链完全带起时的风速约为22.9m/s
%则当v1=12m/s时，有部分锚链未翻起；而v2=24m/s时，锚链全部翻起

%%当v<22.9m/s时进行计算
while abs(res_u0-x_u0)>eps1&&abs(res_d0-x_d0)>eps1
    x_fix=0.5*(x_u0+x_d0);%从中间开始带入x_fix
    x1=[10000,0];
    options = optimoptions('fsolve','Display','none','TolFun',1e-8);
    obj_fun1=@(x)a2016_fun1_1(x,v(1),g,rho_sea,G_level,x_fix);
    M1=fsolve(obj_fun1,x1,options);
    T1=M1(1);
    theta1=M1(2);
    % disp(M1)

    x2=[10000,0,0];
    obj_fun2=@(x)a2016_fun1_2(x,T1,theta1,G_level);
    M2=fsolve(obj_fun2,x2,options);
    T2=M2(1);
    theta2=M2(2);
    alpha_1(1)=M2(3)*180/pi;
    H1=cos(M2(3));%钢管在竖直方向的投影长度

    x3=[10000,0,0];
    obj_fun3=@(x)a2016_fun1_2(x,T2,theta2,G_level);
    M3=fsolve(obj_fun3,x3,options);
    T3=M3(1);
    theta3=M3(2);
    alpha_1(2)=M3(3)*180/pi;
    H2=cos(M3(3));%钢管在竖直方向的投影长度

    x4=[10000,0,0];
    obj_fun4=@(x)a2016_fun1_2(x,T3,theta3,G_level);
    M4=fsolve(obj_fun4,x4,options);
    T4=M4(1);
    theta4=M4(2);
    alpha_1(3)=M4(3)*180/pi;
    H3=cos(M4(3));%钢管在竖直方向的投影长度

    x5=[10000,0,0];
    obj_fun5=@(x)a2016_fun1_2(x,T4,theta4,G_level);
    M5=fsolve(obj_fun5,x5,options);
    T5=M5(1);
    theta5=M5(2);
    alpha_1(4)=M5(3)*180/pi;
    H4=cos(M5(3));%钢管在竖直方向的投影长度

    obj_fun6=@(x)[T5*sin(theta5)-x(1)*sin(x(2));T5*cos(theta5)-G_level(6)-x(1)*cos(x(2))];
    M6=fsolve(obj_fun6,[10000,1],options);
    T_mao=M6(1);
    theta_mao=M6(2);
    alpha_maofun=@(x)T5*0.5*sin(x-theta5)-T_mao*0.5*sin(theta_mao-x)+G_level(7)*0.5*sin(x);
    alpha_mao=fsolve(alpha_maofun,0,options);
    alpha_1(5)=alpha_mao*180/pi;
    H_tong=cos(alpha_mao);

    x7=ones(630,1);
    for i=1:210
        x7(3*i-2)=1000;
    end
    obj_fun7=@(x)a2016_fun1_4(x,T_mao,theta_mao,G_level);
    M7=fsolve(obj_fun7,x7,options);

    y=zeros(210,1);
    y0=zeros(210,1);
    %寻找M7序列中大于pi/2的后续序列进行截断并分段
    for i= 1:210
        if M7(3*i)>pi/2
            n=i;%n为找到的最后一个alpha小于pi/2的单位锚链
            break
        else
            n=210;
        end
    end
    for i =211-n:210
        y(i)=0.105*cos(M7(3*(n+1)-3*(i+n-210)));
        y0(i)=sum(y);
    end
    Y0=sum(y);

    x=zeros(210,1);
    x0=zeros(210,1);
    for i =211-n:210
        x(i)=0.105*sin(M7(3*(n+1)-3*(i+n-210)));
        x0(i)=0.105*(210-n)+sum(x);
    end
    for i = 1:210-n
        x0(i)=i*0.105;
    end
    X0=0.105*(210-n)+sum(x)+sin(M2(3))+sin(M4(3))+sin(M3(3))+sin(M5(3))+sin(alpha_mao);
    h0=Y0+H_tong+H4+H3+H2+H1+2-H;
    res0=h0;%将h的值赋给中间变量的值res
    if abs(abs(res0-x_fix)-abs(res_d0-x_d0))-abs(abs(res0-x_fix)-abs(res_u0-x_u0))<0
        x_u0=x_fix;
        res_u0=res0;
    elseif abs(abs(res0-x_fix)-abs(res_d0-x_d0))-abs(abs(res0-x_fix)-abs(res_u0-x_u0))>0
        x_d0=x_fix;
        res_d0=res0;
    elseif abs(res0-x_fix)<eps1
        break
    end
end
%%画图象
figure
subplot(1,2,1);
plot(x0,y0,'r')
ylim([0,14])
xlabel('x轴/m');
ylabel('y轴/m');
legend('12m/s下的锚链形状');
disp('12m/s下的各项数据')
disp('钢管与钢桶与竖直方向的倾斜角度为(Rad):')
disp(alpha_1);
disp('浮标水面以下的部分长度h0:')
disp(2-res0);
disp('游动区域半径R:')
disp(X0);

%对h初始值进行估计,假设theta0=0°，由悬链线公式可近似得h在1.3~1.31之间
%由公式得，x=1.3时,h=2.4165;x=1.4时,h=-1.8857
%%利用二分法对h进行查找
x_d=1.3;%下界
res_d=2.4165;%下界的值
x_u=1.31;%上界
res_u=-1.8857;%上界的值
eps2=1e-05;

%%当v>22.9m/s时进行计算
while abs(res_u-x_u)>eps&&abs(res_d-x_d)>eps
    x_fix=0.5*(x_u+x_d);%从中间开始带入x_fix
    x1=[10000,0];
    options = optimoptions('fsolve','Display','none','TolFun',1e-8);
    obj_fun1=@(x)a2016_fun1_1(x,v(2),g,rho_sea,G_level,x_fix);
    M1=fsolve(obj_fun1,x1,options);
    T1=M1(1);
    theta1=M1(2);
    % disp(M1)

    x2=[10000,0,0];
    obj_fun2=@(x)a2016_fun1_2(x,T1,theta1,G_level);
    M2=fsolve(obj_fun2,x2,options);
    T2=M2(1);
    theta2=M2(2);
    alpha_2(1)=M2(3)*180/pi;
    H1=cos(M2(3));%钢管在竖直方向的投影长度

    x3=[10000,0,0];
    obj_fun3=@(x)a2016_fun1_2(x,T2,theta2,G_level);
    M3=fsolve(obj_fun3,x3,options);
    T3=M3(1);
    theta3=M3(2);
    alpha_2(2)=M3(3)*180/pi;
    H2=cos(M3(3));%钢管在竖直方向的投影长度

    x4=[10000,0,0];
    obj_fun4=@(x)a2016_fun1_2(x,T3,theta3,G_level);
    M4=fsolve(obj_fun4,x4,options);
    T4=M4(1);
    theta4=M4(2);
    alpha_2(3)=M4(3)*180/pi;
    H3=cos(M4(3));%钢管在竖直方向的投影长度

    x5=[10000,0,0];
    obj_fun5=@(x)a2016_fun1_2(x,T4,theta4,G_level);
    M5=fsolve(obj_fun5,x5,options);
    T5=M5(1);
    theta5=M5(2);
    alpha_2(4)=M5(3)*180/pi;
    H4=cos(M5(3));%钢管在竖直方向的投影长度

    obj_fun6=@(x)[T5*sin(theta5)-x(1)*sin(x(2));T5*cos(theta5)-G_level(6)-x(1)*cos(x(2))];
    M6=fsolve(obj_fun6,[10000,1],options);
    T_mao=M6(1);
    theta_mao=M6(2);
    alpha_maofun=@(x)T5*0.5*sin(x-theta5)-T_mao*0.5*sin(theta_mao-x)+G_level(7)*0.5*sin(x);
    alpha_mao=fsolve(alpha_maofun,0,options);
    alpha_2(5)=alpha_mao*180/pi;
    H_tong=cos(alpha_mao);

    x7=ones(630,1);
    for i=1:210
        x7(3*i-2)=1000;
    end
    obj_fun7=@(x)a2016_fun1_5(x,T_mao,theta_mao,G_level);
    M7=fsolve(obj_fun7,x7,options);
    T0=M7(end-2);
    theta0=M7(end-1);

    % Y=(s/G_level(8))*(sqrt((G_level(8)+T0*cos(theta0))^2+(T0*sin(theta0))^2)-sqrt((T0*sin(theta0))^2+(T0*cos(theta0))^2));
    y=zeros(210,1);
    y0=zeros(210,1);
    for i =1:210
        y(i)=0.105*cos(M7(633-3*i));
        y0(i)=sum(y);
    end
    Y=sum(y);

    x=zeros(210,1);
    x0=zeros(210,1);
    for i =1:210
        x(i)=0.105*sin(M7(633-3*i));
        x0(i)=sum(x);
    end
    X=sum(x)+sin(M2(3))+sin(M4(3))+sin(M3(3))+sin(M5(3))+sin(alpha_mao);
    h=Y+H_tong+H4+H3+H2+H1+2-H;
    res=h;%将h的值赋给中间变量的值res
    if abs(abs(res-x_fix)-abs(res_d-x_d))-abs(abs(res-x_fix)-abs(res_u-x_u))<0
        x_u=x_fix;
        res_u=res;
    elseif abs(abs(res-x_fix)-abs(res_d-x_d))-abs(abs(res-x_fix)-abs(res_u-x_u))>0
        x_d=x_fix;
        res_d=res;
    elseif abs(res-x_fix)<eps2
        break
    end
end
%%画图象
hold on;
subplot(1,2,2);
plot(x0,y0,'r')
xlabel('x轴/m');
ylabel('y轴/m');
legend('24m/s下的锚链形状');
disp('24m/s下的各项数据')
disp('钢管与钢桶与竖直方向的倾斜角度为(Rad):')
disp(alpha_2);
disp('浮标水面以下的部分长度h:')
disp(2-res);
disp('游动区域半径R:')
disp(X);