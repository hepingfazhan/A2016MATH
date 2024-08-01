%%将第一问程序进行再计算

clc;clear;
rho_sea=1025;%海水密度kg/m3
rho_steel=7850;%钢的密度
m_hc=1200;%重物球的质量kg
v_hc=m_hc/rho_steel;%球的体积
m_fb=1000;%浮标的质量kg
g=9.80665;%重力加速度
H=18;%水面高度

v=36;%风速大小
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

alpha=zeros(1,5);%储存v速度下alpha所得大小

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

%由公式得，x=1.3时,h=2.4165;x=1.4时,h=-1.8857
%%利用二分法对h进行查找%由迭代得,x的范围处于1.27~1.29之间
x_d=1.27;%下界
res_d=2.9582;%下界的值
x_u=1.29;%上界
res_u=-0.9449;%上界的值
eps=1e-05;
% x_fix=1.29;

%%当v>22.9m/s时进行计算
while abs(res_u-x_u)>eps&&abs(res_d-x_d)>eps
    x_fix=0.5*(x_u+x_d);%从中间开始带入x_fix
    x1=[100000,0];
    options = optimoptions('fsolve','Display','none','TolFun',1e-8);
    obj_fun1=@(x)a2016_fun1_1(x,v,g,rho_sea,G_level,x_fix);
    M1=fsolve(obj_fun1,x1,options);
    T1=M1(1);
    theta1=M1(2);
    % disp(M1)

    x2=[10000,0,0];
    obj_fun2=@(x)a2016_fun1_2(x,T1,theta1,G_level);
    M2=fsolve(obj_fun2,x2,options);
    T2=M2(1);
    theta2=M2(2);
    alpha(1)=M2(3)*180/pi;
    H1=cos(M2(3));%钢管在竖直方向的投影长度

    x3=[10000,0,0];
    obj_fun3=@(x)a2016_fun1_2(x,T2,theta2,G_level);
    M3=fsolve(obj_fun3,x3,options);
    T3=M3(1);
    theta3=M3(2);
    alpha(2)=M3(3)*180/pi;
    H2=cos(M3(3));%钢管在竖直方向的投影长度

    x4=[10000,0,0];
    obj_fun4=@(x)a2016_fun1_2(x,T3,theta3,G_level);
    M4=fsolve(obj_fun4,x4,options);
    T4=M4(1);
    theta4=M4(2);
    alpha(3)=M4(3)*180/pi;
    H3=cos(M4(3));%钢管在竖直方向的投影长度

    x5=[10000,0,0];
    obj_fun5=@(x)a2016_fun1_2(x,T4,theta4,G_level);
    M5=fsolve(obj_fun5,x5,options);
    T5=M5(1);
    theta5=M5(2);
    alpha(4)=M5(3)*180/pi;
    H4=cos(M5(3));%钢管在竖直方向的投影长度

    obj_fun6=@(x)[T5*sin(theta5)-x(1)*sin(x(2));T5*cos(theta5)-G_level(6)-x(1)*cos(x(2))];
    M6=fsolve(obj_fun6,[10000,1],options);
    T_mao=M6(1);
    theta_mao=M6(2);
    alpha_maofun=@(x)T5*0.5*sin(x-theta5)-T_mao*0.5*sin(theta_mao-x)+G_level(7)*0.5*sin(x);
    alpha_mao=fsolve(alpha_maofun,0,options);
    alpha(5)=alpha_mao*180/pi;
    H_tong=cos(alpha_mao);

    x7=ones(630,1);
    for i=1:210
        x7(3*i-2)=1000;
    end
    obj_fun7=@(x)a2016_fun1_5(x,T_mao,theta_mao,G_level);
    M7=fsolve(obj_fun7,x7,options);
    T0=M7(end-2);
    theta0=M7(end);

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
    res0=h;%将h的值赋给中间变量的值res
    if abs(abs(res0-x_fix)-abs(res_d-x_d))-abs(abs(res0-x_fix)-abs(res_u-x_u))<0
        x_u=x_fix;
        res_u=res0;
    elseif abs(abs(res0-x_fix)-abs(res_d-x_d))-abs(abs(res0-x_fix)-abs(res_u-x_u))>0
        x_d=x_fix;
        res_d=res0;
    elseif abs(res0-x_fix)<eps
        break
    end
end
%%画图象
hold on;
plot(x0,y0,'k')
xlabel('x轴/m');
ylabel('y轴/m');
disp('钢桶与竖直方向的夹角:');
disp(alpha_mao*180/pi);%输出钢桶与竖直方向的夹角判断
disp('锚链与锚连接处与水平面的夹角:')
disp(90-theta0*180/pi);%输出锚链与锚连接处与水平面的夹角进行判断

%%优化问题
%质量最小且满足alpha_mao<5°，theta0<16°
%1.20195,1.12411,1.0461,0.968,0.88975
%1.6,-0.6
%2,-0.2788
%2.42,0.06
%2.875,0.0422
%3.3566,0.817
h_x=zeros(4,1);
h_x(1)=res0;
for j=1:5
    m_hc=1200+300*j;
    Gc=m_hc*g;%重物的重力
    v_hc=m_hc/rho_steel;
    Fc=rho_sea*g*v_hc;
    G_level(6)=(G5-F5)+(Gc-Fc);
    G_level(7)=Gc-Fc;
    x_u1=[1.21,1.13,1.05,0.97,0.89];
    res_u1=[-0.6,-0.2788,0.06,0.0422,0.817];
    x_d1=[1.2,1.12,1.04,0.965,0.88];
    res_d1=[1.6,2,2.42,2.875,3.3566];
    while abs(res_u1(j)-x_u1(j))>eps&&abs(res_d1(j)-x_d1(j))>eps
        x_fix=0.5*(x_u1(j)+x_d1(j));%从中间开始带入x_fix
        % disp(x_fix);
        x1=[100000,0];
        options = optimoptions('fsolve','Display','none','TolFun',1e-8);
        obj_fun1=@(x)a2016_fun1_1(x,v,g,rho_sea,G_level,x_fix);
        M1=fsolve(obj_fun1,x1,options);
        T1=M1(1);
        theta1=M1(2);
        % disp(M1)

        x2=[10000,0,0];
        obj_fun2=@(x)a2016_fun1_2(x,T1,theta1,G_level);
        M2=fsolve(obj_fun2,x2,options);
        T2=M2(1);
        theta2=M2(2);
        alpha(1)=M2(3)*180/pi;
        H1=cos(M2(3));%钢管在竖直方向的投影长度

        x3=[10000,0,0];
        obj_fun3=@(x)a2016_fun1_2(x,T2,theta2,G_level);
        M3=fsolve(obj_fun3,x3,options);
        T3=M3(1);
        theta3=M3(2);
        alpha(2)=M3(3)*180/pi;
        H2=cos(M3(3));%钢管在竖直方向的投影长度

        x4=[10000,0,0];
        obj_fun4=@(x)a2016_fun1_2(x,T3,theta3,G_level);
        M4=fsolve(obj_fun4,x4,options);
        T4=M4(1);
        theta4=M4(2);
        alpha(3)=M4(3)*180/pi;
        H3=cos(M4(3));%钢管在竖直方向的投影长度

        x5=[10000,0,0];
        obj_fun5=@(x)a2016_fun1_2(x,T4,theta4,G_level);
        M5=fsolve(obj_fun5,x5,options);
        T5=M5(1);
        theta5=M5(2);
        alpha(4)=M5(3)*180/pi;
        H4=cos(M5(3));%钢管在竖直方向的投影长度

        obj_fun6=@(x)[T5*sin(theta5)-x(1)*sin(x(2));T5*cos(theta5)-G_level(6)-x(1)*cos(x(2))];
        M6=fsolve(obj_fun6,[10000,1],options);
        T_mao=M6(1);
        theta_mao=M6(2);
        alpha_maofun=@(x)T5*0.5*sin(x-theta5)-T_mao*0.5*sin(theta_mao-x)+G_level(7)*0.5*sin(x);
        alpha_mao=fsolve(alpha_maofun,0,options);
        alpha(5)=alpha_mao*180/pi;
        H_tong=cos(alpha_mao);

        x7=ones(630,1);
        for i=1:210
            x7(3*i-2)=1000;
        end
        obj_fun7=@(x)a2016_fun1_5(x,T_mao,theta_mao,G_level);
        M7=fsolve(obj_fun7,x7,options);
        T0=M7(end-2);
        theta0=M7(end);

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
        if abs(abs(res-x_fix)-abs(res_d1(j)-x_d1(j)))-abs(abs(res-x_fix)-abs(res_u1(j)-x_u1(j)))<0
            x_u1(j)=x_fix;
            res_u1(j)=res;
        elseif abs(abs(res-x_fix)-abs(res_d1(j)-x_d1(j)))-abs(abs(res-x_fix)-abs(res_u1(j)-x_u1(j)))>0
            x_d1(j)=x_fix;
            res_d1(j)=res;
        elseif abs(res-x_fix)<eps
            break
        end
    end
    h_x(j+1)=res;
    color=['b','y','g','r','m'];
    hold on;
    plot(x0,y0,color(j));
    xlabel('x轴/m');
    ylabel('y轴/m');
    print0=sprintf('36m/s,以300kg为间隔的重物,Ⅱ型锚链下的锚链形状');
    title(print0);
    % disp(res);
    if (90-theta0*180/pi)<=16&&(alpha_mao*180/pi)<=5%跳出循环条件
        print1=sprintf('重物球为%dkg时的钢桶倾斜角度(Deg)',m_hc);
        disp(print1)
        disp(alpha_mao*180/pi);
        print2=sprintf('重物球为%dkg时锚链与海平面的夹角(Deg)',m_hc);
        disp(print2);
        disp(90-theta0*180/pi);
        break
    end
end
%%由上述程序可知，m_hc在区间2100~2400kg之间
legend('36m/s,1200kg重物,Ⅱ型锚链下的锚链形状','36m/s,1500kg重物,Ⅱ型锚链下的锚链形状','36m/s,1800kg重物,Ⅱ型锚链下的锚链形状','36m/s,2100kg重物,Ⅱ型锚链下的锚链形状','36m/s,2400kg重物,Ⅱ型锚链下的锚链形状');
figure;
plot(1200:300:2400,h_x,'b^');
%进行线性插值表示h与m的关系
xi=linspace(1200,2400,10000);
yi=interp1(1200:300:2400,h_x,xi,'linear');
hold on;
plot(xi,yi,'r-');
title('重物球质量m与超出水面高度h的函数关系');
xlabel('重物球质量m(kg)');
ylabel('超出水面高度h(m)');


m_u=2400;
m_d=2100;
while true
    m_hc=(m_u+m_d)*0.5;
    Gc=m_hc*g;%重物的重力
    v_hc=m_hc/rho_steel;
    Fc=rho_sea*g*v_hc;
    G_level(6)=(G5-F5)+(Gc-Fc);
    G_level(7)=Gc-Fc;
    h0=interp1(1200:300:2400,h_x,m_hc,'linear');
    % disp(h0);
    x1=[100000,0];
    options = optimoptions('fsolve','Display','none','TolFun',1e-8);
    obj_fun1=@(x)a2016_fun1_1(x,v,g,rho_sea,G_level,h0);
    M1=fsolve(obj_fun1,x1,options);
    T1=M1(1);
    theta1=M1(2);

    x2=[10000,0,0];
    obj_fun2=@(x)a2016_fun1_2(x,T1,theta1,G_level);
    M2=fsolve(obj_fun2,x2,options);
    T2=M2(1);
    theta2=M2(2);
    alpha(1)=M2(3)*180/pi;
    H1=cos(M2(3));%钢管在竖直方向的投影长度

    x3=[10000,0,0];
    obj_fun3=@(x)a2016_fun1_2(x,T2,theta2,G_level);
    M3=fsolve(obj_fun3,x3,options);
    T3=M3(1);
    theta3=M3(2);
    alpha(2)=M3(3)*180/pi;
    H2=cos(M3(3));%钢管在竖直方向的投影长度

    x4=[10000,0,0];
    obj_fun4=@(x)a2016_fun1_2(x,T3,theta3,G_level);
    M4=fsolve(obj_fun4,x4,options);
    T4=M4(1);
    theta4=M4(2);
    alpha(3)=M4(3)*180/pi;
    H3=cos(M4(3));%钢管在竖直方向的投影长度

    x5=[10000,0,0];
    obj_fun5=@(x)a2016_fun1_2(x,T4,theta4,G_level);
    M5=fsolve(obj_fun5,x5,options);
    T5=M5(1);
    theta5=M5(2);
    alpha(4)=M5(3)*180/pi;
    H4=cos(M5(3));%钢管在竖直方向的投影长度

    obj_fun6=@(x)[T5*sin(theta5)-x(1)*sin(x(2));T5*cos(theta5)-G_level(6)-x(1)*cos(x(2))];
    M6=fsolve(obj_fun6,[10000,1],options);
    T_mao=M6(1);
    theta_mao=M6(2);
    alpha_maofun=@(x)T5*0.5*sin(x-theta5)-T_mao*0.5*sin(theta_mao-x)+G_level(7)*0.5*sin(x);
    alpha_mao=fsolve(alpha_maofun,0,options);
    alpha(5)=alpha_mao*180/pi;
    H_tong=cos(alpha_mao);
    % disp(alpha_mao*180/pi);

    x7=ones(630,1);
    for i=1:210
        x7(3*i-2)=1000;
    end
    obj_fun7=@(x)a2016_fun1_5(x,T_mao,theta_mao,G_level);
    M7=fsolve(obj_fun7,x7,options);
    T0=M7(end-2);
    theta0=M7(end);
    % disp(90-theta0*180/pi);

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
    % disp(h);
    res=h;%将h的值赋给中间变量的值res
    if (90-theta0*180/pi)<=16&&(alpha_mao*180/pi)<=5
        m_u=m_hc;
    else
        m_d=m_hc;
    end
    if abs(m_u-m_d)<eps
        break;
    end
end
disp('符合条件的m的最小值')
disp(m_hc);
disp('此时的钢桶倾斜角度');
disp(alpha_mao*180/pi);
disp('此时的锚链与海平面夹角');
disp(90-theta0*180/pi);
figure;
plot(x0,y0,'m');
title('m最小时的锚链形状');
xlabel('x轴');
ylabel('y轴');

