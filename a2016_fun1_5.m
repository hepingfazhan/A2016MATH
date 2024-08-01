function f=a2016_fun1_5(x,T_mao,theta_mao,G_level)
f(1)=T_mao*sin(theta_mao-x(3))-x(1)*sin(x(3)-x(2));
f(2)=x(1)*sin(x(2))-T_mao*sin(theta_mao);
f(3)=x(1)*cos(x(2))-T_mao*cos(theta_mao)+G_level(9);
for i =2:210
    f(3*i-2)=x(3*i-5)*sin(x(3*i-4)-x(3*i))-x(3*i-2)*sin(x(3*i)-x(3*i-1));
    f(3*i-1)=x(3*i-2)*sin(x(3*i-1))-x(3*i-5)*sin(x(3*i-4));
    f(3*i)=x(3*i-2)*cos(x(3*i-1))-x(3*i-5)*cos(x(3*i-4))+G_level(9);
end
end

