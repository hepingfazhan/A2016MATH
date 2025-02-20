function f=a2016_fun1_1(x,v,g,rho1,G,x0)
%x0为超出水面高度h的估计值
%x(1)为T1[第一根钢管对浮标的作用力]
%x(2)为theta1[T1与垂直面的夹角]
f=zeros(2,1);
f(1)=0.625*(v^2)*2*x0-x(1)*sin(x(2));%风力大小
f(2)=rho1*g*pi*(2-x0)-G(1)-x(1)*cos(x(2));
end

