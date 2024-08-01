function f=a2016_fun1_2(x,T,theta,G)
%x(1)为T2[第一根钢管对第二根钢管]
%x(2)为theta2[T2与垂直面的夹角]
%x(3)为alpha1[第一根钢管与竖直方向的夹角]
f=zeros(3,1);
f(1)=T*theta-x(1)*sin(x(2));
f(2)=x(1)*cos(x(2))-T*cos(theta)+G(2);
f(3)=T*0.5*sin(x(3)-theta)-x(1)*0.5*sin(x(2)-x(3));
%%以此类推
end

