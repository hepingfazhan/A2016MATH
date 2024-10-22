function f=a2016_fun1_3(x,T5,theta5,G)
%x(1)为T锚[下方锚链对其的拉力]
%x(2)为theta锚[下方锚链的拉力T锚与竖直平面的夹角]
%x(3)为alpha锚[钢桶与竖直方向的夹角]
f=zeros(3,1);
f(1)=T5*sin(theta5)-x(1)*sin(x(2));
f(2)=T5*cos(theta5)-G(6)-x(1)*cos(x(2));
f(3)=T5*0.5*sin(x(3)-theta5)-x(1)*0.5*sin(x(2)-x(3))+G(7)*0.5*sin(x(3));
end

