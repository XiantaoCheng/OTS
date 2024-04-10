%{
地址::/home/xiantao/Documents/Projects/OTS/src/funcs/w_ZhaoPaper2.m
+[保存文本](,w_ZhaoPaper2)
%}

function y=w_ZhaoPaper2(x)

y=1-2.*x.*exp(-x.^2).*trapz(linspace(0,x,1000),exp(linspace(0,x,1000).^2))+1i.*sqrt(pi).*x.*exp(-x.^2);

A1=-0.0007745012789756754*0;
A2=0.020428443946704566*0;
A3=-0.15289532104922568*0;

if abs(x)>10
    y=A3-A2*x+A1*x.^2;
end

end


%{
复制变量"A(3)"
%}