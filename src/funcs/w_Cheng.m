%{
+[保存文本]"/home/xiantao/Documents/Projects/OTS/src/funcs/w_Cheng.m"(,w_Cheng)
%}

function f_w=w_Cheng(xi)

x1=100;
N1=10001;

g=@(x) (1)./(sqrt(pi)).*exp(-x.^2);
f=@(x) g(x).*(2.*x.^2+x.*xi.*log((xi-x)./(xi+x)));

x=linspace(0,x1,N1);
f1=f(x);
f1(isnan(f1))=0;
f1(isinf(f1))=0;
f_w=trapz(x,f(x))*2;


end

%{
+[M函数](,验证公式)
%}