%{
+[保存文本]"/home/xiantao/Documents/Projects/OTS/src/funcs/w_Cheng_sg.m"(,w_Cheng_sg)
%}

function f_w=w_Cheng_sg(xi,P)

x1=100;
N1=10001;

g=@(x) (1)./(sqrt(pi)).*exp(-(x.^2).^P);
f=@(x) P*g(x).*x.^(2*P-2).*(2.*x.^2+x.*xi.*log((xi-x)./(xi+x)));

x=linspace(0,x1,N1);
f1=f(x);
g0=g(x);
f1(isnan(f1))=0;
f1(isinf(f1))=0;
f_w=trapz(x,f1)./trapz(x,g0);


end

%{
+[M函数](,验证公式)
%}