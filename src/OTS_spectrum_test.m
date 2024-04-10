%{
地址::/home/xiantao/Documents/Projects/OTS/src/OTS_spectrum_test.m
+[保存文本](,OTS_spectrum_test)
%}

addpath("./funcs");
physics_constant;

k_s=[1,0,0];
s=[1,1,0];
E_i0=[0,0,1];
lm_0=263e-9;
T_e=1.4e3*e/kB;
T_i=0.4e3*e/kB;
P=1;

lm_d=50e-9;
Z_i=1;
%{
m_p
+[M函数](,Zhao论文)
%}

s=s/norm(s);
omega_0=2*pi*c/lm_0;
k_0=2*pi/lm_0*k_s;
v_te=sqrt(k_B.*T_e/m_e);
v_ti=sqrt(k_B.*T_i/m_p);

%{
v_te
+[M函数](,验证公式)
%}
% f_e=@(v) 1/sqrt(2*pi)/v_te*exp(-norm(v)^2/(2*v_te^2));
f_e=@(v) 1/sqrt(2*pi)/v_te*exp(-(norm(v)^2/(2*v_te^2)).^P);
f_i=@(v) 1/sqrt(2*pi)/v_ti*exp(-norm(v)^2/(2*v_ti^2));

w=@(x) w_ZhaoPaper2(x);
% w=@(x) w_Cheng_sg(x,P);
% w=@(x) w_Cheng(x);
chi_e=@(omega,k) (1./k./lm_d).^2.*w(omega./(sqrt(2).*k.*v_te));
chi_i=@(omega,k) (1./k./lm_d).^2.*Z_i.*(T_e/T_i).*w(omega./(sqrt(2).*k.*v_ti));
epsilon=@(omega,k) 1+chi_e(omega,k)+chi_i(omega,k);

%{
+[M函数](,Zhao论文)
%}

n_e=@(k,omega) 2*pi/k*abs(1-chi_e(omega,k)./epsilon(omega,k))^2*f_e(omega./k)+2*pi/k*abs(chi_e(omega,k)./epsilon(omega,k))^2*f_i(omega./k);

P_s=@(omega) (1)./(16.*pi.^2).*(e.^4)./(m_e.^2.*c.^3).*norm(cross(s,(cross(s,E_i0)))).^2.*abs(n_e(norm(k_0-(s)./(c).*(omega-omega_0)),omega-omega_0)).^2;


%{
+[M函数](,Zhao论文)
%}

lms=linspace(220e-9,310e-9,10000);
I=zeros(size(lms));
ne=zeros(size(lms));
fe=zeros(size(lms));
chis=zeros(size(lms));
chis_ep=zeros(size(lms));

for i=1:length(lms)
    omega=2*pi*c/lms(i);
    Dk=k_0-(s)./(c).*(omega-omega_0);
    Domega=omega-omega_0;
    omega=2*pi*c/lms(i);
    I(i)=P_s(omega);
    chis(i)=chi_e(Domega,norm(Dk));
    chiis(i)=chi_i(Domega,norm(Dk));
    if isnan(chiis(i))
        chiis(i)=0;
    end
    chis_ep(i)=abs(1-chis(i)/(1+chis(i)+chiis(i)))^2;
    chis_ip(i)=abs(chis(i)/(1+chis(i)+chiis(i)))^2;
    ne(i)=n_e(norm(Dk),Domega);
    fe(i)=f_e(Domega./norm(Dk));
end

chis=chis/max(chis);
chiis=chiis/max(chiis);
chiis(isnan(chiis))=0;
chis_ep=chis_ep/max(chis_ep);
chis_ip=chis_ip/max(chis_ip);
fe=fe/max(fe);
I=I/max(I);

subplot(1,3,1)
plot(lms/1e-9,I)
xlabel("Wavelength [nm]")
ylabel("Intensity [a.u.]")
axis([220,310,0,1.1])

subplot(1,3,2)
plot(lms/1e-9,I)
xlabel("Wavelength [nm]")
ylabel("Intensity [a.u.]")
axis([262,264,0,1.1])

subplot(1,3,3)
plot(lms/1e-9,I)
xlabel("Wavelength [nm]")
ylabel("Intensity [a.u.]")
axis([220,310,0,0.01])

%{
omega=omega_0
+[M函数](,验证公式)
%}