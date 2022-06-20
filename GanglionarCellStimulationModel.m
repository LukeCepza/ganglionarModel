%%Proyecto Final
% By Luis Kevin Cepeda Zapata 
h = 0.01;
N = length(0:h:20);

%constantes
C_M         = 1;        %uF/cm2
g_Na        = 50;       %ms/cm2
g_K         = 12;       %ms/cm2
g_Ca        = 2.2 ;     %ms/cm2
g_AK        = 36 ;      %ms/cm2
g_CaK       = 0.05;     %mS /cm2
Ca2_diss    = 1;        
E_Na        = 35;       % mV
E_K         = -75;      % mV

m    = zeros(1,N+1); H   = zeros(1,N+1);  Ca_n       = zeros(1,N+1); 
A    = zeros(1,N+1); h_A = zeros(1,N+1);  Ca2_in  = zeros(1,N+1); 
V_M_G  = zeros(1,N+1);K_n       = zeros(1,N+1);
E_Ca = zeros(1,N+1); G_K_g = zeros(1,N+1); G_Na_g = zeros(1,N+1);
G_AK_g = zeros(1,N+1);G_Ca_g = zeros(1,N+1);G_CaK = zeros(1,N+1);
I_Ca_g = zeros(1,N+1); I_Na_g = zeros(1,N+1);I_ION_g = zeros(1,N+1);
I_GK_g = zeros(1,N+1);I_G_AK_g  = zeros(1,N+1); I_G_CaK_g= zeros(1,N+1);

%E_Calcio
Ca2_ext = 1800; %mol
R = 8.31;      %CONSTANTE DE GAS
T = 293;       %°K
F = 96500;     %Constante faradhay

%Condiciones iniciales
m(1) = 0.02794;
V_M_G(1) = -70;
H(1) = 0.8871;
Ca_n(1) = 0.003019;
K_n(1) = 0.107809;
A(1) = 0.070467;
h_A(1) = 0.300453;
Ca2_in(1) = 0.1;    

r = 10; %radio de la celula uM
tau_ca = 50;
Ca_r = 0.1; % Calcio reposo

t = (1:N+1)*h;

for k = 1:3       
    J_E = [10,20,40];
    J_E = J_E(k);
    for i = 1:N
        E = V_M_G(i);

        G_Na_g(i)    = g_Na * m(i)^3 * H(i);
        G_Ca_g(i)    = g_Ca * Ca_n(i)^3;
        G_K_g(i)     = g_K  * K_n(i)^4;
        G_AK_g(i)    = g_AK * A(i)^3 * h_A(i);
        G_CaK(i)   = g_CaK.*(Ca2_in(i)/Ca2_diss)^2/(1+(Ca2_in(i)/Ca2_diss)^2); %Value of G_Cak at t = i+1

        E_Ca(i) = 1000*R*T/(2*F)*log(Ca2_ext/Ca2_in(i));
        I_Ca_g(i) = G_Ca_g(i)*(V_M_G(i)-E_Ca(i));

        Ca2_in(i+1) = h*(-3*I_Ca_g(i)/(2*F*r) - (Ca2_in(i)-Ca_r)/tau_ca) + Ca2_in(i); % Here we solve for Ca2(i+1)

        I_Na_g(i) = G_Na_g(i)*(V_M_G(i)-E_Na);
        I_GK_g(i) = G_K_g(i)*(V_M_G(i)-E_K);
        I_G_AK_g(i)  = G_AK_g(i)*(V_M_G(i)-E_K);
        I_G_CaK_g(i) = G_CaK(i)*(V_M_G(i)-E_K);

        I_ION_g(i) = I_Na_g(i)+I_Ca_g(i)+I_GK_g(i)+ I_G_AK_g(i) +I_G_CaK_g(i);

        m(i+1) = gateconstant(h,fun_m(E),m(i));
        H(i+1) = gateconstant(h,fun_H(E),H(i));
        Ca_n(i+1) = gateconstant(h,fun_c(E),Ca_n(i));
        K_n(i+1) = gateconstant(h,fun_n(E),K_n(i));
        A(i+1) = gateconstant(h,fun_A(E),A(i));
        h_A(i+1) = gateconstant(h,fun_h_A(E),h_A(i));

        V_M_G(i+1) = h*(J_E-I_ION_g(i))/C_M + V_M_G(i);
    end
    subplot(3,1,k)
    plot(t,V_M_G, "LineWidth",2);
    title("Membrane Voltage With Constant current J_E = "+J_E+"pA")
    xlim([0 20])
end

subplot(3,1,2)
ylabel("Voltage (mV) // Current uA/Cm^2");
subplot(3,1,1)
legend("Membrane Voltage (V_m)", "Location","northeast");
subplot(3,1,3)
xlabel("Time (ms)");

xlim([0,20])

function a_b = fun_m(E)
    alpha = (-0.6*(E+30))/(exp(-0.1*(E+30))-1);
    quadbeta = 20*exp(-(E+55)/18);
    a_b = [alpha, quadbeta];
end
function a_b = fun_H(E)
    alpha     =   0.4*exp(-(E+50)/20);
    quadbeta  =   6/(exp(-0.1*(E+20))+1);
    a_b = [alpha, quadbeta];

end
function a_b = fun_c(E)
    alpha     =  -0.3*(E+13)/(exp(1)^(-0.1*(E+13))-1);
    quadbeta  =   10*exp(1)^(-(E+38)/18);
    a_b = [alpha, quadbeta];

end
function a_b = fun_n(E)
    alpha     =  -0.02*(E+40)/(exp(1)^(-0.1*(E+40))-1);
    quadbeta  =   0.4*exp(1)^(-(E+50)/80);
    a_b = [alpha, quadbeta];

end
function a_b = fun_A(E)
    alpha     =   (-0.006*(E+90))/(exp(1)^(-0.1*(E+90))-1);
    quadbeta  =   0.1*exp(1)^(-(E+30)/10);
    a_b = [alpha, quadbeta];

end 
function a_b = fun_h_A(E)
    alpha   =   0.04*exp(1)^(-(E+70)/20);
    quadbeta=   (0.6)/(exp(1)^(-0.1*(E+40))+1);
    a_b = [alpha, quadbeta];

end
function final = gateconstant(h,a_b,inicial)
    final = h*(a_b(1)*(1-inicial)-a_b(2)*inicial)+inicial;
end