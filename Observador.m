%Calculo del observador
%Alumno: Herman Raymundo Salvatierra
%TP 2 item 2
%Profesor : Pucheta
%clear all, close all, clc
%clc
%variables
%Laa = 366e-6;  J = 5e-9;Ra = 55.6 ; B = 0;Ki = 6.49e-3 ; Km = 6.53e-3;

Laa= 5e-3;
J= 0.004;
Ra= 0.2;
B= 0.005;
Ki= 6.5e-5;
Km= 0.055;

%El Torque es dato
Tl=1.15e-3;%va a ir variando para pi/2 y -pi/2

%Planteo matrices A B C D de estado a partir de las ecuaciones del modelo
%siendo
% x1 = ia    :la corriente que circula
% x2 = wr    :la velocidad angular
% x3 = theta :el angulo posicion del motor

Mat_A=[-Ra/Laa -Km/Laa 0;Ki/J -B/J 0 ;0 1 0 ];
Mat_B=[1/Laa;0;0];
Mat_C=[0 0 1]; 
Mat_D=[0];

%Planteo las Matrices ampliadas
Mat_Aa=[Mat_A zeros(3,1); -Mat_C 0];
Mat_Ba=[Mat_B; 0];
Mat_Cc=[Mat_C 0];

%Diseño del controlador por LQR
Q=diag([1 1/10000 100 10000000000000]);
R=10.0;
%Q=diag([1 1/10000 100 1000000]);R=0.1; 
%Si achico R tiende a ser mas rapida la respuesta, o tambien voy aumentando Q

%Construcción del Hamiltoniano para el cálculo del controlador
Ha=[Mat_Aa -Mat_Ba*inv(R)*Mat_Ba'; -Q -Mat_Aa'];
[n,va]=size(Ha);
[V,D]=eig(Ha);Mx1x2=[];
for ii=1:n
 if real(D(ii,ii))<0
 Mx1x2=[Mx1x2 V(:,ii)];
 end
end
Mx1=Mx1x2(1:n/2,:); Mx2=Mx1x2(n/2+1:end,:);

P=real(Mx2*inv(Mx1));%tomo la parte real por un tema de residuos
Ka=inv(R)*Mat_Ba'*P;
K=Ka(1:3); KI=Ka(4);    
eig(Mat_Aa-Mat_Ba*Ka);%verifico polos con parte real negativa
%Fin cálculo del controlador-

    %Cálculo del Observador---------------------------------------------------
Mat_A_d=Mat_A';
Mat_B_d=Mat_C';
Mat_C_d=Mat_B';

%Fui probando con distintos valores de Q y R hasta dar alguno que me diera
%un comportamiento razonable
Qo = diag([1 1/10000 200000]);
Ro = 10;



%Construcción del Hamiltoniano para el cálculo del observador
Ho=[Mat_A_d -Mat_B_d*inv(Ro)*Mat_B_d'; -Qo -Mat_A_d'];
[no,va]=size(Ho);
[V,D]=eig(Ho);Mx1x2=[];
for ii=1:no
 if real(D(ii,ii))<0
 Mx1x2=[Mx1x2 V(:,ii)];
 end
end
Mx1=Mx1x2(1:no/2,:); Mx2=Mx1x2(no/2+1:end,:);
Po=real(Mx2*inv(Mx1));
Ko=(inv(Ro)*Mat_B_d'*Po)';%controlador con observador
eig(Mat_A-Ko*Mat_C);
%FIN cálculo del Observador

J_(1)=0; V_(1)=0; 
psi(1)=0;
x_hat=[0;0;0];
%inicializo el observador

%para la simulacion
delta_t=1e-3;
tiempo=10;
pasos=round(tiempo/delta_t);
Ci=[0 0 0 0];
t=0:delta_t:(tiempo-delta_t);

%referencia
ref=(pi/2)*square(2*pi*t/4);

%Torque
Tll=(Tl/2)+(Tl/2)*square(2*pi*t/4);%funcion para ir variando el torque 

%condiciones iniciales
x=zeros(4,pasos);
x(1,1)=Ci(1);%ia
x(2,1)=Ci(2);%w
x(3,1)=Ci(3);%theta
x(4,1)=Ci(4);

ua(1)=0;

for i=2:1:pasos
    estado=[x(1,i-1);x(2,i-1);x(3,i-1);x(4,i-1)];%guardo el estado
    integracion=x(4,i-1)+delta_t*(ref(1,i-1)-Mat_Cc*estado);
    
%     u_actual=-K*estado(1:3)-integracion*KI;color='r';%sin observador
    u_actual=-K*x_hat-integracion*KI;color='g';%Con observador
    
    ua=[ua u_actual];
    
    %ecuaciones del motor
    ia_p=(-Ra/Laa)*estado(1)-(Km/Laa)*estado(2)+(1/Laa)*u_actual;
    w_p=(Ki/J)*estado(1)-(B/J)*estado(2)-(1/J)*Tll(i-1);
    theta_p=estado(2);
    
    xp_actual=[ia_p; w_p; theta_p];
    
    xsig=estado(1:3)+delta_t*xp_actual;
    x(1,i)=xsig(1);
    x(2,i)=xsig(2);
    x(3,i)=xsig(3);
    x(4,i)=integracion;
    
    %___OBSERVADOR___
    y_sal_O(i)=Mat_C*x_hat;
    y_sal(i)=Mat_Cc*estado;
    x_hatp=Mat_A*x_hat+Mat_B*u_actual+Ko*(y_sal(i)-y_sal_O(i));
    x_hat=x_hat+delta_t*x_hatp;
    
end

figure(1)
plot(t,x(3,:),color);title('Angulo Theta y referencia');hold on;
plot(t,ref);hold on;
figure(2)
plot(t,x(2,:),color);title('Velocidad angular, w');hold on
figure(3)
plot(t,x(1,:),color);title('Corriente Ia');hold on;
 %figure(4)
% plot(t,ua,color);title('Accion de control');hold on;
% figure(5)
% plot(t,Tll,color);title('Torque');hold on;
