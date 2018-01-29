function UNIQUAC
%------Para cambiar los modelos, se deben activar los par�metros----------%
%------marcados con &&& en los comentarios, y desactivar aquellos---------%
%------modelos que no se usar�n-------------------------------------------%
 
%composiciones l�quidas de agua%
x(1) = 0.018;
x(2) = 0.052;
x(3) = 0.112;
x(4) = 0.172;
x(5) = 0.208;
x(6) = 0.251;
x(7) = 0.314;
x(8) = 0.47;
x(9) = 0.993;
x(10) = 0.994;
x(11) = 0.996;
x(12) = 0.997;
x(13) = 0.998;
 
%---- para todo el codigo el sub�ndice 1 corresponde al agua ----%
%---- el sub�ndice 2 corresponde al ciclohexanol -----%
 
% par�metros de volumen, (Gmehling, Onken y Arlt, 1984). %
r1 = 0.92; 
r2 = 4.3489;
% par�metros de area, (Gmehling, Onken y Arlt, 1984). %
q1 = 1.4;
q2 = 3.512;
% n�mero de coordinaci�n (Gmehling, Onken y Arlt, 1984). %
z=10;
 
l1 =(z/2)*(r1-q1)-(r1-1);
l2 =(z/2)*(r2-q2)-(r2-1);
 
 
% constante universal R y temperatura del sistema (T=90�C)
R = 1.98721;  %[=] Cal/mol*K
T = 363.15; %grad Kelvin
 
 
% temperatura de saturaci�n, ecuaci�n de Antoine
Ta= 90; %[=]�C
% Presi�n de saturaci�n de agua, Antoine (Gmehling, Onken y Arlt, 1984). %
pas = 10^(8.07131 - (1730.63/(Ta+233.426))); %[=] mmHg 
% Presi�n de saturaci�n del ciclohexanol, Antoine (Gmehling, Onken y Arlt, 1984).%
pbs = 10^(8.35237 - (2258.56/(Ta+251.624))); %[=] mmHg 
 
 
 
%Definici�n de taos (Gmehling, Onken y Arlt, 1984). %
T12= exp(-255.3078/(R*T));  %-- &&&
T21 = exp(-316.5399/(R*T)); %-- &&&
 
%Definici�n de taos (Lazz�s, 2010) %
%T12= exp(((919.49-233.22))/(R*T)); %-- &&& 
%T21 = exp((-(919.49-233.22))/(R*T)); %-- &&&
 
for i=1:13
 
%xb representa la fracci�n liquida de ciclohexanol%
%x(i) representa la fracci�n liquida del agua%
xb= 1-x(i); 
 
%fracci�n de �rea de los componentes%
theta1 = q1*x(i)/(q2*xb+q1*x(i));
theta2 = q2*xb/(q1*x(i)+ q2*xb);
%fracci�n de volumen de los compenentes%
phi1 = r1*x(i)/(r2*xb + r1*x(i));
phi2 = r2*xb/(r1*x(i)+ r2*xb);
 
%coeficiente de actividad del agua%
Ac1 = exp( log(phi1/x(i)) + (z/2)*q1*log(theta1/phi1) +phi2*(l1 -(r1/r2)*l2)...
    -q1*log(theta1+theta2*T21) + theta2*q1*( (T21/(theta1+theta2*T21)) ...
    - (T12/(theta1*T12 +theta2))));
 
%coeficiente de actividad del ciclohexanol%
Ac2 = exp(  log(phi2/xb) +(z/2)*q2*log(theta2/phi2) +phi1*(l2 - (r2/r1)*l1)...
    -q2*log(theta1*T12+theta2) + theta1*q2*((T12/(theta1*T12+theta2))...
    - (T21/(theta1+theta2*T21))));
 
%presi�n de burbuja, Raoult modificado%
Pb= x(i)*Ac1*pas + xb*Ac2*pbs;
 
%fracci�n de vapor del agua%
ya= pas*Ac1*x(i)/Pb;
%fracci�n de vapor del ciclohexanol%
yb= 1-ya;
 
%matriz que almacena todos los valores%
 mat(i,1) = Pb;
 mat(i,2) = x(i);
 mat(i,3) = xb;
 mat(i,4) = ya;
 mat(i,5) = yb;
 mat(i,6) = Ac1;
 mat(i,7) = Ac2; 
end
 
mat
 
end
