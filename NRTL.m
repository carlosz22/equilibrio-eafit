function NRTL
%------Para cambiar los modelos, se deben activar los parámetros----------%
%------marcados con &&& en los comentarios, y desactivar aquellos---------%
%------modelos que no se usarán-------------------------------------------%
%parámetro alfa para la mezcla (Gmehling, Onken y Arlt, 1984).%
alf=0.3458; %-- &&&
%Parámetro alfa (Steyer y Sundmacher, 2004)%
%alf=0.26904; %-- &&&
%Parámetro alfa (Lazzús, 2010)%
%alf=0.287;%-- &&&
%(Gmehling, Onken y Arlt, 1984).%
g12=3102.6349;
g21=390.4777;
 
 
 
 
%composiciones líquidas del agua%
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
 
% temperatura de saturación, ecuación de Antoine
Ta= 90; %[=]°C
% Presión de saturación de agua, Antoine (Gmehling, Onken y Arlt, 1984). %
pas = 10^(8.07131 - (1730.63/(Ta+233.426))); %[=] mmHg 
% Presión de saturación del ciclohexanol, Antoine (Gmehling, Onken y Arlt, 1984).%
pbs = 10^(8.35237 - (2258.56/(Ta+251.624))); %[=] mmHg 
 
 
 
g11=0; %sustancia pura%
g22=0; %sustancia pura%
 
% constante universal R y temperatura del sistema (T=90°C)
R = 1.98721;  %[=] Cal/mol*K
T = 363.15; %grad Kelvin
R2 = 8.314; %[=] J/molK , únicamente para información de Steyer
 
%definición de taos (Gmehling, Onken y Arlt, 1984).%
T21 = (g21-g11)/(R*T); %Tao21 %-- &&&
T12 = (g12-g22)/(R*T); %Tao12  %-- &&&
 
 
%definición de taos (Lazzús, 2010)%
%T21 = (99.55-1029.1)/(R*T); %Tao21 %-- &&&
%T12 = (1029.1-99.55)/(R*T); %Tao12 %-- &&&
 
 
%definición de taos (Steyer y Sundmacher, 2004)%
%T21 = (-569.408)/(R2*T); %Tao21 %-- &&&
%T12 = (12237.7)/(R2*T); %Tao12 %-- &&&
 
 
G12=exp(-alf*T12); 
G21=exp(-alf*T21);
for i= 1:13
    xb= 1 - x(i);
 
b = (G21/(x(i)+xb*G21));
d = (G12)/((xb+x(i)*G12)^2);
 
Ac1 = exp((xb^2)*(T21*(b^2)+T12*d));%coeficiente de actividad 1
 
ba =(G21/(x(i)+xb*G21)^2);
da = (G12)/(xb+x(i)*G12);
 
Ac2 = exp((x(i)^2)*(T12*(da^2) + T21*ba));%coeficiente de actividad 2
 
%presión de burbuja, Raoult modificado%
Pb= x(i)*Ac1*pas + xb*Ac2*pbs;
    
%fracción de vapor del agua%
ya= pas*Ac1*x(i)/Pb;
%fracción de vapor del ciclohexanol%
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
