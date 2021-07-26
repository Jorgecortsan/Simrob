clear all
clear %Limpiar ventana de comandos
clc   %Limpiar la pantalla
cla

clear a;
a=arduino('COM5','Uno');
servo1=servo(a,'D3');
servo2=servo(a,'D9');
servo3=servo(a,'D10');

% Posición deseada espacial (X,Y,Z)
n=100;

%Trayectoria lineal
X=linspace(-188,-130,n);
Y=linspace(0,130,n);
% Z=linspace(144,84,n); %Altura máxima: 144 Altura mínima: 84

Z=130+10*sin(linspace(-2*pi,2*pi,n))

%Trayectoria circular [X,Y]=[xc,yc]
% r=20;xc=0;yc=160;
% X=xc+r*cos(linspace(-pi,pi,n));
% Y=yc+r*sin(linspace(-pi,pi,n));
% Z=linspace(144,84,n); %%Altura máxima: 144 Altura mínima: 84

disp('Cinemática inversa de un robot configuración SCARA (3 GDL)')
disp(' ')

for trayect=1:n
px=X(trayect);
py=Y(trayect);
%cla
% Longitudes de los eslabones un robot de 3GDL
l1=98;
l2=90;
% Calcula el valor angular de theta 2
C2=(px^2+py^2-l1^2-l2^2)/(2*l1*l2);
S2=sqrt(1-C2^2);
t2=[atan2(S2,C2) atan2(-S2,C2)]*(180/pi); %S2 se evalúa con +/-
t1=[atan2((-(px*l2*S2)+((py)*(l1+(l2*C2)))),(((px)*(l1+(l2*C2)))+(py*l2*S2))) atan2((-(px*l2*-S2)+((py)*(l1+(l2*C2)))),(((px)*(l1+(l2*C2)))+(py*l2*-S2)))]*(180/pi);
Qd=[t1(:) t2(:)];
qd=Qd*(pi/180); % Coordenada generalizada "Ángulos/desplazamientos"

for j=1 %:size(qd,2)
    q1=qd(j,1);    % Coordenada generalizada 1 "theta 1"
    q2=qd(j,2);    % Coordenada generalizada 2 "theta 2"
    q3=Z(trayect); %Coordenada generalizada 3 "d3"
    lg=2; l=hsv;k=fix(64/4);s=200;view(170,135);
    % Construcción del sistema de referencia
    line([0 s],[0 0],[0 0],'Marker','.','Color','r','LineWidth',1);
    line([0 0],[0 s],[0 0],'Marker','.','Color','g','LineWidth',1);
    line([0 0],[0 0],[0 s],'Marker','.','Color','b','LineWidth',1);
    grid on;
    pox=0;poy=0;poz=0;
    % Construcción de los eslabones del robot de 3GDL "[l1 l2 l3]"
    linkx=[l1*cos(q1) l1*cos(q1)+l2*cos(q1+q2)];
    linky=[l1*sin(q1) l1*sin(q1)+l2*sin(q1+q2)];
    linkz=[110 190];
    % Gráfico de cada eslabón
    for i=1:lg
        line([pox linkx(i)],[poy linky(i)],[linkz(i) linkz(i)],'MarkerEdgeColor',[0 0 0],'Marker','o','Color',l(i*k,:),'LineWidth',4);
        pox=linkx(i);poy=linky(i);
    end
    line([linkx(1,1) linkx(1,1)],[linky(1,1) linky(1,1)],[linkz(1,1) linkz(1,2)],'MarkerEdgeColor',[0 0 0],'Marker','o','Color',l(3*k,:),'LineWidth',4);
    line([linkx(1,2) linkx(1,2)],[linky(1,2) linky(1,2)],[190 (q3)],'MarkerEdgeColor',[0 0 0],'Marker','o','Color','y','LineWidth',4);
    line([linkx(1,2) (((30*cos(q1+q2)))+(linkx(1,2)))],[linky(1,2) (((30*sin(q1+q2)))+(linky(1,2)))],[(q3) (q3)],'Marker','.','Color','r','LineWidth',2);
    line([linkx(1,2) (((30*cos(q1+q2-pi/2)))+(linkx(1,2)))],[linky(1,2) (((30*sin(q1+q2-pi/2)))+(linky(1,2)))],[(q3) (q3)],'Marker','.','Color','g','LineWidth',2);
    line([linkx(1,2) linkx(1,2)],[linky(1,2) linky(1,2)],[(q3) (q3-30)],'Marker','.','Color','b','LineWidth',2);
end
teta1=(((q1)*(180/pi)))/180;
teta2=(((q2)*(180/pi))+90)/180;
teta3=((180)+(((-180)/(60))*(144-(q3))))/(180); %al mandar un 1 al servo se queda arriba, al mandar un 0 se baja
writePosition(servo1,teta1);
writePosition(servo2,teta2);
writePosition(servo3,teta3);
pause(0.1)
end

disp('Coordenadas generalizadas para configuración codo abajo: ')
q1=Qd(1,1)
q2=Qd(1,2)
disp('Coordenadas generalizadas para configuración codo arriba: ')
q1=Qd(2,1)
q2=Qd(2,2)
disp('Posición en X,Y y Z de efector final')
X=px
Y=py
Z=q3