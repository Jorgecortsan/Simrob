clc     %Limpiar ventana de comandos
clear   %Limpiar pantalla
cla     %Limpiar gráfico
disp('Cinemática de Robot OTC Daihen FD-NV05')
disp(' ')
%%%syms a2%%%
dh=[0 525 270 185 1 0;
    0 0 90 0 1 0;
    0 760 270 0 1 0;
    270 0 270 150 1 0;
    0 750 90 0 1 0;
    90 0 270 0 1 0;
    0 115 0 0 1 0];  %Matriz D-H 

ap=dh(:,3)*(pi/180); a=dh(:,4); dis=dh(:,2); grados=dh(:,1)*(pi/180); 
a=sym(a); ap=sym(ap); dis=sym(dis); grados=sym(grados); %Volvemos simbólicos nuestros valores para eliminar el problema de la variable de precisión aritmética en las funciones seno y coseno
T=eye(4,4);   %Matriz identidad de inicio
for i=1:length(dh(:,1))  %Creación de variables simbólicas "t" y "d"
     eval([' t(',num2str(i),')=[sym(''t',num2str(i),''')];'])
     eval([' d(',num2str(i),')=[sym(''d',num2str(i),''')];']) 
end

for i=1:length(dh(:,1)) %Ciclo para la multiplicación de matrices desde i=1 hasta i=n
    if dh(i,5)==1 %¿Es variable de rotación? R: Sí
    t(1,i)=t(1,i)+grados(i,1); %Agregamos el offset (si existe) a la variable simbólica 
    elseif dh(i,5)==0 %¿Es variable de rotación? R: No
    t(1,i)=grados(i,1); %Sustituimos la variable simbólica por su valor numérico para que lo evalúe en la matriz de transformación homogénea como valor numérico y no como valor simbólico
    end
    if dh(i,6)==1 %¿Es variable de prismática? R: Sí
    d(1,i)=d(1,i)+dis(i,1); %Agregamos el offset (si existe) a la variable simbólica 
    elseif dh(i,6)==0 %¿Es variable de prismática? R: No
    d(1,i)=dis(i,1); %Sustituimos la variable simbólica por su valor numérico para que lo evalúe en la matriz de transformación homogénea como valor numérico y no como valor simbólico
    end
    eval(['A' num2str(i-1),num2str(i) ' = [cos(t(i)),-cos(ap(i))*sin(t(i)), sin(ap(i))*sin(t(i)),a(i)*cos(t(i));'...
                                          'sin(t(i)), cos(ap(i))*cos(t(i)),-sin(ap(i))*cos(t(i)),a(i)*sin(t(i));'...
                                          '0,sin(ap(i)),cos(ap(i)),d(i);0 0 0 1];']);  %Evaluamos la matriz de transformación homogénea para ese par cinemático
                            eval(['T' num2str(0),num2str(i) ' =T*A' num2str(i-1),num2str(i) ';' ]); %Guaradmos el valor de esa matriz en T que va desde 0 hasta i
                            eval(['T = T' num2str(0),num2str(i) ';']); %Guardamos el valor de la T anterior en esta variable para multiplicarla después
                            fix(T);
end
disp('Matriz de posición y orientación final:')
Tf=simplify(T) %Reducimos la matriz aplicando todas las funciones trigonométricas posibles
disp(' ')

Q=[29;-35;-35;91;91;61;-35]; % Condición inicial para Tetha 1 y Tetha 2
q=Q*(pi/180);  % Condición inicial en radianes
epsilon=1e-3; % Tolerancia de desviación
max_iter=100;  % Máximo número de iteraciones
n=10;
X=linspace(233.96,500,n);
Y=linspace(86.98,100,n);
Z=linspace(704.86,900,n);

Xd=[X;Y;Z];

for tr=1:size(Xd,2)
    xd=Xd(:,tr);
% Iteraciones del método de Newton
for i=1:max_iter
    t1=q(1);
    t2=q(2);
    t3=q(3);
    t4=q(4);
    t5=q(5);
    t6=q(6);
    t7=q(7);
    J=Jacobian2(t1,t2,t3,t4,t5,t6,t7);
    Jinversa=pinv(J);
    f=[185*cos(t1) + 760*cos(t1)*sin(t2) + 115*sin(t6)*(cos(t1)*sin(t2)*sin(t4) + cos(t4)*sin(t1)*sin(t3) - cos(t1)*cos(t2)*cos(t3)*cos(t4)) - 750*cos(t4)*(sin(t1)*sin(t3) - cos(t1)*cos(t2)*cos(t3)) - 150*sin(t4)*(sin(t1)*sin(t3) - cos(t1)*cos(t2)*cos(t3)) - 115*cos(t6)*(cos(t5)*(cos(t1)*cos(t4)*sin(t2) - sin(t1)*sin(t3)*sin(t4) + cos(t1)*cos(t2)*cos(t3)*sin(t4)) + sin(t5)*(cos(t3)*sin(t1) + cos(t1)*cos(t2)*sin(t3))) - 750*cos(t1)*sin(t2)*sin(t4) + 150*cos(t1)*cos(t4)*sin(t2);185*sin(t1) + 760*sin(t1)*sin(t2) - 115*sin(t6)*(cos(t1)*cos(t4)*sin(t3) - sin(t1)*sin(t2)*sin(t4) + cos(t2)*cos(t3)*cos(t4)*sin(t1)) + 750*cos(t4)*(cos(t1)*sin(t3) + cos(t2)*cos(t3)*sin(t1)) + 150*sin(t4)*(cos(t1)*sin(t3) + cos(t2)*cos(t3)*sin(t1)) - 115*cos(t6)*(cos(t5)*(cos(t4)*sin(t1)*sin(t2) + cos(t1)*sin(t3)*sin(t4) + cos(t2)*cos(t3)*sin(t1)*sin(t4)) - sin(t5)*(cos(t1)*cos(t3) - cos(t2)*sin(t1)*sin(t3))) + 150*cos(t4)*sin(t1)*sin(t2) - 750*sin(t1)*sin(t2)*sin(t4);760*cos(t2) + 150*cos(t2)*cos(t4) - 750*cos(t2)*sin(t4) - 150*cos(t3)*sin(t2)*sin(t4) + 115*cos(t2)*sin(t4)*sin(t6) - 750*cos(t3)*cos(t4)*sin(t2) - 115*cos(t2)*cos(t4)*cos(t5)*cos(t6) + 115*cos(t3)*cos(t4)*sin(t2)*sin(t6) + 115*cos(t6)*sin(t2)*sin(t3)*sin(t5) + 115*cos(t3)*cos(t5)*cos(t6)*sin(t2)*sin(t4) + 525]; % Modelo cinemático directo
    e=xd-f; % Error valor deseado y el valor actual numérico
    q=q+Jinversa*e; % Ecuación de Newton
    % Condición de término o criterio de paro
    if (norm(e)<epsilon)
        break;
    end
end

t1=q(1);
t2=q(2);
t3=q(3);
t4=q(4);
t5=q(5);
t6=q(6);
t7=q(7);
%disp('Matriz de transformación final con los ángulos propuestos:')
U=eval(Tf);
ET01=eval(T01);
ET02=eval(T02);
ET03=eval(T03);
ET04=eval(T04);
ET05=eval(T05);
ET06=eval(T06);
ET07=eval(T07);
grid on
view(135,135)
axis equal
T0=eye(4);
%cla
%line([T0(1,4) ],[linky(1,2) (((30*sin(q1+q2-pi/2)))+(linky(1,2)))],[(137-q3) (137-q3)],'Marker','.','Color','g','LineWidth',2);
line([T0(1,4) ET01(1,4)],[T0(2,4) ET01(2,4)],[525 525],'Color','r','linewidth',3)
line([T0(1,4) T0(1,4)],[T0(2,4) T0(2,4)],[0 525],'Color','k','linewidth',3)
%line([ET07(1,4) (ET07(1,4)+70)],[ET07(2,4) ET07(2,4)],[ET07(3,4) ET07(3,4)],'Marker','.','Color','r','LineWidth',2);
%line([ET07(1,4) ET07(1,4)],[ET07(2,4) (ET07(2,4)-70)],[ET07(3,4) ET07(3,4)],'Marker','.','Color','g','LineWidth',2);
%line([ET07(1,4) ET07(1,4)],[ET07(2,4) ET07(2,4)],[ET07(3,4) (ET07(3,4)-70)],'Marker','.','Color','b','LineWidth',2);

line([0 500],[0 0],[0 0],'Marker','.','Color','r','LineWidth',1);%%%%%sistema base
line([0 0],[0 500],[0 0],'Marker','.','Color','g','LineWidth',1);
line([0 0],[0 0],[0 700],'Marker','.','Color','b','LineWidth',1);

line([ET01(1,4) ET01(1,4)],[ET01(2,4) ET01(2,4)],[525 ET01(3,4)],'Color','k','linewidth',3)
line([ET01(1,4) ET02(1,4)],[ET01(2,4) ET02(2,4)],[ET01(3,4) ET02(3,4)],'Color','k','linewidth',3)
line([ET02(1,4) ET03(1,4)],[ET02(2,4) ET03(2,4)],[ET02(3,4) ET03(3,4)],'Color','b','linewidth',3)
line([ET03(1,4) ET04(1,4)],[ET03(2,4) ET04(2,4)],[ET03(3,4) ET04(3,4)],'Color','g','linewidth',3)
line([ET04(1,4) ET05(1,4)],[ET04(2,4) ET05(2,4)],[ET04(3,4) ET05(3,4)],'Color','y','linewidth',3)
line([ET05(1,4) ET06(1,4)],[ET05(2,4) ET06(2,4)],[ET05(3,4) ET06(3,4)],'Color','g','linewidth',3)
line([ET06(1,4) ET07(1,4)],[ET06(2,4) ET07(2,4)],[ET06(3,4) ET07(3,4)],'Color','r','linewidth',3)
drawnow

end
