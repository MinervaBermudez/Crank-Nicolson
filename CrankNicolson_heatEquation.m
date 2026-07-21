clc;
clear all;

% Minerva Bermúdez Ferrer
% Doble Grado en Física y Matemáticas

%-------------------------------------------------------------------------------
% Ejercicio 3. Usando el método de diferencias finitas de tipo Crank-Nicolson visto en
% clase, vamos a resolver numéricamente la ecuación del calor.

% El método de Crank-Nicolson consiste en calcular la iteración
% u_m+1 = u_m +k/2(Lu_m+1 + Lu_m) + g,
% con una matriz L y un vector g adecuados.
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% 1. Escribe un programa que calcule y represente la solución aproximada desde tiempo
% t = 0 hasta tiempo t = 5.
%-------------------------------------------------------------------------------

% Defino la función que me dará u(0,x):

  function u0 = u0_function(x)
      u0 = exp(-20. * (x - 0.5).^2);
  endfunction

% Defino las variables

  global L = 1;
  global T = 5;

  N = 50; % nodos espaciales
  M = 100; % nodos temporales

  h = L / N;       % paso espacial
  k = T / M;       % paso temporal
  x = (L * [0:N]) / N;  % N+1 puntos (de 0 a 1)
  t = (T * [0:M]) / M;  % M+1 puntos (de 0 a 5)

% Del enunciado se tiene que, para t=0:

  u0 = u0_function(x(2:N))';

% Relleno la matriz del método.

  clear A;
  A = zeros(N-1);

  for i = [1:N-1]
    A(i,i) = -2 / h^2; % Término diagonal

   if (i>1)
     A(i, i-1) = 1 / h^2; % Término subdiagonal
   endif

   if (i<N-1)
     A(i, i+1) = 1 / h^2; % Término superdiagonal
   endif

  endfor

  A(1,1)=-1/h^2;

% Definimos g

  g = zeros(N-1,1);
  g(1) = k/h;

% Resolver

 u = u0;
 u_save = u0;

% Matrices del sistema (Crank-Nicolson)

  I = eye(N-1);
  Izquierda = I - (k/2)*A;
  Derecha = I + (k/2)*A;

% Iteración temporal

for m = 1:M

    b = Derecha * u + g;
    u = Izquierda \ b;  % resolvemos sistema lineal

    u_save = [u_save, u];

    figure(1)
    plot(x, [u(1)+h;u;0]);
    axis([0, 1, 0, 1]);
    xlabel('t'); ylabel('x');
    title('Solución numérica de la ecuación del calor (Crank-Nicolson)');
    drawnow();

endfor

% Representación
figure(2)
[X, Tgrid] = meshgrid(x(2:N), t);  % malla para graficar
surf(X, Tgrid, u_save')          % transponemos u_save para surf
xlabel('x'); ylabel('t'); zlabel('u(x,t)');
title('Solución numérica de la ecuación del calor (Crank-Nicolson)');
shading interp;            % suaviza la superficie
colormap jet;                  % mapa de color bonito
colorbar                   % barra de colores
view(45, 45);              % vista 3D en ángulo


%-------------------------------------------------------------------------------
% 2. Calcula el único equilibrio posible de la EDP, y represéntalo gráficamente.
%-------------------------------------------------------------------------------

% El equilibrio de la EDP es tal que ∂^2u/∂x^2=0. Por tanto, u(t,x)=Ax+B.

% Definimos la función del equilibrio:

  function u = u0_function2(x)
    u = -x.+1;
  endfunction

% Introducimos los valores en un vector

  u_equilibrio = u0_function2(x);

% Representamos el equilibrio

  figure(3)
  plot(x, [u_equilibrio]);
  xlabel('x'); ylabel('Temperatura');
  title('Equilibrio térmico');
  grid on;

% Podemos observar que si el sistema se encuentra en el equilibrio,
% efectivamente no cambia en el tiempo:

  u0 = u_equilibrio(2:N)'
  u_saveeq = u0;

  for m = 1:M

    b = Derecha * u + g;
    u = Izquierda \ b;  % resolvemos sistema lineal

    u_saveeq = [u_saveeq, u];

  endfor

  % Representación
  figure(3)
  [X, Tgrid] = meshgrid(x(2:N), t);  % malla para graficar
  surf(X, Tgrid, u_saveeq')          % transponemos u_save para surf
  xlabel('x'); ylabel('t'); zlabel('u(x,t)');
  title('Solución numérica de la ecuación del calor (Crank-Nicolson)');
  shading interp;            % suaviza la superficie
  colormap jet;                  % mapa de color bonito
  colorbar                   % barra de colores
  view(45, 45);              % vista 3D en ángulo

%-------------------------------------------------------------------------------
% 3. Haz una gráfica de la distancia ℓ^2 entre el equilibrio y la solución, dependiendo
% del tiempo.
%-------------------------------------------------------------------------------

dist_eucl = zeros(M+1,1);

% La distancia entre el equilibrio y la solución para cada tiempo será

 for j=[1:M+1] % bucle de tiempo
    for i=[1:N-1] % bucle de espacio
      dist_eucl(j) += h^2*(u_equilibrio(i+1)-u_save(i,j))^2; % i+1 porque u_equilibrio incluye el borde
    endfor
    dist_eucl(j) = sqrt(dist_eucl(j));
  endfor

  figure(5)
  plot(t, dist_eucl);
  xlabel('t'); ylabel('Distancia');
  title('Distancia euclídea al equilibrio en función del tiempo');
  grid on;

  printf("Como podemos observar, la solución se acerca al equilibrio conforme pasa el tiempo.");


 % Puesto que no se especifica la norma a utilizar, veámoslo también con la norma infinito

  for j = [1:M+1]
  max_dif = 0;
  for i = [1:N-1]
    dif = abs(u_equilibrio(i+1) - u_save(i, j));  % i+1 porque u_equilibrio incluye el borde
    if dif > max_dif
      max_dif = dif;
    endif
  endfor
  dist_inf(j) = max_dif;
  endfor

  figure(6)
  plot(t, dist_inf);
  xlabel('t'); ylabel('Distancia');
  title('Distancia con la norma infinito');
  grid on;

  printf("Como podemos observar, la solución se acerca al equilibrio conforme pasa el tiempo.");



