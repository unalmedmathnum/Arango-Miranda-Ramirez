% Metodo de las potencias
 function M = power_iteration(A, num_iterations, tol)
     
     format long % Muestra números con hasta 15 decimales para valores en punto flotante de doble precisión.
     zise=size(A, 2); % Dice el # de columbnas que tiene A.
     M= zeros(num_iterations+1,zise+2); %Crea una matriz de ceros de tamaño (num_iterations+1)x(zise+2), es esta matrix 
     % se guardará en la primera columna las aproximaciones de lambda. En
     % la segunda columna se guardara la tolerancias unicamente en la fila
     % donde se encuentra la proximacion de lambda que la diferencia con respecto a la anterior aprox está por debajo
     % de la tolerancia.
     % El resto de columnas estan designadas para guardar la aproximación
     % del vector propio en la i_esima iteración.
     
     b_k = rand(zise, 1); % Genera un vector aleatorio inicial.  

     % Actualiza la matriz M con la información de b_k.
     for i=1:zise
         M(1,i+2)=b_k(i,1);
     end
     % Acá se itera la matrix A en el vector b_k.
     for i = 1:num_iterations            
         b_k1 = A * b_k;   % Calcula el producto matriz-por-vector Ab. 
         lambda_k1=(( b_k1)'*b_k)/(b_k'*b_k); %Aproximas el eigenvalor dominate lambda.
         b_k1 = b_k1 / norm(b_k1); % Normaliza el vector.         
         M(i+1,1)=lambda_k1;% Actualiza la matriz M con la información de lambda_k1 en la i_esima iteración, sin embargo la información se guarda en la fila i+1 de M porque la primera fila guarda el vector inicial elegido.
         for j=1:zise % Actualiza la matriz M con la información de b_k1 en la i_esima iteración, en la fila i+1 de M.
         M(i+1,j+2)=b_k1(j,1);
         end 
         distancia = pdist2(b_k1', b_k', 'chebychev'); % calcula la norma infinito de b_k1-b_k
         if distancia<tol  % Examina si  distancia esta por debajo de tol   
            M(i+1,2)=tol; % Modifica la fila i+1 de la  matix M cuando la diferencia de la i_esima aproximación de lambda con la anterior aprox está por debajo de la tolerancia pedida.
             return
         end        
         % Renobra el vector b_k        
         b_k = b_k1;                           
     end     
 end
