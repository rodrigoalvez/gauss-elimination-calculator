.data 
mensaje_fallo: .string "El sistema no coincide con la dimension ingresada\n"
matriz: .float 1 , 1 , 1,
		1 , 3 , 2,
		1 , 1 , 3	#Se define la matriz A como un vector de floats dados por fila
fin_matriz: .space 4		#Espacio nulo que denota el fin de la matriz
vectorB: .float 1, 1, 1		#Vector B de terminos independientes dado como vector columna
fin_vector: .space 4		#Espacio nulo como fin del vector
dimension: .word 3 		#Dimension de la matriz A
resultado:			#Espacio para almacenar el vector solucion X

.text

#Primero preparamos los argumentos para invocar la subrutina
	la a0, vectorB		
	la a1, dimension	
	la a2,matriz	
	la a3, fin_vector
	la a4,fin_matriz
	la a6, resultado
	jal gauss

	
#Se muestran los resultados del método de Gauss	
	beq t2,zero,fin		
	la a6, resultado	
	lw t1, (a1)			
	add t0, zero, zero
loop_resultado:		
	li  a7, 2          	
    	flw fa0, (a6)  		
    	ecall
	
	li  a7, 11          	
    	addi a0, zero, '\n'  	
    	ecall
    	
    	addi t0, t0, 1
    	addi a6, a6, 4
	bne t0,t1, loop_resultado
	
	j fin
	
	
#------------SubRutinas------------------------------------------------------------------------#
#-----------------------Gauss------------------------------------------------------------------#
gauss:
#------------VERIFICACION DE DIMENSION---------------------------------------------------------#
		add t5,zero,a0			#Guardamos en t5 y t4 las direcciones inicialies de la matriz y del vector
		add t4,zero,a2			
		
		lw t0,(a1)			#Cargo en t0 la dimension de A
		add t1,zero,zero		#Inicializo el contador de elementos
		
		addi t3, zero, 4		#Define un salto de fila
		mul t3,t3,t0
#Loop de verificación
do:		
		beq t1,t0,verifica_dim_final
		
		addi t1,t1,1			#Sino continuo: Sumo 1 al contador de elementos
		addi a0,a0,4			#Aumento a la dirección del siguiente elemento
		add a2,a2,t3
		j do				#Vuelvo al comienzo del Loop de verificacion
#Fin de Loop de Verificación

verifica_dim_final:
		add t2, zero, zero
		
		bne a0,a3, fin_verifica_dim	#Si la direccion es no igual a la del fin del vector salto por falso
		bne a2,a4, fin_verifica_dim
		
		addi t2, zero, 1
	
fin_verifica_dim:
		beq t2, zero, mensaje_error	#Si verifica_dim retorna 0, dim(A) != dim(B) y finaliza el programa	
		
#----------------------------------Reduccion Matricial------------------------------------------------------------#
	add a0,zero,t5		#Restablecemos los valores de a0 y a2
	add a2,zero,t4

	lw s0,(a1)     		#Cargar la dimension
	addi t5,zero,4		#Carga la cantidad de byte por palabra (util)
	
	add t0,zero,zero 	#Contador de loop_columna
	addi s1,s0,-1		#Corte de loop_columna y loop_operacion
	
loop_columnas:
		addi t0,t0,1		#Aumentamos contador loop_columna
		
		add t1,zero,zero	#Inicializamos contador loop_fila
		sub s2,s0,t0		#Corte de loop_fila
		add a3,a2,t3 		#Calculo de a3, primer elemento a modificar por fila
		
		loop_fila:
			addi t1,t1,1		#Aumentamos contador loop_fila
			
			flw ft0,(a2)		#Cargamos el valor flotante del pivote
			flw ft1,(a3)		#Cargamos el valor flotante del valor a reducir
			fdiv.s ft2,ft1,ft0	#Cargamos el alfa
			
			add a4,a2,zero		#Auxiliar de direccion de columnas
			add a5,a3,zero		#Auxiliar de direccion de filas
			addi t4,t0,-1		#Contador de loop_operacion
			sw zero,(a5)		#Ahorramos una operacion seteando en 0 el elemento debajo del pivot
			loop_operacion:
				addi a4,a4,4			#Aumentar la direccion de la pos del pivote por columna
				addi a5,a5,4			#Aumentar la direccion de la pos del elemento a modificar
				addi t4,t4,1			#Aumenta el contador del loop_operacion
				
				flw ft0,(a4)			#Cargar los valores de las filas a operar
				flw ft1,(a5)			#
				fmul.s ft0,ft0,ft2		#Calculo Y -= alfa*X
				fsub.s ft1,ft1,ft0		#
				fsw ft1,(a5)
				
				bne t4,s1,loop_operacion	#Fin de loop
				
				addi a4,t0,-1			#Calculo de las direcciones del vectorB
				mul a4,a4,t5			#
				add a4,a4,a0			#
				mul a5,t1,t5			#
				add a5,a5,a4			#
				
				flw ft0,(a4)			#Cargar los valores de vectorB
				flw ft1,(a5)			#
				fmul.s ft0,ft0,ft2		#Calculo Y -= alfa*X
				fsub.s ft1,ft1,ft0		#
				fsw ft1,(a5)
				
			
			add a3,a3,t3		#Salto de fila
			bne t1,s2,loop_fila	#Fin de loop
		#fin_loop_fila
		
		add a2,a2,t3			#Salto de fila
		addi a2,a2,4			#Salto de columna
		bne t0,s1,loop_columnas		#Fin de loop
		
		
		
#--------------------Sustitucion hacia atras----------------------

		add t0,zero,zero		#Reiniciamos el valor de t0 para usarlo
		mul t2, s1, t5 			#Salto a la ultima posicion del vectorX (n-1)*4
		add a6, a6, t2 			#Salto a la ultima posicion del vectorX (n-1)*4
loop_fila_sust:	
		add t1,zero,zero		#Reiniciamos el valor de t1 para usarlo
		fcvt.s.w ft2,t1			#Inicializamos ft2 en 0
		flw ft0,(a5)			#Guarda el valor del vectorB
		
		loop_operacion_sust:
				
				beq t1,t0,fin_loop_operacion_sust	#Fin del loop
				
				
				flw ft3,(a6)				#Cargar el X
				flw ft1,(a2)				#Cargamos el elem
				fmul.s ft1,ft1,ft3			#Elem * X
				fadd.s ft2,ft2,ft1			#Acumulador
				
				addi a6,a6,-4				#Uno atras
				addi a2,a2,-4
				addi t1,t1,1
				j loop_operacion_sust
		
		
fin_loop_operacion_sust:
		fsub.s ft0,ft0,ft2		#Restar al valor la acumulacion de restas
		flw ft1,(a2)			#Cargo el valor del pivot
		fdiv.s ft0,ft0,ft1		#Calculo Xn
    		
		fsw ft0,(a6)			#Almacenamos en a6
		
		addi a5,a5,-4			#Salto al elem superior del vectorB
		mul t2,t0,t5			#Salto al utlimo elem de la fila sup
		sub t2,t2,t3			#Salto al utlimo elem de la fila sup
		add a2,a2,t2			#Salto al utlimo elem de la fila sup
		
		mul t2,t5,t1			#Vovler a a6 al final
		add a6,a6,t2			#Vovler a a6 al final
		
		addi t0,t0,1			#aumentamos el contador
		bne t0,s0,loop_fila_sust
		
		j fin_gauss
mensaje_error:
	li  a7, 4
   	la a0,mensaje_fallo
    	ecall
    	ret
fin_gauss:	
		addi t2,zero,1			# Se realizaron las operaciones y se pueden mostrar los resultados
		ret
#---------------------------------------------------------------------------------------------	

fin:
	
		
