// ozisikradiation5.cpp : Este archivo contiene la funciÃƒÂ³n "main". La ejecuciÃƒÂ³n del programa comienza y termina ahÃƒÂ­.
//

#include <iostream>
#include <fstream>
using namespace std;
#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
#include<time.h>
#include<fcntl.h>
#include<math.h>






int main()


{   //variables
    int n = 200;
    int i;  //filas
    int j;  //columnas

    double L = 0.75; // bar length in meters
    double d = 0.002; // diameter
    double r = 0.001; // radio
    double e = 2.71828; 
    double pi= 3.141592654;
    double deltax = (L / n);
    double Tamb = 293.15; // Absolute temperature in Kelvin
    double Tcal = 373.15; // Absolute temperature in Kelvin
    double h = 1.; // film coefficient in W/(m^2 K)
    double kcon = 390.; // thermal conductivity W / (m K)
    double epsilon = 0.1; // emisividad del cobre
    double epsilonrapson = 1.0e-12; // valor minimo
    double fmax = 100.; // valor de error inicial
    double sigma = 5.67e-008; // Stephan Boltzmann constant in W/(m^ K^4)
    double x[n+1];
    double R[n+1];
    double P[n+1]; // perimeter in m
    double A[n+1]; // transversal area in m^2
    double Ap[n+1]; // derivative of area
    double T[n+2];
    double f[n+1];
    double Tpp[n+2];
    double Tp[n+2];
    double Tnueva[n+2];
    double deltaT[n+1];
    double jacobiano[n+1][n+1];
    double jacINV[n+1][n+1];
    double power[n];
    double Qcon[n];
    double Qconv[n+1];
    double Qrad[n+1];
    double sumconv;
    double sumrad;
    double sumcon;
    double Qtotal;
    double Largototal;
	int newra=0;
	int w=0;
	int jj;
	
//	fstream archivo2;
//	archivo2.open("Q en funcion a L", ios::out);			/// CALCULO DE Q EN FUNCION A L
//	Largototal=L/10;
	//for (jj=0;jj<=10;jj++){
		
//		L = Largototal*jj;
		
		
	for (w=0; w<1;w++)
	{
	//epsilon=1.0; // emissivity
	
	
	//para rutina inversion
	int Nmax=n+1;
	int  IAUXI [Nmax], IAUXJ[Nmax]; 
	int k;
	double hold,biga;
	double Aux[Nmax][Nmax]; 

    //funcion f

	for (i = 0; i <= n ; i ++)
        x[i] = 1.0*L*i/n; //+ (0.01-0.02)*i/n;

    for (i = 0; i <= n ; i ++)
        R[i] = 0.005; // radio fijo
        //R[i] = 0.02 + (0.001-0.02)*i/n; // radio variable lineal

	//perimetro variable 
	for (i=0; i <= n; i++)
		//P[i]= 2*pi*R[i];
		P[i]= 2*pi*(0.0175/2);
		
	//area variable
	for (i=0; i <= n; i++)
		//A[i] = pi*pow(R[i],2);
		A[i] = pi*(pow((0.0175/2),2)-pow((0.0155/2),2));
		
	//derivada del area
	for (i=1; i < n; i++)
		Ap[i] = (A[i+1]-A[i-1])/(2*deltax);		
		//Ap[n]=0;


	for (i = 0; i <= n+1 ; i ++)
        T[i] = Tcal - (Tcal - Tamb) * i * 0.9999/ n; //T supuesta lineal
        //Tnueva[i] = T[i];


	newra = 1;
	
//for (newra = 1; newra <= 20; newra ++)
	while (fmax>1.0e-7) // con fmax mÃ¡s chico a veces no converge
    {/////////////////////////////////////////////////////////////////////////////////////////////////////
	cout << "iteracion " << newra << " error maximo " << fmax << "\n";
 	fmax = 0.;
 	newra = newra+1;
    //Derivada central de primer orden
    for (i = 1; i <= n; i++)
        Tp[i] = (T[i + 1] - T[i - 1]) / (2 * deltax);
				//Tp[n]=0;
				
    //Derivada central de segundo orden
    for (i = 1; i <= n; i++)
        Tpp[i] = (T[i + 1] - 2 * T[i] + T[i - 1]) / (pow(deltax, 2));
        
    
								////REVISAR WHILE EPSILON
    for (i = 1; i <= n; i ++)
    //while(f[i-1]<epsilonrapson) {	
        {
		f[i-1] = h * P[i] * (T[i] - Tamb) + epsilon * sigma * P[i] * (pow(T[i], 4) -
							pow(Tamb, 4)) - kcon * (A[i] * Tpp[i] + Ap[i] * Tp[i]); //funcion f
	// el indice esta corrido !!!!
		if (fabs(f[i-1]) > fmax) {fmax = fabs(f[i-1]);}
		}
		
	// funcion que se debe anular para cumplir la condicion de contorno
    f[n] = h * (T[n] - Tamb) + epsilon * sigma * (pow(T[n], 4) - 
						pow(Tamb, 4)) + kcon * Tp[n];
						
	if (fabs(f[n]) >fmax) {fmax = fabs(f[n]);}
   
   //}
	//for (i = 0; i <= n; i++)
	//	cout << i << " " << P[i] << " " << A[i] << " " << "\n";
		
	//for (i = 0; i <= n+1; i++)
	//	cout << "Temp " << i << " " << T[i] << "\n";		
		
		
	//for (i = 1; i <= n; i++)
	//	cout << "Derivadas" << i << " " << " " << Tp[i] << " " << Tpp[i] << "\n";
				

	//for (i = 0; i <= n; i++)
	//	cout << i << " " << " f: " << f[i] << "\n";


  //matriz jacobiana de (N+1) x (N+1)
 //RECORDAR QUE LOS INDICES DE LA MATRIZ JACOBIANA VAN DE 0 A N, ES DECIR es de N+1 x N+1 
 
	for (i = 0; i <= n; i++)
	    for (j = 0; j <= n; j++)jacobiano[i][j]=0;
	for (i = 0; i <= n; i++)
	    for (j = 0; j <= n; j++)jacINV[i][j]=0;	
			
		// elementos fila 1 correspondiente a T_1 (nodo 1)
		jacobiano[0][0] = h * P[1] + epsilon * sigma * P[1] * 4 * pow(T[1], 3) - 
			kcon * (A[1] * (-2 / pow(deltax, 2)));
		jacobiano[0][1] = -kcon * (A[1] / pow(deltax, 2) + Ap[1] / (2 * deltax));
	
	
		// elementos filas 2 a N, (nodos 2 a N)
	    for (i = 2; i <= n; i++)
			{
		      jacobiano[i-1][i-1] = h * P[i] + epsilon * sigma * P[i] * 4 * pow(T[i], 3) - 
						kcon * (A[i] * (-2 / pow(deltax, 2)));
	        jacobiano[i-1][i - 2] = - kcon * (A[i] / pow(deltax, 2) + Ap[i]/(2 * deltax));
	        jacobiano[i-1][i] = - kcon * (A[i] / pow(deltax, 2) - Ap[i]/(2 * deltax));
		  }
	 
	 	// fila N+1 que corresponde a la condiciÃ¯Â¿Â½n de contorno
		jacobiano[n][n-2] = - kcon / (2 * deltax);
		jacobiano[n][n-1] = h + epsilon * sigma * 4 * pow(T[n], 3);
		jacobiano[n][n] = kcon / (2 * deltax);


/* invierto */

	for (i = 0; i < n+1; i++)
	    for (j = 0; j < n+1; j++)Aux[i][j]=jacobiano[i][j];
	
		
	for(k=0;k<Nmax;k++)
		{     
	                                                 
				biga=Aux[k][k];
				IAUXI[k]=k;
				IAUXJ[k]=k;
	                                                        
	         	for(j=k;j<Nmax;j++)
				{
	         	    for(i=k;i<Nmax;i++)
				    {
						if( abs(biga) < abs(Aux[i][j]) ) 
						{
					        biga=Aux[i][j];                                                      
	                IAUXI[k]=i;                                                       
	                IAUXJ[k]=j;                                                      
						}
					}
				}
	                                                     
				i=IAUXI[k];
				if(i > k)
				{
					for(j=0;j<Nmax;j++)
					{
	
						hold=-Aux[k][j];
						Aux[k][j]=Aux[i][j];
						Aux[i][j]=hold;
					}
				}
	            j=IAUXJ[k];
				if(j > k)
				{
					for(i=0;i<Nmax;i++)
					{
	
						hold=-Aux[i][k];
						Aux[i][k]=Aux[i][j];
						Aux[i][j]=hold;
					}
				}
				if(biga == 0.0)
				{
				    printf("matrix singular!!!\n"); 
					abort(); 	
				}
	                                                         
				for(i=0;i<Nmax;i++)
				{
				    if(i!=k)
						Aux[i][k]=Aux[i][k]/(-biga);
				}
	
				for(i=0;i<Nmax;i++)
				{
				    hold=Aux[i][k];
	     			for(j=0;j<Nmax;j++)
		    		{
						if(i!=k && j!=k)
							Aux[i][j]=hold*Aux[k][j]+Aux[i][j];
					}
				}
	
				for(j=0;j<Nmax;j++)
				{
					 if(j!=k)
						 Aux[k][j]=Aux[k][j]/biga;
				}
				Aux[k][k]=1.0/biga;
	    
		}
		
	
		for(k=Nmax-1;k>=0;k--)
		{
	
	                                               
	        i=IAUXI[k];
			if(i>k)
			{
				for(j=0;j<Nmax;j++)
				{
					hold=Aux[j][k];
					Aux[j][k]=-Aux[j][i];
					Aux[j][i]=hold;
				}
			}
	        j=IAUXJ[k];
			if(j>k)
			{
				for(i=0;i<Nmax;i++)
				{
					hold=Aux[k][i];
					Aux[k][i]=-Aux[j][i];
					Aux[j][i]=hold;
				}
			}
		}

//cargo el jacobiano invertido
	for (i = 0; i < Nmax; i++)
	    for (j = 0; j < Nmax; j++)jacINV[i][j]=Aux[i][j];

///*****************************/
//double sum=0;
//	for (i = 0; i <= n; i++) {
//		for (j = 0; j <= n; j++)
//		{
//			sum=0;
//			for (k = 0; k <= n; k++) sum=sum+jacINV[i][k]*jacobiano[k][j];
//				cout << i << " " << j << " " << sum<< "\n";
//		}
//	}
//
//// filas N y N+1 incluyen al nodo ficticio T_(N+1)


/*****************************/
//multiplicacion Jacobiano inverso por vecto f
	double sum=0;
	for (i = 0; i <= n; i++) {
		sum=0;
		for (k = 0; k <= n; k++) sum=sum+jacINV[i][k]*f[k];
				//cout << i << " " << j << " " << sum<< "\n";
			deltaT[i]=-sum;
			//Tnueva[i+1]=T[i+1]-sum;
		//	cout << "delta T" << i << " " << deltaT[i] << "\n";
	}
	
	for (i = 1; i <= n+1; i++)
		Tnueva[i]=T[i]+deltaT[i-1];      //Calculo del item 5)
	
	
	for (i = 1; i <= n+1; i++) 
		{
		T[i]=Tnueva[i];					//seteo variable para Newton-Rhapson
		cout << i << "Temp " << T[i] << "  " << " f " << f[i] << "\n";	
		}

	cout << " hola " << "\n";


}

///////////////////////////////////


	for (i = 0; i <= n+1; i++)
		cout << "T nueva" << i << " " << Tnueva[i] << "\n";

  
    fstream archivo;
    archivo.open("resultados.dat", ios::out);

    for (i =0; i <= n; i++)
        archivo << "nodo" << i << "Temperatura: " << T[i] << " Tp: " << Tp[i] << " Tpp: " << Tpp[i] << " f(i): " << f[i] << "\n";

	archivo.close();
	

	fstream archi;				/// GENERAR ARCHIVO
	archi.open("Tx.dat", ios::out);

    for (i =0; i <= n; i++)
        archi << x[i] << " " << T[i]-273.15 << "\n"; // en ªC

	archi.close();	
	
	
	
//    for (i = 0; i <= n; i++) cout << T[i] << "  " << "T nueva=  " << Tnueva[i]<< "   " << Tp[i] << "  " << Tpp[i] << " f("<< i<< " ):  "<< f[i] << " " << R[i] << " " << P[i] << " " << A[i] << " " << Ap[i] << "\n";
//	cout << "Se repitio unas"<< "   "  << newra <<"\n";
//    cout << "Hello World\n";

	/*for (i = 0; i <= n; i++) {
		for (j = 0; j <= n; j++) {
			cout << i << " " << j << " " << jacINV[i][j] << "\n";
		}
	}*/

    Qcon[n]=-kcon*A[n]*Tp[n];	//poder extraido de la fuente por conduccion
    for (i=1;i<n;i++) {
    	Qcon[i]=-kcon*A[i]*Tp[i];
    	}
    	sumcon=0;
    for (i=0;i<=n;i++) {
	
	sumcon=sumcon+Qcon[i];
	}
	
	Qcon[1]=-kcon*A[0]*(T[1]-T[0])/deltax;
	cout << "extraido desde la fuente  " << Qcon[1] << "flujo conduccion en  L  "  << Qcon[n] << "\n";
   
   
   
   
    for (i=0;i<=n;i++) {				//calcular el poder extrido mediante conveccion
	
    	Qconv[i]= P[i]*deltax*h*(T[i]-Tamb);
	    }		 
	Qconv[n+1]=A[n]*h*(T[n]-Tamb);  			//perdida por conveccion en extremo derecho
	
	
	sumconv=Qconv[0]/2.;
	for (i=1;i<=n-1;i++){
		sumconv=sumconv+Qconv[i];
		
	}
	sumconv=sumconv+Qconv[n]/2.+Qconv[n+1]; // incluyendo tapa
	//cout << "suma de lo convectivo  " << sumconv << "extraido por conveccion en L  " << Qconv[n+1] << "\n";
	
	
	
	for (i=0;i<=n;i++){								//Suma del poder extraido por radiacion
		Qrad[i]=P[i]*deltax*sigma*epsilon*(pow(T[i],4)-pow(Tamb,4));
		}
	Qrad[n+1]=A[n]*epsilon*sigma*(pow(T[n],4)-pow(Tamb,4)); 
	sumrad=Qrad[0]/2.;
	for (i=1;i<=n-1;i++){
		sumrad=sumrad+Qrad[i];
		
	}
	sumrad=sumrad+Qrad[n]/2.+Qrad[n+1]; // incluyendo tapa
	Qtotal=sumrad+sumconv;
	
	cout << "  convectivo  " << sumconv << "  radiativo  " << sumrad << "\n";
	
	//cout << "suma de lo radiativo  " << sumrad << "extraido por rad en L  " << Qrad[n+1] << "\n";
	
	
	
	cout << "total disipado " << Qtotal << " porcentaje radiativo/total " << sumrad*100/Qtotal <<"\n";
	
	
	
	
	
	}
	
 
	fstream file; ////	Q EN FUNCION A EPSILON
	file.open("Flujo en funcion a emisividad", ios::out);
	
	file << epsilon << "  " <<   sumconv << "  "<< sumrad <<"    "   << Qtotal <<   "\n";
	
	file.close();


	
			
	//archivo2 << L << "   " << Qtotal <<   "\n";


//	}

//	archivo2.close();

    return 0;
}
