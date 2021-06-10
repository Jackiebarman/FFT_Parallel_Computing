
#include <stdio.h>
#include <mpi.h> //To use MPI
#include <complex.h> //to use complex numbers
#include <math.h>	//for cos() and sin()
#include<stdlib.h> //to use rand()

#define PI 3.14159265
#define howmanytimesavg 3 //How many times to do avg

int powerOfTwo(int n) 
{
   int i = 0, m = 1;
   while (m < n) 
   {
     m <<= 1;
     ++i;
   }
   return i;
}
void serial(int bigN)
{
	int x = powerOfTwo(bigN);
    bigN = 1 << x;

    double table[bigN][3];
	double avgtime = 0;
	int h;
	FILE *outfile1;
	outfile1 = fopen("serial1.txt", "w"); //open or create a file in current directory

	for(h = 0;h < howmanytimesavg; h++ )  //loop how many times you want to avg over
	{
		//For time
		double start=MPI_Wtime(); //start the timer
											
							  
		double complex evenpart[bigN/2]; //array to save the data for EVENHALF
		double complex oddpart[bigN/2]; //array to save the data for ODDHALF
		double storeKsumreal[bigN]; //store the K real variable so we can abuse symmerty
		double storeKsumimag[bigN]; //store the K imaginary variable so we can abuse symmerty
		int k, i ,j; 
		  
		
		for(i = 0; i < bigN; i++)
		{
			table[i][0] = i;
			for(j = 1; j < 3;j++)
			{
				table[i][j] = rand()%10; 
			}
		}
		

		for (k = 0; k < bigN / 2; k++ ) //K loop
		{	
			/* Variables used for the computation */
			double sumrealeven = 0.0; //sum of real numbers for even
			double sumimageven = 0.0; //sum of imaginary numbers for even
			double sumrealodd = 0.0; //sum of real numbers for odd
			double sumimagodd = 0.0; //sum of imaginary numbers for odd
			
			for (i = 0; i <= (bigN/2 - 1); i++) //loop for series 0->N/2 -1
			{
				/* -------- EVEN PART -------- */
				double realeven = table[2*i][1]; //Access table for real number at spot 2i
				double complex imaginaryeven = table[2*i][2]; //Access table for imaginary number at spot 2i
				double complex componeeven = (realeven + imaginaryeven * I); //Create the first component from table

				double factoreven = ((2*PI)*((2*i)*k))/bigN; //Calculates the even factor for Cos() and Sin()
				
				double complex comptwoeven = (cos(factoreven) - (sin(factoreven)*I)); //Create the second component

				evenpart[i] = (componeeven * comptwoeven); //store in the evenpart array
				
				/* -------- ODD PART -------- */
				double realodd = table[2*i + 1][1]; //Access table for real number at spot 2i+1
				double complex imaginaryodd = table[2*i + 1][2]; //Access table for imaginary number at spot 2i+1
				double complex componeodd = (realodd + imaginaryodd * I); //Create the first component from table
				
				double factorodd = ((2*PI)*((2*i+1)*k))/bigN;//Calculates the odd factor for Cos() and Sin()
															
				double complex comptwoodd = (cos(factorodd) - (sin(factorodd)*I));//Create the second component

				oddpart[i] = (componeodd * comptwoodd); //store in the oddpart array
				
			}
			
			for(i = 0; i < bigN/2; i++) //loop to sum the EVEN and ODD parts
			{
				sumrealeven += creal(evenpart[i]); //sums the realpart of the even half
				sumimageven += cimag(evenpart[i]); //sums the imaginarypart of the even half
				
				sumrealodd += creal(oddpart[i]); //sums the realpart of the odd half
				sumimagodd += cimag(oddpart[i]); //sums the imaginary part of the odd half
			}
			
			storeKsumreal[k] = sumrealeven + sumrealodd; //add the calculated reals from even and odd
			storeKsumimag[k] = sumimageven + sumimagodd; //add the calculated imaginary from even and odd
			
			storeKsumreal[k + bigN/2] = sumrealeven - sumrealodd; //ABUSE symmetry Xkreal + N/2 = Evenk - OddK
			storeKsumimag[k + bigN/2] = sumimageven - sumimagodd; //ABUSE symmetry Xkimag + N/2 = Evenk - OddK

			if(k <= 10) //for demonstartion print the first 10 values
			{
				if(k == 0)
				{
					fprintf(outfile1," \n\nTOTAL PROCESSED SAMPLES : %d\n\n\n",bigN);
				}
				fprintf(outfile1,"================================\n");
				fprintf(outfile1,"XR[%d]: %.4f XI[%d]: %.4f \n",k,storeKsumreal[k],k,storeKsumimag[k]);
				fprintf(outfile1,"================================\n");
			}
		}
		double finish=MPI_Wtime();//stop timer
		double timeElapsed = finish-start; //Time for that iteration
		avgtime = avgtime + timeElapsed; //AVG the time 
		fprintf(outfile1,"Time taken on Iteration %d: %f Seconds\n", (h+1),timeElapsed);
	}

	avgtime = avgtime / howmanytimesavg;
	printf("\n %d point DFT calculated.\n",bigN);
	fprintf(outfile1,"\nAverage Time taken by serial: %f Seconds", avgtime);
	printf("\nAverage Time taken by serial:%f seconds\n",avgtime );
	fclose(outfile1); //close file
}
void parallel(int argc,char** argv,int bigN)
{
	int x = powerOfTwo(bigN);
    bigN = 1 << x;

    double table[bigN][3];

	int my_rank,comm_sz;

	MPI_Init(&argc,&argv); //start MPI
	MPI_Comm_size(MPI_COMM_WORLD,&comm_sz);  
	MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);   

	double start,finish;
	double avgtime = 0;

	FILE *outfile;

	int h;

	if(my_rank == 0) 
	{
		outfile = fopen("parallel_mpi.txt", "w"); 
	}
	for(h = 0; h < howmanytimesavg; h++) 
	{
		if(my_rank == 0) 
		{	
			start = MPI_Wtime();
		}
		int i,k,n,j; 

		double complex evenpart[(bigN / comm_sz / 2)]; //array to save the data for EVENHALF
		double complex oddpart[(bigN / comm_sz / 2)]; //array to save the data for ODDHALF
		double complex evenpartmaster[ (bigN / comm_sz / 2) * comm_sz]; //array to save the data for EVENHALF
		double complex oddpartmaster[ (bigN / comm_sz / 2) * comm_sz]; //array to save the data for ODDHALF
		double storeKsumreal[bigN]; //store the K real variable so we can abuse symmerty
		double storeKsumimag[bigN]; //store the K imaginary variable so we can abuse symmerty
		
		double subtable[(bigN / comm_sz)][3]; //Each process owns a subtable from the table below 
		
		
			
		for(i = 0; i < bigN; i++)
			{
				table[i][0] = i;
				for(j = 1; j < 3;j++)
				{
					table[i][j] = rand()%10; 
				}
			}
			
		int sendandrecvct = (bigN / comm_sz) * 3; //count

		MPI_Scatter(table,sendandrecvct,MPI_DOUBLE,subtable,sendandrecvct,MPI_DOUBLE,0,MPI_COMM_WORLD); //scatter the table to subtables

		for (k = 0; k < bigN / 2; k++) 
		{
					
			double sumrealeven = 0.0; 
			double sumimageven = 0.0; 
			double sumrealodd = 0.0; 
			double sumimagodd = 0.0; 
			
			for(i = 0; i < (bigN/comm_sz)/2; i++) 
			{
				double factoreven , factorodd = 0.0;
				int shiftevenonnonzeroP = my_rank * subtable[2*i][0]; 
				int shiftoddonnonzeroP = my_rank * subtable[2*i + 1][0]; 
				
				/* -------- EVEN PART -------- */
				double realeven = subtable[2*i][1]; 
				double complex imaginaryeven = subtable[2*i][2];
				double complex componeeven = (realeven + imaginaryeven * I); 
				if(my_rank == 0) 
				{
					factoreven = ((2*PI)*((2*i)*k))/bigN; 										
								
				}
				else 
				{
					factoreven = ((2*PI)*((shiftevenonnonzeroP)*k))/bigN; 
								
				}
				double complex comptwoeven = (cos(factoreven) - (sin(factoreven)*I)); 
				evenpart[i] = (componeeven * comptwoeven); 
				
				/* -------- ODD PART -------- */
				double realodd = subtable[2*i + 1][1]; 
				double complex imaginaryodd = subtable[2*i + 1][2]; 
				double complex componeodd = (realodd + imaginaryodd * I);
				if (my_rank == 0)
				{
					factorodd = ((2*PI)*((2*i+1)*k))/bigN;										
								
				}
				else 
				{
					factorodd = ((2*PI)*((shiftoddonnonzeroP)*k))/bigN;
								
				}
							
				double complex comptwoodd = (cos(factorodd) - (sin(factorodd)*I));

				oddpart[i] = (componeodd * comptwoodd); 
				
			}

			/*Process ZERO gathers the even and odd part arrays and creates a evenpartmaster and oddpartmaster array*/
			MPI_Gather(evenpart,(bigN / comm_sz / 2),MPI_DOUBLE_COMPLEX,evenpartmaster,(bigN / comm_sz / 2), MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD);
			MPI_Gather(oddpart,(bigN / comm_sz / 2),MPI_DOUBLE_COMPLEX,oddpartmaster,(bigN / comm_sz / 2), MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD);

			if(my_rank == 0)
			{
				for(i = 0; i < (bigN / comm_sz / 2) * comm_sz; i++) 
				{
					sumrealeven += creal(evenpartmaster[i]); 
					sumimageven += cimag(evenpartmaster[i]); 
					sumrealodd += creal(oddpartmaster[i]); 
					sumimagodd += cimag(oddpartmaster[i]); 
				}
				storeKsumreal[k] = sumrealeven + sumrealodd; 
				storeKsumimag[k]  = sumimageven + sumimagodd;
				storeKsumreal[k + bigN/2] = sumrealeven - sumrealodd;
				storeKsumimag[k + bigN/2] = sumimageven - sumimagodd;
				if(k <= 10) 
				{
					if(k == 0)
					{
						fprintf(outfile," \n\n TOTAL PROCESSED SAMPLES : %d\n",bigN);
					}
					fprintf(outfile,"================================\n");
					fprintf(outfile,"XR[%d]: %.4f XI[%d]: %.4f \n",k,storeKsumreal[k],k,storeKsumimag[k]);
					fprintf(outfile,"================================\n");
				}
			}
		}
		if(my_rank == 0)
		{
			finish=MPI_Wtime(); 
			double timeElapsed = finish-start; 
			avgtime = avgtime + timeElapsed; 
			fprintf(outfile,"Time taken on Iteration %d: %f Seconds\n", (h+1),timeElapsed);
		}
	}
	if(my_rank == 0)
	{
		avgtime = avgtime / howmanytimesavg; 
		printf("\n %d point DFT calculated.\n",bigN);
		fprintf(outfile,"\nAverage Time taken by parallel mpi: %f Seconds", avgtime);
		printf("\nAverage Time taken by parallel mpi:%f seconds\n",avgtime );
		fclose(outfile); 
		double start_time = MPI_Wtime();
       serial(bigN);
       double time_serial = MPI_Wtime() - start_time;
		FILE *fp;
       char *filename = "time1.txt";
       fp = fopen(filename, "a");
       fprintf(fp,"%d\n",bigN);
       fprintf(fp, "%f\n", avgtime);
       fprintf(fp, "%f\n", time_serial);
       fclose(fp);
	}
	MPI_Barrier(MPI_COMM_WORLD); 
	MPI_Finalize(); 
} 
int main(int argc,char** argv)
{
	int bigN=atoi(argv[1]);
	parallel(argc,argv,bigN);

	
	return 0;
}
