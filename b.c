
#include <stdio.h>
#include <omp.h> //To use MPI
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
	outfile1 = fopen("serial2.txt", "w"); //open or create a file in current directory

	for(h = 0;h < howmanytimesavg; h++ )  //loop how many times you want to avg over
	{
		//For time
		double start=omp_get_wtime(); //start the timer
											
							  
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
		double finish=omp_get_wtime();//stop timer
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
void parallel(int bigN,int n)
{
	int no_thread=n;
	int x = powerOfTwo(bigN);
    bigN = 1 << x;

    double table[bigN][3];
	double avgtime = 0;
	int h;
	FILE *outfile;
	outfile = fopen("parallel_openmp.txt", "w"); 

	for(h = 0;h < howmanytimesavg; h++ )  
	{
		double start,finish; 
		start=omp_get_wtime(); 
											
							  
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
		
		#pragma omp parallel num_threads(no_thread)
		#pragma omp for schedule(guided,10)
		for (k = 0; k < bigN / 2; k++ ) 
		{	
			
			double sumrealeven = 0.0; 
			double sumimageven = 0.0; 
			double sumrealodd = 0.0; 
			double sumimagodd = 0.0; 
			
			#pragma task for schedule(guided,10)
			for (i = 0; i <= (bigN/2 - 1); i++) 
			{
				/* -------- EVEN PART -------- */
				double realeven = table[2*i][1]; 
				double complex imaginaryeven = table[2*i][2]; 
				double complex componeeven = (realeven + imaginaryeven * I); 

				double factoreven = ((2*PI)*((2*i)*k))/bigN; 
				
				double complex comptwoeven = (cos(factoreven) - (sin(factoreven)*I)); 

				evenpart[i] = (componeeven * comptwoeven); 
				
				/* -------- ODD PART -------- */
				double realodd = table[2*i + 1][1]; 
				double complex imaginaryodd = table[2*i + 1][2]; 
				double complex componeodd = (realodd + imaginaryodd * I); 
				
				double factorodd = ((2*PI)*((2*i+1)*k))/bigN;
															
				double complex comptwoodd = (cos(factorodd) - (sin(factorodd)*I));

				oddpart[i] = (componeodd * comptwoodd); 
				
			}
			
			#pragma task for schedule(guided,10)
			for(i = 0; i < bigN/2; i++) 
			{
				sumrealeven += creal(evenpart[i]); 
				sumimageven += cimag(evenpart[i]); 
				
				sumrealodd += creal(oddpart[i]); 
				sumimagodd += cimag(oddpart[i]); 
			}
			
			storeKsumreal[k] = sumrealeven + sumrealodd; 
			storeKsumimag[k] = sumimageven + sumimagodd; 
			
			storeKsumreal[k + bigN/2] = sumrealeven - sumrealodd; 
			storeKsumimag[k + bigN/2] = sumimageven - sumimagodd; 

			if(k <= 10) 
			{
				if(k == 0)
				{
					fprintf(outfile," \n\nTOTAL PROCESSED SAMPLES : %d\n\n\n",bigN);
				}
				fprintf(outfile,"================================\n");
				fprintf(outfile,"XR[%d]: %.4f XI[%d]: %.4f \n",k,storeKsumreal[k],k,storeKsumimag[k]);
				fprintf(outfile,"================================\n");
			}
		}
		finish=omp_get_wtime(); 
		double timeElapsed = finish-start; 
		avgtime = avgtime + timeElapsed; 
		fprintf(outfile,"Time taken on Iteration %d: %f Seconds\n", (h+1),timeElapsed);
	}

	avgtime = avgtime / howmanytimesavg;
	printf("\n %d point DFT calculated.\n",bigN);
	fprintf(outfile,"\nAverage Time taken by parallel openmp: %f Seconds", avgtime);
	printf("\nAverage Time taken by parallel openmp:%f seconds\n",avgtime );
	fclose(outfile);
	double start_time = omp_get_wtime();
       serial(bigN);
       double time_serial = omp_get_wtime() - start_time;
		FILE *fp;
       char *filename = "time2.txt";
       fp = fopen(filename, "a");
       fprintf(fp,"%d\n",bigN);
       fprintf(fp, "%f\n", avgtime);
       fprintf(fp, "%f\n", time_serial);
       fclose(fp); 
} 
int main(int argc,char** argv)
{
	int bigN=atoi(argv[1]);
	parallel(bigN,4);

	
	return 0;
}
