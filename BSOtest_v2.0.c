/*
 ============================================================================
 Name        : BSOtest.c
 Author      : 511
 Version     : 2.0
 Copyright   : www.DigVan.com
 Description : 
 			V1.0:Sequential implementation of Brain Storm Optimization Algorithm 
 			in C based on "Y. Shi, 'Brain storm optimization algorithm', ICSI 2011"
			V2.0:output a csv file show different results via different paraments
 ============================================================================
 */

#include<stdio.h>
#include<math.h>
#include<time.h>
#include<stdlib.h>

#define	PI 3.141592653589
#define File_PATH "file\\OptiResuDiffPara.csv"
#define File_title "Population_Size,Cluster_Number,Dimension,Noise_Range,Iteration,Value,Time\n"

#define P_Size 100          //population size
#define DIM 10         		  //solution dimension
#define Iteration 100000       //iteration number
#define Clu_Num 10           //cluster number
#define NoiseRange 0.02       //the range a noise will be generated
#define Range_L -10.0         //left range
#define Range_R 10.0          //right range
#define P_NewCenter 0.05    //probability to choose a cluster and replace the center
#define P_OneClu 0.85       //probability to generate new solution from one cluster
#define P_One_Center 0.60   //probability to genereta new solution from center of one cluster
#define P_Two_Center 0.60   //probability to generate new solution from center of two clusters

#define Nrand() (((double)rand()/(RAND_MAX+1.0))*NoiseRange)							//generate noise
#define Srand() (Range_L+(((double)rand()/(RAND_MAX+1.0))*(Range_R-Range_L)))	//generate new solution dimensionly
#define Prand() ((double)rand()/(RAND_MAX+1.0))									//generate probability
#define Crand() ((int)(rand()/(RAND_MAX+1.0)*Clu_Num))							//generate a cluster index
#define SIrand() ((int)(rand()/(RAND_MAX+1.0)*(P_Size/Clu_Num)))				//generate solution index

double initSolution[P_Size][DIM];        			//initial solution population
double solution[Clu_Num][P_Size/Clu_Num][DIM];   	//solution population after cluster and sorting

//evaluation functions: Rastrigin :global min = 0 in (0,0,0,0,0,...)
double Rastrigin(double x[]){
	int i;
	double y=0.0;
	for(i=0;i<DIM;i++){
		y+=x[i]*x[i]-10.0*cos(2*PI*x[i])+10.0;
	}
	return y;
}
//evaluation functions: Rosenbrock :global min = 0 in (1,1,1,1,1,1,...)
double Rosenbrock(double a[]){

	int i;
	double sum=0.0;
	for(i=0;i<DIM-1; i++){
        sum+= 100*(a[i+1]-a[i]*a[i])*(a[i+1]-a[i]*a[i])+(a[i]-1)*(a[i]-1);
	}
  	    return sum;
}
//evaluation functions: Sphere :global min = 0 in (0,0,0,0,0,...)
double Sphere(double a[]){

	int i;
	double sum=0.0;
	for(i=0; i<DIM; i++){
		sum+=a[i]*a[i];
	}
	return sum;
}
//call the evaluation function and return fitness value
double evaluation(double a[]){
	return Sphere(a);
	//  return Rastrigin(a);
	//  return Rosenbrock(a);
}
//initial solution
void init(){
	int i,j;
	// srand((unsigned)time(NULL));
	for(i=0; i<P_Size; i++)
		for(j=0; j<DIM; j++)
			initSolution[i][j] = Srand();
}
//cluster solution
void cluster(){
	int i,j,k,count;
	for(i = 0, count = 0;i < Clu_Num && count < P_Size; i++)
		for(j = 0;j < P_Size/Clu_Num;j++){
			for(k = 0;k < DIM;k++)
				solution[i][j][k] = initSolution[count][k];
			count++;
		}
}
//sorting solution:te best solution get lowest index in each cluster
void sorting(){
	int i,j,k,m;
	double temp[DIM];
	for(i = 0;i < Clu_Num;i++)
		for(j = 0; j < P_Size/Clu_Num;j++){
			for(k = j + 1; k < P_Size/Clu_Num;k++){
				if(evaluation(solution[i][k]) < evaluation(solution[i][j])){
					for(m = 0;m < DIM;m ++){
						temp[m] = solution[i][j][m];
						solution[i][j][m] = solution[i][k][m];
						solution[i][k][m] = temp[m];
					}
				}
			}
		}
}
//to show the best solution yet
double test(){
	sorting();
	int i,m;
	double temp[DIM];
	for(m = 0;m < DIM;m ++)
		temp[m] = solution[0][0][m];

	//find and store the best solution in all clusters yet
	for(i = 0;i < Clu_Num;i++)
		if(evaluation(solution[i][0]) < evaluation(temp))
			for(m = 0;m < DIM;m ++)
				temp[m] = solution[i][0][m];

	return evaluation(temp);

}
//generate new solution
void refresh(){
	/**
	 * cluIndex0X: randomly generated cluster index
	 * soluIndex0X: randomly generated solution index in a cluster
	 * m: count dimension for each solution
	 */
	int cluIndex,cluIndex01,cluIndex02,soluIndex,soluIndex01,soluIndex02,m;
	double newIndivi[DIM];

	//refresh center in a randomly cluster
	if(Prand() < P_NewCenter){
		//randomly choose a cluster
		cluIndex = Crand();
		//replace this cluster center
		for(m = 0;m < DIM;m ++)
			solution[cluIndex][0][m] = Srand();
	}

	//generate solutions from one cluster
	if(Prand() < P_OneClu){
		//randomly choose a cluster
		cluIndex = Crand();

		if(Prand() < P_One_Center){
			//add noise to center to generate new solution
			for(m = 0;m < DIM;m ++)
				newIndivi[m] = solution[cluIndex][0][m] + Nrand();
			//keep the better one
			if(evaluation(newIndivi) < evaluation(solution[cluIndex][0]))
				for(m = 0;m < DIM;m ++)
					solution[cluIndex][0][m] = newIndivi[m];
		}else{
			//randomly choose a solution in this cluster
			soluIndex = SIrand();
			//add noise to this solution to generate new solution
			for(m = 0;m < DIM;m ++)
				newIndivi[m] = solution[cluIndex][soluIndex][m] + Nrand();
			//keep the better one
			if(evaluation(newIndivi) < evaluation(solution[cluIndex][soluIndex]))
				for(m = 0;m < DIM;m ++)
					solution[cluIndex][soluIndex][m] = newIndivi[m];
		}
	}
	//generate solutions from two clusters
	else{
		//randomly choose two clusters
		cluIndex01 = Crand();
		cluIndex02 = Crand();

		if(Prand() < P_Two_Center){
			//generate new solution from two centers
			for(m = 0;m < DIM;m ++)
				newIndivi[m] = (solution[cluIndex01][0][m]+solution[cluIndex02][0][m])/2 + Nrand();
			//keep the better one
			if(evaluation(newIndivi) < evaluation(solution[cluIndex01][0]))
				for(m = 0;m < DIM;m ++)
					solution[cluIndex01][0][m] = newIndivi[m];
			if(evaluation(newIndivi) < evaluation(solution[cluIndex02][0]))
				for(m = 0;m < DIM;m ++)
					solution[cluIndex02][0][m] = newIndivi[m];
		}else{
			//randomly choose two solutions in clusters
			soluIndex01 = SIrand();
			soluIndex02 = SIrand();

			//generate new solution from these two solutions
			for(m = 0;m < DIM;m ++)
				newIndivi[m] = (solution[cluIndex01][soluIndex01][m]+solution[cluIndex02][soluIndex02][m])/2 + Nrand();
			//keep the better one
			if(evaluation(newIndivi) < evaluation(solution[cluIndex01][soluIndex01]))
				for(m = 0;m < DIM;m ++)
					solution[cluIndex01][soluIndex01][m] = newIndivi[m];
			if(evaluation(newIndivi) < evaluation(solution[cluIndex02][soluIndex02]))
				for(m = 0;m < DIM;m ++)
					solution[cluIndex02][soluIndex02][m] = newIndivi[m];
		}
	}
}
//main function to test
int main(){
	//to count the time
	clock_t start,finish;
 	start=clock();
 	float duration;

	int i,j;
	double bestYet;			//to store the best solution
	FILE *fp = NULL;		//pointer to output file

	init();

/*
	//test init function
	for(int i = 0;i < P_Size; i++){
		for(int j = 0;j < DIM;j++)
			printf("%f\t",initSolution[i][j]);
		printf("\n");
	}
	printf("\n this is initial solution matrix.\n\n");
*/

	cluster();

/*
	//test cluster function
	for(i = 0;i < Clu_Num; i++){
		for(j = 0;j < P_Size/Clu_Num;j++){
			for(k = 0;k < DIM;k++)
				printf("%f\t",solution[i][j][k]);
			printf("\n");
		}
		printf("\n this is a cluster.\n");
	}
	printf("\n this is clustered solution matrix.\n\n");
*/

	sorting();

/*
	//test sorting function
	for(i = 0;i < Clu_Num; i++){
		for(j = 0;j < P_Size/Clu_Num;j++){
			for(k = 0;k < DIM;k++)
				printf("%f\t",solution[i][j][k]);
			printf("\n");
		}
		printf("\n this is a cluster.\n");
	}
	printf("\n this is sorted solution matrix.\n\n");
*/

	bestYet = evaluation(solution[0][0]);
	for(i = 0; i < Iteration;i++){
		for(j = 0;j < P_Size;j++){
			refresh();
		}
		if(test() < bestYet)
			bestYet = test();
	}

	//finish counting
	finish=clock();
	duration=(float)(finish-start)/1000;//1000 in windows; 1000000 in linux
	fp = fopen(File_PATH, "a+");
	if(NULL == fp){
		return -1;
	}
//	fprintf(fp, File_title);  //write title,only run at the first time

	//write data
	fprintf(fp, "%d,",P_Size);
	fprintf(fp, "%d,",Clu_Num);
	fprintf(fp, "%d,",DIM);
	fprintf(fp, "%.4f,",NoiseRange);
	fprintf(fp, "%d,",Iteration);
	fprintf(fp, "%.6f,%f\n",bestYet,duration);

	fclose(fp);
	fp = NULL; //free the pointer
	return 0;
}
