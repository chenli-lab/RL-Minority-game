#include <iostream>
#include <math.h>
#include <fstream>
#include <stdlib.h>
#include <time.h>
#include <queue>
#include <stack>
#include <iomanip>
#define PI 3.1415926

using namespace std;

//-- array creation and destruction --//
//1D
void make1Darray(int *& array,int rows)
{
    array=new int [rows];
}

void delete1Darray(int *array)
{
    delete []array;
    array=0;
}

void make1Darray(double *& array,int rows)
{
    array=new double [rows];
}

void delete1Darray(double *array)
{
    delete []array;
    array=0;
}

//2D
void make2Darray(int **& array,int rows,int cols)
{
    array=new int *[rows];
    for(int i=0;i<rows;i++) 
        array[i]=new int [cols];
}

void delete2Darray(int **array,int rows)
{
    for(int i=0;i<rows;i++)
        delete []array[i];
    delete []array;
    array=0;
}

void make2Darray(double **& array,int rows,int cols)
{
    array=new double *[rows];
    for(int i=0;i<rows;i++)  
        array[i]=new double [cols];
}

void delete2Darray(double **array,int rows)
{
    for(int i=0;i<rows;i++)
        delete []array[i];
    delete []array;
    array=0;
}

//3D
void make3Darray(int *** &array,int high,int row,int col)
{
    int i,j;
    array=new int **[high];
    for(i=0;i<high;i++)
        array[i]=new int*[row];
    
    for(i=0;i<high;i++)
        for(j=0;j<row;j++)
            array[i][j]=new int[col];
}

void delete3Darray(int ***array,int hight, int row)
{
    int i,j;
    for(i=0;i<hight;i++)
        for(j=0;j<row;j++)
            delete []array[i][j] ;
    
    for(i=0;i<hight;i++)
        delete []array[i] ;
    delete []array;
    array=0 ;
}

void make3Darray(double ***&array,int high,int row,int col)
{
    int i,j;
    array=new double **[high];
    for(i=0;i<high;i++)
        array[i]=new double* [row];
    
    for(i=0;i<high;i++)
        for(j=0;j<row;j++)
            array[i][j]=new double[col];
}

void delete3Darray(double ***array,int hight, int row)
{
    int i,j;
    for(i=0;i<hight;i++)
        for(j=0;j<row;j++)
            delete []array[i][j] ;
    for(i=0;i<hight;i++)
        delete []array[i] ;
    delete []array;
    array=0 ;
}

//--------------Define a function with arguments max------------//
double max(double x, double y)
{
	double z;
	z=x>y? x:y;
	return(z);
}


#define size 101							//Define the total number of people
#define tmax 55000                          // Maximum number of iterations  


int main()
{
	int i,j,k,m,realization;
	int success,u;
	double value,va;
    int it,iNN,N,dir;
    int num_state,num_action;
    double ***Qtable,*Reward,Q_max,Q_max_next;
    double molecule[2],pro[2],pro_new[2];
    double tau,denominator;
    int *Action,*State_new,*State,state_temp;
    int A,winning_side;
    double r,alpha,gamma;
    double epsilon,variance;
    double A_norm,capacity;
    int n_ave;
    double a,b,M;
    
    N=size;                                   // Set the number of agents
    num_action=2;                             // Number of possible actions: go to the bar or not go
    num_state=102;                            // All possible states: The number of people who went to the bar last time
    
  //  tau = 0.011;
    alpha = 0.1;							 // Learning rate
    gamma = 0.9;                             // Discount factor
    capacity=N/2.0;    						 // Maximum bar capacity defined as half of N
    
    srand48((long)time(NULL));
 //   srand((unsigned)time(NULL));
    
    make1Darray(State,N);                                   
    make1Darray(State_new,N);                             
    
    make1Darray(Action,N);                                  //action for each player
    make1Darray(Reward,N);                                  //the Reward obtained
    make3Darray(Qtable,N,num_state,num_action);             //Q table for each individual
    
    ofstream outfile("pdf_tauxiao.txt");
    if(!outfile){
        cout<<"You can't open the pointed fme!"<<endl;
        exit(1);
    }
    for(u=0;u<120;u++)
	{
		tau=(u+1)*0.0001;
	    
		success = 0;   
		for(realization=0;realization<100;realization++)
		{ 
		 //---------------------- System Initialization -----------------//
	        for(i=0;i<N;i++)
	        {
	            for(j=0;j<num_state;j++)
	            {
	                for(k=0;k<num_action;k++)
	                {
	                    Qtable[i][j][k] = 0;                  //all iterms in Q-table are set to be 0 (default)
	                }
	            }
	            Action[i] = rand()%2;                         //set to be 1 or 0 randomly
	            Reward[i] = 0;                                //No reward at the beginning
	        }
	    
	        A=0;
	        for(i=0;i<N;i++)
	        {
	            A += Action[i];	      //Add up the number of people who choose to go to the bar
	        }
	        
	    	state_temp = A;
	        for(i=0;i<N;i++)
	        {
	            State[i] = state_temp;        //The state of each individual is obtained
	        }
	        
	        variance=0;                               //measure the effeciency
	        n_ave=0;                                  //for doing average of variance 
	    //--------- The evolution of the game by Q-learning -----------//
	        for(it=0;it<tmax;it++)
	        {
	        //--------Choose an action------------//
	            for(i=0;i<N;i++)
	            {
	            	a = Qtable[i][State[i]][0]/tau;
		            b = Qtable[i][State[i]][1]/tau;
		            	
			       	M = max(a,b);
		            	
		            denominator  = 0;
		            for(k=0;k<num_action;k++)
					{
						molecule[k] = exp(Qtable[i][State[i]][k]/tau-M); 
						denominator  += exp(Qtable[i][State[i]][k]/tau-M);	
					} 
					
					for(k=0;k<num_action;k++)
					{
						pro[k] = molecule[k]/denominator;
					}
					
	                r = drand48();
	                if(r<pro[0])
	                {
	                	Action[i] = 0;
					}
					else
					{
						Action[i] = 1;
					}
	            }
	        //--------Perform action and measure the reward------------//
	            A=0;
	            for(i=0;i<N;i++)
	            {
	                A += Action[i];
	            }
	        
	            if(2*A > N)
	                winning_side = 0;           //the minority side is 0
	            else
	                winning_side = 1;           //otherwise is 1
	        
	            for(i=0;i<N;i++)
	            {
	                if(Action[i] == winning_side)
	                    Reward[i] = 1;
	                else
	                    Reward[i] = 0;
	            }
	        
	         	state_temp = A;
	            for(i=0;i<N;i++)
	            {
	                State_new[i] = state_temp;
	            }
	        
	            for(i=0;i<N;i++)
	            {
	           		a = Qtable[i][State_new[i]][0]/tau;
		            b = Qtable[i][State_new[i]][1]/tau;
		            	
		            M = max(a,b);
		        
		            denominator  = 0;
		            for(k=0;k<num_action;k++)
					{
						molecule[k] = exp(Qtable[i][State_new[i]][k]/tau-M); 
						denominator  += exp(Qtable[i][State_new[i]][k]/tau-M);	
					} 
					
					for(k=0;k<num_action;k++)
					{
						pro_new[k] = molecule[k]/denominator;
					} 
				//--------Update the Q table---------//	
	                Q_max_next = 0;
	                for(k=0;k<num_action;k++)                                   //go through all actions
	                {
	                	Q_max_next += pro_new[k]*Qtable[i][State_new[i]][k];
	                }
	            
	                Qtable[i][State[i]][Action[i]] = (1-alpha)*Qtable[i][State[i]][Action[i]] + alpha*(Reward[i] + gamma*Q_max_next);
	            }
	            
	            for(i=0;i<N;i++)
	            {
	                State[i]=State_new[i];
	            }
	            
	        //----------------The end of an evolutionary round------------//
	   			if(it > 50000)      // Calculate the variance at this temperature
	            {
	                variance += (A-capacity)*(A-capacity);
	                n_ave++;
	            }
	  			
	        }
	        variance = variance/(n_ave*N);
	        value = (51-capacity)*(51-capacity);
		    va = value/N;             // Calculate the variance coordinated to an optimal state
		    
		    if(variance == va)  //Compared with the variance of optimal coordination, if equal, success++
		    {
		    	success ++;
			}
    	}
    	outfile<<tau<<"   "<<double(success)/100<<"  "<<success<<endl;
	}
        
        delete1Darray(State);
	    delete1Darray(State_new);
	    
	    delete1Darray(Action);
	    delete1Darray(Reward);
	    delete3Darray(Qtable,N,num_state);
        
	return 0;
}