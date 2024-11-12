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

/*
 * File: MingorityGame_for_typical_time_series.cpp
 * Description: This program provides foundational code to reproduce
 * the typical time series of the number of people going to the bar 
 * and the volatility with the temperature tau=0.1 from the article 
 * "Optimal coordination of resource: A solution from reinforcement learning."
 * modifying this value allows for the exploration of different results.
 * 
 * Author: [Guozhong Zheng]
 * Date: [2024.10.31]
 * 
 * Usage: Compile and run the program. The output will be written to 
 * one text file:
 * - time_series_tau=0.1.txt: Contains time series data of actions 
 *   and variance.
 */

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
    array=new int *[rows];//make a point for every row
    for(int i=0;i<rows;i++)  //make a point for every col
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

//   Define a function with arguments max------------//
double max(double x, double y)
{
	double z;
	z=x>y? x:y;
	return(z);
}

#define size 101                 //Define the total number of people  
#define tmax 100000              // Maximum number of iterations           

int main()
{
	int i,j,k;
    int it,N;
    int num_state,num_action;
    double ***Qtable,*Reward,Q_max,Q_max_next;
    double molecule[2],pro[2],pro_new[2];
    double pr[size],prb[size];
    double tau,denominator;
    int *Action,*State_new,*State,state_temp;
    int A,winning_side;
    double r,alpha,gamma;
    double epsilon,variance;
    double A_norm,capacity;
    int n_ave;
    double a,b,M;
    
    N=size;            // Set the number of agents                       
    num_action=2;      // Number of possible actions: go to the bar or not go                       
    num_state=102;     // All possible states: The number of people who went to the bar last time                       
    
    tau = 0.1;         // Temperature parameter for the softmax function
    alpha = 0.1;       // Learning rate
    gamma = 0.9;       // Discount factor
    capacity=N/2.0;    // Maximum bar capacity defined as half of N
    
    
    
    // Random seed initialization (optional)  
	srand48((long)time(NULL));
 	srand((unsigned)time(NULL));
    
    // Create arrays for state, action, reward, and Q-table
    make1Darray(State,N);                                  
    make1Darray(State_new,N);                                   
    
    make1Darray(Action,N);                                  
    make1Darray(Reward,N);                                 
    make3Darray(Qtable,N,num_state,num_action);             
    
    // Open output files for results
    ofstream outfile("time_series_tau=0.1.txt");
    if(!outfile){
        cout<<"You can't open the pointed fme!"<<endl;
        exit(1);
    }
	 
	 //---------------------- System Initialization -----------------//
        for(i=0;i<N;i++)
        {
            for(j=0;j<num_state;j++)
            {
                for(k=0;k<num_action;k++)
                {
                    Qtable[i][j][k] = 0; // Initialize Q-table entries to zero            
                }
            }
            Action[i] = rand()%2;     // Randomly initialize actions  0 is not go, 1 is go                 
            Reward[i] = 0;            // Initialize rewards to zero                  
        }
    
        A=0;     
        for(i=0;i<N;i++)
        {
            A += Action[i];   //Add up the number of people who choose to go to the bar
        }
        
    	state_temp = A;
        for(i=0;i<N;i++)
        {
            State[i] = state_temp;     //The state of each individual is obtained
        }
        
        //--------- The evolution of the game by Q-learning -----------//
        for(it=0;it<tmax;it++)
        {
        //--------Choose an action------------//
            for(i=0;i<N;i++)
            {
          	
          	//-------Prevents downward overflows of values----------//
            	a = Qtable[i][State[i]][0]/tau;
            	b = Qtable[i][State[i]][1]/tau;
            	
            	M = max(a,b);
            	
			//---------the softmax manner is employed to balance the trail-and-error exploration and the exploitation of Q-table---------//	
            	denominator  = 0;
            	for(k=0;k<num_action;k++)
				{
					molecule[k] = exp(Qtable[i][State[i]][k]/tau-M); 
					denominator  += exp(Qtable[i][State[i]][k]/tau-M);	
				} 
				for(k=0;k<num_action;k++)
				{
					pro[k] = molecule[k]/denominator; // the probability of the two actions being selected is calculated 
				} 
				
                r = drand48();  // Generate random number
                if(r<pro[0])
                {
                	Action[i] = 0;  // Choose action 0
				}
				else
				{
					Action[i] = 1;  // Choose action 1
				}
            }
        //--------Perform action and measure the reward------------//
            A=0;
            for(i=0;i<N;i++)
            {
                A += Action[i];  //Add up the number of people who choose to go to the bar
            }
        
            if(2*A > N)
                winning_side = 0;      // Side 0 wins     
            else
                winning_side = 1;      // Side 1 wins  
        
            for(i=0;i<N;i++)
            {
                if(Action[i] == winning_side)
                    Reward[i] = 1;    // Reward for winning action
                else
                    Reward[i] = 0;    // No reward for losing action
            }
        //--------Update the Q table---------//
         	state_temp = A;
            for(i=0;i<N;i++)
            {
                State_new[i] = state_temp;    // Update new state
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
					pro_new[k] = molecule[k]/denominator;    // Calculate new action probabilities
				} 
				
                Q_max_next = 0;
                for(k=0;k<num_action;k++)                                   //go through all actions
                {
                	Q_max_next += pro_new[k]*Qtable[i][State_new[i]][k];
                }
            	// Update Q-table using the Q-learning formula
                Qtable[i][State[i]][Action[i]] = (1-alpha)*Qtable[i][State[i]][Action[i]] + alpha*(Reward[i] + gamma*Q_max_next);
            }
        
            for(i=0;i<N;i++)
            {
                State[i]=State_new[i];  // Update states for next iteration
            }  
        //----------------The end of an evolutionary round------------//
        	variance = (A-capacity)*(A-capacity);
        //-----------------output----------------------//
  			outfile<<it<<"   "<<A<<"  "<<variance/N<<endl;
        }
        
        delete1Darray(State);
	    delete1Darray(State_new);
	    
	    delete1Darray(Action);
	    delete1Darray(Reward);
	    delete3Darray(Qtable,N,num_state);
        
	return 0;
}
