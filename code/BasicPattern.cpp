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
 * File: MingorityGame_for_basic_pattern.cpp
 * Description: This file provides the basic code for spatio-temporal patterns.
 * Running the program can produce the results for Figure 4 in the article "Optimal coordination of resource: A solution from reinforcement learning."
 * One can modify the temperature and time window in this program to obtain other spatio-temporal patterns from the article.
 * Author: [Guozhong Zheng]
 * Date: [2024.11.3]
 * 
 * Usage: Compile and run the program. The output will be written to 
 * two text files:
 * - Pattern50C.txt and Pattern51C.txt. 
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
    array=new double *[rows];//make a point for every row
    for(int i=0;i<rows;i++)  //make a point for every col
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


#define size 101
#define tmax 10000                          // max nr. of time steps


int main()
{
	int i,j,k,m,e;
    int it,iNN,N,dir;
    int a[size][200],b[200][size];
    int c[size][200],d[200][size];
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
    
    N=size;                                   //the Population size
    num_action=2;                             //0 or +1
    num_state=102;                         
    
    tau = 0.1;
    alpha = 0.1;
    gamma = 0.9;
    capacity=N/2.0;
    
 //   srand48((long)time(NULL));
 //   srand((unsigned)time(NULL));
    
    make1Darray(State,N);                                   //the number of players with +1 state
    make1Darray(State_new,N);                                   //the number of players with +1 state
    
    make1Darray(Action,N);                                  //action for each player
    make1Darray(Reward,N);                                  //the Reward obtained
    make3Darray(Qtable,N,num_state,num_action);             //Q table for each individual
    
    ofstream outfile7("Pattern50C.txt");
    if(!outfile7){
        cout<<"You can't open the pointed fme!"<<endl;
        exit(1);
    }
     ofstream outfile8("Pattern51C.txt");
    if(!outfile8){
        cout<<"You can't open the pointed fme!"<<endl;
        exit(1);
    }

	    for(i=0;i<N;i++)
		{
			for(j=0;j<200;j++)
			{
				a[i][j]=0;
				c[i][j]=0;
			}
		}
		for(i=0;i<200;i++)
		{	
			for(j=0;j<N;j++)
			{
				b[i][j]=0;
				d[i][j]=0;
			}
		}
	 //---------------------- System Initialization -----------------//
        for(i=0;i<N;i++)
        {
            for(j=0;j<num_state;j++)
            {
                for(k=0;k<num_action;k++)
                {
                    Qtable[i][j][k] = drand48()*0.02-0.01;                  
                }
            }
            Action[i] = rand()%2;                         
            Reward[i] = 0;                                
        }
    
        A=0;
        for(i=0;i<N;i++)
        {
            A += Action[i];        //recording how many guys in the neighborhood went to the bar (i.e. with value 1)
        }
        
    	state_temp = A;
        for(i=0;i<N;i++)
        {
            State[i] = state_temp;                                     
        }
        
    	//--------- The evolution of the game by Q-learning -----------//
        for(it=0;it<tmax;it++)
        {
        //--------Choose an action------------//
            for(i=0;i<N;i++)
            {
            	denominator  = 0;
            	for(k=0;k<num_action;k++)
				{
					molecule[k] = exp(Qtable[i][State[i]][k]/tau); 
					denominator  += exp(Qtable[i][State[i]][k]/tau);	
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
			if(it>9799&&it<10000)
		    {
		    	if(A==50)
		    	{
		    		for(e=0;e<N;e++)
	        		{
						b[it-9800][e]=Action[e];
				
					}
				}
				else
				{
					for(e=0;e<N;e++)
	        		{
						b[it-9800][e]=2;
				
					}
				}
				
				if(A==51)
		    	{
		    		for(e=0;e<N;e++)
	        		{
						d[it-9800][e]=Action[e];
				
					}
				}
				else
				{
					for(e=0;e<N;e++)
	        		{
						d[it-9800][e]=2;
				
					}
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
            	denominator  = 0;
            	for(k=0;k<num_action;k++)
				{
					molecule[k] = exp(Qtable[i][State_new[i]][k]/tau); 
					denominator  += exp(Qtable[i][State_new[i]][k]/tau);	
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
        	variance = (A-capacity)*(A-capacity);
        }
        for(i=0;i<200;i++)
		{
			for(j=0;j<N;j++)
			{
				a[j][i]=b[i][j];
				c[j][i]=d[i][j];
			}
		}
		//-------------------output the results-------------------// 
		for(i=0;i<N;i++)
		{
			for(j=0;j<200;j++)
			{
				outfile7<<a[i][j]<<"  ";
				outfile8<<c[i][j]<<"  ";
			}
			outfile7<<endl;
			outfile8<<endl;
		}
        
        delete1Darray(State);
	    delete1Darray(State_new);
	    
	    delete1Darray(Action);
	    delete1Darray(Reward);
	    delete3Darray(Qtable,N,num_state);
        
	return 0;
}
