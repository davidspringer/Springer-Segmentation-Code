/* Many people have requested a simple example on how to create a C
 * MEX-file.  In response to this request, the following C MEX-file,
 * named mexample, is provided as an introduction to cmex
 * programming. mexample is a commented program which describes how to
 * use the following MEX-functions:
 *
 * mexErrMsgTxt
 * mxCreateDoubleMatrix
 * mxGetM
 * mxGetN
 * mxGetPr
 * mxIsComplex
 * mxIsSparse
 * mxIsChar
 *
 * In MATLAB, mexample accepts two inputs and returns one output. The
 * inputs are a 2x2 array denoted as ARRAY_IN and a 2x1 vector denoted as
 * VECTOR_IN.  The function calculates the determinant of ARRAY_IN,
 * multiplies each element of VECTOR_IN by the determinant, and returns
 * this as the output, denoted by VECTOR_OUT.  All inputs and outputs to
 * this function are assumed to be real (not complex). */

/*   First, include some basic header files.  The header file
 * "mex.h" is required for a MEX-file.  Add any other header
 * files that your function may need here. */

#include "mex.h"
#include <limits.h>
#include <float.h>
#include <math.h>       /* log */
/*   A C MEX-file generally consists of two sections.  The first
 * section is a function or set of functions which performs
 * the actual mathematical calculation that the MEX-function
 * is to carry out.  In this example, the function is called
 * workFcn().  The second section is a gateway between MATLAB
 * and the first section, and consists of a function called
 * mexFunction.  The gateway is responsible for several tasks,
 * including:
 *
 * I)  error checking,
 * II)  allocating memory for return arguments,
 * III)  converting data from MATLAB into a format that
 * the workFcn function can use, and vice versa.
 *
 * The first function to be written in this example, then, is
 * workFcn:
 *
 * Since C and MATLAB handle two-dimensional arrays
 * differently, we will explicitly declare the dimension of
 * the variable theArray.  The variables, theVector and
 * theResult, are both one-dimensional arrays, and therefore
 * do not need such rigid typing. */


void viterbi(
        int N,
        int T,
        double a_matrix[4][4],
        int max_duration_D,
        double *delta,
        double *observation_probs,
        double duration_probs [4][150],
        double *psi,
        double *psi_duration_out
        )
        
{
    
    int i;
    int j;
    int t;
    int d;
    
    
    
    for (t = 1; t<T+ max_duration_D-1;t++){
        
        /*For each state */
        for (j = 0; j<N;j++){
            double emission_probs = 0;
            
            /*        max_duration_D*/
            for (d = 1; d<max_duration_D; d++){
                int start; int max_index = 0;
                int end_t = 0;
                int duration = 0;
                double probs = 0;
                double duration_sum = 0;
                double delta_temp = 0;
                double max_delta = -1*DBL_MAX ;
                
                
                
                /*  Get the maximum value for delta at this t, and record the state where it was found:
                 * This is the first half of the expression of equation 33a from Rabiner:*/
                
                /*
                 * %The start of the analysis window, which is the current time
                 * %step, minus d (the time horizon we are currently looking back),
                 * %plus 1. The analysis window can be seen to be starting one
                 * %step back each time the variable d is increased.
                 * % This is clamped to 1 if extending past the start of the
                 * % record, and T-1 is extending past the end of the record:
                 */
                
                start = t - d;
                
                if(start < 0){
                    start = 0;
                }
                
                if(start > T-1){
                    start = T-1;
                }
                
                /*
                 * %The end of the analysis window, which is the current time
                 * %step, unless the time has gone past T, the end of the record, in
                 * %which case it is truncated to T. This allows the analysis
                 * %window to extend past the end of the record, so that the
                 * %timing durations of the states do not have to "end" at the end
                 * %of the record.
                 */
                
                end_t = t;
                if(end_t>T){
                    end_t = T;
                }
                
                duration = end_t -start;
                
                
                
                
                for(i = 0; i<N; i++)
                {
                    double temp = delta[t-(d) +i*T] + log(a_matrix[i][j]);
                    
                    
                    
                    if(temp > max_delta){
                        max_delta = temp;
                        max_index = i;
                    }
                }
                
                
                
                
                /*//Find the normaliser for the observations at the start of the
                 * //analysis window. The probability of seeing all the
                 * //observations in the analysis window in state j is updated each
                 * //time d is incrememented two lines below, so we only need to
                 * //find the observation probabilities for one time step, each
                 * //time d is updated:*/
                
                for(i = 0; i<N; i++)
                {
                    if(observation_probs[start +i*T] == 0){
                        observation_probs[start +i*T] = DBL_MIN;
                    }
                }
                
                probs = 0;
                for(i = start; i<end_t; i++){
                    probs = probs + log(observation_probs[i +j*T]);
                }
                
                
                if(probs == 0){
                    probs = DBL_MIN;
                }
                
                
                emission_probs = (probs);
                
                
//                     if(t<5){
//                         mexPrintf ("emission_probs:%f \n",emission_probs);
//                     }
                
                
                
                /*Find the normaliser for the duration probablities, being the sum across all allowed times in each state:*/
                for(i = 0; i<max_duration_D; i++)
                {
                    if (isnan(duration_probs[j][i])) {
                        
                    }
                    else{
                        duration_sum = duration_sum + duration_probs[j][i];
                    }
                    
                    
                }
                
                
                
                /*Find the total probability of transitioning from the last
                 * //state to this one, with the observations and being in the same
                 * //state for the analysis window. This is the duration-dependant
                 * //variation of equation 33a from Rabiner:*/
                delta_temp = max_delta + (emission_probs)+ (log((duration_probs[j][d]/duration_sum)));
                
                
//                     if(t<5){
//                         mexPrintf ("max_delta: %f \n",max_delta);
//
//                         mexPrintf ("emission_probs: %f \n",emission_probs);
//                         mexPrintf ("duration_probs[j][d]: %f \n",duration_probs[j][d]);
//                         mexPrintf ("duration_sum: %f \n",duration_sum);
//
//                         mexPrintf ("log((duration_probs[j][d]/duration_sum)): %f \n",log((duration_probs[j][d]/duration_sum)));
//
//                         mexPrintf ("delta_temp: %f \n",delta_temp);
//
//                     }
                
                
                
                /*
                 * Unlike equation 33a from Rabiner, the maximum delta could come
                 * from multiple d values, or from multiple size of the analysis
                 * window. Therefore, only keep the maximum delta value over the
                 * entire analysis window:
                 * If this probability is greater than the last greatest,
                 * update the delta matrix and the time duration variable:
                 */
                
                
                
                
//                     mexPrintf ("delta_temp: %f \n",delta_temp);
//
//                     mexPrintf ("delta_temp>delta?: %d \n",delta_temp>delta[t+j*T]);
                
                if(delta_temp>delta[t+j*T]){
                    
                    delta[t+j*T] = delta_temp;
                    psi[t+j*T] = max_index+1;
                    
                    psi_duration_out[t + j*T] = d;
                    
                }
            }
        }
    }
    
}

/*   Now, define the gateway function, i.e., mexFunction.Below
 * is the standard, predeclared header to mexFunction.  nlhs
 * and nrhs are the number of left-hand and right-hand side
 * arguments that mexample was called with from within MATLAB.
 * In this example, nlhs equals 1 and nrhs should equal 2.  If
 * not, then the user has called mexample the wrong way and
 * should be informed of this.  plhs and prhs are arrays which
 * contain the pointers to the MATLAB arrays, which are
 * stored in a C struct called an Array.  prhs is an array of
 * length rhs,and its pointers point to valid input data.
 * plhs is an array of length nlhs, and its pointers point to
 * invalid data (i.e., garbage).  It is the job of mexFunction
 * to fill plhs with valid data.
 *
 * First, define the following values.  This makes it much
 * easier to change the order of inputs to mexample, should we
 * want to change the function later.  In addition, it makes
 * the code easier to read. */

#define N prhs[0]
#define T prhs[1]
#define a_matrix prhs[2]
#define max_duration_D prhs[3]
#define delta  prhs[4]
#define observation_probs prhs[5]
#define duration_probs prhs[6]
#define psi prhs[7]


#define delta_out plhs[0]
#define psi_out plhs[1]
#define psi_duration plhs[2]


void mexFunction(
        int     nlhs,
        mxArray  *plhs[],
        int     nrhs,
        const mxArray  *prhs[]
        )
{
    double a_matrix_in[4][4];/* 2 dimensional C array to pass to workFcn() */
    double *delta_in_matrix;/* 2 dimensional C array to pass to workFcn() */
    double *observation_probs_matrix;/* 2 dimensional C array to pass to workFcn() */
    double *psi_matrix;/* 2 dimensional C array to pass to workFcn() */
    
    double duration_probs_matrix[4][150];/* 2 dimensional C array to pass to workFcn() */
    
    int actual_T;
    int fake_T_extended;
    int actual_N;
    int max_duration_D_val;
    
    int    row,col;        /* loop indices */
    int    m,n;            /* temporary array size holders */
    
    /*   Step 1: Error Checking Step 1a: is nlhs 1?  If not,
     * generate an error message and exit mexample (mexErrMsgTxt
     * does this for us!) */
    if (nlhs!=3)
        mexErrMsgTxt("mexample requires 3 output argument.");
    
    /*   Step 1b: is nrhs 2? */
    if (nrhs!=8)
        mexErrMsgTxt("mexample requires 8 input arguments");
    
    
    actual_T = mxGetM(observation_probs);
    actual_N = mxGetN(observation_probs);
    
    max_duration_D_val = mxGetScalar(max_duration_D);

        
    /*   Step 2:  Allocate memory for return argument(s) */
    delta_out = mxCreateDoubleMatrix((actual_T+max_duration_D_val-1), actual_N, mxREAL);
    psi_out = mxCreateDoubleMatrix((actual_T+max_duration_D_val-1), actual_N, mxREAL);
    psi_duration = mxCreateDoubleMatrix((actual_T+max_duration_D_val-1), actual_N, mxREAL);
    
    /*   Step 3:  Convert ARRAY_IN to a 2x2 C array
     * MATLAB stores a two-dimensional matrix in memory as a one-
     * dimensional array.  If the matrix is size MxN, then the
     * first M elements of the one-dimensional array correspond to
     * the first column of the matrix, and the next M elements
     * correspond to the second column, etc. The following loop
     * converts from MATLAB format to C format: */
    
    for (col=0; col < mxGetN(a_matrix); col++){
        for (row=0; row < mxGetM(a_matrix); row++){
            a_matrix_in[row][col] =(mxGetPr(a_matrix))[row+col*mxGetM(a_matrix)];
        }
    }
    
    
    delta_in_matrix = mxGetPr(delta);
    observation_probs_matrix = mxGetPr(observation_probs);
    psi_matrix = mxGetPr(psi);
    
    /*     for (col=0; col < mxGetN(delta); col++){
     * //         for (row=0; row < mxGetM(delta); row++){
     * //
     * //
     * //             observation_probs_matrix[row][col] =(mxGetPr(observation_probs))[row+col*mxGetM(observation_probs)];
     * //             psi_matrix[row][col] =(mxGetPr(psi))[row+col*mxGetM(psi)];
     * //         }
     * //     }*/
    
    
    for (col=0; col < mxGetN(duration_probs); col++){
        for (row=0; row < mxGetM(duration_probs); row++){
            duration_probs_matrix[row][col] =(mxGetPr(duration_probs))[row+col*mxGetM(duration_probs)];
            
        }
    }
    
    
    
    
    /*   mxGetPr returns a pointer to the real part of the array
     * ARRAY_IN.  In the line above, it is treated as the one-
     * dimensional array mentioned in the previous comment.  */
    
    /*   Step 4:  Call workFcn function */
    viterbi(actual_N,actual_T,a_matrix_in,max_duration_D_val,delta_in_matrix,observation_probs_matrix,duration_probs_matrix,psi_matrix,mxGetPr(psi_duration));
    memcpy ( mxGetPr(delta_out), delta_in_matrix, actual_N*actual_T*8);
    memcpy ( mxGetPr(psi_out), psi_matrix, actual_N*actual_T*8);
    
}