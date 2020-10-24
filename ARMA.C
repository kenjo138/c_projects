#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "rtdspc.h"

float iir_biquad(float input,float *a,float *b);
float iir_adapt_filter(float input,float d,float *a,float *b);

#define LEN 7000

void main()
{
    int i;
    float y;
    static float d[LEN];
    float b[2] = { 0.7477891445, -0.2722149193 };
    float a[3] = { 0.187218, 0.149990698, 0.187218 };

/* set the random seed */
    srand(1);
    for(i = 0 ; i < LEN ; i++) d[i] = iir_biquad(gaussian(),a,b);

/* clear all coefficients to zero at start of adaptive process */
    b[0] = b[1] = 0.0;
    a[0] = a[1] = a[2] = 0.0;

/* reset the random seed to re-generate the same random sequence */
    srand(1);
    for(i = 0 ; i < LEN ; i++) {
        y = iir_adapt_filter(gaussian(),d[i],a,b);
        printf("\n%f %f %f %f %f %f %f",
            d[i],y,a[0],a[1],a[2],b[0],b[1]);
    }
}

/* 2 poles (2 b coefs) and 2 zeros (3 a coefs) iir filter single biquad */

float iir_biquad(float input,float *a,float *b)
{
    static float out_hist1,out_hist2;
    static float in_hist1,in_hist2;
    float output;

    output = out_hist1 * b[0];
    output += out_hist2 * b[1];             /* poles */

    output += input * a[0];
    output += in_hist1 * a[1];            /* zeros */
    output += in_hist2 * a[2];

    out_hist2 = out_hist1;          /* update history */
    out_hist1 = output;

    in_hist2 = in_hist1;
    in_hist1 = input;

    return(output);
}

/* 2 poles (2 b coefs) and 2 zeros (3 a coefs) adaptive iir biquad filter */

float iir_adapt_filter(float input,float d,float *a,float *b)
{
    int i;
    static float out_hist1,out_hist2;
    static float beta[2],beta_h1[2],beta_h2[2];
    static float alpha[3],alpha_h1[3],alpha_h2[3];
    static float in_hist[3];
    float output,e;

    output = out_hist1 * b[0];
    output += out_hist2 * b[1];             /* poles */

    in_hist[0] = input;
    for(i = 0 ; i < 3 ; i++)
        output += in_hist[i] * a[i];            /* zeros */

/* calclulate alpha and beta update coefficients */
    for(i = 0 ; i < 3 ; i++)
        alpha[i] = in_hist[i] + b[0]*alpha_h1[i] + b[1]*alpha_h2[i];

    beta[0] = out_hist1 + b[0]*beta_h1[0] + b[1]*beta_h2[0];
    beta[1] = out_hist2 + b[0]*beta_h1[1] + b[1]*beta_h2[1];

/* error calculation */
    e = d - output;
/* update coefficients */
    a[0] += e*0.2*alpha[0];
    a[1] += e*0.1*alpha[1];
    a[2] += e*0.06*alpha[2];

    b[0] += e*0.04*beta[0];
    b[1] += e*0.02*beta[1];

/* update history for alpha */
    for(i = 0 ; i < 3 ; i++) {
        alpha_h2[i] = alpha_h1[i];
        alpha_h1[i] = alpha[i];
    }

/* update history for beta */
    for(i = 0 ; i < 2 ; i++) {
        beta_h2[i] = beta_h1[i];
        beta_h1[i] = beta[i];
    }

/* update input/output history */
    out_hist2 = out_hist1;
    out_hist1 = output;

    in_hist[2] = in_hist[1];
    in_hist[1] = input;

    return(output);
}
