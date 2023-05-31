#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <time.h>

// This code was kindly provided by Dr. fredrik Jansson.
// It was been slightly modified to read input and save 
// output from/to a file.

int fastDerridaSmall (int N, double *vl, double *vr, double *V, double *D);

int main(void)
{
  int count = 0;
  char c;

  FILE *fpin;

  if ((fpin = fopen("rates.txt", "r")) == NULL)
  { 
    perror("fopen rates.txt");
    return 1;
  }

  for (c = getc(fpin); c != EOF; c = getc(fpin))
  {
    if (c == '\n') 
    {
      count = count + 1;
    }
  }

  int N = count;
  int i; 
  double forward_rates[N];
  double reverse_rates[N];
  double vl[N+1];
  double vr[N];

  rewind(fpin);

  for(i = 0; i < N; i++)
  {
    fscanf(fpin, "%lf,%lf", &forward_rates[i], &reverse_rates[i]);
    vl[i] =  forward_rates[i];
    if(i != 0) {vr[i] = reverse_rates[i-1];}
    if(i == N-1) {vr[0] = reverse_rates[N-1];}
  }

  printf("The are %d electron transfer steps \n", N);
  for(i = 0; i < N; i++) {printf("%E %E\n", vl[i], vr[i]);}

  fclose(fpin);

  double V, D;    // variables for the velocity and diffucion results
  int inaccurate; // flag for inaccurate result
  
  inaccurate = fastDerridaSmall(N, vl, vr, &V, &D);

  if(inaccurate)
    printf("Inaccurate result. Try mirroring the chain so that the drift is to the right.\n");
  else
    {
      printf("V: %f\n", V);
      printf("D: %e\n", D);
    }
}

void /*@noreturn@*/ die (char *s)
{
  perror (s);
  exit (EXIT_FAILURE);
}

/*
 Derrida's formula calculated in linear time.
 Version that uses less memory.
 Returns true if precision is bad.
 Aborts calculations and leave *V and *D unchanged in this case.
 The routine works better if the average velocity is to the right, positive.
 Inaccurate results can occur if the chain is long and the avrage 
 velocity is negative. Mirror the chain in that case.

  vl MUST have one extra element, used for optimizing boundary conditions.
 */
int fastDerridaSmall (int N, double *vl, double *vr, double *V, double *D)
{
  FILE *fpout;
  fpout = fopen("D.txt", "w");
  int n, i;
  double s, t, x;
  double u_n;
  double rsum = 0;

  double G;  //product of all g
  double H;  //product of all h
  double *r;
  //double *g;  //the fractions vl[n]/vr[n]   
  //double *h;  //the fractions vl[n+1]/vr[n]
  
  double Gn;
  double Hn;
  int inaccurate = 0;

  r = malloc((N+2)*sizeof(r[0]));      //one extra for PBC, the second one is read but not used
  if (r == NULL) 
    die ("fastDerridaSmall, malloc r");

  vl[N] = vl[0];  /*!!! vl[N] must be defined !!!*/

#define g(n) (vl[n]/vr[n])
#define h(n) (vl[(n)+1]/vr[n])

  //Eq. 50
  Gn = 1;
  t = 1;
  n = N-1;
  for (i = 1; i <= N-1; i++)
    {
      t *= g(i-1); // Really g[n+i], but n=N-1 =-1
      Gn += t;
//    printf ("i %f n %e vl[i-1] %e vr[i-1] %e g(i-1) %e Gn %e \n", (double)i, (double)n, (double)vl[i-1], (double)vr[i-1], (double)g(i-1), (double)Gn);
    }
  G = t*g(n);  //product of all g
     //printf ("t %e g(n) %e G %e \n", (double)t, (double)g(n), (double)G);

  if (G > 1e12 || !isfinite(G))
    {
      fprintf (stderr, "G=%e Result may be inaccurate! Giving up.\n",(double)G);
      inaccurate = 1;
      goto quit;
    }

  s = 0;

  for (n = N-1; n >= 0; n--)
    {
      r[n] = Gn / vr[n];
      //printf ("n %e Gn %e vr[n] %e r[n] %e \n", (double)n, (double)Gn, (double)vr[n], (double)r[n]);
      Gn = g(n) * Gn - G + 1;
      rsum += r[n];
      s += n * r[n];       //initial "triangle"-sum S_0
    }
  
  s += N * r[0];

  //Eq. 46
  Hn = 1;
  t = 1;
  n = 0;
  for (i = 1; i <= N-1; i++)
    {
      t *= h(N+n-i); //n is 0, bound OK
      Hn += t;
      //printf ("N %e n %e i %e h %e t %e Hn %e \n", (double)N, (double)n, (double)i, (double)h(N+n-i), (double)t, (double)Hn);
    }
  H = t*h(n);  //product of all h (= G, product of all g)
  
  //Eq. 49
  *V = N/rsum * (1-H);
  //printf ("rsum %e H %e \n", (double)rsum, (double)H);
  
  //Eq. 47
  t = 0;
  x = 0;
  
  //  for (i = 1; i <= N; i++)             
  // s += i * r[i];        

  for (n = 0; n < N; n++)
    {
      u_n = Hn / vr[n];
      Hn = h(n) * Hn - H + 1;
      //usum += u_n;

      t += u_n * s;
      s += N*r[n+1] - rsum;

      x += vr[n] * u_n * r[n];
    //printf ("x %e vr[n] %e u_n %e r[n] %e \n", (double)x, (double)vr[n], (double)u_n, (double)r[n]);
    }

  //rsum = usum according to Eq 48
  //  printf ("rsum: %10.15g\n usum: %10.15g\n", (double)rsum, (double)usum);
  //  printf ("rsum %e \n", (double)rsum);

  t *= *V;
  x *= N;  

//  printf ("rsum %e t %e x %e *V %e D %e \n", (double)rsum, (double)t, (double)x, (double)*V, (double)*D);
  
    *D = ( 1 / (rsum*rsum) * (t + x) - *V * (N+2.0)/2.0 ) * 2.5e-15;

    fprintf(fpout, "Diffusion constant = %E cm^2/S", *D);
    fclose(fpout);

 quit:
  free (r);

  return inaccurate;
}
