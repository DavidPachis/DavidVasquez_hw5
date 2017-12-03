// el siguiente ejercicio se baso en el matrial proporcionado en clase,de  Bayesian parameter estimation
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define bb  0.2497 // en kiloparecs, 1 parsec = 206265 ua = 3,2616 años luz = 3,0857 × 10^16 m, mil parsecs, "3262 años luz".
#define bd  5.16   //mil parsecs, "3262 años luz".
#define ad  0.3105 //mil parsecs, "3262 años luz".
#define ah  64.3   //mil parsecs, "3262 años luz".
#define n 2000 // numero de trotaditas (iteraciones en la caminata aleatoria).
#define escala 100000000
double likelihood(double *Vc,double *Vr,double chifinal);
void mi_model(double mb,double md,double mh, double *R,double *Vc);
void trotadita (double *Vr,double *Vc, double chifinal,double *R);

int main(void)
{ 
 double MB, MD, MH;
 double *R,*Vr,*Vc,*Vfinal; 
 R  = malloc(300 * sizeof(double));
 Vr = malloc(300 * sizeof(double));// velocidad medida
 Vc = malloc(300 * sizeof(double));

 double chifinal;
 MB= fabs((200000*drand48())-100000);
 MD= fabs((200000*drand48())-100000);
 MH= fabs((200000*drand48())-100000);
 int i;
 FILE *in;
 float var;
 char filename[100]="RadialVelocities.dat";

 in = fopen(filename, "r");
 if(!in)
 {
  printf("la habeis cagado abriendo el archivo tio %s\n", filename);
  exit(1);
 }
 //printf("ahora si, como asi que como asi que como fue \n");
 char greeting1[60];
 char greeting2[60];
 fscanf(in, "%s %s \n", &greeting1,&greeting2);
 for (i=0;i<300;i++)
 {
   fscanf(in, "%lf %lf\n", &R[i],&Vr[i]);
   //R[i] = var;
   //fscanf(in, "%f\n", &var);
   //Vr[i] = var;
  //printf("%lf \n",Vr[i]);
 }
  //printf("archivos leidos y copiados en arreglos\n");
 mi_model(MB,MD,MH,R,Vc);
 likelihood(Vc,Vr,chifinal);
 trotadita (Vr,Vc,chifinal,R);
 
                                                     
 //printf("FIN FIN FIN\n");
 FILE *in2;
 char filename2[100]="new_data.dat";
 int z;
 in2 = fopen(filename2, "w");
 if(!in2)
 {
  printf("la habeis cagado abriendo el archivo tio %s\n", filename2);
  exit(1);
 }
 for (z=0;z<300;z++)
  {
    fprintf(in2, "%f %f\n",Vr[z],Vc[z]);
  }
 fclose(in2);
 printf ("FIN \n");
}
double likelihood(double *Vc,double *Vr,double chifinal)
{ 
 double *sum,*chis;
 int i;
 double cuenta;
 double total=0;
 sum = malloc(300 * sizeof(double));
 chis = malloc(300 * sizeof(double));
 double chi_squared;
 for (i=0;i<300;i++)
 {
   sum[i]= (Vr[i]-Vc[i])*(Vr[i]-Vc[i]);
 } 
 for (i=0;i<300;i++)
 {
  chi_squared = 0.5*sum[i];
  chis[i]= chi_squared;
  cuenta = cuenta+chi_squared;
  total = total+1;
 } 
 chifinal = exp(-( cuenta/total));
 return chifinal;
 //printf("FIN like\n");
}

void mi_model(double mb,double md,double mh, double *R,double *Vc)
{ 
 int i;
 double Mb =pow(mb,0.5);
 double Md =pow(md,0.5);
 double Mh =pow(mh,0.5);
 for (i=0; i<300;i++)// verificar si es desde el paso uno 
 { 
  Vc[i] = ((Mb*R[i])/pow((R[i]*R[i]+bb*bb),0.75))+((Md*R[i])/pow((R[i]*R[i]+(bd+ad)*(bd+ad)),0.75))+((Mh)/pow((R[i]*R[i]+(ah*ah)),0.25));
 } 
 //printf("FIN mi_model\n");
}

void trotadita (double *Vr,double *Vc, double chifinal,double *R)
{
 double alpha;
 double beta;
 double *mb_run,*md_run,*mh_run,*l_run;
 mb_run = malloc(n * sizeof(double));
 md_run = malloc(n * sizeof(double));
 mh_run = malloc(n * sizeof(double));
 l_run  = malloc(n * sizeof(double));
 double mb_prima;  
 double md_prima; 
 double mh_prima; 
 double Vc_incial;
 double Vc_prima;
 double l_prima;
 double l_inicial;
 int i; 
 double grande;
 grande = pow(10.0,7.0);
 // generar los primeros pasos
 mb_run[0]= fabs((200000*drand48())-100000);
 md_run[0]= fabs((200000*drand48())-100000);
 mh_run[0]= fabs((200000*drand48())-100000);
 //printf("%f %f %f\n",mb_run[0],md_run[0],mh_run[0]);
 mi_model(mb_run[0],md_run[0],mh_run[0], R, Vc);
 Vc_incial = Vc[0];
 
 l_run[0] = likelihood(Vc,Vr,chifinal); // Vr es velocidad medida y Vc ls velocidad del modelo 
 
 for (i=0; i<n;i++)
 { 
  mb_prima = fabs((200000*drand48())-100000); 
  //printf("%f \n",mb_prima);
  md_prima = fabs((200000*drand48())-100000);
  mh_prima = fabs((200000*drand48())-100000);
  mi_model(mb_run[i],md_run[i],mh_run[i], R, Vc);
  Vc_incial = Vc[i];
  l_inicial = likelihood(Vc,Vr,chifinal);
  mi_model(mb_prima,md_prima,mh_prima, R, Vc);
  Vc_prima  = Vc[i]; 
  //printf("%f \n",Vc[i]);
  l_prima   = likelihood(Vc,Vr,chifinal) ;
  alpha =  l_prima/l_inicial;
  if (alpha>=1.0)
  { 
   mb_run[i]= mb_prima;
   md_run[i]= md_prima;
   mh_run[i]= mh_prima;
   l_run[i] = l_prima;
   
  }
  else
  {
    beta = drand48();
    if (beta<=alpha)
    {
      mb_run[i]= mb_prima;
      md_run[i]= md_prima;
      mh_run[i]= mh_prima;
      l_run[i] = l_prima;
    }
   else
    {
      mb_run[i]= mb_run[i-1];
      md_run[i]= md_run[i-1];
      mh_run[i]= mh_run[i-1];
      l_run[i] = l_run[i-1];
     
    }
  }
 }
  double max;
  int c;
  int lugar = 0;
  max = l_run[0];
 
  for (c = 1; c < n; c++)
  {
    if (l_run[c] > max)
    {
       max = l_run[c];
       lugar = c;
    }
  }
  mi_model(mb_run[lugar],md_run[lugar],mh_run[lugar],R,Vc);
 printf("%f %f %f\n",mb_run[lugar],md_run[lugar],mh_run[lugar]);

}









