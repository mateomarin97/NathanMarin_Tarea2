#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define L 128.0
#define N 1000
#define gamma 1.4
#define E0 1.0
void inicializar(float *u,float *rho,float *P,float *u2,float *u3,float *F2,float *F3,float *c,float *uc);
void LaxWendroff(float *u,float *rho,float *P,float *u2,float *u3,float *F2,float *F3,float *c,float *uc);
float max(float *A);
float min(float *A);
int wherefrente(float *A);
int wheremax(float *A);

int main(){
  /*Declaro las listas que necesitare*/
  float *u;
  float *rho;
  float *P;
  float *u2;
  float *u3;
  float *F2;
  float *F3;
  float *c;
  float *uc;
  /*Inicializo las listas*/
  u=malloc(N*sizeof(float));
  rho=malloc(N*sizeof(float));
  P=malloc(N*sizeof(float));
  u2=malloc(N*sizeof(float));
  u3=malloc(N*sizeof(float));
  F2=malloc(N*sizeof(float));
  F3=malloc(N*sizeof(float));
  c=malloc(N*sizeof(float));
  uc=malloc(N*sizeof(float));
  inicializar(u,rho,P,u2,u3,F2,F3,c,uc);
  /*Ejecutamos LaxWendroff*/
  LaxWendroff(u,rho,P,u2,u3,F2,F3,c,uc);
  /*Imprimimos los datos*/
  int i;
  FILE *archivo;
  archivo=fopen("u3sedov.dat","w");
  fclose(archivo);
  archivo=fopen("u3sedov.dat","a");
  for(i=0;i<N;i++){
    fprintf(archivo,"%f ",u[i]);
  }
  fclose(archivo);

  archivo=fopen("p3sedov.dat","w");
  fclose(archivo);
  archivo=fopen("p3sedov.dat","a");
  for(i=0;i<N;i++){
    fprintf(archivo,"%f ",P[i]);
  }
  fclose(archivo);

  archivo=fopen("rho3sedov.dat","w");
  fclose(archivo);
  archivo=fopen("rho3sedov.dat","a");
  for(i=0;i<N;i++){
    fprintf(archivo,"%f ",rho[i]);
  }
  fclose(archivo);







  return 0;
}


void inicializar(float *u,float *rho,float *P,float *u2,float *u3,float *F2,float *F3,float *c,float *uc){
  int i;
  int dL = L/N;
  for(i=0;i<N;i++){
    u[i]=0.0;
    rho[i]=1.0;
    P[i]=0.8403;


    u2[i]=rho[i]*u[i];
    u3[i]=(P[i]/(gamma-1.0))+0.5*rho[i]*u[i]*u[i];
    F2[i]=(rho[i]*u[i]*u[i])+P[i];
    F3[i]=u[i]*((gamma*P[i]/(gamma-1.0))+0.5*rho[i]*u[i]*u[i]);
    c[i]=sqrt(((gamma+1.0)/rho[i])*((gamma*P[i]/(gamma-1.0))+0.5*rho[i]*u[i]*u[i]));
    uc[i]=u[i]+c[i];
  }
  // P[0] = (E0*(gamma-1.0)/(dL*dL*dL))/(101000.0);
  P[0] = 5.0;
  u2[0]=rho[0]*u[0];
  u3[0]=(P[0]/(gamma-1.0))+0.5*rho[0]*u[0]*u[0];
  F2[0]=(rho[0]*u[0]*u[0])+P[0];
  F3[0]=u[0]*((gamma*P[0]/(gamma-1.0))+0.5*rho[0]*u[0]*u[0]);
  c[0]=sqrt(((gamma+1.0)/rho[0])*((gamma*P[0]/(gamma-1.0))+0.5*rho[0]*u[0]*u[0]));
  uc[0]=u[0]+c[0];

}


void LaxWendroff(float *u,float *rho,float *P,float *u2,float *u3,float *F2,float *F3,float *c,float *uc){

  /*Declaro variables y constantes importantes*/
  float dx;
  int ic;
  int icc;
  int iccc;
  float tiempo;
  dx= L/N;
  ic=(int)(120.0/dx);
  ic= ic+1;
  icc=(int)(60.0/dx);
  icc= icc+1;
  iccc=(int)(10.0/dx);
  iccc= iccc+1;
  int i;
  float dt;
  dt=2.0*dx/max(uc);
  tiempo=dt;
  /*Defino las listas intermedias*/
  float *rhoip;
  float *rhoin;
  float *u2ip;
  float *u2in;
  float *u3ip;
  float *u3in;
  float *uip;
  float *uin;
  float *Pip;
  float *Pin;
  float *F2ip;
  float *F2in;
  float *F3ip;
  float *F3in;
  FILE *archivo;

  rhoip=malloc(N*sizeof(float));
  rhoin=malloc(N*sizeof(float));
  u2ip=malloc(N*sizeof(float));
  u2in=malloc(N*sizeof(float));
  u3ip=malloc(N*sizeof(float));
  u3in=malloc(N*sizeof(float));
  uip=malloc(N*sizeof(float));
  uin=malloc(N*sizeof(float));
  Pip=malloc(N*sizeof(float));
  Pin=malloc(N*sizeof(float));
  F2ip=malloc(N*sizeof(float));
  F2in=malloc(N*sizeof(float));
  F3ip=malloc(N*sizeof(float));
  F3in=malloc(N*sizeof(float));


  int j;
  j=0;
  int con;
  con=0;
  int j2;
  j2=0;
  int j3;
  j3=0;


  /*Empieza el metodo*/
  while(j==0 && con<100000){
    con++;
    /*Primer paso*/
    i=0;
    rhoip[i]=0.5*(rho[i+1]+rho[i])-(0.5*dt/dx)*(u2[i+1]-u2[i]);
    rhoin[i]=rho[i]-(dt/dx)*u2[i];
    u2ip[i]=0.5*(u2[i+1]+u2[i])-(0.5*dt/dx)*(F2[i+1]-F2[i]);
    u2in[i]=0.0;
    u3ip[i]=0.5*(u3[i+1]+u3[i])-(0.5*dt/dx)*(F3[i+1]-F3[i]);
    u3in[i]=u3[i]-(dt/dx)*F3[i];

    for(i=1;i<N-1;i++){
      rhoip[i]=0.5*(rho[i+1]+rho[i])-(0.5*dt/dx)*(u2[i+1]-u2[i]);
      rhoin[i]=0.5*(rho[i]+rho[i-1])-(0.5*dt/dx)*(u2[i]-u2[i-1]);
      u2ip[i]=0.5*(u2[i+1]+u2[i])-(0.5*dt/dx)*(F2[i+1]-F2[i]);
      u2in[i]=0.5*(u2[i]+u2[i-1])-(0.5*dt/dx)*(F2[i]-F2[i-1]);
      u3ip[i]=0.5*(u3[i+1]+u3[i])-(0.5*dt/dx)*(F3[i+1]-F3[i]);
      u3in[i]=0.5*(u3[i]+u3[i-1])-(0.5*dt/dx)*(F3[i]-F3[i-1]);
    }
    i=N-1;
    rhoip[i]=rho[i]+(dt/dx)*u2[i];
    rhoin[i]=0.5*(rho[i]+rho[i-1])-(0.5*dt/dx)*(u2[i]-u2[i-1]);
    u2ip[i]=0.0;
    u2in[i]=0.5*(u2[i]+u2[i-1])-(0.5*dt/dx)*(F2[i]-F2[i-1]);
    u3ip[i]=u3[i]+(dt/dx)*F3[i];
    u3in[i]=0.5*(u3[i]+u3[i-1])-(0.5*dt/dx)*(F3[i]-F3[i-1]);

    for(i=0;i<N;i++){
      /*Rece porque rho no sea cero nunca*/
      uip[i]=u2ip[i]/rhoip[i];
      uin[i]=u2in[i]/rhoin[i];
      Pip[i]=(gamma-1.0)*(u3ip[i]-0.5*u2ip[i]*uip[i]);
      Pin[i]=(gamma-1.0)*(u3in[i]-0.5*u2in[i]*uin[i]);
      F2ip[i]=u2ip[i]*uip[i]+Pip[i];
      F2in[i]=u2in[i]*uin[i]+Pin[i];
      F3ip[i]=uip[i]*((gamma*Pip[i]/(gamma-1.0))+0.5*u2ip[i]*uip[i]);
      F3in[i]=uin[i]*((gamma*Pin[i]/(gamma-1.0))+0.5*u2in[i]*uin[i]);
      /*Segundo paso*/
      rho[i]=rho[i]-(dt/dx)*(u2ip[i]-u2in[i]);
      u2[i]=u2[i]-(dt/dx)*(F2ip[i]-F2in[i]);
      u3[i]=u3[i]-(dt/dx)*(F3ip[i]-F3in[i]);
      /*Esperemos que rho no sea cero*/
      u[i]=u2[i]/rho[i];
      P[i]=(gamma-1.0)*(u3[i]-0.5*u2[i]*u[i]);
      F2[i]=u2[i]*u[i]+P[i];
      F3[i]=u[i]*((gamma*P[i]/(gamma-1.0))+0.5*u2[i]*u[i]);
      c[i]=sqrt(((gamma+1.0)/rho[i])*((gamma*P[i]/(gamma-1.0))+0.5*rho[i]*u[i]*u[i]));
      uc[i]=u[i]+c[i];

    }
    dt=2.0*dx/max(uc);
    tiempo=tiempo+dt;
    if(wherefrente(rho)>=ic){
      j=1;
    }

    if(wherefrente(rho)>=icc && j2==0){
      j2=1;
      archivo=fopen("t2sedov.dat","w");
      fclose(archivo);
      archivo=fopen("t2sedov.dat","a");
      fprintf(archivo,"%f ",tiempo);
      fclose(archivo);
      archivo=fopen("u2sedov.dat","w");
      fclose(archivo);
      archivo=fopen("u2sedov.dat","a");
      for(i=0;i<N;i++){
	fprintf(archivo,"%f ",u[i]);
      }
      fclose(archivo);

      archivo=fopen("p2sedov.dat","w");
      fclose(archivo);
      archivo=fopen("p2sedov.dat","a");
      for(i=0;i<N;i++){
	fprintf(archivo,"%f ",P[i]);
      }
      fclose(archivo);

      archivo=fopen("rho2sedov.dat","w");
      fclose(archivo);
      archivo=fopen("rho2sedov.dat","a");
      for(i=0;i<N;i++){
	fprintf(archivo,"%f ",rho[i]);
      }
      fclose(archivo);

    }



    if(wherefrente(rho)>=iccc && j3==0){
      j3=1;
      archivo=fopen("t1sedov.dat","w");
      fclose(archivo);
      archivo=fopen("t1sedov.dat","a");
      fprintf(archivo,"%f ",tiempo);
      fclose(archivo);
      archivo=fopen("u1sedov.dat","w");
      fclose(archivo);
      archivo=fopen("u1sedov.dat","a");
      for(i=0;i<N;i++){
	fprintf(archivo,"%f ",u[i]);
      }
      fclose(archivo);

      archivo=fopen("p1sedov.dat","w");
      fclose(archivo);
      archivo=fopen("p1sedov.dat","a");
      for(i=0;i<N;i++){
	fprintf(archivo,"%f ",P[i]);
      }
      fclose(archivo);

      archivo=fopen("rho1sedov.dat","w");
      fclose(archivo);
      archivo=fopen("rho1sedov.dat","a");
      for(i=0;i<N;i++){
	fprintf(archivo,"%f ",rho[i]);
      }
      fclose(archivo);

    }
  



    
  }
  
  archivo=fopen("t3sedov.dat","w");
  fclose(archivo);
  archivo=fopen("t3sedov.dat","a");
  fprintf(archivo,"%f ",tiempo);
  fclose(archivo);

}


float max(float *A){
  float maximo;
  maximo=A[0];
  int i;
  for(i=1;i<N;i++){
    if(A[i]>maximo){
      maximo=A[i];
    }
  }

  return maximo;
}



float min(float *A){
  float minimo;
  minimo=A[0];
  int i;
  for(i=1;i<N;i++){
    if(A[i]<minimo){
      minimo=A[i];
    }
  }

  return minimo;
}

int wherefrente(float *A){
  int i;
  int f;
  f=0;
  for(i=wheremax(A);i<N;i++){
    if(A[i]<A[N-1]+0.0001 && f==0){
      f=i;
    }
  }
  return f;
}

int wheremax(float *A){
  float maximo;
  maximo=A[0];
  int where;
  where=0;
  int i;
  for(i=1;i<N;i++){
    if(A[i]>maximo){
      maximo=A[i];
      where=i;
    }
  }

  return where;

}
