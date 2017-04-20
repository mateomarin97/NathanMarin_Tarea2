#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define L 1.0
#define N 1000
#define gamma 1.4
void inicializar(float *u,float *rho,float *P,float *u2,float *u3,float *F2,float *F3,float *c,float *uc);
void LaxWendroff(float *u,float *rho,float *P,float *u2,float *u3,float *F2,float *F3,float *c,float *uc);
float max(float *A);
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
  archivo=fopen("ushock.dat","w");
  fclose(archivo);
  archivo=fopen("ushock.dat","a");
  for(i=0;i<N;i++){
    fprintf(archivo,"%f ",u[i]);
  }
  fclose(archivo);

  archivo=fopen("pshock.dat","w");
  fclose(archivo);
  archivo=fopen("pshock.dat","a");
  for(i=0;i<N;i++){
    fprintf(archivo,"%f ",P[i]);
  }
  fclose(archivo);

  archivo=fopen("rhoshock.dat","w");
  fclose(archivo);
  archivo=fopen("rhoshock.dat","a");
  for(i=0;i<N;i++){
    fprintf(archivo,"%f ",rho[i]);
  }
  fclose(archivo);


  
  
 
  

  return 0;
}


void inicializar(float *u,float *rho,float *P,float *u2,float *u3,float *F2,float *F3,float *c,float *uc){
  int i;
  for(i=0;i<N;i++){
    u[i]=0.0;

    if(i<N/2){
      rho[i]=1.0;
      P[i]=1.0;
    }

    else {
      rho[i]=0.1;
      P[i]=0.1;
    }

    u2[i]=rho[i]*u[i];
    u3[i]=(P[i]/(gamma-1.0))+0.5*rho[i]*u[i]*u[i];
    F2[i]=(rho[i]*u[i]*u[i])+P[i];
    F3[i]=u[i]*((gamma*P[i]/(gamma-1.0))+0.5*rho[i]*u[i]*u[i]);
    c[i]=343.2;
    uc[i]=u[i]+c[i];
  }


}


void LaxWendroff(float *u,float *rho,float *P,float *u2,float *u3,float *F2,float *F3,float *c,float *uc){
 
  /*Declaro variables y constantes importantes*/
  float dx;
  int ic;
  dx= L/N;
  ic=(int)(0.9*N-1);
  ic= ic+1;
  int i;
  float dt;
  dt=2.0*dx/max(uc);
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

  /*Empieza el metodo*/
  while(j==0){

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
      uc[i]=u[i]+c[i];
     
    }
    dt=2.0*dx/max(uc);
    
    if(wherefrente(u)>=ic){
      j=1;
    }
  }
 
  
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

int wherefrente(float *A){
  int i;
  int f;
  f=0;
  for(i=wheremax(A);i<N;i++){
    if(A[i]<0.01 && f==0){
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
