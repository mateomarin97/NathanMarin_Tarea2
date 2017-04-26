#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define dL 2.0
#define L 256.0
#define gamma 1.4
#define E0 10000000000.0
float sgnx(int i,int Nf);
float sgny(int i,int Np);
float sgnz(int i,int N);
void inicializar(int N,int Nf,int Np,int ic, float *rho,float *c,float *cu,float *cv,float *cw,float *r,float *P,float *u,float *v,float *w);
void sedov(int N,int Nf,int Np,int ic, float *rho,float *c,float *cu,float *cv,float *cw,float *r,float *P,float *u,float *v,float *w);
float max(float *A, int N);
int wherefrente(float *A, int N);
int wheremax(float *A, int N);

int main(){
  /*Declaro unas constantes importantes*/
  int N;
  N= (int)pow(L/dL,3.0);
  int Nf;
  Nf= (int)(L/dL);
  int Np;
  Np= (int)pow(L/dL,2.0);
  int ic;
  ic= (int)(0.5*(N+Np+Nf));
  /*Declaro las listas*/
  float *rho;
  float *c;
  float *cu;
  float *cv;
  float *cw;
  float *r;
  float *P;
  float *u;
  float *v;
  float *w;
  
  rho=malloc(N*sizeof(float));
  c=malloc(N*sizeof(float));
  cu=malloc(N*sizeof(float));
  cv=malloc(N*sizeof(float));
  cw=malloc(N*sizeof(float));
  r=malloc(N*sizeof(float));
  P=malloc(N*sizeof(float));
  u=malloc(N*sizeof(float));
  v=malloc(N*sizeof(float));
  w=malloc(N*sizeof(float));
 
  /*Ahora las inicilaizo*/
  inicializar(N,Nf,Np,ic,rho,c,cu,cv,cw,r,P,u,v,w);
  /*Ahora hago volumes finitos*/
  sedov(N,Nf,Np,ic,rho,c,cu,cv,cw,r,P,u,v,w);



  return 0;
}



float sgnx(int i,int Nf){
  if(i%Nf>=Nf/2){
    return 1.0;
  }

  else{
    return -1.0;
  }

}



float sgny(int i,int Np){
  
  if(i%Np>=Np/2){
    return 1.0;
  }

  else{
    return -1.0;
  }

}



float sgnz(int i,int N){
  if(i>=N/2){
    return 1.0;
  }

  else{
    return -1.0;
  }

}


void inicializar(int N,int Nf,int Np,int ic, float *rho,float *c,float *cu,float *cv,float *cw,float *r,float *P,float *u,float *v,float *w){
  int i;
  for(i=0;i<N;i++){
    rho[i]=1.0;
    P[i]=1.0;
    c[i]=pow(gamma*P[i]/rho[i],0.5);
    u[i]=0.0;
    v[i]=0.0;
    w[i]=0.0;
    cu[i]=c[i]+u[i];
    cv[i]=c[i]+v[i];
    cw[i]=c[i]+w[i];
  }

  float xc;
  float yc;
  float zc;
  xc=dL*(0.5+(ic%Nf));
  yc=dL*(0.5+(int)((ic%Np)/Nf));
  zc=dL*(0.5+(int)(ic/Np));
  for(i=0;i<N;i++){
    r[i]=pow((xc-(dL*(0.5+(i%Nf))))*(xc-(dL*(0.5+(i%Nf))))+(yc-(dL*(0.5+(int)((i%Np)/Nf))))*(yc-(dL*(0.5+(int)((i%Np)/Nf))))+(zc-(dL*(0.5+(int)(i/Np))))*(zc-(dL*(0.5+(int)(i/Np)))),0.5);
  }
    rho[ic]=1.0;
    P[ic]=(E0*(gamma-1.0)/(dL*dL*dL))/(101000.0);
    c[ic]=pow(gamma*P[ic]/rho[ic],0.5);
    u[ic]=0.0;
    v[ic]=0.0;
    w[ic]=0.0;
    cu[ic]=c[ic]+u[ic];
    cv[ic]=c[ic]+v[ic];
    cw[ic]=c[ic]+w[ic];
    

}



float max(float *A, int N){
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




int wherefrente(float *A, int N){
  int i;
  int f;
  f=0;
  for(i=wheremax(A,N);i<N;i++){
    if(A[i]<0.01 && f==0){
      f=i;
    }
  }
  return f;
}


int wheremax(float *A, int N){
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





void sedov(int N,int Nf,int Np,int ic, float *rho,float *c,float *cu,float *cv,float *cw,float *r,float *P,float *u,float *v,float *w){

  float tiempo;
  float dt;
  float *caux;
  caux=malloc(3*sizeof(float));
  caux[0]=max(cu,N);
  caux[1]=max(cv,N);
  caux[2]=max(cw,N);
  dt=2.0*dL/max(caux,N);
  int cont;
  cont=0;
  int i;
  /*Defino las listas intermedias*/
  float u2ipx;
  float Fx2ipx;
  float Fx3ipx;
  float Fx4ipx;
  float Fx5ipx;
  float u2inx;
  float Fx2inx;
  float Fx3inx;
  float Fx4inx;
  float Fx5inx;
  float u3ipy;
  float Fx3ipy;
  float Fy3ipy;
  float Fy4ipy;
  float Fy5ipy;
  float u3iny;
  float Fx3iny;
  float Fy3iny;
  float Fy4iny;
  float Fy5iny;
  float u4ipz;
  float Fx4ipz;
  float Fy4ipz;
  float Fz4ipz;
  float Fz5ipz;
  float u4inz;
  float Fx4inz;
  float Fy4inz;
  float Fz4inz;
  float Fz5inz;

  float rhoa;
  float ua;
  float va;
  float wa;

  
  
  while(cont<100){
    for(i=0;i<N;i++){

      if(i%Nf!=Nf-1){
	u2ipx=0.5*(rho[i+1]*u[i+1]+rho[i]*u[i]);
	Fx2ipx=0.5*((rho[i+1]*u[i+1]*u[i+1])+P[i+1]+P[i]+(rho[i]*u[i]*u[i]));
	Fx3ipx=0.5*(rho[i+1]*u[i+1]*v[i+1]+rho[i]*u[i]*v[i]);
	Fx4ipx=0.5*(rho[i+1]*u[i+1]*w[i+1]+rho[i]*u[i]*w[i]);
	Fx5ipx=0.5*((u[i+1]*(P[i+1]+(P[i+1]/(gamma-1.0))+0.5*rho[i+1]*(u[i+1]*u[i+1]+v[i+1]*v[i+1]+w[i+1]*w[i+1])))+(u[i]*(P[i]+(P[i]/(gamma-1.0))+0.5*rho[i]*(u[i]*u[i]+v[i]*v[i]+w[i]*w[i]))));
      }
      else {
	u2ipx=rho[i]*u[i];
	Fx2ipx=rho[i]*u[i]*u[i]+P[i];
	Fx3ipx=rho[i]*u[i]*v[i];
	Fx4ipx=rho[i]*u[i]*w[i];
	Fx5ipx=u[i]*(P[i]+(P[i]/(gamma-1.0))+0.5*rho[i]*(u[i]*u[i]+v[i]*v[i]+w[i]*w[i]));
      }

      if(i%Nf!=0){
	u2inx=0.5*(rho[i-1]*u[i-1]+rho[i]*u[i]);
	Fx2inx=0.5*((rho[i-1]*u[i-1]*u[i-1])+P[i-1]+P[i]+(rho[i]*u[i]*u[i]));
	Fx3inx=0.5*(rho[i-1]*u[i-1]*v[i-1]+rho[i]*u[i]*v[i]);
	Fx4inx=0.5*(rho[i-1]*u[i-1]*w[i-1]+rho[i]*u[i]*w[i]);
	Fx5inx=0.5*((u[i-1]*(P[i-1]+(P[i-1]/(gamma-1.0))+0.5*rho[i-1]*(u[i-1]*u[i-1]+v[i-1]*v[i-1]+w[i-1]*w[i-1])))+(u[i]*(P[i]+(P[i]/(gamma-1.0))+0.5*rho[i]*(u[i]*u[i]+v[i]*v[i]+w[i]*w[i]))));
      }
      else{
	u2inx=rho[i]*u[i];
	Fx2inx=rho[i]*u[i]*u[i]+P[i];
	Fx3inx=rho[i]*u[i]*v[i];
	Fx4inx=rho[i]*u[i]*w[i];
	Fx5inx=u[i]*(P[i]+(P[i]/(gamma-1.0))+0.5*rho[i]*(u[i]*u[i]+v[i]*v[i]+w[i]*w[i]));
      }


      if(i%Np<Np-Nf){
	u3ipy=0.5*(rho[i+Nf]*v[i+Nf]+rho[i]*v[i]);
	Fx3ipy=0.5*(rho[i+Nf]*u[i+Nf]*v[i+Nf]+rho[i]*u[i]*v[i]);
	Fy3ipy=0.5*(rho[i+Nf]*v[i+Nf]*v[i+Nf]+P[i+Nf]+rho[i]*v[i]*v[i]+P[i]);
	Fy4ipy=0.5*(rho[i+Nf]*v[i+Nf]*w[i+Nf]+rho[i]*v[i]*w[i]);
       	Fy5ipy=0.5*((v[i+Nf]*(P[i+Nf]+(P[i+Nf]/(gamma-1.0))+0.5*rho[i+Nf]*(u[i+Nf]*u[i+Nf]+v[i+Nf]*v[i+Nf]+w[i+Nf]*w[i+Nf])))+(v[i]*(P[i]+(P[i]/(gamma-1.0))+0.5*rho[i]*(u[i]*u[i]+v[i]*v[i]+w[i]*w[i]))));
      }
      else{
	u3ipy=rho[i]*v[i];
	Fx3ipy=rho[i]*u[i]*v[i];
	Fy3ipy=rho[i]*v[i]*v[i]+P[i];
	Fy4ipy=rho[i]*v[i]*w[i];
       	Fy5ipy=v[i]*(P[i]+(P[i]/(gamma-1.0))+0.5*rho[i]*(u[i]*u[i]+v[i]*v[i]+w[i]*w[i]));
      }


      if(i%Np>=Nf){
       	u3iny=0.5*(rho[i-Nf]*v[i-Nf]+rho[i]*v[i]);
	Fx3iny=0.5*(rho[i-Nf]*u[i-Nf]*v[i-Nf]+rho[i]*u[i]*v[i]);
	Fy3iny=0.5*(rho[i-Nf]*v[i-Nf]*v[i-Nf]+P[i-Nf]+rho[i]*v[i]*v[i]+P[i]);
	Fy4iny=0.5*(rho[i-Nf]*v[i-Nf]*w[i-Nf]+rho[i]*v[i]*w[i]);
       	Fy5iny=0.5*((v[i-Nf]*(P[i-Nf]+(P[i-Nf]/(gamma-1.0))+0.5*rho[i-Nf]*(u[i-Nf]*u[i-Nf]+v[i-Nf]*v[i-Nf]+w[i-Nf]*w[i-Nf])))+(v[i]*(P[i]+(P[i]/(gamma-1.0))+0.5*rho[i]*(u[i]*u[i]+v[i]*v[i]+w[i]*w[i]))));
      }

      else{
	u3iny=rho[i]*v[i];
	Fx3iny=rho[i]*u[i]*v[i];
	Fy3iny=rho[i]*v[i]*v[i]+P[i];
	Fy4iny=rho[i]*v[i]*w[i];
       	Fy5iny=v[i]*(P[i]+(P[i]/(gamma-1.0))+0.5*rho[i]*(u[i]*u[i]+v[i]*v[i]+w[i]*w[i]));
      }


      if(i<N-Np){
	u4ipz=0.5*(rho[i+Np]*w[i+Np]+rho[i]*w[i]);
	Fx4ipz=0.5*(rho[i+Np]*u[i+Np]*w[i+Np]+rho[i]*u[i]*w[i]);
	Fy4ipz=0.5*(rho[i+Np]*v[i+Np]*w[i+Np]+rho[i]*v[i]*w[i]);
	Fz4ipz=0.5*(rho[i+Np]*w[i+Np]*w[i+Np]+P[i+Np]+rho[i]*w[i]*w[i]+P[i]);
	Fz5ipz=0.5*((w[i+Np]*(P[i+Np]+(P[i+Np]/(gamma-1.0))+0.5*rho[i+Np]*(u[i+Np]*u[i+Np]+v[i+Np]*v[i+Np]+w[i+Np]*w[i+Np])))+(w[i]*(P[i]+(P[i]/(gamma-1.0))+0.5*rho[i]*(u[i]*u[i]+v[i]*v[i]+w[i]*w[i]))));
      }

      else{
	u4ipz=rho[i]*w[i];
	Fx4ipz=rho[i]*u[i]*w[i];
	Fy4ipz=rho[i]*v[i]*w[i];
	Fz4ipz=rho[i]*w[i]*w[i]+P[i];
	Fz5ipz=w[i]*(P[i]+(P[i]/(gamma-1.0))+0.5*rho[i]*(u[i]*u[i]+v[i]*v[i]+w[i]*w[i]));
      }

      

      if(i>=Np){
       	u4inz=0.5*(rho[i-Np]*w[i-Np]+rho[i]*w[i]);
	Fx4inz=0.5*(rho[i-Np]*u[i-Np]*w[i-Np]+rho[i]*u[i]*w[i]);
	Fy4inz=0.5*(rho[i-Np]*v[i-Np]*w[i-Np]+rho[i]*v[i]*w[i]);
	Fz4inz=0.5*(rho[i-Np]*w[i-Np]*w[i-Np]+P[i-Np]+rho[i]*w[i]*w[i]+P[i]);
	Fz5inz=0.5*((w[i-Np]*(P[i-Np]+(P[i-Np]/(gamma-1.0))+0.5*rho[i-Np]*(u[i-Np]*u[i-Np]+v[i-Np]*v[i-Np]+w[i-Np]*w[i-Np])))+(w[i]*(P[i]+(P[i]/(gamma-1.0))+0.5*rho[i]*(u[i]*u[i]+v[i]*v[i]+w[i]*w[i]))));
      }

       else{
	u4inz=rho[i]*w[i];
	Fx4inz=rho[i]*u[i]*w[i];
	Fy4inz=rho[i]*v[i]*w[i];
	Fz4inz=rho[i]*w[i]*w[i]+P[i];
	Fz5inz=w[i]*(P[i]+(P[i]/(gamma-1.0))+0.5*rho[i]*(u[i]*u[i]+v[i]*v[i]+w[i]*w[i]));
      }


      rhoa=rho[i];
      ua=u[i];
      va=v[i];
      wa=w[i];
      
      rho[i]=rho[i]+sgnx(i,Nf)*(dt/dL)*(u2inx-u2ipx)+sgny(i,Np)*(dt/dL)*(u3iny-u3ipy)+sgnz(i,N)*(dt/dL)*(u4inz-u4ipz);
      u[i]=(u[i]*rhoa/rho[i])+((sgnx(i,Nf)*(dt/dL)*(Fx2inx-Fx2ipx)+sgny(i,Np)*(dt/dL)*(Fx3iny-Fx3ipy)+sgnz(i,N)*(dt/dL)*(Fx4inz-Fx4ipz))/rho[i]);
      v[i]=(v[i]*rhoa/rho[i])+((sgnx(i,Nf)*(dt/dL)*(Fx3inx-Fx3ipx)+sgny(i,Np)*(dt/dL)*(Fy3iny-Fy3ipy)+sgnz(i,N)*(dt/dL)*(Fy4inz-Fy4ipz))/rho[i]);
      w[i]=(w[i]*rhoa/rho[i])+((sgnx(i,Nf)*(dt/dL)*(Fx4inx-Fx4ipx)+sgny(i,Np)*(dt/dL)*(Fy4iny-Fy4ipy)+sgnz(i,N)*(dt/dL)*(Fz4inz-Fz4ipz))/rho[i]);
      P[i]=((P[i]/(gamma-1.0))+0.5*rhoa*(ua*ua+va*va+wa*wa)+sgnx(i,Nf)*(dt/dL)*(Fx5inx-Fx5ipx)+sgny(i,Np)*(dt/dL)*(Fy5iny-Fy5ipy)+sgnz(i,N)*(dt/dL)*(Fz5inz-Fz5ipz)-rho[i]*0.5*(u[i]*u[i]+v[i]*v[i]+w[i]*w[i]))*(gamma-1.0);
      c[i]=pow(gamma*P[i]/rho[i],0.5);
      cu[i]=c[i]+u[i];
      cv[i]=c[i]+v[i];
      cw[i]=c[i]+w[i];
    }
   
    caux[0]=max(cu,N);
    caux[1]=max(cv,N);
    caux[2]=max(cw,N);
    dt=2.0*dL/max(caux,N);
    cont++;
  }

  /*
  FILE *archivo;
  archivo=fopen("rhosedov.dat","w");
  fclose(archivo);
  archivo=fopen("rhosedov.dat","a");
  for(i=ic;i<ic+0.5*Nf;i++){
    fprintf(archivo,"%f ",rho[i]);
  }
  fclose(archivo);
  */


}






