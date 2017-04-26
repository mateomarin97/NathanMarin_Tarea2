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
void inicializar(int N,int Nf,int Np,int ic, float *rho, float *u2,float *u3,float *u4,float *u5,float *Fx2,float *Fx3,float *Fx4,float *Fx5,float *Fy3,float *Fy4,float *Fy5,float *Fz4,float *Fz5,float *c,float *cu,float *cv,float *cw,float *x,float *y,float *z,float *r,float *P,float *u,float *v,float *w);
void sedov(int N,int Nf,int Np,int ic, float *rho, float *u2,float *u3,float *u4,float *u5,float *Fx2,float *Fx3,float *Fx4,float *Fx5,float *Fy3,float *Fy4,float *Fy5,float *Fz4,float *Fz5,float *c,float *cu,float *cv,float *cw,float *x,float *y,float *z,float *r,float *P,float *u,float *v,float *w);
float max(float *A);
int wherefrente(float *A);
int wheremax(float *A);

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
  float *u2;
  float *u3;
  float *u4;
  float *u5;
  float *Fx2;
  float *Fx3;
  float *Fx4;
  float *Fx5;
  float *Fy3;
  float *Fy4;
  float *Fy5;
  float *Fz4;
  float *Fz5;
  float *c;
  float *cu;
  float *cv;
  float *cw;
  float *x;
  float *y;
  float *z;
  float *r;
  float *P;
  float *u;
  float *v;
  float *w;
  
  rho=malloc(N*sizeof(float));
  u2=malloc(N*sizeof(float));
  u3=malloc(N*sizeof(float));
  u4=malloc(N*sizeof(float));
  u5=malloc(N*sizeof(float));
  Fx2=malloc(N*sizeof(float));
  Fx3=malloc(N*sizeof(float));
  Fx4=malloc(N*sizeof(float));
  Fx5=malloc(N*sizeof(float));
  Fy3=malloc(N*sizeof(float));
  Fy4=malloc(N*sizeof(float));
  Fy5=malloc(N*sizeof(float));
  Fz4=malloc(N*sizeof(float));
  Fz5=malloc(N*sizeof(float));
  c=malloc(N*sizeof(float));
  cu=malloc(N*sizeof(float));
  cv=malloc(N*sizeof(float));
  cw=malloc(N*sizeof(float));
  x=malloc(N*sizeof(float));
  y=malloc(N*sizeof(float));
  z=malloc(N*sizeof(float));
  r=malloc(N*sizeof(float));
  P=malloc(N*sizeof(float));
  u=malloc(N*sizeof(float));
  v=malloc(N*sizeof(float));
  w=malloc(N*sizeof(float));
 
  /*Ahora las inicilaizo*/
  inicializar(N,Nf,Np,ic,rho,u2,u3,u4,u5,Fx2,Fx3,Fx4,Fx5,Fy3,Fy4,Fy5,Fz4,Fz5,c,cu,cv,cw,x,y,z,r,P,u,v,w);
  /*Ahora hago volumes finitos*/
  sedov(N,Nf,Np,ic,rho,u2,u3,u4,u5,Fx2,Fx3,Fx4,Fx5,Fy3,Fy4,Fy5,Fz4,Fz5,c,cu,cv,cw,x,y,z,r,P,u,v,w);



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
    return -1.0
  }

}


void inicializar(int N,int Nf,int Np,int ic, float *rho, float *u2,float *u3,float *u4,float *u5,float *Fx2,float *Fx3,float *Fx4,float *Fx5,float *Fy3,float *Fy4,float *Fy5,float *Fz4,float *Fz5,float *c,float *cu,float *cv,float *cw,float *x,float *y,float *z,float *r,float *P,float *u,float *v,float *w){
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
    x[i]=dl*(0.5+(i%Nf));
    y[i]=dl*(0.5+(int)((i%Np)/Nf));
    z[i]=dl*(0.5+(int)(i/Np));
    u2[i]=rho[i]*u[i];
    u3[i]=rho[i]*v[i];
    u4[i]=rho[i]*w[i];
    u5[i]=(P[i]/(gamma-1.0))+0.5*rho[i]*(u[i]*u[i]+v[i]*v[i]+w[i]*w[i]);
    Fx2[i]=u2[i]*u[i]+P[i];
    Fx3[i]=u2[i]*v[i];
    Fx4[i]=u2[i]*w[i];
    Fx5[i]=u[i]*(P[i]+u5[i]);
    Fy3[i]=u3[i]*v[i]+P[i];
    Fy4[i]=u3[i]*w[i];
    Fy5[i]=v[i]*(P[i]+u5[i]);
    Fz4[i]=u4[i]*w[i]+P[i];
    Fz5[i]=w[i]*(P[i]+u5[i]);
  }

  float xc;
  float yc;
  float zc;
  xc=x[ic];
  yc=y[ic];
  zc=z[ic];
  for(i=0;i<N;i++){
    r[i]=pow((xc-x[i])*(xc-x[i])+(yc-y[i])*(yc-y[i])+(zc-z[i])*(zc-z[i]),0.5);
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
    x[ic]=dl*(0.5+(ic%Nf));
    y[ic]=dl*(0.5+(int)((ic%Np)/Nf));
    z[ic]=dl*(0.5+(int)(ic/Np));
    u2[ic]=rho[ic]*u[ic];
    u3[ic]=rho[ic]*v[ic];
    u4[ic]=rho[ic]*w[ic];
    u5[ic]=(P[ic]/(gamma-1.0))+0.5*rho[ic]*(u[ic]*u[ic]+v[ic]*v[ic]+w[ic]*w[ic]);
    Fx2[ic]=u2[ic]*u[ic]+P[ic];
    Fx3[ic]=u2[ic]*v[ic];
    Fx4[ic]=u2[ic]*w[ic];
    Fx5[ic]=u[ic]*(P[ic]+u5[ic]);
    Fy3[ic]=u3[ic]*v[ic]+P[ic];
    Fy4[ic]=u3[ic]*w[ic];
    Fy5[ic]=v[ic]*(P[ic]+u5[ic]);
    Fz4[ic]=u4[ic]*w[ic]+P[ic];
    Fz5[ic]=w[ic]*(P[ic]+u5[ic]);

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





void sedov(int N,int Nf,int Np,int ic, float *rho, float *u2,float *u3,float *u4,float *u5,float *Fx2,float *Fx3,float *Fx4,float *Fx5,float *Fy3,float *Fy4,float *Fy5,float *Fz4,float *Fz5,float *c,float *cu,float *cv,float *cw,float *x,float *y,float *z,float *r,float *P,float *u,float *v,float *w){

  float tiempo;
  float dt;
  float *caux;
  caux=malloc(3.0*sizeof(float));
  caux[0]=max(cu);
  caux[1]=max(cv);
  caux[2]=max(cw);
  dt=2.0*dL/max(caux);
  int cont;
  cont=0;
  int i;
  /*Defino las listas intermedias*/
  float *u2ipx;
  float *Fx2ipx;
  float *Fx3ipx;
  float *Fx4ipx;
  float *Fx5ipx;
  float *u2inx;
  float *Fx2inx;
  float *Fx3inx;
  float *Fx4inx;
  float *Fx5inx;
  float *u3ipy;
  float *Fx3ipy;
  float *Fy3ipy;
  float *Fy4ipy;
  float *Fy5ipy;
  float *u3iny;
  float *Fx3iny;
  float *Fy3iny;
  float *Fy4iny;
  float *Fy5iny;
  float *u4ipz;
  float *Fx4ipz;
  float *Fy4ipz;
  float *Fz4ipz;
  float *Fz5ipz;
  float *u4inz;
  float *Fx4inz;
  float *Fy4inz;
  float *Fz4inz;
  float *Fz5inz;


  u2ipx=malloc(N*sizeof(float));
  Fx2ipx=malloc(N*sizeof(float));
  Fx3ipx=malloc(N*sizeof(float));
  Fx4ipx=malloc(N*sizeof(float));
  Fx5ipx=malloc(N*sizeof(float));
  u2inx=malloc(N*sizeof(float));
  Fx2inx=malloc(N*sizeof(float));
  Fx3inx=malloc(N*sizeof(float));
  Fx4inx=malloc(N*sizeof(float));
  Fx5inx=malloc(N*sizeof(float));
  u3ipy=malloc(N*sizeof(float));
  Fx3ipy=malloc(N*sizeof(float));
  Fy3ipy=malloc(N*sizeof(float));
  Fy4ipy=malloc(N*sizeof(float));
  Fy5ipy=malloc(N*sizeof(float));
  u3iny=malloc(N*sizeof(float));
  Fx3iny=malloc(N*sizeof(float));
  Fy3iny=malloc(N*sizeof(float));
  Fy4iny=malloc(N*sizeof(float));
  Fy5iny=malloc(N*sizeof(float));
  u4ipz=malloc(N*sizeof(float));
  Fx4ipz=malloc(N*sizeof(float));
  Fy4ipz=malloc(N*sizeof(float));
  Fz4ipz=malloc(N*sizeof(float));
  Fz5ipz=malloc(N*sizeof(float));
  u4inz=malloc(N*sizeof(float));
  Fx4inz=malloc(N*sizeof(float));
  Fy4inz=malloc(N*sizeof(float));
  Fz4inz=malloc(N*sizeof(float));
  Fz5inz=malloc(N*sizeof(float));




  
  
  while(cont<100){
    for(i=0;i<N;i++){

      if(i%Nf!=Nf-1){
	u2ipx[i]=0.5*(u2[i+1]+u2[i]);
	Fx2ipx[i]=0.5*(Fx2[i+1]+Fx2[i]);
	Fx3ipx[i]=0.5*(Fx3[i+1]+Fx3[i]);
	Fx4ipx[i]=0.5*(Fx4[i+1]+Fx4[i]);
	Fx5ipx[i]=0.5*(Fx5[i+1]+Fx5[i]);
      }
      else {
	u2ipx[i]=u2[i];
	Fx2ipx[i]=Fx2[i];
	Fx3ipx[i]=Fx3[i];
	Fx4ipx[i]=Fx4[i];
	Fx5ipx[i]=Fx5[i];
      }

      if(i%Nf!=0){
	u2inx[i]=0.5*(u2[i-1]+u2[i]);
	Fx2inx[i]=0.5*(Fx2[i-1]+Fx2[i]);
	Fx3inx[i]=0.5*(Fx3[i-1]+Fx3[i]);
	Fx4inx[i]=0.5*(Fx4[i-1]+Fx4[i]);
	Fx5inx[i]=0.5*(Fx5[i-1]+Fx5[i]);
      }
      else{
	u2inx[i]=u2[i];
	Fx2inx[i]=Fx2[i];
	Fx3inx[i]=Fx3[i];
	Fx4inx[i]=Fx4[i];
	Fx5inx[i]=Fx5[i];
      }


      if(i%Np<Np-Nf){
	u3ipy[i]=0.5*(u3[i+Nf]+u3[i]);
	Fx3ipy[i]=0.5*(Fx3[i+Nf]+Fx3[i]);
	Fy3ipy[i]=0.5*(Fy3[i+Nf]+Fy3[i]);
	Fy4ipy[i]=0.5*(Fy4[i+Nf]+Fy4[i]);
       	Fy5ipy[i]=0.5*(Fy5[i+Nf]+Fy5[i]);
      }
      else{
	u3ipy[i]=u3[i];
	Fx3ipy[i]=Fx3[i];
	Fy3ipy[i]=Fy3[i];
	Fy4ipy[i]=Fy4[i];
       	Fy5ipy[i]=Fy5[i];
      }


      if(i%Np>=Nf){
       	u3iny[i]=0.5*(u3[i-Nf]+u3[i]);
	Fx3iny[i]=0.5*(Fx3[i-Nf]+Fx3[i]);
	Fy3iny[i]=0.5*(Fy3[i-Nf]+Fy3[i]);
	Fy4iny[i]=0.5*(Fy4[i-Nf]+Fy4[i]);
       	Fy5iny[i]=0.5*(Fy5[i-Nf]+Fy5[i]);
      }

      else{
	u3iny[i]=u3[i];
	Fx3iny[i]=Fx3[i];
	Fy3iny[i]=Fy3[i];
	Fy4iny[i]=Fy4[i];
       	Fy5iny[i]=Fy5[i];
      }


      if(i<N-Np){
	u4ipz[i]=0.5*(u4[i+Np]+u4[i]);
	Fx4ipz[i]=0.5*(Fx4[i+Np]+Fx4[i]);
	Fy4ipz[i]=0.5*(Fy4[i+Np]+Fy4[i]);
	Fz4ipz[i]=0.5*(Fz4[i+Np]+Fz4[i]);
	Fz5ipz[i]=0.5*(Fz5[i+Np]+Fz5[i]);
      }

      else{
	u4ipz[i]=u4[i];
	Fx4ipz[i]=Fx4[i];
	Fy4ipz[i]=Fy4[i];
	Fz4ipz[i]=Fz4[i];
	Fz5ipz[i]=Fz5[i];
      }

      

      if(i>=Np){
       	u4inz[i]=0.5*(u4[i-Np]+u4[i]);
	Fx4inz[i]=0.5*(Fx4[i-Np]+Fx4[i]);
	Fy4inz[i]=0.5*(Fy4[i-Np]+Fy4[i]);
	Fz4inz[i]=0.5*(Fz4[i-Np]+Fz4[i]);
	Fz5inz[i]=0.5*(Fz5[i-Np]+Fz5[i]);
      }

       else{
	u4inz[i]=u4[i];
	Fx4inz[i]=Fx4[i];
	Fy4inz[i]=Fy4[i];
	Fz4inz[i]=Fz4[i];
	Fz5inz[i]=Fz5[i];
      }

      rho[i]=rho[i]+sgnx(i,Nf)*(dt/dL)*(u2inx[i]-u2ipx[i])+sgny(i,Np)*(dt/dL)*(u3iny[i]-u3ipy[i])+sgnz(i,N)*(dt/dL)*(u4inz[i]-u4ipz[i]);
      u2[i]=u2[i]+sgnx(i,Nf)*(dt/dL)*(Fx2inx[i]-Fx2ipx[i])+sgny(i,Np)*(dt/dL)*(Fx3iny[i]-Fx3ipy[i])+sgnz(i,N)*(dt/dL)*(Fx4inz[i]-Fx4ipz[i]);
      u3[i]=u3[i]+sgnx(i,Nf)*(dt/dL)*(Fx3inx[i]-Fx3ipx[i])+sgny(i,Np)*(dt/dL)*(Fy3iny[i]-Fy3ipy[i])+sgnz(i,N)*(dt/dL)*(Fy4inz[i]-Fy4ipz[i]);
      u4[i]=u4[i]+sgnx(i,Nf)*(dt/dL)*(Fx4inx[i]-Fx4ipx[i])+sgny(i,Np)*(dt/dL)*(Fy4iny[i]-Fy4ipy[i])+sgnz(i,N)*(dt/dL)*(Fz4inz[i]-Fz4ipz[i]);
      u5[i]=u5[i]+sgnx(i,Nf)*(dt/dL)*(Fx5inx[i]-Fx5ipx[i])+sgny(i,Np)*(dt/dL)*(Fy5iny[i]-Fy5ipy[i])+sgnz(i,N)*(dt/dL)*(Fz5inz[i]-Fz5ipz[i]);

      u[i]=u2[i]/rho[i];
      v[i]=u3[i]/rho[i];
      w[i]=u4[i]/rho[i];
      P[i]=(gamma-1.0)*(u5[i]-0.5*rho[i]*(u[i]*u[i]+v[i]*v[i]+w[i]*w[i]));
      c[i]=pow(gamma*P[i]/rho[i],0.5);
      cu[i]=c[i]+u[i];
      cv[i]=c[i]+v[i];
      cw[i]=c[i]+w[i];
      Fx2[i]=u2[i]*u[i]+P[i];
      Fx3[i]=u2[i]*v[i];
      Fx4[i]=u2[i]*w[i];
      Fx5[i]=u[i]*(P[i]+u5[i]);
      Fy3[i]=u3[i]*v[i]+P[i];
      Fy4[i]=u3[i]*w[i];
      Fy5[i]=v[i]*(P[i]+u5[i]);
      Fz4[i]=u4[i]*w[i]+P[i];
      Fz5[i]=w[i]*(P[i]+u5[i]);
    }
   
    caux[0]=max(cu);
    caux[1]=max(cv);
    caux[2]=max(cw);
    dt=2.0*dL/max(caux);
    cont++;
  }
  FILE *archivo;
  archivo=fopen("rhosedov.dat","w");
  fclose(archivo);
  archivo=fopen("rhosedov.dat","a");
  for(i=ic;i<ic+0.5*Nf;i++){
    fprintf(archivo,"%f ",rho[i]);
  }
  fclose(archivo);
  


}






