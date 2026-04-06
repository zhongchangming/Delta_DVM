/*limiter style

*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


#define N 51
#define N1 (N+1)
#define Qx 201
#define imin 1
#define imax (N-1)
#define R 0.5
#define sgn(x) (((x)<0.0?-1.0:1.0))
#define CFL (1000.0e0)
//#define lim(x) ((x)=((x+fabs(x))/(1.0+fabs(x)))

double gamma1; //ratio of specific heat
double omega; //temperature dependence index in HS/VHS/VSS model
double Pr;     //Prandtl number
double mu_ref; //viscosity coefficient in reference state
int cL, cK, cLK,center,xL,xR;  //CL=2 (1D): rest freedom of velocity space; cK: internal degree of freedom; cLK=cL+cK
double lambda; //reference mean-free-path
double rhoL, uL, TL, rhoR, uR, TR; //upstream and downstream variables
double L; // length of the domain
double Vmin, Vmax; //min and max discrete velocities
double dx, dt;
double X[N-1];
double g[N1];       //distribution function
double h[N1];       //distribution function
double dg[N1],dh[N1];

double rho[N1],ux[N1],T[N1],E[N1],MM[N1]; //density, xmoment, and total energy in cell center
double ex[Qx], tp[Qx];       //ex is the discete velocity and tp is the weight of the discrete velocity
double Co_X[Qx],Co_W[Qx]; // Cotes points and weights
double PI;   //pi=3.1415926
double qflux[N1],q1flux[N1],stress[N1];
double rho1[N1], ux1[N1], T1[N1],E1[N1],MM1[N1];
double uxold[N1], Told[N1];

double lam1,lam2,ma1,ma2,T0,Tlimit,cp,error1=0.0,cc,error2=0.0;
int num, m,ss;
double sg[N1],sh[N1];
double s,s1,s2,rhs;
double kk;


void initial();
void datadeal();
void updatemac();
void update();
double tauf(double heat,double density);
double feq(int kx, double RHO, double Ux, double TT);
void Cotes(double Vmin, double Vmax, int np);
void parameter();
double gshakhov(int kx, double RHO, double Ux, double TT, double HQ);
double hshakhov(int kx, double RHO, double Ux, double TT, double HQ);
double derror[10000],dtime[10000];
double aa[N1],gg[N1],hh[N1];

double A[N1][3][3],W1[N1][3],B[N1][3],W2[N1][3];
double A1[N1][3][3],A2[N1][3][3];   //1 is plus; 2 is minus
double r[N1];
double U1[N1],U2[N1],U3[N1];

void Error()
{


  double temp1,temp2;
  int i;
  temp1=0;
  temp2=0;
  for(i=1;i<N;i++)
  {
    temp1+=(T[i]-Told[i])*(T[i]-Told[i])  ;
    temp2+=(qflux[i]*qflux[i]);
  }

  temp1=sqrt(temp1/(N-1));
  temp2=sqrt(temp2);

  error1=temp1;
 //error1 = temp1;


}


double gshakhov(int kx, double RHO, double Ux, double TT, double HQ)
{

  double FM, g_Pr, g_Sh, cq, cn, fac;

  FM=feq(kx,RHO, Ux, TT);
  cq=(ex[kx]-Ux)*HQ/(5*RHO*R*R*TT*TT);
  cn=(ex[kx]-Ux)*(ex[kx]-Ux)/(R*TT)-5;
  fac=(1-Pr)*FM*cq;
  g_Pr=fac*(cn+cL);
  g_Sh=FM+g_Pr;
  return g_Sh;
}


double hshakhov(int kx, double RHO, double Ux, double TT, double HQ)
{
  double FM1, h_Sh, h_Pr, cq, cn, fac;
  FM1=feq(kx,RHO, Ux, TT);
  cq=(ex[kx]-Ux)*HQ/(5*RHO*R*R*TT*TT);
  cn=(ex[kx]-Ux)*(ex[kx]-Ux)/(R*TT)-5;
  fac=(1-Pr)*FM1*cq;
  h_Pr=fac*(cn*cLK+cL*(cLK+2));
  h_Sh=(cLK*FM1+h_Pr)*R*TT;
  return h_Sh;
}

double tauf(double heat,double density)
{
  double mu1;
  mu1=mu_ref*exp(omega*log(heat/TL));
  mu1=mu1/(density*R*heat);
  return mu1;
}


void parameter()
{
  int i;
  printf("the density velocity temperature is\n");
  for(i=0;i<=N;i++)
  {
    printf("%6.3f\t%6.3f\t%6.3f\n",rho[i],ux[i],T[i]);
  }

}

int main()
{
  clock_t start,finish;

  m=1;
  num=1;  // Delta 内迭代数
  ss=0;
  initial();
  start=clock();
  update();
  Error();
 // error2=error1;
  while( m<= 1e10 )
  {
    //updatemac();
    update();


    if(m%1==0)
    {
     error2=error1;
      Error();
      //derror[ss]=error1;

      printf("\nerror1 is %10.8e\n",error1);

      double err=fabs(error1-error2);

      printf("abs(error1-error2) is %10.8e\n",err);
      finish=clock();
    //  dtime[ss]=(double)(finish-start)/CLOCKS_PER_SEC;
    //  datadeal();
      printf("rhs=%e  iterations= %d \n",rhs, m);
      ss++;
    }

    if(error1<=1.0e-6)
    //if(fabs(error1-error2)<=3.0e-8)
    {
        break;
    }


    m++;
  }
  finish=clock();
 // parameter();
  printf("\nrunning nums is %d \n",m);
  printf("error1 is %10.8e\n",error1);
  printf("rhoL is %10.8e, and rhoR is  %10.8e\n",rhoL,rhoR);
  printf("TL   is %10.8e, and TR   is  %10.8e\n",TL,TR);
  datadeal();
  printf("total running time is %e seconds\n",(double)(finish-start)/CLOCKS_PER_SEC);



}


// use the code from Guo laoshi

void Cotes(double Vmin, double Vmax, int np)
{
  int k,n;
  double dh, b, a;

  a=Vmin; b=Vmax;
  n=(np-1)/4;
  dh=(b-a)/n;
//  n1=np-1;

  for(k=0;k<np;k++)
  {
    Co_X[k]=Vmin+k*dh/4;//Co_X[k]=(k-0.5*n1)/n1*(b-a);
  }

  for(k=0;k<n;k++)
  {
    Co_W[4*k]=14.0;
    Co_W[4*k+1]=32.0;
    Co_W[4*k+2]=12.0;
    Co_W[4*k+3]=32.0;
  }
  Co_W[0]=7.0;
  Co_W[np-1]=7.0;

  for(k=0;k<np;k++) Co_W[k]*=dh/90;
}



double feq(int kx, double RHO, double Ux, double TT)
{
  double eu2, x,ee;
  eu2=(ex[kx]-Ux)*(ex[kx]-Ux);
  ee=sqrt(2*PI*R*TT);
  x=tp[kx]*RHO*exp(-eu2/(2*R*TT))/ee;
  return x;
}

void initial()
{


  int i,kx;
  double u_local,Umax;
  PI=4.0*atan(1.0);
  center=N/2;
  Pr =2.0/3;
  cL=2; cK = 0; cLK=cL+cK;
  gamma1 = (cL+ cK+3.0)/(cL+cK+1.0);
  omega = 0.5;   // omega数

  for(i=0;i<100;i++)
  {
    derror[i]=0;
    dtime[i]=0;
  }

   // upstream
   ma1=3.0; //进口马赫数
   rhoL=1.0; TL=1.0;
   uL=ma1*sqrt(gamma1*R*TL);
   xL=imin-1;
   rho[xL] = rhoL;  ux[xL]=uL; T[xL]=TL;
   cp=gamma1/(gamma1-1)*R;

   T0=(cp*TL+0.5*uL*uL)/cp;
   Tlimit=2.0/(gamma1+1)*T0;
   lam1=uL/sqrt(gamma1*R*Tlimit);
   lam2=1.0/lam1;
   uR=lam2*sqrt(gamma1*R*Tlimit);
   TR=(cp*T0-0.5*uR*uR)/cp;
   rhoR=rhoL*uL/uR;
   ma2=uR/sqrt(TR*R*gamma1);



//downstream
   xR=imax+1;
   rho[xR]=rhoR;   ux[xR]=uR; T[xR]=TR;
   /*

   lambda=2*mu_ref*(7-2*omega)*(5-2*omega)/(15*rhoL*sqrt(2*PI*R*TL));

   dx=1000.0*lambda;      //cell size
   L=dx*(imax-imin+1);
   */

   //discrete velocities and weights

   Vmin=-15.0;
   Vmax=15.0;

  //Vmin=(ma1-15.0)*sqrt(gamma1*R);
 // Vmax=(ma1+15.0)*sqrt(gamma1*R); //min and max discrete velocities
  Cotes(Vmin, Vmax, Qx);
  for(kx=0;kx<Qx;kx++)
  {

    ex[kx]=Co_X[kx];
    tp[kx]=Co_W[kx];

  }
   L=1.0;
   kk=0.01;   // Kn数的倒数
   dx=L/(imax-imin+1);
   lambda=dx/kk;
   mu_ref=lambda*(15*rhoL*sqrt(2*PI*R*TL))/(2*(7-2*omega)*(5-2*omega));
   printf("mu_ref=%e\n",mu_ref);
   //initial macroscopic variables
   for(i=0;i<=N;i++)
   {
    if(i<=center)
    {
      rho[i]=rho[xL];
      ux[i]=ux[xL];
      T[i]=T[xL];
    }
    else
    {
      rho[i]=rho[xR];
      ux[i]=ux[xR];
      T[i]=T[xR];
    }
    E[i]=rho[i]*(ux[i]*ux[i]+3*R*T[i])/2;
    MM[i]=rho[i]*ux[i];
  }
  for(i=0;i<=N;i++)
  {
    q1flux[i]=0;
    stress[i]=0;
    qflux[i]=0;
    dg[i]=0;
    dh[i]=0;


  }



  for(kx=0;kx<Qx;kx++)
  {
      for(i=0;i<=N;i++)
      {
          g[i]=feq(kx,rho[i],ux[i],T[i]);
          h[i]=2*R*T[i]*g[i];

      }

    for(i=0;i<=N;i++)
    {
      q1flux[i]+=0.5*(ex[kx]-ux[i])*((ex[kx]-ux[i])*(ex[kx]-ux[i])*g[i]+h[i]);
      stress[i]+=(ex[kx]-ux[i])*(ex[kx]-ux[i])*g[i];
    }
  }
  Umax=0.0;

  for(i=0;i<=N;i++)
  {
    u_local=Vmax+sqrt(gamma1*R*T[i]);
    Umax=((Umax>u_local)?Umax:u_local);
  }
  dt=dx/Umax*CFL;
  for(i=0;i<=N;i++)
  {
    U3[i]=rho[i]*(ux[i]*ux[i]+3.0*R*T[i])/2.0;
    U2[i]=rho[i]*ux[i];
    U1[i]=rho[i];
  }



}

void update()
{


  int kx,i,j;
  for(i=1;i<N;i++)
  {

    qflux[i]=q1flux[i];
    q1flux[i]=0;
    stress[i]=0;
    rho1[i]=0;
    ux1[i]=0;
    T1[i]=0;
    MM1[i]=0;
    E1[i]=0;
    B[i][0]=0;
    B[i][1]=0;
    B[i][2]=0;

    uxold[i] = ux[i];
    Told[i] = T[i];

  }
  for(i=0;i<=N;i++)
  {
    aa[i]=1.0;
  }
  dg[0]=0;
  dh[0]=0;
  dg[N]=0;
  dh[N]=0;
  sg[0]=0;
  sh[0]=0;
  sg[N]=0;
  sh[N]=0;

  for(kx=0;kx<Qx;kx++)
  {
    for(i=0;i<=N;i++)
    {
      gg[i]=gshakhov(kx,rho[i],ux[i],T[i],qflux[i]);
      hh[i]=hshakhov(kx,rho[i],ux[i],T[i],qflux[i]);
      g[i]=gg[i];
      h[i]=hh[i];
    }
    for(j=0;j<num;j++)
    {
      for(i=1;i<N;i++)
      {
        s1=(g[i+1]-g[i])/dx;
        s2=(g[i]-g[i-1])/dx;
        sg[i]=(sgn(s1)+sgn(s2))*fabs(s1)*fabs(s2)/(fabs(s1)+fabs(s2)+1.0e-10);
        s1=(h[i+1]-h[i])/dx;
        s2=(h[i]-h[i-1])/dx;
        sh[i]=(sgn(s1)+sgn(s2))*fabs(s1)*fabs(s2)/(fabs(s1)+fabs(s2)+1.0e-10);
      }

      if(ex[kx]>=0)
      {
        for(i=1;i<N;i++)
        {
          cc=ex[kx]/dx*tauf(T[i],rho[i]);
          dg[i]=((gg[i]-g[i])*aa[i]-cc*(g[i]-g[i-1]+dx/2*(sg[i]-sg[i-1]))+cc*dg[i-1])/(cc+aa[i]);
          dh[i]=((hh[i]-h[i])*aa[i]-cc*(h[i]-h[i-1]+dx/2*(sh[i]-sh[i-1]))+cc*dh[i-1])/(cc+aa[i]);

        }
      }
      if(ex[kx]<0)
      {
        for(i=N-1;i>=1;i--)
        {
          cc=ex[kx]/dx*tauf(T[i],rho[i]);
          dg[i]=((gg[i]-g[i])*aa[i]-cc*(g[i+1]-g[i]-dx/2*(sg[i+1]-sg[i]))-cc*dg[i+1])/(-cc+aa[i]);
          dh[i]=((hh[i]-h[i])*aa[i]-cc*(h[i+1]-h[i]-dx/2*(sh[i+1]-sh[i]))-cc*dh[i+1])/(-cc+aa[i]);

        }
      }
      for(i=1;i<N;i++)
      {
        g[i]+=dg[i];
        h[i]+=dh[i];

      }
    }
    for(i=1;i<N;i++)
    {
      s1=(g[i+1]-g[i])/dx;
      s2=(g[i]-g[i-1])/dx;
      sg[i]=(sgn(s1)+sgn(s2))*fabs(s1)*fabs(s2)/(fabs(s1)+fabs(s2)+1.0e-10);
      s1=(h[i+1]-h[i])/dx;
      s2=(h[i]-h[i-1])/dx;
      sh[i]=(sgn(s1)+sgn(s2))*fabs(s1)*fabs(s2)/(fabs(s1)+fabs(s2)+1.0e-10);
    }
    for(i=1;i<N;i++)
    {

      if(ex[kx]>=0)
      {
        B[i][0]-=ex[kx]*(g[i]-g[i-1]+dx/2*(sg[i]-sg[i-1]));
        B[i][1]-=ex[kx]*(g[i]-g[i-1]+dx/2*(sg[i]-sg[i-1]))*ex[kx];
        B[i][2]-=ex[kx]*0.5*(ex[kx]*ex[kx]*(g[i]-g[i-1]+dx/2*(sg[i]-sg[i-1]))
          +(h[i]-h[i-1]+dx/2*(sh[i]-sh[i-1])));
      }
      if(ex[kx]<0)
      {
        B[i][0]-=ex[kx]*(g[i+1]-g[i]-dx/2*(sg[i+1]-sg[i]));
        B[i][1]-=ex[kx]*(g[i+1]-g[i]-dx/2*(sg[i+1]-sg[i]))*ex[kx];
        B[i][2]-=ex[kx]*0.5*(ex[kx]*ex[kx]*(g[i+1]-g[i]-dx/2*(sg[i+1]-sg[i]))
          +(h[i+1]-h[i]-dx/2*(sh[i+1]-sh[i])));
      }
    }
    for(i=1;i<N;i++)
    {
      rho1[i]+=g[i];
      MM1[i]+=ex[kx]*g[i];
      E1[i]+=0.5*(ex[kx]*ex[kx]*g[i]+h[i]);
    }
    for(i=1;i<N;i++)
    {
      q1flux[i]+=0.5*(ex[kx]-ux[i])*((ex[kx]-ux[i])*(ex[kx]-ux[i])*g[i]+h[i]);
      stress[i]+=(ex[kx]-ux[i])*(ex[kx]-ux[i])*g[i];
    }

  }
  for(i=1;i<N;i++)
  {
    rho[i]=rho1[i];
    ux[i]=MM1[i]/rho1[i];
    T[i]=(E1[i]-MM1[i]*MM1[i]/(rho1[i]*2))/(1.5*R*rho1[i]);
    stress[i]-=rho[i]*R*T[i];
  }


}
void updatemac()
{
  int i,j,m,k;
  double s1,s2,s3,s4,s5,s6;
  rhs=0;
  for(i=1;i<N;i++)
  {
    for(m=0;m<3;m++)
    {
      if(rhs<B[i][m])
      {
        rhs=B[i][m];
      }
    }
  }


  for(i=0;i<=N;i++)
  {
    r[i]=sqrt(gamma1*T[i]*R)+fabs(ux[i]);

  }
  /*
  for(m=0;m<=N;m++)
  {

    A[m][0][0]=0;
    A[m][0][1]=1.0;
    A[m][0][2]=0;
    A[m][1][0]=-(3-gamma1)*U2[m]*U2[m]/(U1[m]*U1[m]*2);
    A[m][1][1]=(3-gamma1)*U2[m]/U1[m];
    A[m][1][2]=gamma1-1;
    A[m][2][0]=-(gamma1-1)/2.0*U2[m]*U2[m]/(U1[m]*U1[m])-U2[m]*U3[m]/(U1[m]*U1[m])*gamma1;
    A[m][2][1]=(gamma1-1)*U2[m]/U1[m]+U3[m]/U1[m]*gamma1;
    A[m][2][2]=gamma1*U2[m]/U1[m];
  }
  for(m=0;m<N1;m++)
  {
    for(i=0;i<3;i++)
    {
      for(j=0;j<3;j++)
      {
        if(i==j)
        {
          A1[m][i][j]=0.5*(A[m][i][j]+r[m]);
          A2[m][i][j]=0.5*(A[m][i][j]-r[m]);

        }
        else
        {
          A1[m][i][j]=0.5*A[m][i][j];
          A2[m][i][j]=0.5*A[m][i][j];

        }
      }
    }
  }

  */


  for(i=1;i<N;i++)
  {
    if(i>1)
    {
      s1=W1[i-1][0]+U1[i-1];
      s2=W1[i-1][1]+U2[i-1];
      s3=W1[i-1][2]+U3[i-1];

      s4=s2-U2[i-1];
      s5=s2*s2/s1+(gamma1-1)*(s3-s2*s2/2/s1);
      s5-=(U2[i-1]*U2[i-1]/U1[i-1]+(gamma1-1)*(U3[i-1]-U2[i-1]*U2[i-1]/2/U1[i-1]));

      s6=s2/s1*(s3+(gamma1-1)*(s3-s2*s2/2/s1));
      s6-=(U2[i-1]/U1[i-1]*((gamma1-1)*(U3[i-1]-U2[i-1]*U2[i-1]/2/U1[i-1])+U3[i-1]));


      B[i][0]+=0.5*(s4+r[i]*W1[i-1][0]);
      B[i][1]+=0.5*(s5+r[i]*W1[i-1][1]);
      B[i][2]+=0.5*(s6+r[i]*W1[i-1][2]);

    }
    for(k=0;k<3;k++)
    {
      W1[i][k]=B[i][k]/(r[i]+dx/dt);
    }
  }

  for(i=(N-1);i>=1;i--)
  {
    for(j=2;j>=0;j--)
    {
      B[i][j]=(r[i]+dx/dt)*W1[i][j];
    }
  }
  for(i=(N-1);i>=1;i--)
  {
    if(i<(N-1))
    {
      s1=W2[i+1][0]+U1[i+1];
      s2=W2[i+1][1]+U2[i+1];
      s3=W2[i+1][2]+U3[i+1];

      s4=s2-U2[i+1];
      s5=s2*s2/s1+(gamma1-1)*(s3-s2*s2/2/s1);
      s5-=(U2[i+1]*U2[i+1]/U1[i+1]+(gamma1-1)*(U3[i+1]-U2[i+1]*U2[i+1]/2/U1[i+1]));

      s6=s2/s1*(s3+(gamma1-1)*(s3-s2*s2/2/s1));
      s6-=(U2[i+1]/U1[i+1]*((gamma1-1)*(U3[i+1]-U2[i+1]*U2[i+1]/2/U1[i+1])+U3[i+1]));

      B[i][0]-=0.5*(s4-r[i]*W2[i+1][0]);
      B[i][1]-=0.5*(s5-r[i]*W2[i+1][1]);
      B[i][2]-=0.5*(s6-r[i]*W2[i+1][2]);


    }
    for(j=2;j>=0;j--)
    {
      W2[i][j]=B[i][j]/(dx/dt+r[i]);
    }
  }
  for(i=1;i<N;i++)
  {
    U1[i]+=W2[i][0];
    U2[i]+=W2[i][1];
    U3[i]+=W2[i][2];
  }


  /*
  for(i=1;i<N;i++)
  {
    if(i>1)
    {
      for(j=0;j<=2;j++)
      {
        B[i][0]+=A1[i-1][0][j]*W1[i-1][j];
        B[i][1]+=A1[i-1][1][j]*W1[i-1][j];
        B[i][2]+=A1[i-1][2][j]*W1[i-1][j];
      }
    }
    for(k=0;k<3;k++)
    {
      W1[i][k]=B[i][k]/(r[i]+dx/dt);
    }
  }


  for(i=(N-1);i>=1;i--)
  {
    for(j=2;j>=0;j--)
    {
      B[i][j]=(r[i]+dx/dt)*W1[i][j];
    }
  }
  for(i=(N-1);i>=1;i--)
  {
    if(i<(N-1))
    {
      for(j=0;j<3;j++)
      {
        B[i][0]-=A2[i+1][0][j]*W2[i+1][j];
        B[i][1]-=A2[i+1][1][j]*W2[i+1][j];
        B[i][2]-=A2[i+1][2][j]*W2[i+1][j];
      }
    }
    for(j=2;j>=0;j--)
    {
      W2[i][j]=B[i][j]/(dx/dt+r[i]);
    }
  }
  for(i=1;i<N;i++)
  {
    U1[i]+=W2[i][0];
    U2[i]+=W2[i][1];
    U3[i]+=W2[i][2];
  }

  */
  for(i=1;i<N;i++)
  {
    rho[i]=U1[i];
    ux[i]=U2[i]/U1[i];
    T[i]=(U3[i]-0.5*U2[i]*U2[i]/U1[i])*(gamma1-1)/(R*U1[i]);
  }

}
void  datadeal()
{
     int i;
     double xmid;
     char data[1000];
     FILE *fp;
     sprintf(data,"ma%.1fdx%.1ft%d.dat",ma1,kk,ss);

     xmid=0.5;
     fp=fopen(data,"w");
    // fp=fopen("xdm.dat","w");
     for (i=imin; i<=imax; i++) fprintf(fp,"%e ", (i-xmid)*dx);
     fprintf(fp,"\n");
     for (i=imin; i<=imax; i++) fprintf(fp,"%e ",rho[i]);
    fprintf(fp,"\n");
     for (i=imin; i<=imax; i++) fprintf(fp,"%e ",ux[i]);
    fprintf(fp,"\n");
     for (i=imin; i<=imax; i++) fprintf(fp,"%e ",T[i]);
     fprintf(fp,"\n");
      for (i=imin; i<=imax; i++) fprintf(fp,"%e ",q1flux[i]);
     fprintf(fp,"\n");
      for (i=imin; i<=imax; i++) fprintf(fp,"%e ",stress[i]);
     fprintf(fp,"\n");
      for (i=0; i<100; i++) fprintf(fp,"%e ",derror[i]);
     fprintf(fp,"\n");
     for (i=0; i<100; i++) fprintf(fp,"%e ",dtime[i]);
     fprintf(fp,"\n");
      fclose(fp);


     fp=fopen("rhodm.dat","w");
     for (i=imin; i<=imax; i++)
     {
         fprintf(fp,"%e\n ",rho[i]);
     }
     fclose(fp);

    fp=fopen("dx.dat","w");
    for (i=imin; i<=imax; i++)
    {
        fprintf(fp,"%e\n ",(i-xmid)*dx);
    }
    fclose(fp);


     fp=fopen("udm.dat","w");
     for (i=imin; i<=imax; i++) fprintf(fp,"%e\n ",ux[i]);
     fclose(fp);

     fp=fopen("Tdm.dat","w");
     for (i=imin; i<=imax; i++) fprintf(fp,"%e\n ",T[i]);
     fclose(fp);
     fp=fopen("qfluxdm.dat","w");
     for (i=imin; i<=imax; i++) fprintf(fp,"%e\n ",q1flux[i]);
     fclose(fp);
     fp=fopen("stressdm.dat","w");
     for (i=imin; i<=imax; i++) fprintf(fp,"%e\n ",stress[i]);
     fclose(fp);
     fp=fopen("derrordm.dat","w");
     for (i=0; i<100; i++) fprintf(fp,"%e\n ",derror[i]);
     fclose(fp);
     fp=fopen("dtimedm.dat","w");
     for (i=0; i<100; i++) fprintf(fp,"%e\n ",dtime[i]);
     fclose(fp);



}

