//---------------------------------------------------------------------------
#include <vcl.h>
#pragma hdrstop
#include <iostream.h>
#include <strstrea.h>
#include <fstream.h>
#include <math.h>
#include <stdio.h>
#include <conio.h>
#include "Unit1.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"
TForm1 *Form1;


// ****************
// INPUT PARAMETERS
// ****************

double const
a=1.0,

// FPU potential parameters
k2=1.0,
k3=-0.0,
k4=3.0,

AtomMass=1.0,            // Mass of an atom

// Strain components for unit cell
exx=0.0;

int const
Nx = 2048,   // number of primitive superstructure cells along X
Nx2=Nx/2,
N=Nx+1,

NumberOfTimeStepsMD = 400000000; // Number of MD time steps

const char
ReadDisplFrom[] = "Displ_a.txt", // File to read the Initial Conditions
WriteDisplTo[]  = "Displ_b.txt"; // File to write the Initial Conditions

double const
Magnif = 1.0, // Magnification for displacements
MagnifDispl=4000.0,
Zoom = 3.0;

// Name of file for output data
const char FileName[] ="Output.txt";

// *******************
// ACCURACY PARAMETERS
// *******************
const int
NumberOfNeighbors=1;                 // Number of neighbours  Lennard-Jones

double const
tau = 0.001, t2 = tau*tau,
pi=3.1415926535897932385;

// *******************
// Variables *********
// *******************
int nit;

double
SKappax,t,
w1x,a1x,A11,A12,A22,LocalizParam,
F,Fp,FxD,EEE,Rx,Rx2,RLen,RLen2,RLen3,
shift,kappax,Skapx,M11,TotEn,
Pot,Kin,Volume,Vol,Cosinus,Sinus,
MaxForce,r,fric,FFx,mult,sxx,p1;

// *******************
// Arrays ************
// *******************
double
Dx[Nx],Dxp[Nx],Dxpp[Nx],Dxppp[Nx],Fx[Nx],Tot[Nx];

// Arrays
double *fun,*funinv,*fi1,*fi2,*fi3,*gRe,*gIm,*GRe,*GIm,*gAmp;

// *******************
//Functions **********
// *******************

double sqr (double r) {return r*r;}
double FPU(double r);
double DFPU(double r);
double DDFPU(double r);
void Geometry();
void ForceMD();
void Stress();
void Modulus();
void PotEn();
void KinEn();
void TotEnAtoms();
void Localization();
void ZeroInitialCond();
void ReadDisplacements();
void WriteDisplacements();
void ZeroIntoForces();
void PrintForces(int nx, int TopMargin, int LeftMargine);
double MaxF();
void TimeStepMD();
void ScreenView(int Freq);
void ShowAtoms();
void FourierTransform (double x[]);
void FFT1D (int n, bool inverse, double x[], double y[]);
void FourierInverceTransform (double x1[], double x2[]);

//---------------------------------------------------------------------------
__fastcall TForm1::TForm1(TComponent* Owner) : TForm(Owner)
{
}
//---------------------------------------------------------------------------

void __fastcall TForm1::Button1Click(TObject *Sender){

fun = new double [Nx*8]; funinv = new double [Nx*8];
fi1= new double [(Nx2+1)*8];
fi2= new double [(Nx2+1)*8];
fi3= new double [(Nx2+1)*8];
gRe = new double [Nx*8]; gIm = new double [Nx*8];
GRe = new double [Nx*8]; GIm = new double [Nx*8];
gAmp = new double [Nx*8];

ofstream OFSO("Dx.txt",ios::out);
ofstream OFSOe("Energy.txt",ios::out);
OFSOe.precision(10);
ofstream OFSOea("EnergyAver.txt",ios::out);
OFSOea.precision(10);
ofstream OFSOp("Potential.txt",ios::out);
OFSOe.precision(10);
ofstream OFSOf("Fourier.txt",ios::out);
OFSOf.precision(10);
ofstream OFSOd("Displ.txt",ios::out);
OFSOd.precision(10);
ofstream OFSOdbL("DB_num_ener.txt",ios::out);
OFSOdbL.precision(10);

// Description of geometry of a unit cell and supercell
Geometry();

double R=0.1;
while (R<3.0){
   OFSOp << R;        OFSOp << "   "; // 1
   OFSOp << FPU(R);   OFSOp << "   "; // 2
   OFSOp << DFPU(R);  OFSOp << "  ";  // 3
   OFSOp << DDFPU(R); OFSOp << "\n";  // 4
   R+=0.01;
}

Form1->Canvas->TextOut(10,10,"k2 = ");
Form1->Canvas->TextOut(70,10,k2);
Form1->Canvas->TextOut(10,25,"k3=");
Form1->Canvas->TextOut(70,25,k3);
Form1->Canvas->TextOut(10,40,"k4 = ");
Form1->Canvas->TextOut(70,40,k4);

Form1->Canvas->TextOut(10,55,"a = ");
Form1->Canvas->TextOut(70,55,a);
Form1->Canvas->TextOut(10,70,"w1x=");
Form1->Canvas->TextOut(70,70,w1x);

// Initial conditions (Zeroes)
ZeroInitialCond();

ShowAtoms();

// Read Initial Conditions from file
//ReadDisplacements();

Stress();
Form1->Canvas->TextOut(10,85,"sxx =");
Form1->Canvas->TextOut(70,85,sxx);

Modulus();
Form1->Canvas->TextOut(210,100,"M11 =");
Form1->Canvas->TextOut(270,100,M11);

//************************************************
//***** MD for the whole computational cell ******
//************************************************

// Staggered mode
for (int i=0; i<Nx; i++){
   double A=0.07;
   if (i%2==0) Dx[i]+=A;
   else        Dx[i]-=A;
   Dxp[i]=Dx[i]; Dxpp[i]=Dx[i]; Dxppp[i]=Dx[i];
}

/*
int Kcos=Nx2-0;
int Ksin=Nx2-0;
double Wcos=2.0*sqrt(k2)*sin(pi*Kcos/(1.0*Nx));
double Wsin=2.0*sqrt(k2)*sin(pi*Ksin/(1.0*Nx));
double Acos=0.01;
double Asin=0*0.005;
double FIcos=0*random(6283)/1000.0;
double FIsin=0*random(6283)/1000.0;
for (int i=0; i<Nx; i++){
    Dx[i]=Acos*cos(2.0*pi*i*Kcos/(1.0*Nx)+FIcos)+Asin*sin(2.0*pi*(i)*Ksin/(1.0*Nx)+FIsin);
    Dxp[i]=Dx[i]; Dxpp[i]=Dx[i]; Dxppp[i]=Dx[i];
}

// Adjust energy
double WantedEnergy=0.005;
PotEn();
for (int i=0; i<Nx; i++){
    Dx[i]=Dx[i]*sqrt(WantedEnergy/Pot);
    Dxp[i]=Dx[i]; Dxpp[i]=Dx[i]; Dxppp[i]=Dx[i];
}


Form1->Canvas->TextOut(410,40,"Ecos = ");
Form1->Canvas->TextOut(470,40,Acos*Acos*Wcos*Wcos);
Form1->Canvas->TextOut(410,55,"Wcos = ");
Form1->Canvas->TextOut(470,55,Wcos);
Form1->Canvas->TextOut(410,70,"Wsin=");
Form1->Canvas->TextOut(470,70,Wsin);
Form1->Canvas->TextOut(410,85,"Tcos = ");
Form1->Canvas->TextOut(470,85,2.0*pi/Wcos);
Form1->Canvas->TextOut(410,100,"Tsin=");
Form1->Canvas->TextOut(470,100,2.0*pi/Wsin);

FourierTransform(Dx);
*/
for (int i=0; i<Nx; i++){
  OFSO << i;      OFSO << "  "; // 1
  OFSO << Dx[i]; OFSO << "\n"; // 2
}

for (int i=0; i<=Nx2; i++){
  OFSOf << i;      OFSOf << "  "; // 1
  OFSOf << fi1[i]; OFSOf << "  "; // 2
  OFSOf << fi2[i]; OFSOf << "  "; // 3
  OFSOf << fi3[i]; OFSOf << "\n"; // 4
}


PotEn();
Form1->Canvas->TextOut(10,100," Pot = ");
Form1->Canvas->TextOut(70,100,(ceil(Pot*1000000)/1000000.0));

double KinAv=0.0, PotAv=0.0, sxxAv=0.0, M11Av=0.0, LocalizParamAv=0.0;
int CountAv=0;
int NumAver;
NumAver=random(10000)+5000;

t=0.0;
for (nit=0; nit<NumberOfTimeStepsMD; nit++){
  // Zeroes into forces
  ZeroIntoForces();
  // Calculation of Short range forces
  ForceMD();

  // Calculations of displacements at t=t+tau
  TimeStepMD();

  if ( fmod(nit,15000) == 0 ){
     Form1->Canvas->TextOut(5,125,"nit =                                             ");
     Form1->Canvas->TextOut(130,125,nit);
     ShowAtoms();
  }

  if ( fmod(nit,10) == 0 ){
     KinEn();
     PotEn();
     Stress();
     Localization();
     Modulus();
//     OFSOe << t;            OFSOe << "   "; // 1
//     OFSOe << Kin;          OFSOe << "   "; // 2
//     OFSOe << Pot;          OFSOe << "   "; // 3
//     OFSOe << Kin+Pot;      OFSOe << "  ";  // 4
//     OFSOe << sxx;          OFSOe << "  ";  // 5
//     OFSOe << M11;          OFSOe << "  ";  // 6
//     OFSOe << LocalizParam; OFSOe << "\n";  // 7
       if ( fmod(nit,15000) == 0 ){
          Form1->Canvas->TextOut(210,10,"Pot = ");
          Form1->Canvas->TextOut(270,10,Pot);
          Form1->Canvas->TextOut(210,25,"Kin=");
          Form1->Canvas->TextOut(270,25,Kin);
          Form1->Canvas->TextOut(210,40,"Tot = ");
          Form1->Canvas->TextOut(270,40,Kin+Pot);
       }
     if (CountAv<NumAver){
        KinAv+=Kin; PotAv+=Pot; sxxAv+=sxx; M11Av+=M11; LocalizParamAv+=LocalizParam;
        CountAv+=1;
     }
     else {
        KinAv=KinAv/CountAv;
        PotAv=PotAv/CountAv;
        sxxAv=sxxAv/CountAv;
        M11Av=M11Av/CountAv;
        LocalizParamAv=LocalizParamAv/CountAv;
        OFSOea << t;              OFSOea << "  "; // 1
        OFSOea << KinAv;          OFSOea << "  "; // 2
        OFSOea << PotAv;          OFSOea << "  "; // 3
        OFSOea << sxxAv;          OFSOea << "  "; // 4
        OFSOea << M11Av;          OFSOea << "  "; // 5
        OFSOea << LocalizParamAv; OFSOea << "\n"; // 6
        KinAv=0.0, PotAv=0.0, sxxAv=0.0, M11Av=0.0, LocalizParamAv=0.0;
        CountAv=0;
        NumAver=random(10000)+5000;
     }
  }

// Count DB
  if ( fmod(nit,10000) == 0 ){
        PotEn();
        KinEn();
        TotEn=Kin+Pot;
        TotEnAtoms();
        int DBN=0;
        double DBen[1000];
        for(int i=0; i<1000; i++) DBen[i]=0.0;
        for(int i=10; i<Nx-10; i++){
           if ((Tot[i]>Tot[i-1])&&(Tot[i]>Tot[i+1])){
              if (Tot[i]>10.0*TotEn){
                 // count DB energy
                 DBen[DBN]=Tot[i];
                 int k=1;
                 while (Tot[i+k-1]>Tot[i+k]){
                    DBen[DBN]+=Tot[i+k];
                    k+=1;
                 }
                 k=-1;
                 while (Tot[i+k+1]>Tot[i+k]){
                    DBen[DBN]+=Tot[i+k];
                    k-=1;
                 }
                 DBN+=1;
              }
           }
        }// for(int i=1; i<Nat-1; i++)
        // find average DB energy
        double DBenerAv=0.0;
        for(int i=0; i<DBN; i++) DBenerAv+=DBen[i];
        DBenerAv=DBenerAv/(DBN+1.0);
        OFSOdbL.width(12); OFSOdbL << t;  OFSOdbL << " ";  // 1
        OFSOdbL.width(12); OFSOdbL << DBN;  OFSOdbL << " ";  // 2
        OFSOdbL.width(12); OFSOdbL << DBenerAv; OFSOdbL << "\n"; // 3
  }
//  if (nit>100){
//     if ((Dxp[0]>Dx[0])&&(Dxp[0]>Dxpp[0])){
//        Form1->Canvas->TextOut(410,115,"Tcos=");
//        Form1->Canvas->TextOut(470,115,t);
//     }
//  }

//  if ( fmod(nit,5) == 0 ){
//     OFSOd << t;     OFSOd << "  ";    // 1
//     OFSOd << Dx[0]; OFSOd << "\n";     // 5
//  }

  t+=tau;
}// for (nit=0; nit<NumberOfTimeStepsMD; nit++)

OFSO.close();


} // for void __fastcall TForm1::Button1Click(TObject *Sender)
//---------------------------------------------------------------------------


//************************************
//   PROCEDURES
//************************************

// Morse potential
double FPU(double r){
     double A=sqr(r-a);
     return 0.5*k2*A+k3*A*(r-a)/3.0+0.25*k4*A*A;
}

// Derivative of Morse potential
double DFPU(double r){
     double A=sqr(r-a);
     return k2*(r-a)+k3*A+k4*A*(r-a);
}

// Second derivative of Morse potential
double DDFPU(double r){
     return k2+2.0*k3*(r-a)+3.0*k4*(r-a)*(r-a);
}

void Geometry(){
    // Direct lattice basis rotated to make w1y=0
    w1x= (1.0 + exx)*a;
    // Volume of lattice cell and SUPERCELL
    Vol=w1x;
    Volume=w1x*Nx;
    // Direct basis for SUPERCELL
    a1x=w1x;
}// for void Geometry()


void ShowAtoms(){
   TotEnAtoms();
   for (int i=0; i<Nx; i++){
      double rx = i*a1x;
      int x = 30  + floor(rx*Zoom);
//      int y = 400 + MagnifDispl*Dx[i];
      int y = 400 - MagnifDispl*Tot[i];
      int size;
      size=4;
      Form1->Canvas->Brush->Color = clGreen;
      Form1->Canvas->Ellipse(x-size,y-size,x+size,y+size);
   }
   Form1->Canvas->Brush->Color = clBtnFace;
} // void ShowAtoms()

void ZeroInitialCond(){
    for (int i=0; i<Nx; i++){
       Dx[i]=0.0; Dxp[i]=0.0; Dxpp[i]=0.0; Dxppp[i]=0.0;
    }
}// for void ZeroInitialCond()

void ZeroIntoForces(){
    for (int i=0; i<Nx; i++){
       Fx[i]=0.0;
    }
}// for void ZeroIntoForces()

void PotEn(){
  Pot = 0.0;
  for (int i=0; i<Nx-1; i++){
      Pot += FPU(a1x+Dx[i+1]-Dx[i]);
  } // for i
  Pot += FPU(a1x+Dx[0]-Dx[Nx-1]);
  Pot=Pot/Nx;
}// for void PotEn();


void KinEn(){
  double Vx;
  Kin=0.0;
  for (int m=0; m<Nx; m++){
      Vx=(11.0*Dx[m]-18.0*Dxp[m]+9.0*Dxpp[m]-2.0*Dxppp[m])/(6.0*tau);
      Kin+=Vx*Vx;
  } //  for i
  Kin=0.5*AtomMass*Kin/Nx;
} // void KinEn()

void TotEnAtoms(){
  double Vx;
  for (int i=0; i<Nx; i++){
      Tot[i] = 0.0;
  } // for i
  for (int i=1; i<Nx-1; i++){
      Tot[i] = 0.5*FPU(a1x+Dx[i]-Dx[i-1]) + 0.5*FPU(a1x+Dx[i+1]-Dx[i]);
      Vx=(11.0*Dx[i]-18.0*Dxp[i]+9.0*Dxpp[i]-2.0*Dxppp[i])/(6.0*tau);
      Tot[i]+=0.5*AtomMass*Vx*Vx;
  } // for i
  Tot[0] = 0.5*FPU(a1x+Dx[0]-Dx[Nx-1]) + 0.5*FPU(a1x+Dx[1]-Dx[0]);
  Vx=(11.0*Dx[0]-18.0*Dxp[0]+9.0*Dxpp[0]-2.0*Dxppp[0])/(6.0*tau);
  Tot[0]+=0.5*AtomMass*Vx*Vx;
  Tot[Nx-1] = 0.5*FPU(a1x+Dx[Nx-1]-Dx[Nx-2]) + 0.5*FPU(a1x+Dx[0]-Dx[Nx-1]);
  Vx=(11.0*Dx[Nx-1]-18.0*Dxp[Nx-1]+9.0*Dxpp[Nx-1]-2.0*Dxppp[Nx-1])/(6.0*tau);
  Tot[Nx-1]+=0.5*AtomMass*Vx*Vx;
}// for void TotEnAtoms();

void Localization(){
  double Sum=0.0,Sum2=0.0;
  for (int i=0; i<Nx; i++){
      Sum+=Tot[i];
      Sum2+=sqr(Tot[i]);
  } // for i
  LocalizParam=Nx*Sum2/sqr(Sum);
}// for void Localization();

void ForceMD(){
  double Rx;
  for (int i=0; i<Nx-1; i++){
      Rx=a1x+Dx[i+1]-Dx[i];
      Fx[i]   += DFPU(Rx);
      Fx[i+1] -= DFPU(Rx);
  } // for i
  Rx=a1x+Dx[0]-Dx[Nx-1];
  Fx[Nx-1] += DFPU(Rx);
  Fx[0]    -= DFPU(Rx);
}// for void ForceMD();


double MaxF(){
  double MaxForce = 0.0, r;
  for (int nx=0; nx<Nx; nx++){
     r = fabs(Fx[nx]);
     if (r>MaxForce) MaxForce = r;
  }
  return MaxForce;
}// for double MaxF()


void TimeStepMD(){
  for (int i=0; i<Nx; i++){
      p1=(20.0*Dx[i]-6.0*Dxp[i]-4.0*Dxpp[i]+Dxppp[i]+Fx[i]*12.0*t2/AtomMass)/11.0;
      Dxppp[i]=Dxpp[i]; Dxpp[i]=Dxp[i]; Dxp[i]=Dx[i]; Dx[i]=p1;
  }
}// for void TimeStepMD()


void Stress(){
  double Rx;
  sxx=0.0;
  for (int i=0; i<Nx-1; i++){
      Rx=a1x+Dx[i+1]-Dx[i];
      sxx += DFPU(Rx);
  } // for i
  Rx=a1x+Dx[0]-Dx[Nx-1];
  sxx += DFPU(Rx);
  sxx=sxx/Volume;
}// for void Stress();



void FFT1D (int n, bool inverse, double x[], double y[]){
    //Calculate m=log_2(n)
    int m = 0;
    int p = 1;
    while(p < n){
        p *= 2;
        m++;
    }
    //Bit reversal
    GRe[n - 1] = gRe[n - 1];
    GIm[n - 1] = gIm[n - 1];
    int j = 0;
    for(int i = 0; i < n - 1; i++){
        GRe[i] = gRe[j];
        GIm[i] = gIm[j];
        int k = n / 2;
        while(k <= j){
            j -= k;
            k /= 2;
        }
        j += k;
    }
    //Calculate the FFT
    double ca = -1.0; 
    double sa = 0.0;
    int l1 = 1, l2 = 1;
    for(int l = 0; l < m; l++){
        l1 = l2;
        l2 *= 2;
        double u1 = 1.0;
        double u2 = 0.0;
        for(int j = 0; j < l1; j++){
            for(int i = j; i < n; i += l2){
                int i1 = i + l1;
                double t1 = u1 * GRe[i1] - u2 * GIm[i1];
                double t2 = u1 * GIm[i1] + u2 * GRe[i1];
                GRe[i1] = GRe[i] - t1;
                GIm[i1] = GIm[i] - t2;
                GRe[i] += t1;
                GIm[i] += t2;
            }
            double z =  u1 * ca - u2 * sa;
            u2 = u1 * sa + u2 * ca;
            u1 = z;
        }
        sa = sqrt((1.0 - ca) / 2.0);
        if(!inverse) sa =- sa;
        ca = sqrt((1.0 + ca) / 2.0);
    }
    //Divide through n if it isn't the IDFT
    if(!inverse)
    for(int i = 0; i < n; i++){
        GRe[i] /= n;
        GIm[i] /= n;
    }
    // Calculate the amplitudes
    for(int x = 0; x < n; x++){
        gAmp[x] = sqrt(GRe[x]*GRe[x] + GIm[x]*GIm[x]);
    }
}// void FFT1D (double x[], double y[])


void FourierTransform (double x[]){
  // Koeff. calculation
  for (int i=0; i<=Nx2; i++){
      fi1[i]=0.0; fi2[i]=0.0;
  }
  for (int k=0; k<=Nx2; k++){
      for (int i=0; i<Nx; i++){
          fi1[k]=fi1[k]+x[i]*cos(2.0*k*pi*i/(1.0*Nx));
      }// for i
  }// for k
  for (int k=1; k<Nx2; k++){
      for (int i=1; i<Nx; i++){
          fi2[k]=fi2[k]+x[i]*sin(2.0*k*pi*i/(1.0*Nx));
      }// for i
  }// for k
  // calculate amplitudes
  fi3[0]=sqrt(sqr(fi1[0]));
  for (int k=1; k<Nx2; k++) fi3[k]=sqrt(sqr(fi1[k])+sqr(fi2[k]));
  fi3[Nx2]=sqrt(sqr(fi1[Nx2]));
}// for void FourierTransform (double x[])


void FourierInverceTransform (double x1[], double x2[]){
  // Inverse transformation
  for (int j=0; j<Nx; j++) funinv[j]=0.0;
  for (int j=0; j<Nx; j++){
      for (int k=1; k<Nx2; k++)
      funinv[j]+=x1[k]*cos(2.0*pi*k*j/(1.0*Nx))+x2[k]*sin(2.0*pi*k*j/(1.0*Nx));
      funinv[j]+=0.5*x1[0];
      funinv[j]+=0.5*x1[Nx2]*cos(pi*j);
  }
  for (int j=0; j<Nx; j++) funinv[j]=funinv[j]*(2.0/(1.0*Nx));
}// for void FourierInverceTransform (double x1[], double x2[])


void Modulus(){
   int mm,m;
   double Rx,r,rr,rrr,Rx4;
   M11=0.0;
   for (int i=0; i<Nx; i++){
      m=i+1;
      mm=m;
      while (mm>=Nx) mm=mm-Nx;
      while (mm<0)  mm=mm+Nx;
      Rx=(m-i)*a1x-Dx[i]+Dx[mm];
      rr=Rx*Rx;
      r=sqrt(rr);
      rrr=rr*r;
      Rx4=Rx*Rx*Rx*Rx;
      M11+=DDFPU(r)*Rx4/rr+DFPU(r)*Rx4/rrr;
   }
   M11=M11/Nx;
}



