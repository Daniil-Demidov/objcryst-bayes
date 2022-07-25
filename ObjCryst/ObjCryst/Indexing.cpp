/*  ObjCryst++ Object-Oriented Crystallographic Library
    (c) 2006- Vincent Favre-Nicolin vincefn@users.sourceforge.net

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
/* 
   modified by Daniil Demidov
*/
/*
*  source file for Indexing classes & functions
*
*/
#include <algorithm>
#include <iomanip>
#include <fstream>
#include <cmath>

// mine
#include <sstream>
// for debug
#include <cassert>

#include <numeric>

#include <boost/integer/common_factor_rt.hpp>
#include <boost/math/special_functions/expm1.hpp>

#include <boost/algorithm/cxx11/all_of.hpp>

#include "ObjCryst/ObjCryst/Indexing.h"
#include "ObjCryst/Quirks/VFNDebug.h"
#include "ObjCryst/Quirks/VFNStreamFormat.h"
#include "ObjCryst/Quirks/Chronometer.h"

using namespace std;

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

#ifndef DEG2RAD
#define DEG2RAD (M_PI/180.)
#endif
#ifndef RAD2DEG
#define RAD2DEG (180./M_PI)
#endif

namespace ObjCryst
{

//debug globals
extern double global_p = 0.0;
extern double global_pmax = 0.0;
extern double global_pkmax = 0.0;
extern double global_kmax = 0.0; // not overall

extern double global_pk_integralmax = 0.0;


extern long long global_check_inside_count = 0;
extern long long global_is_inside_count = 0;

extern bool global_skip_all = false;

// just a variable to increment by anything so that compiler optimization is not possible
extern float global_s = 0;

// clean soon
// extern float global_min_ratio = 1000;
// extern float global_max_ratio = 0;

vector<unsigned long> repeated_calls_per_level(20);

vector<unsigned long> calls_per_level(20);
vector<double> sum_PK_per_level(20);
vector<double> sum_passed_PK_per_level(20);
// end debug globals

float EstimateCellVolume(const float dmin, const float dmax, const float nbrefl,
                         const CrystalSystem system,const CrystalCentering centering,const float kappa)
{
   const float q1=dmin*dmin*dmin-dmax*dmax*dmax;
   const float q2=dmin*dmin     -dmax*dmax;
   float D0,C0;
   if(system==TRICLINIC)
   {
      C0=2.095;
      return nbrefl/(C0*kappa*q1);
   }
   if(system==CUBIC)
   {
      if(centering==LATTICE_P) D0=0.862;
      if(centering==LATTICE_I) D0=0.475;
      if(centering==LATTICE_F) D0=0.354;
      return pow(nbrefl/(D0*kappa*q2),(float)1.5);
   }
   // "*.85" means using D0_min rather than D0
   if(system==MONOCLINIC) {C0=1.047;D0=0.786*.85;}
   if(system==ORTHOROMBIC){C0=0.524;D0=1.36 *.85;}
   if(system==HEXAGONAL)  {C0=0.150;D0=1.04 *.85;}
   if(system==RHOMBOHEDRAL){C0=0.230;D0=1.04 *.85;}
   if(system==TETRAGONAL) {C0=0.214;D0=1.25 *.85;}
   if((centering==LATTICE_I)||(centering==LATTICE_A)||(centering==LATTICE_B)||(centering==LATTICE_C)) {C0/=2;D0/=2;}
   if(centering==LATTICE_F){C0/=4;D0/=4;}
   const float alpha=D0*q2/(3*C0*q1), beta=nbrefl/(2*kappa*C0*q1);
   const float eta=beta-alpha*alpha*alpha,gamma=sqrt(beta*beta-2*beta*alpha*alpha*alpha);
   const float v=pow(pow(eta+gamma,(float)(1/3.))+pow(fabs(eta-gamma),(float)(1/3.))-alpha,(float)3);
   return v;
}

/** light-weight class storing the reciprocal space unitcell
*/

RecUnitCell::RecUnitCell(const float zero,const float p0,const float p1,const float p2,
                         const float p3,const float p4,const float p5,CrystalSystem lattice,
                         const CrystalCentering cent, const unsigned int nbspurious):
mlattice(lattice),mCentering(cent),mNbSpurious(nbspurious), P(0)
{
   this->par[0]=zero;
   this->par[1]=p0;
   this->par[2]=p1;
   this->par[3]=p2;
   this->par[4]=p3;
   this->par[5]=p4;
   this->par[6]=p5;
}

RecUnitCell::RecUnitCell(const RecUnitCell &old)
{
   *this=old;
}

void RecUnitCell::operator=(const RecUnitCell &rhs)
{
   for(unsigned int i=0;i<7;++i) par[i]=rhs.par[i];
   mlattice=rhs.mlattice;
   mCentering=rhs.mCentering;
   mNbSpurious=rhs.mNbSpurious;
   P = rhs.P;
}

float RecUnitCell::hkl2d(const float h,const float k,const float l,REAL *derivpar,const unsigned int derivhkl) const
{
   if((derivpar==NULL)&&(derivhkl==0))
   {
      switch(mlattice)
      {
         case TRICLINIC:
         {
            return par[0]+par[1]*h*h + par[2]*k*k + par[3]*l*l + par[4]*h*k + par[5]*k*l + par[6]*h*l;
            break;
         }
         case MONOCLINIC:
         {
            return par[0]+par[1]*par[1]*h*h + par[2]*par[2]*k*k + par[3]*par[3]*l*l + 2*par[1]*par[3]*par[4]*h*l;
            break;
         }
         case ORTHOROMBIC:
         {
            return par[0]+par[1]*par[1]*h*h + par[2]*par[2]*k*k + par[3]*par[3]*l*l;
            break;
         }
         case HEXAGONAL:
         {
            return par[0]+par[1]*par[1]*(h*h + k*k + h*k)+ par[2]*par[2]*l*l ;
            break;
         }
         case RHOMBOHEDRAL:
         {
            return par[0]+par[1]*par[1]*(h*h + k*k + l*l + 2*par[2]*(h*k + k*l + h*l));
            break;
         }
         case TETRAGONAL:
         {
            return par[0]+par[1]*par[1]*(h*h + k*k) + par[2]*par[2]*l*l;
            break;
         }
         case CUBIC:
         {
            return par[0]+par[1]*par[1]*(h*h+k*k+l*l);
            break;
         }
      }
   }
   if(derivhkl==1)
   {
      switch(mlattice)
      {
         case TRICLINIC:
         {
            return 2*par[1]*h + par[4]*k + par[6]*l;
            break;
         }
         case MONOCLINIC:
         {
            return 2*par[1]*par[1]*h + 2*par[1]*par[3]*par[4]*l;
            break;
         }
         case ORTHOROMBIC:
         {
            return 2*par[1]*par[1]*h;
            break;
         }
         case HEXAGONAL:
         {
            return par[1]*par[1]*(2*h + k);
            break;
         }
         case RHOMBOHEDRAL:
         {
            return par[1]*par[1]*(2*h + 2*par[2]*(k + l));
            break;
         }
         case TETRAGONAL:
         {
            return 2*par[1]*par[1]*h;
            break;
         }
         case CUBIC:
         {
            return 2*par[1]*par[1]*h;
            break;
         }
      }
   }
   if(derivhkl==2)
   {
      switch(mlattice)
      {
         case TRICLINIC:
         {
            return 2*par[2]*k + par[4]*h + par[5]*l;
            break;
         }
         case MONOCLINIC:
         {
            return 2*par[2]*par[2]*k;
            break;
         }
         case ORTHOROMBIC:
         {
            return 2*par[2]*par[2]*k;
            break;
         }
         case HEXAGONAL:
         {
            return par[1]*par[1]*(2*k + h);
            break;
         }
         case RHOMBOHEDRAL:
         {
            return par[1]*par[1]*(2*k + l*l + 2*par[2]*(h + l));
            break;
         }
         case TETRAGONAL:
         {
            return 2*par[1]*par[1]*k;
            break;
         }
         case CUBIC:
         {
            return 2*par[1]*par[1]*k;
            break;
         }
      }
   }
   if(derivhkl==3)
   {
      switch(mlattice)
      {
         case TRICLINIC:
         {
            return 2*par[3]*l + par[5]*k + par[6]*h;
            break;
         }
         case MONOCLINIC:
         {
            return 2*par[3]*par[3]*l + 2*par[1]*par[3]*par[4]*h;
            break;
         }
         case ORTHOROMBIC:
         {
            return 2*par[3]*par[3]*l;
            break;
         }
         case HEXAGONAL:
         {
            return 2*par[2]*par[2]*l;
            break;
         }
         case RHOMBOHEDRAL:
         {
            return par[1]*par[1]*(2*l + 2*par[2]*(k + h));
            break;
         }
         case TETRAGONAL:
         {
            return 2*par[2]*par[2]*l;
            break;
         }
         case CUBIC:
         {
            return 2*par[1]*par[1]*l;
            break;
         }
      }
   }

   if(derivpar==&par[0]) return 1.0;
   if(derivpar==(par+1))
   {
      switch(mlattice)
      {
         case TRICLINIC:
         {
            return h*h;
            break;
         }
         case MONOCLINIC:
         {
            return 2*par[1]*h*h + 2*par[3]*par[4]*h*l;
            break;
         }
         case ORTHOROMBIC:
         {
            return 2*par[1]*h*h;
            break;
         }
         case HEXAGONAL:
         {
            return 2*par[1]*(h*h + k*k + h*k);
            break;
         }
         case RHOMBOHEDRAL:
         {
            return 2*par[1]*(h*h + k*k + l*l + 2*par[2]*(h*k + k*l + h*l));
            break;
         }
         case TETRAGONAL:
         {
            return 2*par[1]*(h*h + k*k);
            break;
         }
         case CUBIC:
         {
            return 2*par[1]*(h*h+k*k+l*l);
            break;
         }
      }
   }
   if(derivpar==(par+2))
   {
      switch(mlattice)
      {
         case TRICLINIC:
         {
            return k*k;
            break;
         }
         case MONOCLINIC:
         {
            return 2*par[2]*k*k;
            break;
         }
         case ORTHOROMBIC:
         {
            return 2*par[2]*k*k;
            break;
         }
         case HEXAGONAL:
         {
            return 2*par[2]*l*l ;
            break;
         }
         case RHOMBOHEDRAL:
         {
            return par[1]*par[1]*(h*h + k*k + l*l + 2*(h*k + k*l + h*l));
            break;
         }
         case TETRAGONAL:
         {
            return 2*par[2]*l*l;
            break;
         }
         case CUBIC:
         {
            throw 0;
            break;
         }
      }
   }
   if(derivpar==(par+3))
   {
      switch(mlattice)
      {
         case TRICLINIC:
         {
            return l*l;
            break;
         }
         case MONOCLINIC:
         {
            return 2*par[3]*l*l + 2*par[1]*par[4]*h*l;
            break;
         }
         case ORTHOROMBIC:
         {
            return 2*par[3]*l*l;
            break;
         }
         case HEXAGONAL:
         {
            throw 0;
            break;
         }
         case RHOMBOHEDRAL:
         {
            throw 0;
            break;
         }
         case TETRAGONAL:
         {
            throw 0;
            break;
         }
         case CUBIC:
         {
            throw 0;
            break;
         }
      }
   }
   if(derivpar==(par+4))
   {
      switch(mlattice)
      {
         case TRICLINIC:
         {
            return h*k;
            break;
         }
         case MONOCLINIC:
         {
            return 2*par[1]*par[3]*h*l;
            break;
         }
         default:
         {
            throw 0;
            break;
         }
      }
   }
   if(derivpar==(par+5))
   {
      switch(mlattice)
      {
         case TRICLINIC:
         {
            return k*l;
            break;
         }
         default:
         {
            throw 0;
            break;
         }
      }
   }
   if(derivpar==(par+6))
   {
      switch(mlattice)
      {
         case TRICLINIC:
         {
            return h*l;
            break;
         }
         default:
         {
            throw 0;
            break;
         }
      }
   }
   throw 0;
   return 0.0;
}

void RecUnitCell::hkl2d_delta(const float h,const float k,const float l,
                              const RecUnitCell &delta, float & dmin, float &dmax) const
{
   const float p0m=par[0]-delta.par[0] , p0p=par[0]+delta.par[0],
               p1m=par[1]-delta.par[1] , p1p=par[1]+delta.par[1],
               p2m=par[2]-delta.par[2] , p2p=par[2]+delta.par[2],
               p3m=par[3]-delta.par[3] , p3p=par[3]+delta.par[3],
               p4m=par[4]-delta.par[4] , p4p=par[4]+delta.par[4],
               p5m=par[5]-delta.par[5] , p5p=par[5]+delta.par[5],
               p6m=par[6]-delta.par[6] , p6p=par[6]+delta.par[6];
   switch(mlattice)
   {
      case TRICLINIC:
      {//par[0]+par[1]*h*h + par[2]*k*k + par[3]*l*l + par[4]*h*k + par[5]*k*l + par[6]*h*l;
         float p4mm,p5mm,p6mm,p4pp,p5pp,p6pp;
         if((h*k)>0){p4mm=p4m;p4pp=p4p;}else{p4mm=p4p;p4pp=p4m;}
         if((k*l)>0){p5mm=p5m;p5pp=p5p;}else{p5mm=p5p;p5pp=p5m;}
         if((h*l)>0){p6mm=p6m;p6pp=p6p;}else{p6mm=p6p;p6pp=p6m;}
         dmin=p0m+p1m*h*h+p2m*k*k+p3m*l*l+p4mm*h*k+p5mm*k*l+p6mm*h*l;
         dmax=p0p+p1p*h*h+p2p*k*k+p3p*l*l+p4pp*h*k+p5pp*k*l+p6pp*h*l;
         /*
         if(dmin<0)
         {
            cout<<"hkl2d_delta: dmin<0 ! "<<int(h)<<","<<int(k)<<","<<int(l)<<endl;
            for(unsigned int i=0;i<7;++i) cout<<par[i]<<" +/-"<<delta.par[i]<<endl;
            exit(0);
         }*/
         return;
      }
      case MONOCLINIC: //OK
      {
         if((h*l)>0)
         {
            dmin = p0m + p1m*p1m*h*h + p2m*p2m*k*k + p3m*p3m*l*l + 2*p1m*p3m*p4m*h*l;
            dmax = p0p + p1p*p1p*h*h + p2p*p2p*k*k + p3p*p3p*l*l + 2*p1p*p3p*p4p*h*l;
            return;
         }
         else
         {
            const bool b1=(h*(par[1]*h+par[3]*par[4]*l))>0;// d(d*^2)/dp1
            const bool b3=(l*(par[3]*l+par[1]*par[4]*h))>0;// d(d*^2)/dp2
            if(b1 && b3)
            {
               dmin = p0m + p1m*p1m*h*h + p2m*p2m*k*k + p3m*p3m*l*l + 2*p1m*p3m*p4p*h*l;
               dmax = p0p + p1p*p1p*h*h + p2p*p2p*k*k + p3p*p3p*l*l + 2*p1p*p3p*p4m*h*l;
               return;
            }
            else if(b1 && (!b3))
               {
                  dmin = p0m + p1m*p1m*h*h + p2m*p2m*k*k + p3p*p3p*l*l + 2*p1m*p3p*p4p*h*l;
                  dmax = p0p + p1p*p1p*h*h + p2p*p2p*k*k + p3m*p3m*l*l + 2*p1p*p3m*p4m*h*l;
                  return;
               }
               else if((!b1) && b3)
                  {
                     dmin = p0m + p1p*p1p*h*h + p2m*p2m*k*k + p3m*p3m*l*l + 2*p1p*p3m*p4p*h*l;
                     dmax = p0p + p1m*p1m*h*h + p2p*p2p*k*k + p3p*p3p*l*l + 2*p1m*p3p*p4m*h*l;
                     return;
                  }
                  else
                  {
                     dmin = p0m + p1p*p1p*h*h + p2m*p2m*k*k + p3p*p3p*l*l + 2*p1p*p3p*p4p*h*l;
                     dmax = p0p + p1m*p1m*h*h + p2p*p2p*k*k + p3m*p3m*l*l + 2*p1m*p3m*p4m*h*l;
                     return;
                  }
            }
      }
      case ORTHOROMBIC: //OK
      {
         dmin= p0m + p1m*p1m*h*h + p2m*p2m*k*k + p3m*p3m*l*l;
         dmax= p0p + p1p*p1p*h*h + p2p*p2p*k*k + p3p*p3p*l*l;
         return;
      }
      case HEXAGONAL: //OK
      {
         dmin=p0m+p1m*p1m*(h*h + k*k + h*k)+ p2m*p2m*l*l ;
         dmax=p0p+p1p*p1p*(h*h + k*k + h*k)+ p2p*p2p*l*l ;
         return;
      }
      case RHOMBOHEDRAL:
      {
         if((h*k + k*l + h*l)>=0)
         {
            dmin= p0m+p1m*p1m*(h*h + k*k + l*l + 2*p2m*(h*k + k*l + h*l));
            dmax= p0p+p1p*p1p*(h*h + k*k + l*l + 2*p2p*(h*k + k*l + h*l));
         }
         else
         {
            dmin= p0m+p1m*p1m*(h*h + k*k + l*l + 2*p2p*(h*k + k*l + h*l));
            dmax= p0p+p1p*p1p*(h*h + k*k + l*l + 2*p2m*(h*k + k*l + h*l));
         }
         return;
      }
      case TETRAGONAL: //OK
      {
         dmin= p0m + p1m*p1m*(h*h + k*k) + p2m*p2m*l*l;
         dmax= p0p + p1p*p1p*(h*h + k*k) + p2p*p2p*l*l;
         return;
      }
      case CUBIC: //OK
      {
         dmin=p0m + p1m*p1m*(h*h+k*k+l*l);
         dmax=p0p + p1p*p1p*(h*h+k*k+l*l);
         return;
      }
   }
}

vector<float> RecUnitCell::DirectUnitCell(const bool equiv)const
{
   // reciprocal unit cell parameters
   float aa,bb,cc,calphaa,cbetaa,cgammaa,salphaa,sbetaa,sgammaa;
   switch(mlattice)
   {
      case TRICLINIC:
      {
         aa=sqrt(par[1]);
         bb=sqrt(par[2]);
         cc=sqrt(par[3]);
         calphaa=par[5]/(2*bb*cc);
         cbetaa =par[6]/(2*aa*cc);
         cgammaa=par[4]/(2*aa*bb);
         salphaa=sqrt(abs(1-calphaa*calphaa));
         sbetaa =sqrt(abs(1-cbetaa *cbetaa));
         sgammaa=sqrt(abs(1-cgammaa*cgammaa));
         break;
      }
      case MONOCLINIC:
      {
         aa=par[1];
         bb=par[2];
         cc=par[3];
         calphaa=0;
         cbetaa=par[4];
         cgammaa=0;
         salphaa=1;
         sbetaa =sqrt(abs(1-cbetaa *cbetaa));
         sgammaa=1;
         break;
      }
      case ORTHOROMBIC:
      {
         aa=par[1];
         bb=par[2];
         cc=par[3];
         calphaa=0;
         cbetaa =0;
         cgammaa=0;
         salphaa=1;
         sbetaa =1;
         sgammaa=1;
         break;
      }
      case HEXAGONAL:
      {
         aa=par[1];
         bb=par[1];
         cc=par[2];
         calphaa=0;
         cbetaa =0;
         cgammaa=0.5;
         salphaa=1;
         sbetaa =1;
         sgammaa=0.8660254037844386;
         break;
      }
      case RHOMBOHEDRAL:
      {
         aa=par[1];
         bb=par[1];
         cc=par[1];
         calphaa=par[4];
         cbetaa =par[4];
         cgammaa=par[4];
         salphaa=sqrt(abs(1-calphaa *calphaa));
         sbetaa =salphaa;
         sgammaa=salphaa;
         break;
      }
      case TETRAGONAL:
      {
         aa=par[1];
         bb=par[1];
         cc=par[2];
         calphaa=0;
         cbetaa =0;
         cgammaa=0;
         salphaa=1;
         sbetaa =1;
         sgammaa=1;
         break;
      }
      case CUBIC:
      {
         aa=par[1];
         bb=par[1];
         cc=par[1];
         calphaa=0;
         cbetaa =0;
         cgammaa=0;
         salphaa=1;
         sbetaa =1;
         sgammaa=1;
         break;
      }
      // This should never happen.  Avoid using unitialized cell parameters.
      default:
         throw 0;
   }
   // Volume of reciprocal unit cell
   const float vv=sqrt(abs(1-calphaa*calphaa-cbetaa*cbetaa-cgammaa*cgammaa+2*calphaa*cbetaa*cgammaa));

   const float a=salphaa/(aa*vv);
   const float b=sbetaa /(bb*vv);
   const float c=sgammaa/(cc*vv);

   const float calpha=(cbetaa *cgammaa-calphaa)/(sbetaa *sgammaa);
   const float cbeta =(calphaa*cgammaa-cbetaa )/(salphaa*sgammaa);
   const float cgamma=(calphaa*cbetaa -cgammaa)/(salphaa*sbetaa );

   const float alpha=acos( calpha );
   const float beta =acos( cbeta );
   const float gamma=acos( cgamma );

   const float v=a*b*c*sqrt(1-calpha*calpha-cbeta*cbeta-cgamma*cgamma+2*calpha*cbeta*cgamma);

   vector<float> par(7);
   par[0]=a;
   par[1]=b;
   par[2]=c;
   par[3]=alpha;
   par[4]=beta;
   par[5]=gamma;
   par[6]=v;
   return par;
}

///////////////////////////////////////////////// FUNCTIONTABLE /////////////////////
FunctionTable::FunctionTable()
: scale(1.0), xmax(0.0)
{}

FunctionTable::FunctionTable(vector<float> & data, float scale, float xmax, bool rob)
: scale(scale), xmax(xmax)
{
    if (rob)
    {
        this->data.swap(data);
    }
    else
    {
        this->data = data;
    }
}

float FunctionTable::operator[](float x) const
{
   // temporary
   if (x < 0) return (float)0.0;

   if (x > xmax) return (float)0.0;
   return data[scale*x];
}

///////////////////////////////////////////////// PEAKLIST:INTEGRAL /////////////////////
PeakList::integral::integral(const float d2, const float Rho, const float F, const float Rho_ext)
: d2(d2), Rho(Rho), F(F), Rho_ext(Rho_ext)
{}

///////////////////////////////////////////////// PEAKLIST:THKL /////////////////////
PeakList::thkl::thkl(int h, int k, int l, int count, float d2min, float d2max, double K, double dRho, double P, bool is_disabled)
: h(h), k(k), l(l), count(count), d2min(d2min), d2max(d2max), K(K), dRho(dRho), P(P), is_disabled(is_disabled)
{}

bool operator<(const PeakList::thkl & thkl1, const PeakList::thkl & thkl2)
{
   if (thkl1.h == thkl2.h)
   {
      if (thkl1.k == thkl2.k)
      {
         return thkl1.l < thkl2.l;
      }
      return thkl1.k < thkl2.k;
   }
   return thkl1.h < thkl2.h;
}

bool operator>(const PeakList::thkl & thkl1, const PeakList::thkl & thkl2)
{
   return (thkl2 < thkl1);
}

bool operator<=(const PeakList::thkl & thkl1, const PeakList::thkl & thkl2)
{
   return !(thkl1 > thkl2);
}

bool operator>=(const PeakList::thkl & thkl1, const PeakList::thkl & thkl2)
{
   return !(thkl1 < thkl2);
}

///////////////////////////////////////////////// PEAKLIST:HKL0 /////////////////////
PeakList::hkl0::hkl0(const int h0,const int k0, const int l0):
h(h0),k(k0),l(l0)
{}

bool operator<(const PeakList::hkl0 & hkl1, const PeakList::hkl0 & hkl2)
{
   if (hkl1.h == hkl2.h)
   {
      if (hkl1.k == hkl2.k)
      {
         return hkl1.l < hkl2.l;
      }
      return hkl1.k < hkl2.k;
   }
   return hkl1.h < hkl2.h;
}

bool operator>(const PeakList::hkl0 & hkl1, const PeakList::hkl0 & hkl2)
{
   return (hkl2 < hkl1);
}

bool operator<=(const PeakList::hkl0 & hkl1, const PeakList::hkl0 & hkl2)
{
   return !(hkl2 < hkl1);
}

bool operator>=(const PeakList::hkl0 & hkl1, const PeakList::hkl0 & hkl2)
{
   return !(hkl1 < hkl2);
}

///////////////////////////////////////////////// PEAKLIST:HKL /////////////////////
// changed from d2obsmin((d-ds/2)*(d-ds/2)),d2obsmax((d+ds/2)*(d+ds/2))
PeakList::hkl::hkl(const float d,const float i,const float ds,const float is,
                   const int h0,const int k0, const int l0,const float dc0):
dobs(d),dobssigma(ds),d2obs(d*d),d2obsmin((d-ds)*(d-ds)),d2obsmax((d+ds)*(d+ds)),iobs(i),iobssigma(is),
h(h0),k(k0),l(l0),isIndexed(false),isSpurious(false),stats(0),
d2calc(dc0),d2diff(0)
{}

bool compareHKL_d(const PeakList::hkl &d1, const PeakList::hkl &d2)
{
   return d1.dobs < d2.dobs;
}

///////////////////////////////////////////////// PEAKLIST:AX /////////////////////
PeakList::Ax::Ax(vector<thkl *> const & vthkl, double P, double N_with1, double N_without1)
: vthkl(vthkl), P(P), N_with1(N_with1), N_without1(N_without1)
{}

///////////////////////////////////////////////// PEAKLIST /////////////////////

PeakList::PeakList()
{}

PeakList::PeakList(const PeakList &old)
{
   *this=old;
}

void PeakList::operator=(const PeakList &rhs)
{
   VFN_DEBUG_ENTRY("PeakList::operator=(PeakList &old)",10);
   mvHKL=rhs.mvHKL;
   integrals = rhs.integrals;
   vthkl = rhs.vthkl;
   vax = rhs.vax;
   mvcount = rhs.mvcount;
   mvvcriterion = rhs.mvvcriterion;
   mvvPK = rhs.mvvPK;
   mvvis_good = rhs.mvvis_good;
   mvPK_integral = rhs.mvPK_integral;
   mis_exta_run = rhs.mis_exta_run;
   mPK_solution_critical = rhs.mPK_solution_critical;
   active_depth = rhs.active_depth;
   critical_depth = rhs.critical_depth;
   critical_passability = rhs.critical_passability;

   VFN_DEBUG_EXIT("PeakList::operator=(PeakList &old)",10);
}

PeakList::~PeakList()
{}

void PeakList::ImportDhklDSigmaIntensity(istream &is,float defaultsigma)
{
   float d,sigma,iobs;
   while(true)
   {// :TODO: use readline to make sure when the end is reached
      is >>d;
      cout<<__FILE__<<":"<<__LINE__<<"  "<<mvHKL.size()<<":d="<<d;
      if(is.good()==false) break;
      is>>sigma;
      if(is.good()==false) break;
      is>>iobs;
      if(sigma<=0) sigma=d*defaultsigma;
      if(iobs<=0) iobs=1.0;
      mvHKL.push_back(hkl(1/d,iobs,1/(d-sigma/2)-1/(d+sigma/2)));
      cout<<"+/-"<<sigma<<", I="<<iobs<<" 1/d="<<1/d<<endl;
      if(is.good()==false) break;
   }
   sort(mvHKL.begin(),mvHKL.end(),compareHKL_d);
   cout<<"Imported "<<mvHKL.size()<<" observed reflection positions."<<endl;
}

void PeakList::ImportDhklIntensity(istream &is)
{
   float d,iobs;
   while(true)
   {// :TODO: use readline to make sure when the end is reached
      is >>d;
      if(is.eof()) break;
      is>>iobs;
      mvHKL.push_back(hkl(1/d,iobs));
      cout<<__FILE__<<":"<<__LINE__<<"  "<<mvHKL.size()<<":d="<<d<<", I="<<iobs<<" 1/d="<<1/d<<endl;
      if(is.eof()) break;
   }
   sort(mvHKL.begin(),mvHKL.end(),compareHKL_d);
   cout<<"Imported "<<mvHKL.size()<<" observed reflection positions."<<endl;
}

void PeakList::ImportDhkl(istream &is)
{
   std::vector<std::pair<float,float> > v;
   float d;
   while(true)
   {// :TODO: use readline to make sure when the end is reached
      is >>d;
      if(is.eof()) break;
      mvHKL.push_back(hkl(1/d));
      cout<<__FILE__<<":"<<__LINE__<<"  "<<mvHKL.size()<<":d="<<d<<" 1/d="<<1/d<<endl;
      if(is.eof()) break;
   }
   sort(mvHKL.begin(),mvHKL.end(),compareHKL_d);
   cout<<"Imported "<<v.size()<<" observed reflection positions."<<endl;
}

void PeakList::ImportDhklDSigma(istream &is)
{
   char sharp_char;
   is >> sharp_char;
   assert (sharp_char == '#');
   is.ignore(100, '\n');

   float d, sigma, iobs;
   iobs = 1;
   while(true)
   {// :TODO: use readline to make sure when the end is reached
      is >> d;
      cout<<__FILE__<<":"<<__LINE__<<"  "<<mvHKL.size()<<":d="<<d;
      if(is.good()==false) break;
      is >> sigma;
      if(is.good()==false) break;
      mvHKL.push_back(hkl(1/d,iobs,1/(d-sigma/2)-1/(d+sigma/2)));
      cout<<"+/-"<<sigma<<", I="<<iobs<<" 1/d="<<1/d<<endl;
      if(is.good()==false) break;
   }
   sort(mvHKL.begin(),mvHKL.end(),compareHKL_d);
   cout<<"Imported "<<mvHKL.size()<<" observed reflection positions."<<endl;
}

void PeakList::ImportIntegrals(istream &is)
{
   char sharp_char;
   is >> sharp_char;
   assert (sharp_char == '#');
   is.ignore(100, '\n');
   float d2, Rho, F, Rho_ext;
   while(true)
   {// :TODO: use readline to make sure when the end is reached
      is >> d2;
      if(is.good()==false) break;
      is >> Rho;
      if(is.good()==false) break;
      is >> F;
      if(is.good()==false) break;
      is >> Rho_ext;
      integrals.push_back(integral(d2, Rho, F, Rho_ext));
      if(is.good()==false) break;
   }
   coef = 1/integrals[1].d2;
   d2max = integrals.back().d2*0.999;
   cout<<"Imported "<<integrals.size()<<" integral values."<<endl;
}

template<class T,class U> bool comparePairFirst(std::pair<T,U> &p1, std::pair<T,U> &p2)
{
   return p1.first < p2.first;
}

void PeakList::Import2ThetaIntensity(istream &is, const float wavelength)
{
   std::list<std::pair<float,float> > v;
   float d,iobs;
   while(true)
   {// :TODO: use readline to make sure when the end is reached
      is >>d;
      if(is.eof()) break;
      is>>iobs;
      d=2*sin(d/2*DEG2RAD)/wavelength;
      mvHKL.push_back(hkl(1/d,iobs));
      cout<<__FILE__<<":"<<__LINE__<<"  "<<mvHKL.size()<<":d="<<1/d<<", I="<<iobs<<" 1/d="<<d<<endl;
      if((is.eof())||v.size()>=20) break;
   }
   sort(mvHKL.begin(),mvHKL.end(),compareHKL_d);
   cout<<"Imported "<<v.size()<<" observed reflection positions."<<endl;
}
float PeakList::Simulate(float zero, float a, float b, float c,
              float alpha, float beta, float gamma,
              bool deg, unsigned int nb, unsigned int nbspurious,
              float sigma,float percentMissing, const bool verbose)
{
   if(deg){alpha*=DEG2RAD;beta*=DEG2RAD;gamma*=DEG2RAD;}
   const float v=sqrt(1-cos(alpha)*cos(alpha)-cos(beta)*cos(beta)-cos(gamma)*cos(gamma)
            +2*cos(alpha)*cos(beta)*cos(gamma));

   const float aa=sin(alpha)/a/v;
   const float bb=sin(beta )/b/v;
   const float cc=sin(gamma)/c/v;

   const float alphaa=acos( (cos(beta )*cos(gamma)-cos(alpha))/sin(beta )/sin(gamma) );
   const float betaa =acos( (cos(alpha)*cos(gamma)-cos(beta ))/sin(alpha)/sin(gamma) );
   const float gammaa=acos( (cos(alpha)*cos(beta )-cos(gamma))/sin(alpha)/sin(beta ) );

   RecUnitCell ruc(zero,aa*aa,bb*bb,cc*cc,2*aa*bb*cos(gammaa),2*bb*cc*cos(alphaa),2*aa*cc*cos(betaa),TRICLINIC,LATTICE_P);
   std::list<float> vd2;

   for(int h=0;h<=20;h++)
      for(int k=-20;k<=20;k++)
      {
         if((h==0) && (k<0)) k=0;
         for(int l=-20;l<=20;l++)
         {
           if((h==0) && (k==0) && (l<=0)) l=1;
           vd2.push_back(sqrt(ruc.hkl2d(h,k,l)));
         }
      }
   //
   std::list<float>::iterator pos=vd2.begin();
   if(percentMissing>0.90) percentMissing=0.90;
   for(;pos!=vd2.end();++pos)
   {
      if((rand()/float(RAND_MAX))<percentMissing) *pos=1e10;
   }
   vd2.sort();
   pos=vd2.begin();
   const float dmin=*pos/2;
   for(unsigned int i=0;i<nb;++i)pos++;
   const float dmax=*pos;

   for(unsigned int i=0;i<nbspurious;++i)
   {
      const unsigned int idx=1+i*nb/nbspurious+(rand()%nbspurious);
      pos=vd2.begin();
      for(unsigned int j=0;j<idx;++j) pos++;
      *pos=dmin+rand()/float(RAND_MAX)*(dmax-dmin);
   }

   pos=vd2.begin();
   for(unsigned int i=0;i<nb;++i)
   {
      float d=*pos++;
      const float ds=d*sigma;
      float d1=d+ds*(rand()/float(RAND_MAX)*2-1);
      //cout<<d<<"  "<<ds<<"  "<<d1<<"   "<<sigma<<endl;
      mvHKL.push_back(hkl(d1,1.0,ds));
   }

   if(verbose)
   {
      char buf[200];
      sprintf(buf,"a=%6.3f b=%6.3f c=%6.3f alpha=%6.2f beta=%6.2f gamma=%6.2f V=%8.2f",a,b,c,alpha*RAD2DEG,beta*RAD2DEG,gamma*RAD2DEG,v*a*b*c);
      cout<<"PeakList::Simulate():"<<buf<<endl;
   }
   sort(mvHKL.begin(),mvHKL.end(),compareHKL_d);
   return v*a*b*c;
}

void PeakList::ExportDhklDSigmaIntensity(std::ostream &os)const
{
   for(vector<PeakList::hkl>::const_iterator pos=mvHKL.begin();pos!=mvHKL.end();++pos)
   {
      const float sigma=1/(pos->dobs-pos->dobssigma/2)-1/(pos->dobs+pos->dobssigma/2);
      os<< std::fixed << setw(6) << setprecision(3) << 1/pos->dobs <<" "<< sigma <<" "<< std::scientific << pos->iobs <<endl;
   }
}

void PeakList::AddPeak(const float d, const float iobs,const float dobssigma,const float iobssigma,
                       const int h,const int k, const int l,const float d2calc)
{
   if(dobssigma<=0)
   {// Manually added peak ? Use other reflection's sigmas to evaluate sigma for this reflection
      float s=0;
      for(vector<hkl>::const_iterator pos=mvHKL.begin();pos!=mvHKL.end();++pos)
         s+= pos->dobssigma;
      s/=mvHKL.size();
      if(s>0) mvHKL.push_back(hkl(d,iobs,s,iobssigma,h,k,l,d2calc));
      else mvHKL.push_back(hkl(d,iobs,1e-4,iobssigma,h,k,l,d2calc));
   }
   else mvHKL.push_back(hkl(d,iobs,dobssigma,iobssigma,h,k,l,d2calc));
   sort(mvHKL.begin(),mvHKL.end(),compareHKL_d);
   //this->Print(cout);
}

void PeakList::RemovePeak(unsigned int idx)
{
   for(unsigned int i=idx;i<(mvHKL.size()-1);++i) mvHKL[i]=mvHKL[i+1];
   mvHKL.pop_back();
}

void PeakList::Print(std::ostream &os) const
{
   unsigned int i=0;
   char buf[200];
   os<<"PeakList, with "<<mvHKL.size()<<" peaks"<<endl;
   for(vector<PeakList::hkl>::const_iterator pos=mvHKL.begin();pos!=mvHKL.end();++pos)
   {
      const float sigma=1/(pos->dobs-pos->dobssigma/2)-1/(pos->dobs+pos->dobssigma/2);
      if(pos->isIndexed)
         sprintf(buf,"#%3d d=%6.3f+/-%7.4f dcalc=%6.3f, diff=%7.4f, iobs=%6.3f HKL=%2d %2d %2d Spurious=%1d stats=%6lu",
                 i++,1/pos->dobs,sigma,
                 1/sqrt(abs(pos->d2calc)),1/sqrt(abs(pos->d2calc))-1/pos->dobs,
                 pos->iobs,pos->h,pos->k,pos->l,pos->isSpurious,pos->stats);
      else
         sprintf(buf,"#%3d d=%6.3f+/-%6.3f              iobs=%6.3f  UNINDEXED   Spurious=%1d stats=%6lu",
                 i++,1/pos->dobs,1/(pos->dobs-pos->dobssigma/2)-1/(pos->dobs+pos->dobssigma/2),
                 pos->iobs,pos->isSpurious,pos->stats);
      os<<buf<<endl;
   }
}

vector<PeakList::hkl> & PeakList::GetPeakList(){return mvHKL;}

const vector<PeakList::hkl> & PeakList::GetPeakList()const {return mvHKL;}

void PeakList::clear_locals() const
{
   vector<PeakList::thkl>::const_iterator pthkl, vthkl_end;
   vthkl_end = vthkl.end();
   for (pthkl = vthkl.begin(); pthkl != vthkl_end; ++pthkl)
   {
      pthkl->count = 0;
      pthkl->d2min = pthkl->d2max = 1;
      pthkl->K = 1;
      pthkl->dRho = -1;
      pthkl->P = -1;
      pthkl->is_disabled = false;
   }
   
   for (vector<PeakList::Ax>::const_iterator pax = vax.begin(); pax != vax.end(); ++pax)
   {
      pax->N_with1 = 0;
      pax->N_without1 = 0;
      pax->P = 1;
   }

   for (vector<PeakList::hkl>::const_iterator pos = mvHKL.begin(); pos != mvHKL.end(); ++pos)
   {
      pos->h = pos->k = pos->l = 0;
      pos->isIndexed = pos->isSpurious = false;
      pos->stats = 0;
      pos->d2calc = pos->d2diff = 0;
   }
}

/////////////////////////////////////////////////////// SCORE ///////////////////////////////////////

float Score(const PeakList &dhkl, const RecUnitCell &rpar, const unsigned int nbSpurious,
            const bool verbose,const bool storehkl,const bool storePredictedHKL)
{
   const bool autozero=false;
   vector<PeakList::hkl>::const_iterator pos,first,last;
   for(pos=dhkl.GetPeakList().begin();pos!=dhkl.GetPeakList().end();++pos)
   {
      if(storehkl) pos->isIndexed=false;
      pos->d2calc=0;
      pos->d2diff=1000;
   }
   const unsigned long nb=dhkl.GetPeakList().size();
   if(storePredictedHKL) dhkl.mvPredictedHKL.clear();

   unsigned long nbCalc=0;
   int h,k,l;
   float predict_coeff=1;
   if(storePredictedHKL)predict_coeff=2;
   const float dmax=dhkl.mvHKL[nb-1].d2obs*predict_coeff*1.05;
   int sk0,sl0;// do we need >0 *and* <0 indices for k,l ?
   switch(rpar.mlattice)
   {
      case TRICLINIC:
         sk0=-1;sl0=-1;
	 break;
      case MONOCLINIC:
         sk0=1;sl0=-1;
         break;
      case ORTHOROMBIC:
         sk0=1;sl0=1;
         break;
      case HEXAGONAL:
         sk0=-1;sl0=1;
         break;
      case RHOMBOHEDRAL:
         sk0=-1;sl0=-1;
         break;
      case TETRAGONAL:
         sk0=1;sl0=1;
         break;
      case CUBIC:
         sk0=1;sl0=1;
         break;
      // This should never happen.  Avoid using unitialized values.
      default:
         throw 0;
   }
   int stepk,stepl;// steps in k,l to use for centered lattices
   switch(rpar.mCentering)
   {
      case LATTICE_P:stepk=1;stepl=1;break;
      case LATTICE_I:stepk=1;stepl=2;break;
      case LATTICE_A:stepk=1;stepl=2;break;
      case LATTICE_B:stepk=1;stepl=2;break;
      case LATTICE_C:stepk=2;stepl=1;break;
      case LATTICE_F:stepk=2;stepl=2;break;
      // This should never happen.  Avoid using unitialized values.
      default: throw 0;
   }
   first=dhkl.GetPeakList().begin();last=dhkl.GetPeakList().end();
   unsigned long nbCalcH,nbCalcK;// Number of calculated lines below dmax for each h,k
   for(h=0;;++h)
   {
      nbCalcH=0;
      for(int sk=sk0;sk<=1;sk+=2)
      {
         if(h==0) sk=1;// no need to explore 0kl with both sk -1 and 1
         if(stepk==2) k=(h%2);// For LATTICE_C,LATTICE_F: h odd => k odd
         else k=0;
         for(;;k+=stepk)
         {
            nbCalcK=0;
            for(int sl=sl0;sl<=1;sl+=2)
            {
               if((h+k)==0)
               {
                  sl=1;// No need to list 0 0 l with l<0
                  l=1;
               }
               else
               {
                  if(h==0)
                  {
                     if(rpar.mlattice==MONOCLINIC) sl=1;// 0 k l and 0 k -l are equivalent
                     if((sk<0)||(sl<0)) l=1;// Do not list 0 k 0 with k<0
                     else l=0;// h==k==0 already covered
                  }
                  else
                  {
                     if(sl<0) l=1;// Do not list h k 0 twice
                     else l=0;
                  }
               }
               if(stepl==2)
               {
                  if(rpar.mCentering==LATTICE_I) l+=(h+k+l)%2;
                  if(rpar.mCentering==LATTICE_A) l+=(k+l)%2;// Start at hk1 if k odd
                  if(  (rpar.mCentering==LATTICE_B)
                     ||(rpar.mCentering==LATTICE_F)) l+=(h+l)%2;// Start at hk1 if h odd
               }
               for(;;l+=stepl)
               {
                  const float d2=rpar.hkl2d(h,sk*k,sl*l);
                  if(d2>dmax)
                  {
                     //cout<<__FILE__<<":"<<__LINE__<<" hkl: "<<h<<" "<<sk*k<<" "<<sl*l<<":"<<sqrt(d2)<<" deriv="<<sl*rpar.hkl2d(h,sk*k,sl*l,NULL,3)<<"/"<<sqrt(dmax)<<endl;
                     // Only break if d is increasing with l
                     if((sl*rpar.hkl2d(h,sk*k,sl*l,NULL,3))>=0) break;
                     else continue;
                  }
                  nbCalc++;nbCalcK++;nbCalcH++;
                  if(storePredictedHKL)
                  {
                     dhkl.mvPredictedHKL.push_back(PeakList::hkl(0,0,0,0,h,sk*k,sl*l,d2));
                     //continue;
                  }
                  for(pos=first;pos!=last;++pos)
                  {
                     const float tmp=d2-pos->d2obs;
                     if(tmp<.1)
                     {
                        if(tmp<-.1) break;
                        if(fabs(tmp)<fabs(pos->d2diff))
                        {
                           pos->d2diff=tmp;
                           if(storehkl)
                           {
                              pos->h=h;
                              pos->k=sk*k;
                              pos->l=sl*l;
                              pos->isIndexed=true;
                              pos->d2calc=d2;
                           }
                        }
                        /*
                        if((verbose)&&(fabs(tmp)<.01))
                           cout<<__FILE__<<":"<<__LINE__<<"      hkl: "<<h<<" "<<k<<" "<<l
                              <<"#"<<i<<": calc="<<sqrt(d2)<<", obs="<<sqrt(*pd2)<<", min_epsilon="<<*pdq2<<", tmp="<<tmp<<endl;
                        */
                     }
                  }
               }
            }
            if(nbCalcK==0) //d(hk0)>dmax
            {
               //cout<<__FILE__<<":"<<__LINE__<<" hkl: "<<h<<" "<<sk*k<<" "<<0<<" deriv="<<sk*rpar.hkl2d(h,sk*k,0,NULL,2)<<endl;
               if((sk*rpar.hkl2d(h,sk*k,0,NULL,2))>=0) break;
            }
         }
      }
      if(nbCalcH==0) break;//h00 beyond limit
   }
   float epsilon=0.0,zero=0.0;
   if(autozero)
   {
      for(pos=dhkl.GetPeakList().begin();pos!=dhkl.GetPeakList().end();++pos) zero+=pos->d2diff;
      zero/=nb;
   }
   for(pos=dhkl.GetPeakList().begin();pos!=dhkl.GetPeakList().end();++pos)
   {
      epsilon +=fabs(pos->d2diff-zero);
   }
   if(nbSpurious>0)
   {// find worst fitting lines and remove them from epsilon calculation
      list<pair<float,unsigned int> > vdiff_idx;
      unsigned int i=0;
      for(pos=dhkl.GetPeakList().begin();pos!=dhkl.GetPeakList().end();++pos)
         vdiff_idx.push_back(make_pair(fabs(pos->d2diff),i++));
      vdiff_idx.sort(comparePairFirst<float,unsigned int>);
      i=0;
      for(list<pair<float,unsigned int> >::reverse_iterator rpos=vdiff_idx.rbegin();rpos!=vdiff_idx.rend();++rpos)
      {// :TODO: correct zero after removing spurious lines
         epsilon -= fabs(rpos->first-zero);
         if(storehkl) dhkl.GetPeakList()[rpos->second].isIndexed=false;
         dhkl.GetPeakList()[rpos->second].stats++;
         if(++i==nbSpurious) break;
      }
   }
   if(verbose)
   {
      float epstmp=0;
      //unsigned long i=0;
      for(pos=dhkl.GetPeakList().begin();pos!=dhkl.GetPeakList().end();++pos)
      {
         epstmp+=fabs(pos->d2diff-zero);
         //cout<<"Line #"<<i++<<": obs="<<pos->d2obs<<", diff="<<pos->d2diff<<" -> epsilon="<<epstmp<<endl;
      }
      cout<<"epsilon="<<epstmp<<", dmax="<<dmax<<" ,nb="<<nb<<" ,nbcalc="<<nbCalc<<endl;
   }
   /*
   else
   {//Only stat+ the worst
      float max=-1;
      unsigned int worst=0;
      unsigned int i=0;
      for(pos=dhkl.GetPeakList().begin();pos!=dhkl.GetPeakList().end();++pos)
         if(abs(pos->d2diff)>max) {worst=i++;max=abs(pos->d2diff);}
         else i++;
      dhkl.GetPeakList()[worst].stats++;
   }
   */
   if(nbCalc==0) return 0;
   const float score=(dmax/predict_coeff)*nb/(2*epsilon*nbCalc);
   if(verbose)
   {
      dhkl.Print(cout);
      cout<<"Final score:"<<score<<", nbCalc="<<nbCalc<<" ,<epsilon>="<<epsilon<<" nb="<<nb<<" Qn="<<sqrt(dmax)<<endl;
   }
   return score;
}

///////////////////////////////////////////////// PARENTINFO ///////////////////////////////////////
ParentInfo::ParentInfo(int ndhkl, int nax)
: vdhkl_info(ndhkl), vax_info(nax), pbest_cell(0), PK_integral(0)
{}

///////////////////////////////////////////////////AXINFO ///////////////////////////////////////
ParentInfo::AxInfo::AxInfo(int n, double P, double N_with1, double N_without1)
: vK(n, 1), P(P), N_with1(N_with1), N_without1(N_without1)
{}

// ///////////////////////////////////////////////// THKLINFO ///////////////////////////////////////
// ThklInfo::ThklInfo(int count, float d2min, float d2max, float K, float dRho)
// : count(count), d2min(d2min), d2max(d2max), K(K), dRho(dRho)
// {}

/////////////////////////////////////////////////////// CellExplorer ///////////////////////////////////////

CellExplorer::CellExplorer(const PeakList &dhkl, const CrystalSystem lattice, const unsigned int nbSpurious):
//mMaxDicVolDepth changed from 6 for debug
mnpar(3),mpPeakList(&dhkl),
mLengthMin(4),mLengthMax(25),
mAngleMin(M_PI),mAngleMax(2*M_PI/3),
mVolumeMin(0),mVolumeMax(1600),
mZeroShiftMin(0),mZeroShiftMax(0),
mlattice(lattice),mCentering(LATTICE_P),mNbSpurious(nbSpurious), mMaxLevelSize(0),
mObs(0),mCalc(0),mWeight(0),mDeriv(0),mBestScore(0.0),
mMinScoreReport(10),mMaxDicVolDepth(6),mDicVolDepthReport(6),
mNbLSQExcept(0), mdatabase("ICSD")
{
   this->Init();
   this->LoadPrior();
}


void CellExplorer::SetLengthMinMax(const float min,const float max)
{
   mLengthMin=min;
   mLengthMax=max;
}
void CellExplorer::SetAngleMinMax(const float min,const float max)
{
   mAngleMin=min;
   mAngleMax=max;
}
void CellExplorer::SetVolumeMinMax(const float min,const float max)
{
   mVolumeMin=min;
   mVolumeMax=max;
}
void CellExplorer::SetNbSpurious(const unsigned int nb)
{
   mNbSpurious=nb;
}
void CellExplorer::SetMaxLevelSize(const unsigned int nb)
{
   mMaxLevelSize=nb;
}
void CellExplorer::SetMinMaxZeroShift(const float min,const float max)
{
   mZeroShiftMin=min;
   mZeroShiftMax=max;
}

void CellExplorer::SetCrystalSystem(const CrystalSystem system)
{
   mlattice=system;
}

void CellExplorer::SetCrystalCentering(const CrystalCentering cent)
{
   mCentering=cent;
}

void CellExplorer::SetD2Error(const float err){mD2Error=err;}

const string& CellExplorer::GetClassName() const
{
   const static string className="CellExplorer";
   return className;
}
const string& CellExplorer::GetName() const
{
   const static string name="Some CellExplorer Object";
   return name;
}
void CellExplorer::Print() const
{
   this->RefinableObj::Print();
}
unsigned int CellExplorer::GetNbLSQFunction() const
{return 1;}

const CrystVector_REAL& CellExplorer::GetLSQCalc(const unsigned int) const
{
   VFN_DEBUG_ENTRY("CellExplorer::GetLSQCalc()",2)
   unsigned int j=0;
   for(vector<PeakList::hkl>::const_iterator pos=mpPeakList->GetPeakList().begin();pos!=mpPeakList->GetPeakList().end();++pos)
   {
      if(pos->isIndexed)
         mCalc(j++)=mRecUnitCell.hkl2d(pos->h,pos->k,pos->l);
   }
   //cout<<__FILE__<<":"<<__LINE__<<"LSQCalc : Score:"<<Score(*mpPeakList,mRecUnitCell,mNbSpurious,false,true,false)<<endl;
   VFN_DEBUG_EXIT("CellExplorer::GetLSQCalc()",2)
   return mCalc;
}
const CrystVector_REAL& CellExplorer::GetLSQObs(const unsigned int) const
{
   VFN_DEBUG_MESSAGE("CellExplorer::GetLSQObs()",2)
   return mObs;
}
const CrystVector_REAL& CellExplorer::GetLSQWeight(const unsigned int) const
{
   VFN_DEBUG_MESSAGE("CellExplorer::GetLSQWeight()",2)
   //:TODO: exclude the worst points (user-chosen number)
   return mWeight;
}
const CrystVector_REAL& CellExplorer::GetLSQDeriv(const unsigned int, RefinablePar &refpar)
{
   VFN_DEBUG_ENTRY("CellExplorer::GetLSQDeriv()",2)
   REAL *par=NULL;
   if(refpar.GetName()=="Reciprocal unit cell par #0") par=mRecUnitCell.par+1;
   else
      if(refpar.GetName()=="Reciprocal unit cell par #1") par=mRecUnitCell.par+2;
      else
         if(refpar.GetName()=="Reciprocal unit cell par #2") par=mRecUnitCell.par+3;
         else
            if(refpar.GetName()=="Reciprocal unit cell par #3") par=mRecUnitCell.par+4;
            else
               if(refpar.GetName()=="Reciprocal unit cell par #4") par=mRecUnitCell.par+5;
               else
                  if(refpar.GetName()=="Reciprocal unit cell par #5") par=mRecUnitCell.par+6;
                  else
                     if(refpar.GetName()=="Zero") par=mRecUnitCell.par+0;
                     else cout<<__FILE__<<":"<<__LINE__<<":Parameter not found:"<<refpar.GetName()<<endl;
   unsigned int j=0;
   for(vector<PeakList::hkl>::const_iterator pos=mpPeakList->GetPeakList().begin();pos!=mpPeakList->GetPeakList().end();++pos)
   {
      VFN_DEBUG_MESSAGE("CellExplorer::GetLSQDeriv():"<<j<<"/"<<mpPeakList->GetPeakList().size(),2)
      VFN_DEBUG_MESSAGE("CellExplorer::GetLSQDeriv():"<<pos->h<<","<<pos->k<<","<<pos->l,2)
      if(pos->isIndexed)
         mDeriv(j++)=mRecUnitCell.hkl2d(pos->h,pos->k,pos->l,par);
   }
   VFN_DEBUG_EXIT("CellExplorer::GetLSQDeriv()",2)
   return mDeriv;
}

void CellExplorer::BeginOptimization(const bool allowApproximations, const bool enableRestraints)
{
   VFN_DEBUG_ENTRY("CellExplorer::BeginOptimization()",5)
   Score(*mpPeakList,mRecUnitCell,mNbSpurious,false,true,false);
   const unsigned int nb=mpPeakList->GetPeakList().size();
   mCalc.resize(nb-mNbSpurious);
   mObs.resize(nb-mNbSpurious);
   mWeight.resize(nb-mNbSpurious);
   mDeriv.resize(nb-mNbSpurious);
   int j=0;
   float thres=0.0;
   for(vector<PeakList::hkl>::const_iterator pos=mpPeakList->GetPeakList().begin();pos!=mpPeakList->GetPeakList().end();++pos)
      if(thres<pos->iobs) thres=pos->iobs;
   thres/=10;// weight=1 for intensities up to Imax/10

   //cout <<"Beginning optimization with reflexions:"<<endl;
   //char buf[100];
   for(vector<PeakList::hkl>::const_iterator pos=mpPeakList->GetPeakList().begin();pos!=mpPeakList->GetPeakList().end();++pos)
   {
      if(pos->isIndexed)
      {
         mObs(j)=pos->d2obs;
         if(mObs(j)>thres) mWeight(j)=1;
         else mWeight(j)=mObs(j)/thres;
         /*
         sprintf(buf,"#%2d  (%3d %3d %3d) dobs=%6.3f dcalc=%6.3f iobs=%6.3f weight=%6.4f",
                 i,mpPeakList->mvHKL[i].h,mpPeakList->mvHKL[i].k,mpPeakList->mvHKL[i].l,
                 1/mpPeakList->mvdobs[i],1/sqrt(mRecUnitCell.hkl2d(mpPeakList->mvHKL[i].h,mpPeakList->mvHKL[i].k,mpPeakList->mvHKL[i].l)),
                 mpPeakList->mviobs[i],mWeight(j));
         */
         j++;
      }
      /*
      else
      {
         sprintf(buf,"#%2d  (%3d %3d %3d) dobs=%6.3f dcalc=%6.3f iobs=%6.3f               SPURIOUS",
                 i,mpPeakList->mvHKL[i].h,mpPeakList->mvHKL[i].k,mpPeakList->mvHKL[i].l,
                 1/mpPeakList->mvdobs[i],1/sqrt(mRecUnitCell.hkl2d(mpPeakList->mvHKL[i].h,mpPeakList->mvHKL[i].k,mpPeakList->mvHKL[i].l)),
                 mpPeakList->mviobs[i]);
      }
      cout<<buf<<endl;
      */
   }
   this->RefinableObj::BeginOptimization(allowApproximations,enableRestraints);
   VFN_DEBUG_EXIT("CellExplorer::BeginOptimization()",5)
}

void CellExplorer::LSQRefine(int nbCycle, bool useLevenbergMarquardt, const bool silent)
{
   if(mNbLSQExcept>100)
   {
      if(!silent) cout<<"CellExplorer::LSQRefine(): LSQ was disabled, too many (>100) exceptions caught. Weird unit cell parameters ?";
      return;
   }
   VFN_DEBUG_ENTRY("CellExplorer::LSQRefine()",5)
   mLSQObj.SetRefinedObj(*this);
   mLSQObj.PrepareRefParList(true);
   //this->BeginOptimization();
   //cout<<FormatVertVector<REAL>(this->GetLSQObs(0),this->GetLSQCalc(0),this->GetLSQWeight(0),this->GetLSQDeriv(0,this->GetPar((long)0)))<<endl;
   try {mLSQObj.Refine(nbCycle,useLevenbergMarquardt,silent);}
   catch(const ObjCrystException &except)
   {
      if(++mNbLSQExcept>100) cout<<"WARNING CellExplorer::LSQRefine(): LSQ was disabled, too many (>100) exceptions caught. Weird unit cell parameters ?"<<endl ;
   }
   if(!silent) mpPeakList->Print(cout);
   VFN_DEBUG_EXIT("CellExplorer::LSQRefine()",5)
}

/** Number of reflexions found in the intervals calculated between par+dpar and par-dpar
*
* \param useStoredHKL:
*    - if equal to 0, explore all possible hkl values to find possible Miller indices.
*    - if useStoredHKL=1, use the Miller indices already stored in hkl.vDicVolHKL
*    for each observed line as the only possible indices.
*    - if useStoredHKL=2, search all the possible Miller indices for all reflections
*    and store them in hkl.vDicVolHKL for each observed line.
* \param maxNbMissingBelow5: the maximum number of lines that have been calculated before
* the d-value of the 5th observed line, but which have not been observed.
* If nbMissingBelow5>0 and more than nbMissingBelow5 lines have not been observed,
* return false. Recommended to speed up triclinic, with nbMissingBelow5=5
*/

bool check_wanted(const RecUnitCell &par,const RecUnitCell &dpar)
{
   //for debug - sample1. Only [1:5] are relevant
   //optimized_cell
   // REAL wanted_cell[7] = {-3.53719e-06, 0.117223, 0.1352, 0.0968922, 0.0233423, 0, 0};
   //unoptimized best cell from Fox
   // REAL wanted_cell[7] = {0, 0.11726, 0.135243, 0.0967941, 0.0228798, 0, 0};
   //standard cell
   // REAL wanted_cell[7] = {0, 0.11723810495687369, 0.13520463280123687, 0.09685758971154582, 0.023330625962932387, 0, 0};
   // ghost K
   REAL wanted_cell[7] = {0, 0.117239, 0.135195, 0.0968546, 0.0234302, 0, 0};
   bool is_wanted = true;
   for (int i = 1; i < 5; ++i)
   {
      if (abs(par.par[i] - wanted_cell[i]) > dpar.par[i]) return false;
   }
   return true;
}

// clean. useless if  pparent_info_old==0? If so, remove default
ParentInfo make_parent_info(const PeakList &dhkl,  ParentInfo const * pparent_info_old)
{
   const unsigned int ndhkl=dhkl.GetPeakList().size();
   const unsigned int nax=dhkl.vax.size();
   
   ParentInfo parent_info_new(ndhkl, nax);

   if (pparent_info_old)
   {
      for (int idhkl = 0; idhkl != ndhkl; ++idhkl)
      // for(pos=dhkl.GetPeakList().begin(), pdhkl_info=pparent_info->vdhkl_info.begin(); pos!=dhkl.GetPeakList().end(); ++pos, ++pdhkl_info)
      {
         PeakList::hkl const * pos = &dhkl.GetPeakList()[idhkl];
         ParentInfo::DhklInfo const * pdhkl_info_old = &pparent_info_old->vdhkl_info[idhkl];
         ParentInfo::DhklInfo       * pdhkl_info_new = &parent_info_new.vdhkl_info[idhkl];
         // is it really optimization?
         pdhkl_info_new->vpthkl.reserve(pdhkl_info_old->vpthkl.size());
         for (vector<PeakList::thkl *>::const_iterator ppthkl_old=pdhkl_info_old->vpthkl.begin(); ppthkl_old!=pdhkl_info_old->vpthkl.end(); ++ppthkl_old)
         {
            PeakList::thkl * pthkl = *ppthkl_old;
            if((pos->d2obsmax >= pthkl->d2min) && (pthkl->d2max >= pos->d2obsmin))
            {
               pdhkl_info_new->vpthkl.push_back(pthkl);
            }
         }
      }
   }
   
   vector<PeakList::Ax>::const_iterator pax;
   vector<ParentInfo::AxInfo>::iterator pax_info;
   for (pax = dhkl.vax.begin(), pax_info=parent_info_new.vax_info.begin(); pax != dhkl.vax.end(); ++pax, ++pax_info)
   {
      pax_info->N_with1 = pax->N_with1;
      pax_info->N_without1 = pax->N_without1;
      pax_info->P = pax->P;
      pax_info->vK.resize(pax->vthkl.size());
      // clear - rename pthkl as ppthkl
      vector<PeakList::thkl *>::const_iterator pthkl = pax->vthkl.begin();
      vector<double>::iterator pax_info_vK_end = pax_info->vK.end();
      for (vector<double>::iterator pK = pax_info->vK.begin(); pK !=pax_info_vK_end; ++pK, ++pthkl)
      {
         *pK = (*pthkl)->K;
      }

   }

   return parent_info_new;
}

ParentInfo init_block(const PeakList &dhkl, const RecUnitCell &par,const RecUnitCell &dpar ,const bool verbose=false)
{
   const unsigned int ndhkl=dhkl.GetPeakList().size();
   vector<PeakList::hkl>::const_iterator pos,first,last,end;
   // vector<ParentInfo::DhklInfo>::iterator pdhkl_info;
   // vector<PeakList::thkl>::const_iterator pthkl;

   for(pos=dhkl.GetPeakList().begin();pos!=dhkl.GetPeakList().end();++pos)
   {
      pos->isIndexed=false;
   }
   dhkl.vthkl.clear();
   dhkl.vax.clear();

   int h,k,l;
   float dmax=dhkl.GetPeakList()[ndhkl-1].d2obs;
   float dmin=dhkl.GetPeakList()[0   ].d2obs;


   int sk0,sl0;// do we need >0 *and* <0 indices for k,l ?
   switch(par.mlattice)
   {
      case TRICLINIC:
         sk0=-1;sl0=-1;
   break;
      case MONOCLINIC:
      {
         sk0=1;sl0=-1;
         break;
      }
      case ORTHOROMBIC:
         sk0=1;sl0=1;
         break;
      case HEXAGONAL:
         sk0=-1;sl0=1;
         break;
      case RHOMBOHEDRAL:
         sk0=-1;sl0=-1;
         break;
      case TETRAGONAL:
         sk0=1;sl0=1;
         break;
      case CUBIC:
         sk0=1;sl0=1;
         break;
      // This should never happen.  Avoid using unitialized values.
      default:
         throw 0;
   }
   int stepk,stepl;// steps in k,l to use for centered lattices
   switch(par.mCentering)
   {
      case LATTICE_P:stepk=1;stepl=1;break;
      case LATTICE_I:stepk=1;stepl=2;break;
      case LATTICE_A:stepk=1;stepl=2;break;
      case LATTICE_B:stepk=1;stepl=2;break;
      case LATTICE_C:stepk=2;stepl=1;break;
      case LATTICE_F:stepk=2;stepl=2;break;
      // This should never happen.  Avoid using unitialized values.
      default: throw 0;
   }
   //RecUnitCell par0(par),par1(par);
   //for(unsigned int i=0;i<7;++i) {par0.par[i]-=dpar.par[i];par1.par[i]+=dpar.par[i];}

   //currently first & last unindexed dhkl
   first=dhkl.GetPeakList().begin(),last=dhkl.GetPeakList().end(),end=dhkl.GetPeakList().end();

   unsigned long nbCalcH,nbCalcK;// Number of calculated lines below dmax for each h,k
   
   //index search cycles
   for(h=0;;++h)
   {
      if(verbose) cout<<"H="<<h<<endl;
      nbCalcH=0;
      for(int sk=sk0;sk<=1;sk+=2)
      {
         if(h==0) sk=1;
         if(stepk==2) k=(h%2);// For LATTICE_C,LATTICE_F: h odd => k odd
         else k=0;
         for(;;k+=stepk)
         {
            if(verbose) cout<<"K="<<k*sk<<endl;
            nbCalcK=0;
            for(int sl=sl0;sl<=1;sl+=2)
            {
               int l0=0;
               if((h+k)==0)
               {
                  sl=1;// No need to list 0 0 l with l<0
                  l0=1;
               }
               else
               {
                  if(h==0)
                  {
                     if(par.mlattice==MONOCLINIC) sl=1;// 0 k l and 0 k -l are equivalent
                     if((sk<0)||(sl<0)) l0=1;// Do not list 0 k 0 with k<0
                     else l0=0;// h==k==0 already covered
                  }
                  else
                  {
                     if(sl<0) l0=1;// Do not list h k 0 twice
                     else l0=0;
                  }
               }
               if(stepl==2)
               {
                  if(par.mCentering==LATTICE_I) l0+=(h+k+l0)%2;
                  if(par.mCentering==LATTICE_A) l0+=(k+l0)%2;// Start at k+l even
                  if(  (par.mCentering==LATTICE_B)
                     ||(par.mCentering==LATTICE_F)) l0+=(h+l0)%2;// Start at h+l even
               }
               if(verbose) cout<<"SL="<<sl<<", L0="<<l0<<", STEPL="<<stepl<<", Centering="<<par.mCentering<<endl;
               for(l=l0;;l+=stepl)
               {
                  if(verbose) cout<<"L="<<l<<","<<sl<<endl;
                  float d0,d1;
                  par.hkl2d_delta(h,sk*k,sl*l,dpar,d0,d1);
                  if(d0<dmax) {nbCalcH++;nbCalcK++;}
                  if(d0>dmax)
                  {
                     if(par.mlattice==TRICLINIC)
                     {
                        // Must check that d is increasing with l, otherwise we still need to increase it
                        if(verbose) cout<<"L?="<< par.hkl2d(h,sk*k,sl*l,NULL,3)*sl <<", dmax="<<dmax<<endl;
                        if((par.hkl2d(h,sk*k,sl*l,NULL,3)*sl)>0) break;
                     }
                     else break;
                  }
                  //common block
                  dhkl.vthkl.push_back(PeakList::thkl(h, sk*k, sl*l, 0, d0, d1));
               }
            }
            if(nbCalcK==0) break;// && ((par.hkl2d(h,sk*k,0,NULL,2)*sk)>0)) break; // k too large
         }
      }
      if(nbCalcH==0) break;//h too large
   }
   
   //make axes
   sort(dhkl.vthkl.begin(), dhkl.vthkl.end());
   // map hkl0 to corresponding pointer pthkl in vthkl
   map<PeakList::hkl0, vector<PeakList::thkl *> > map_axes;
   for (vector<PeakList::thkl>::iterator pthkl = dhkl.vthkl.begin(); pthkl != dhkl.vthkl.end(); ++pthkl)
   {
      int p = boost::integer::gcd(boost::integer::gcd(pthkl->h, pthkl->k), pthkl->l);
      PeakList::hkl0 base_hkl(pthkl->h/p, pthkl->k/p, pthkl->l/p);
      map_axes[base_hkl].push_back(&*pthkl);
   }

   for (map<PeakList::hkl0, vector<PeakList::thkl *> >::const_iterator pax = map_axes.begin(); pax != map_axes.end(); ++pax)
   {
      // if (true)
      if (pax->second.size() > 1)
      {
         PeakList::Ax ax(pax->second);
         dhkl.vax.push_back(ax);
      }
   }

   const unsigned int nax=dhkl.vax.size();
   
   ParentInfo parent_info(ndhkl, nax);

   vector<PeakList::hkl>::const_iterator dhkl_mvHKL_end = dhkl.GetPeakList().end();
   vector<PeakList::thkl>::iterator dhkl_vthkl_end = dhkl.vthkl.end();
   vector<ParentInfo::DhklInfo>::iterator pdhkl_info;
   for(pos=dhkl.GetPeakList().begin(), pdhkl_info=parent_info.vdhkl_info.begin(); pos!=dhkl_mvHKL_end; ++pos, ++pdhkl_info)
   {
      // for (vector<PeakList::thkl *>::const_iterator ppthkl_old=pdhkl_info_old->vpthkl.begin(); ppthkl_old!=pdhkl_info_old->vpthkl.end(); ++ppthkl_old)
      for(vector<PeakList::thkl>::iterator pthkl=dhkl.vthkl.begin(); pthkl !=dhkl_vthkl_end; ++pthkl)
      {
         if((pos->d2obsmax >= pthkl->d2min) && (pthkl->d2max >= pos->d2obsmin))
         {
            pdhkl_info->vpthkl.push_back(&*pthkl);
         }
      }
   }
   vector<PeakList::Ax>::const_iterator pax = dhkl.vax.begin();
   for (vector<ParentInfo::AxInfo>::iterator pax_info = parent_info.vax_info.begin(); pax_info != parent_info.vax_info.end(); ++pax_info, ++pax)
   {
      pax_info->vK = vector<double>(pax->vthkl.size(), 1);
   }

   return parent_info;
}

float my_expm1(float x)
{
    return x*(1 + x*(0.5 + x*(1.0/6 + x*(1.0/24 + x*(1.0/120 + 1.0/720*x)))));
}

bool DichoIndexed(const PeakList &dhkl, const RecUnitCell &par,const RecUnitCell &dpar,
                  const unsigned int nbUnindexed=0,const bool verbose=false,unsigned int useStoredHKL=0,
                  const unsigned int maxNbMissingBelow5=0, double * pK=0, ParentInfo const * pparent_info=0, bool pseudo_correct = true)
{
   // for debug - sample1.
   bool is_wanted = check_wanted(par, dpar);
   if (is_wanted) cout << "new DichoIndexed" << endl;
   //end debug
   // maximum dRho for a peak that we use in analysis. Corresponds to maximum Pa for an axis to be accounted in K
   double dRho_max = 1;
   const unsigned int nb=dhkl.GetPeakList().size();
   int nbIndexed=nb-nbUnindexed;// Number of reflections we require to be indexed
   float d5=0;
   if(maxNbMissingBelow5>0) d5=dhkl.GetPeakList()[4].d2obs;
   // number of missing reflections calculated below 5th observed line
   unsigned int nbMissingBelow5=0;
   // List of indexed reflections
   vector<PeakList::hkl>::const_iterator pos,first,last,end;
   vector<ParentInfo::DhklInfo>::const_iterator pdhkl_info;
   vector<PeakList::thkl>::const_iterator pthkl;
   // clean - no useStoredHKL
   if(useStoredHKL==1)
   {// We already now possible Miller indices for all reflections
      unsigned int nbUnIx = 0;
      int count_good = 0;
      int count_peaks = 0;
      for (vector<PeakList::thkl>::const_iterator pthkl = dhkl.vthkl.begin(); pthkl != dhkl.vthkl.end(); ++pthkl)
      {
         pthkl->count = 0;
         float d0,d1;
         par.hkl2d_delta(pthkl->h,pthkl->k,pthkl->l,dpar,d0,d1);
         pthkl->d2min = d0;
         pthkl->d2max = d1;
         //common block
         pthkl->K = 1;
         pthkl->dRho = -1;
         if ((0 < d0) && (d1 < dhkl.d2max))
         {
            int index1 = floorf(dhkl.coef*d0);
            int index2 = floorf(dhkl.coef*d1) + 1;
            PeakList::integral const & I1 = dhkl.integrals[index1];
            PeakList::integral const & I2 = dhkl.integrals[index2];
            //not sure if the check is needed
            pthkl->dRho = I2.Rho - I1.Rho;
            // dRho for first peak of an axis is assumed to be < 1 in K accumulation. Need changes if remove this
            ++count_peaks;
            if (pthkl->dRho < dRho_max)
            {
               ++count_good;
               // old formula
               // pthkl->K = (I2.F - I1.F + I2.Rho_ext - I1.Rho_ext) / pthkl->dRho;
               double x = min(I2.F - I1.F, (float)1.0);
               double P_post = x - (1 - x)*my_expm1(I1.Rho_ext - I2.Rho_ext);
               double P_prior = -my_expm1(-pthkl->dRho);
               // float P_post = x - (1 - x)*boost::math::expm1(I2.Rho_ext - I1.Rho_ext);
               // float P_prior = -boost::math::expm1(pthkl->dRho);
               pthkl->K = P_post/P_prior;
               pthkl->P = P_prior;
            }
         }
      }


      for(pos=dhkl.GetPeakList().begin(), pdhkl_info=pparent_info->vdhkl_info.begin(); pos!=dhkl.GetPeakList().end(); ++pos, ++pdhkl_info)
      {
         pos->isIndexed=false;
         // if (is_wanted)
         //    cout << "experimental: " << pos->d2obsmin << endl;
         
         vector<PeakList::thkl *>::const_iterator vpthkl_end = pdhkl_info->vpthkl.end();
         for (vector<PeakList::thkl *>::const_iterator ppthkl=pdhkl_info->vpthkl.begin(); ppthkl!=vpthkl_end; ++ppthkl)
         {
            // if (is_wanted)
            //    cout << (*ppthkl)->h << " " << (*ppthkl)->k << " " << (*ppthkl)->l << endl;
            PeakList::thkl * pthkl = *ppthkl;
            ++global_check_inside_count;
            if((pos->d2obsmax >= pthkl->d2min) && (pthkl->d2max >= pos->d2obsmin))
            {
               ++global_is_inside_count;
               if (!(pos->isIndexed))
                  --nbIndexed;
               pos->isIndexed = true;
               ++(pthkl->count);
            }
         }
      }

   }

   //K -> 1 correction for pseudosymmetry
   if (pseudo_correct)
   {
      for(pos=dhkl.GetPeakList().begin(), pdhkl_info=pparent_info->vdhkl_info.begin(); pos!=dhkl.GetPeakList().end(); ++pos, ++pdhkl_info)
      {
         PeakList::thkl * pbest_thin_thkl = 0;
         vector<PeakList::thkl *>::const_iterator vpthkl_end = pdhkl_info->vpthkl.end();
         for (vector<PeakList::thkl *>::const_iterator ppthkl=pdhkl_info->vpthkl.begin(); ppthkl!=vpthkl_end; ++ppthkl)
         {
            PeakList::thkl * pthkl = *ppthkl;
            if((pos->d2obsmax < pthkl->d2min) || (pthkl->d2max < pos->d2obsmin)) continue;
            if (pthkl->count == 1 && pthkl->K > 1)
            {
               if (pbest_thin_thkl)
               {
                  if (pbest_thin_thkl->K < pthkl->K)
                  {
                     // pbest_thin_thkl->K = 1;
                     pbest_thin_thkl->is_disabled = true;
                     pbest_thin_thkl = pthkl;
                  }
                  else
                  {
                     // Now K = 1 is done in axes. Not done for first axis peak
                     // pthkl->K = 1;
                     pthkl->is_disabled = true;
                  }
               }
               else
                  pbest_thin_thkl = pthkl;
            }
         }
      }
   }
   
   //get overall K
   //need to convert to double. Either everything or just main K and PK calculation
   double K = 1;
   for (vector<PeakList::thkl>::const_iterator pthkl = dhkl.vthkl.begin(); pthkl != dhkl.vthkl.end(); ++pthkl)
   {
      if (!(pthkl->is_disabled)) K *= pthkl->K;
   }

   vector<PeakList::Ax>::iterator pax;
   vector<ParentInfo::AxInfo>::const_iterator pax_info;
   for (pax = dhkl.vax.begin(), pax_info=pparent_info->vax_info.begin(); pax != dhkl.vax.end(); ++pax, ++pax_info)
   {
      double dRho = pax->vthkl.front()->dRho;
      //if the first interval is out of range, same holds for other intervals of that axis
      if (dRho > 0 && dRho < dRho_max)
      {
         // y stands for the parent axis, z for current
         double Pz1 = pax->vthkl.front()->P;
         pax->P = Pz1;
         double Kz1 = pax->vthkl.front()->K;
         double Ky1 = pax_info->vK.front();
         double N_without1_y = pax_info->N_without1;
         double N_with1_y = pax_info->N_with1;
         double K_secondary = 1;
         double K_secondary_y = 1;
         vector<PeakList::thkl *>::const_iterator ppthkl = ++pax->vthkl.begin();
         for (vector<double>::const_iterator pKy = ++pax_info->vK.begin(); ppthkl != pax->vthkl.end(); ++ppthkl, ++pKy)
         {
            if ((*ppthkl)->is_disabled)
            {
               N_with1_y /= *pKy;
               N_without1_y /= *pKy;
               (*ppthkl)->K = 1;
            }
            else
            {
               K_secondary *= (*ppthkl)->K;
               K_secondary_y *= *pKy;
            }
         }
         pax->N_without1 = (pax_info->P - Pz1)*K_secondary_y + N_without1_y;
         pax->N_with1 = (pax_info->P*Ky1 - Pz1*Kz1)*K_secondary_y + N_with1_y;
         double Nz = pax->vthkl.front()->is_disabled ? pax->N_without1 : pax->N_with1;
         // speed up - make more products and less divisions
         double Kz1_selected = pax->vthkl.front()->is_disabled ? 1 : Kz1;
         double Ka_corr = 1/(Pz1*Kz1_selected*K_secondary + Nz);
         double Ka = Kz1_selected*K_secondary*Ka_corr;
         // original
         K *= Ka_corr;
         // we assume Pa < 1 because dRho < 0.5
         if (is_wanted)
         {
            int ax_number = pax - dhkl.vax.begin();
            int num_peaks = pax->vthkl.size();

            if (pax_info->P*Ky1 - Pz1*Kz1 <= -1.0e-15)
            {
               cout << "2479: strange probability:" << endl;
               printf("Py*Ky=%f*%f=%f; Pz*Kz=%f*%f=%f\n", pax_info->P, Ky1, pax_info->P*Ky1, Pz1, Kz1, Pz1*Kz1);
               assert(0);
            }
            // printf(", Pa*Ka=%f*%f=%f\n", Pa, Ka, Pa*Ka);
         }
      }
   }
   if(verbose)
   {
      dhkl.Print(cout);
   }
   if (pK != 0)
      *pK = K;
   return nbIndexed<=0;
}

double CellExplorer::GetBestScore()const{return mBestScore;}
const std::list<std::pair<RecUnitCell,double> >& CellExplorer::GetSolutions()const {return mvSolution;}
std::list<std::pair<RecUnitCell,double> >& CellExplorer::GetSolutions() {return mvSolution;}

double get_cs_probability(CrystalSystem cs)
{
   switch(cs)
   {
   case TRICLINIC:     return 0.0398;
   case MONOCLINIC:    return 0.1679;
   case ORTHOROMBIC:   return 0.2062;
   case HEXAGONAL:     return 0.1102;
   // In my notes its for trigonal. But I probably meant rhombohedral. It's a problem for hexagonal also
   case RHOMBOHEDRAL:  return 0.0949;
   case TETRAGONAL:    return 0.1548;
   case CUBIC:         return 0.2261;
   }
}

unsigned int CellExplorer::RDicVol(RecUnitCell par0,RecUnitCell dpar, unsigned int depth,unsigned long &nbCalc,const float minV,
                                    const float maxV,vector<unsigned int> vdepth, ParentInfo const  * pparent_info)
{
   // not synchronized with analogous constant. Make it CellExplorer member!
   int control_depth = 6;

   mpPeakList->mvcount[depth]++;
   
   bool is_wanted = check_wanted(par0, dpar);

   bool is_final_run = (mpPeakList->active_depth == mMaxDicVolDepth);

   if (is_wanted)
   {
      cout << endl << "new RDicVol, minV: " << minV << " depth: " << depth << endl;
      cout << "number: " << mpPeakList->mvcount[depth] << endl;
   }

   if ((depth != mpPeakList->active_depth) && !mpPeakList->mvvis_good[depth][mpPeakList->mvcount[depth]])
   {  
      if (is_wanted) cout << "not_calculated" << endl;
      if (depth > 0)
      {
         pparent_info->PK_integral += pparent_info->PK/pow(2.0, mnpar - 1);
      }
      return 0;
   }
   if (is_wanted) cout << "calculated" << endl;

   // in case volume check will fail
   // mpPeakList->mvvcriterion[depth].push_back(-20000);
   if (!mpPeakList->mis_exta_run) mpPeakList->mvvPK[depth].push_back(0.0);


   static bool localverbose=false;
   if(mlattice==TRICLINIC)
   {
      const float p1=par0.par[1]    , p2=par0.par[2]    , p3=par0.par[3]    , p4=par0.par[4]    , p5=par0.par[5]    , p6=par0.par[6];
      const float p1m=p1-dpar.par[1], p2m=p2-dpar.par[2], p3m=p3-dpar.par[3], p4m=p4-dpar.par[4], p5m=p5-dpar.par[5], p6m=p6-dpar.par[6];
      const float p1p=p1+dpar.par[1], p2p=p2+dpar.par[2], p3p=p3+dpar.par[3], p4p=p4+dpar.par[4], p5p=p5+dpar.par[5], p6p=p6+dpar.par[6];

      // a*<b*<c*
      if((p1m>p2p)||(p2m>p3p)) return 0;

      // max/min absolute values for p4,p5,p6
      if((p4m>p1p)||(-p4p>p1p)) return 0;//abs(p4)<p1 <=> b* < b*+/-a*
      if((p5m>p2p)||(-p5p>p2p)) return 0;//abs(p5)<p2 <=> c* < c*+/-b*
      if((p6m>p1p)||(-p6p>p1p)) return 0;//abs(p6)<p1 <=> c* < c*+/-a*

      const float max6=p1p+p2p-p4m-p5m;
      if((p6m>max6)||(-p6p>max6)) return 0;//abs(p6)<p1+p2-p4-p5 <=> c* < c*+/-a*+/-b*

      float p6mm,p6pp,p5mm,p5pp,p4mm,p4pp; // p6pp: smaller V*, larger V, etc..
      // derivative of V*^2 with p6
      if((p4*p5-2*p2*p6)>0) {p6pp=p6p;p6mm=p6m;}
      else{p6pp=p6m;p6mm=p6p;}
      // derivative of V*^2 with p5
      if((p4*p6-2*p1*p5)>0) {p5pp=p5p;p5mm=p5m;}
      else{p5pp=p5m;p5mm=p5p;}
      // derivative of V*^2 with p5
      if((p5*p6-2*p3*p4)>0) {p4pp=p4p;p4mm=p4m;}
      else{p4pp=p4m;p4mm=p4p;}

      //const float vmin0=1/sqrt(abs(p1p*p2p*p3p*(1-p5mm*p5mm/(4*p2p*p3p)-p6mm*p6mm/(4*p1p*p3p)-p4mm*p4mm/(4*p1p*p2p)+p4mm*p5mm*p6mm/(4*p1m*p2m*p3m))));
      //const float vmax0=1/sqrt(abs(p1m*p2m*p3m*(1-p5pp*p5pp/(4*p2m*p3m)-p6pp*p6pp/(4*p1m*p3m)-p4pp*p4pp/(4*p1m*p2m)+p4pp*p5pp*p6pp/(4*p1m*p2m*p3m))));
      const float vmin0=1/sqrt(abs(p1p*p2p*p3p*(1-p5mm*p5mm/(4*p2p*p3p)-p6mm*p6mm/(4*p1p*p3p)-p4mm*p4mm/(4*p1p*p2p)+p4mm*p5mm*p6mm/(4*p1m*p2m*p3m))));
      const float vmax0=1/sqrt(abs(p1m*p2m*p3m*(1-p5pp*p5pp/(4*p2m*p3m)-p6pp*p6pp/(4*p1m*p3m)-p4pp*p4pp/(4*p1m*p2m)+p4pp*p5pp*p6pp/(4*p1m*p2m*p3m))));
      const float v0=1/sqrt(abs(p1*p2*p3*(1-p5*p5/(4*p2*p3)-p6*p6/(4*p1*p3)-p4*p4/(4*p1*p2)+p4*p5*p6/(4*p1*p2*p3))));
      if((vmin0>maxV)||(vmax0<minV)) return 0;

      // P recalculation
      // //not used currently
      // const float vav0=1/sqrt(abs(p1*p2*p3*(1-p5*p5/(4*p2*p3)-p6*p6/(4*p1*p3)-p4*p4/(4*p1*p2)+p4*p5*p6/(4*p1*p2*p3))));

      // for p4 negative values are possible.  We need the one smaller in magnitude
      float p4s = (abs(p4m) < abs(p4p)) ? p4m : p4p;
      float ap = vmax0*sqrt(p2p*p3p - float(0.25)*p5m*p5m);
      float bp = vmax0*sqrt(p1p*p3p - float(0.25)*p6m*p6m);
      float cp = vmax0*sqrt(p1p*p2p - float(0.25)*p4s*p4s);
      // float tmp = float(0.25)*vmax0*vmax0*vmax0/(ap*bp*cp);
      // float jacobian = tmp*tmp*tmp/vmax0;
      // float search_volume = 64*dpar.par[1]*dpar.par[2]*dpar.par[3]*dpar.par[4]*dpar.par[5]*dpar.par[6];
      
      // calculation of cos of direct angles
      // copied with slight modifications from DirectUnitCell()
      // a slightly more efficient variant exists (p.21)
      float aa,bb,cc,calphaa,cbetaa,cgammaa,salphaa,sbetaa,sgammaa;
      aa=sqrt(p1);
      bb=sqrt(p2);
      cc=sqrt(p3);
      calphaa=p5/(2*bb*cc);
      cbetaa =p6/(2*aa*cc);
      cgammaa=p4/(2*aa*bb);
      salphaa=sqrt(abs(1-calphaa*calphaa));
      sbetaa =sqrt(abs(1-cbetaa *cbetaa));
      sgammaa=sqrt(abs(1-cgammaa*cgammaa));
      const float vv=sqrt(abs(1-calphaa*calphaa-cbetaa*cbetaa-cgammaa*cgammaa+2*calphaa*cbetaa*cgammaa));
      const float a0=salphaa/(aa*vv);
      const float b0=sbetaa /(bb*vv);
      const float c0=sgammaa/(cc*vv);
      const float calpha=(cbetaa *cgammaa-calphaa)/(sbetaa *sgammaa);
      const float cbeta =(calphaa*cgammaa-cbetaa )/(salphaa*sgammaa);
      const float cgamma=(calphaa*cbetaa -cgammaa)/(salphaa*sbetaa );
      // positive cos used for pdf estimation
      const float cos_alpha = min(abs(calpha), mCosAngMax);
      const float cos_beta = min(abs(cbeta), mCosAngMax);
      const float cos_gamma = min(abs(cgamma), mCosAngMax);
      float tmp = float(0.25)*v0*v0*v0/(a0*b0*c0);
      float jacobian = tmp*tmp*tmp/v0;
      float search_volume = 64*dpar.par[1]*dpar.par[2]*dpar.par[3]*dpar.par[4]*dpar.par[5]*dpar.par[6];
      // multiplier = 3 = 6 permutations of (a, b, c) / 2 signs of cos_gamma
      par0.P = 3*jacobian*search_volume*max(mPrior["a"][ap], mPrior["a"][a0])*max(mPrior["a"][bp], mPrior["a"][b0])*max(mPrior["a"][cp], mPrior["a"][c0])* \
                                    mPrior["cos_alpha"][cos_alpha]*mPrior["cos_alpha"][cos_beta]*mPrior["cos_alpha"][cos_gamma];
      par0.P *= get_cs_probability(mlattice);



   }
   else
      //debug! change back to 2
      if(depth<=2)// test if volume is within range and calculate P from scratch
      {
         RecUnitCell parm=par0,parp=par0;
         for(unsigned int i=0;i<6;++i) {parm.par[i]-=dpar.par[i];parp.par[i]+=dpar.par[i];}
         vector<float> parmd=parm.DirectUnitCell();
         vector<float> parpd=parp.DirectUnitCell();
         if((parpd[6]>maxV)||(parmd[6]<minV))return 0;
         switch (mlattice)
         {
         case MONOCLINIC:
         {
            float jacobian = float(0.25)*(parpd[6] + parmd[6])*(parpd[6] + parmd[6]);
            float search_volume = 16*dpar.par[1]*dpar.par[2]*dpar.par[3]*dpar.par[4];
            float a = parmd[0];
            float b = parmd[1];
            float c = parmd[2];
            float cos_beta = abs(par0.par[4]) - abs(dpar.par[4]);
            cos_beta = max(cos_beta, (float)0.0);


            // multiplier = 2 permutations of (a, c)
            par0.P = 2*jacobian*search_volume*mPrior["a"][a]*mPrior["b"][b]*mPrior["a"][c]*mPrior["cos_beta"][cos_beta];
            break;
         }
         case ORTHOROMBIC:
         {
            float jacobian = parpd[6]*parpd[6];
            float search_volume = 8*dpar.par[1]*dpar.par[2]*dpar.par[3];
            float a = parmd[0];
            float b = parmd[1];
            float c = parmd[2];
            // multiplier = 6 permutations of (a, b, c)
            par0.P = 6*jacobian*search_volume*mPrior["a"][a]*mPrior["a"][b]*mPrior["a"][c];
            break;
         }
         case HEXAGONAL:
         {
            float a = parmd[0];
            float c = parmd[2];
            float jacobian = sqrtf(3)/2*a*a*c*c;
            float search_volume = 4*dpar.par[1]*dpar.par[2];
            // multiplier = 1
            par0.P = jacobian*search_volume*mPrior["a"][a]*mPrior["c"][c];
            break;
         }
         case RHOMBOHEDRAL:
         {
            // r for rhombohedral axes, h - for hexagonal
            vector<float> par0d=par0.DirectUnitCell();
            const float p1 = par0.par[1];
            const float p2 = par0.par[2];
            // const float cos_alphaa = p2/p1/p1;
            // const float cos_alpha = -cos_alphaa/(1 + cos_alphaa);
            const float ar0 = par0d[0];
            // remember cos_alpha is negative. Understand m, p algebraically
            const float cos_alpha0 = cos(par0d[4]);
            // const float cos_alpham = cos(parpd[4]);
            // const float cos_alphap = cos(parmd[4]);
            const float jacobian = 3*sqrtf((float)1.5)*(1 + cos_alpha0)*(1 + cos_alpha0)*ar0*ar0*ar0/(p1*p1);
            const float search_volume = 4*dpar.par[1]*dpar.par[2];
            RecUnitCell big_a_par = par0;
            big_a_par.par[1] += dpar.par[1];
            if (big_a_par.par[2] > dpar.par[2]) big_a_par.par[2] -= dpar.par[2];
            vector<float> big_a_pard = big_a_par.DirectUnitCell();
            const float arp = big_a_pard[0];
            // not necessarily the biggest possible value. but can be bigger than that. Just some elevation above the average to pass the pdf rise
            const float ah0 = sqrt(2)*ar0*sqrt(1 - cos_alpha0);
            const float ahp = sqrt(2)*arp*sqrt(1 - cos_alpha0);
            const float ch0 = sqrt(3)*ar0*sqrt(1 + 2*cos_alpha0);
            const float chp = sqrt(3)*arp*sqrt(1 + 2*cos_alpha0);
            // multiplier = 1
            par0.P = jacobian*search_volume*max(mPrior["a"][ah0], mPrior["a"][ahp])*max(mPrior["c"][ch0], mPrior["c"][chp]);
            break;
         }
         case TETRAGONAL:
         {
            const float a0 = 1/par0.par[1];
            const float c0 = 1/par0.par[2];
            const float jacobian = a0*a0*c0*c0;
            const float search_volume = 4*dpar.par[1]*dpar.par[2];
            assert(dpar.par[1] < par0.par[1]);
            assert(dpar.par[2] < par0.par[2]);
            const float ap = 1/(par0.par[1] - dpar.par[1]);
            const float cp = 1/(par0.par[2] - dpar.par[2]);
            par0.P = jacobian*search_volume*max(mPrior["a"][a0], mPrior["a"][ap])*max(mPrior["c"][c0], mPrior["c"][cp]);
            break;
         }
         case CUBIC:
         {
            const float a0 = 1/par0.par[1];
            const float jacobian = a0*a0;
            const float search_volume = 2*dpar.par[1];
            assert(dpar.par[1] < par0.par[1]);
            const float ap = 1/(par0.par[1] - dpar.par[1]);
            par0.P = jacobian*search_volume*max(mPrior["a"][a0], mPrior["a"][ap]);
            break;
         }
         default:
            throw 0;
         }
         par0.P *= get_cs_probability(mlattice);
      }
   //debug
   if (depth == 0)
   {
      // check par0.P is not NAN
      assert(par0.P == par0.P);
      if (par0.P == par0.P)
      {
         global_p += par0.P;
         global_pmax = max(global_pmax, par0.P);
      }
   }


   unsigned int useStoredHKL=1;//Use already stored hkl
   // if(depth==0) useStoredHKL=2; //Store possible hkl for all observed lines

   unsigned int maxMissingBelow5=0;
   // In the triclinic case, accept a maximum of 5 missing reflections below the 5th observed line
   if(mlattice==TRICLINIC) maxMissingBelow5=5;

   double K = 0;
   ParentInfo this_parent_info;
   if (depth == 0)
   {
      mpPeakList->clear_locals();
      this_parent_info = init_block(*mpPeakList, par0, dpar);
      pparent_info = &this_parent_info;
   }

   bool indexed = false;
   indexed=DichoIndexed(*mpPeakList,par0,dpar,mNbSpurious,localverbose,useStoredHKL,maxMissingBelow5, &K, pparent_info);
   if (!indexed) K = 0.0;
   ++calls_per_level[depth];
   if (depth == mpPeakList->active_depth && !mpPeakList->mis_exta_run)
   {
      
      *mpPeakList->mvvPK[depth].rbegin() = par0.P*K;
      assert(mpPeakList->mvvPK[depth].size() == mpPeakList->mvcount[depth] + 1);
      indexed = false; // not going futher yet
   }
   else
      {
         indexed = true;
         ++repeated_calls_per_level[depth];
      }

   if (is_wanted)
   {
      cout << endl <<"minV: " << minV << " depth: " << depth << endl;
      cout << "main_block\n";
      cout << "P: " << par0.P << endl;
      cout << "K: " << K << endl;
      cout << "PK: " << par0.P*K << endl;
   }



   #if 0
   // If indexation failed but depth>=4, try adding a zero ?
   if( (!indexed) && (depth>=4))
   {//:TODO: Check if this is OK ! Vary value with depth
      dpar.par[0]=.0001;
      indexed=DichoIndexed(*mpPeakList,par0,dpar,mNbSpurious,false,useStoredHKL,maxMissingBelow5);
      //if(indexed) cout<<"Added zero - SUCCESS !"<<endl;
   }
   #endif

   // clean - useStoredHKL != 2
   if((indexed)&&(useStoredHKL==2))
   {
      // Test if two successive lines have been indexed exclusively with the same hkl
      unsigned int nbident=0;
      // for(vector<PeakList::hkl>::const_iterator pos=mpPeakList->GetPeakList().begin();pos!=mpPeakList->GetPeakList().end();)
      vector<PeakList::hkl>::const_iterator pos;
      vector<ParentInfo::DhklInfo>::const_iterator pdhkl_info;
      for(pos=mpPeakList->GetPeakList().begin(), pdhkl_info=pparent_info->vdhkl_info.begin(); pos!=mpPeakList->GetPeakList().end(); ++pos, ++pdhkl_info)
      {
         if(pdhkl_info->vpthkl.size()==1)
         // if(pos->vDicVolTHKL.size()==1)
         {
            // const PeakList::thkl * pthkl0=pos->vDicVolTHKL.front();
            const PeakList::thkl * pthkl0=pdhkl_info->vpthkl.front();
            if(++pos==mpPeakList->GetPeakList().end()) break;
            if(pdhkl_info->vpthkl.size()==1)
            {
               const PeakList::thkl * pthkl1=pdhkl_info->vpthkl.front();
               if (pthkl0 == pthkl1)
                  ++nbident;
               if(nbident>mNbSpurious) {indexed=false;break;}
            }
         }
         else ++pos;
      }
   }

   // if we can zoom in for one parameter directly, we need per-parameter depth
   if(vdepth.size()==0)
   {
      vdepth.resize(mnpar-1);
      for(vector<unsigned int>::iterator pos=vdepth.begin();pos!=vdepth.end();) *pos++=depth;
   }
   else
      for(vector<unsigned int>::iterator pos=vdepth.begin();pos!=vdepth.end();++pos) if(*pos<depth)*pos=depth;

   if(false)//(depth==1)&&(rand()%10==0))
   {
      RecUnitCell parm=par0,parp=par0;
      for(unsigned int i=0;i<4;++i) {parm.par[i]-=dpar.par[i];parp.par[i]+=dpar.par[i];}
      for(unsigned int i=4;i<7;++i) {parm.par[i]+=dpar.par[i];parp.par[i]-=dpar.par[i];}
      vector<float> parmd=parm.DirectUnitCell();
      vector<float> parpd=parp.DirectUnitCell();
      char buf[200];
      sprintf(buf,"a=%5.2f-%5.2f b=%5.2f-%5.2f c=%5.2f-%5.2f alpha=%6.2f-%6.2f beta=%6.2f-%6.2f gamma=%6.2f-%6.2f V=%6.2f-%6.2f",
               parpd[0],parmd[0],parpd[1],parmd[1],parpd[2],parmd[2],parpd[3]*RAD2DEG,parmd[3]*RAD2DEG,
               parpd[4]*RAD2DEG,parmd[4]*RAD2DEG,parpd[5]*RAD2DEG,parmd[5]*RAD2DEG,parpd[6],parmd[6]);
      for(unsigned int i = 0; i < depth; ++i)  cout << " ";
      cout<<buf<<"level="<<depth<<", indexed="<<indexed<<"("<<mvSolution.size()<<" sol.)"<<endl;
   }
  
   nbCalc++;
   // if (nbCalc % 100000 == 0) cout << nbCalc << endl;
   // :TODO: if we failed the dichotomy and reached some depth, try guessing a zero shift from the indexed reflections
   /*
   if((!indexed)&&(depth>=2))
   {
      vector<float> shifts(mpPeakList->GetPeakList().size());
      vector<PeakList::hkl>::const_iterator peakpos=mpPeakList->GetPeakList().begin();
      for(vector<float>::iterator spos=shifts.begin();spos!=shifts.end();)
      {   *spos++ = peakpos->d2diff * (float)(peakpos->isIndexed&&(!peakpos->isSpurious));peakpos++;}
      sort(shifts.begin(),shifts.end());
      par0.par[0]=shifts[mpPeakList->GetPeakList().size()/2];//use median value
      indexed=DichoIndexed(*mpPeakList,par0,dpar,mNbSpurious);
      if(indexed) cout<<"Failed Dicho ? Trying auto-zero shifting :Worked !"<<endl;
   }
   */
   




   global_pkmax = max(par0.P*K, global_pkmax);


   unsigned int deeperSolutions=0;

   ParentInfo parent_info_new;
   // dummy except for control depth
   pair<RecUnitCell, double> best_cell;
   if(indexed)
   {
      deeperSolutions += 1;
      parent_info_new = make_parent_info(*mpPeakList, pparent_info);
      parent_info_new.PK = par0.P*K;
      mpPeakList->clear_locals();
      if (is_final_run)
      {
         if (depth == control_depth)
         {
            best_cell = make_pair(par0, K);
            parent_info_new.pbest_cell = &best_cell;
         }
         if (depth > control_depth)
         {
            parent_info_new.pbest_cell = pparent_info->pbest_cell;
         }
      }

      // seems redundant
      if(depth<mMaxDicVolDepth)
      {
         if(false)//depth>=5)
         {
            RecUnitCell parm=par0,parp=par0;
            for(unsigned int i=0;i<6;++i) {parm.par[i]-=dpar.par[i];parp.par[i]+=dpar.par[i];}
            vector<float> parmd=parm.DirectUnitCell();
            vector<float> parpd=parp.DirectUnitCell();
            char buf[200];
            sprintf(buf,"a=%5.2f-%5.2f b=%5.2f-%5.2f c=%5.2f-%5.2f alpha=%5.2f-%5.2f beta=%5.2f-%5.2f gamma=%5.2f-%5.2f V=%5.2f-%5.2f",
                     parpd[0],parmd[0],parpd[1],parmd[1],parpd[2],parmd[2],parpd[3]*RAD2DEG,parmd[3]*RAD2DEG,
                     parpd[4]*RAD2DEG,parmd[4]*RAD2DEG,parpd[5]*RAD2DEG,parmd[5]*RAD2DEG,parpd[6],parmd[6]);
            for(unsigned int i=0;i<depth;++i) cout<<"  ";
            cout<<"Depth="<<depth<<" "<<buf<<endl;
            for(int i=0;i<=6;++i)cout<<par0.par[i]<<",";
            for(int i=0;i<=6;++i)cout<<dpar.par[i]<<",";
            cout<<endl;
         }
         RecUnitCell par=par0;
         // zero (if used...)
         dpar.par[0]=0.5*dpar.par[0];
         // Divide interval by 2, except if this parameter is already at a higher depth
         // because a main axis has been indexed already.
         for(unsigned int i=1;i<mnpar;++i) dpar.par[i]*=(0.5+0.5*(vdepth[i-1]>depth));
         // same for P
         for(unsigned int i=1;i<mnpar;++i) par.P*=(0.5+0.5*(vdepth[i-1]>depth));

         // clean. There could be something more elegant. However this is more flexible.
         for(int i0=-1;i0<=1;i0+=2)
         {
            //:TODO: dichotomy on zero shift ?
            if(localverbose) cout<<__FILE__<<":"<<__LINE__<<":"<<par.par[3]<<" +/- "<<dpar.par[3]<<" ("<<vdepth[2]<<")"<<endl;
            // Don't change parameter if it is already determined at a higher depth
            if(vdepth[0]==depth) {par.par[1]=par0.par[1]+i0*dpar.par[1];}
            else {i0=2;}// no need to dicho this parameter which is already at higher depth
            if(mnpar==2)
               deeperSolutions+=RDicVol(par,dpar, depth+1,nbCalc,minV,maxV,vdepth, &parent_info_new);
            else
               for(int i1=-1;i1<=1;i1+=2)
               {
                  if(vdepth[1]==depth) {par.par[2]=par0.par[2]+i1*dpar.par[2];}
                  else {i1=2;}// no need to dicho this parameter which is already at higher depth
                  if(mnpar==3)
                     deeperSolutions+=RDicVol(par,dpar, depth+1,nbCalc,minV,maxV,vdepth, &parent_info_new);
                  else
                     for(int i2=-1;i2<=1;i2+=2)
                     {
                        if(vdepth[2]==depth) {par.par[3]=par0.par[3]+i2*dpar.par[3];}
                        else {i2=2;}// no need to dicho this parameter which is already at higher depth
                        if(mnpar==4)
                           deeperSolutions+=RDicVol(par,dpar, depth+1,nbCalc,minV,maxV,vdepth, &parent_info_new);
                        else
                           for(int i3=-1;i3<=1;i3+=2)
                           {
                              if(vdepth[3]==depth)par.par[4]=par0.par[4]+i3*dpar.par[4];
                              else i3=2;
                              if(mnpar==5)
                                 deeperSolutions+=RDicVol(par,dpar, depth+1,nbCalc,minV,maxV,vdepth, &parent_info_new);
                              else
                                 for(int i4=-1;i4<=1;i4+=2)
                                 {
                                    par.par[5]=par0.par[5]+i4*dpar.par[5];
                                    //if(mnpar==7)
                                    //   deeperSolutions+=RDicVol(par,dpar, depth+1,nbCalc,minV,maxV,vdepth);
                                    //else
                                       for(int i5=-1;i5<=1;i5+=2)
                                       {
                                          par.par[6]=par0.par[6]+i5*dpar.par[6];
                                          //if(localverbose) cout<<__FILE__<<":"<<__LINE__<<":"<<par.par[3]<<" +/- "<<dpar.par[3]<<" ("<<vdepth[2]<<")"<<endl;
                                          deeperSolutions+=RDicVol(par,dpar, depth+1,nbCalc,minV,maxV,vdepth, &parent_info_new);
                                       }
                                 }
                           }
                     }
               }
          }
      }
      
      // integration
      double PK_integral_max = max(par0.P*K, parent_info_new.PK_integral);
      pparent_info->PK_integral += PK_integral_max;
      // pparent_info->PK_integral +=  parent_info_new.PK_integral;



   }
   else
   {
      pparent_info->PK_integral += par0.P*K;
   }

   if (is_final_run)
   {
      if (depth == control_depth && indexed)
      {
         if (mpPeakList->mis_exta_run)
         {
            if (parent_info_new.PK_integral >= mpPeakList->mPK_solution_critical)
            {
               mvSolution.push_back(make_pair(parent_info_new.pbest_cell->first, parent_info_new.PK_integral));
            }
         }
         else
         {
            mpPeakList->mvPK_integral.push_back(parent_info_new.PK_integral);
         }
      }
      if (depth > control_depth && K > pparent_info->pbest_cell->second)
      {
         *pparent_info->pbest_cell = make_pair(par0, K);
      }
   }


   if (mpPeakList->mis_exta_run && depth == control_depth)
   {
      if (indexed)
      {
         global_pk_integralmax = max(global_pk_integralmax, parent_info_new.PK_integral);
      }
      if (is_wanted)
      {
         cout << "wanted PK_integral for control depth: " << parent_info_new.PK_integral << endl;
         cout << "wanted best K: " << parent_info_new.pbest_cell->second << endl;
      }
   }

   if (!indexed) mpPeakList->clear_locals();
   return deeperSolutions;
}

vector<float> linspace(float min, float max,unsigned int nb)
{
   vector<float> v(nb);
   for(unsigned int i=0;i<nb;++i) v[i]=min+(max-min)*i/(nb-1);
   return v;
}

void CellExplorer::DicVol(const float minScore,const unsigned int minDepth,const float stopOnScore,const unsigned int stopOnDepth)
{
   mNbLSQExcept=0;
   mDicVolDepthReport=minDepth;
   mMinScoreReport=minScore;
   this->Init();
   if(minDepth>mMaxDicVolDepth) mMaxDicVolDepth=minDepth;
   mvNbSolutionDepth.resize(mMaxDicVolDepth+1);
   for(unsigned int i=0;i<=mMaxDicVolDepth;++i) mvNbSolutionDepth[i]=0;

   float latstep=0.5,
         vstep=(mVolumeMax-mVolumeMin)/(ceil((mVolumeMax-mVolumeMin)/500)-0.0001);
   mCosAngMax=abs(cos(mAngleMax));
   const float cosangstep=mCosAngMax/(ceil(mCosAngMax/.08)-.0001);
   if(((mVolumeMax-mVolumeMin)/vstep)>10) vstep=(mVolumeMax-mVolumeMin)/9.999;
   if(((mLengthMax-mLengthMin)/latstep)>25) latstep=(mLengthMax-mLengthMin)/24.9999;

   cout<<mLengthMin<<"->"<<mLengthMax<<","<<latstep<<","<<(mLengthMax-mLengthMin)/latstep<<endl;
   cout<<mAngleMin<<"->"<<mAngleMax<<","<<cosangstep<<","<<mCosAngMax<<","<<(mAngleMax-mAngleMin)/cosangstep<<endl;
   cout<<mVolumeMin<<"->"<<mVolumeMax<<","<<vstep<<","<<(mVolumeMax-mVolumeMin)/vstep<<endl;
   RecUnitCell par0,dpar;
   par0.mlattice=mlattice;
   dpar.mlattice=mlattice;
   par0.mCentering=mCentering;
   dpar.mCentering=mCentering;
   //Zero shift parameter - not used for dicvol right now ? :TODO:
   par0.par[0]=0.0;
   dpar.par[0]=0.0;
   unsigned long nbCalc=0;
   Chronometer chrono;
   float bestscore=0;
   list<pair<RecUnitCell,double> >::iterator bestpos;
   bool breakDepth=false;

   // now set it here
   mMaxDicVolDepth = 15;
   mpPeakList->mvvcriterion.assign(mMaxDicVolDepth+1, vector<short>());
   mpPeakList->mvvPK.assign(mMaxDicVolDepth+1, vector<double>());
   mpPeakList->mvvis_good.assign(mMaxDicVolDepth+1, vector<bool>());
   mpPeakList->active_depth = -1;
   mpPeakList->critical_depth = -1;
   mpPeakList->critical_passability = 2;
   for (unsigned int dicvol_count = 0; dicvol_count <= mMaxDicVolDepth + 1; ++dicvol_count)
   {
      mpPeakList->active_depth = min(dicvol_count, mMaxDicVolDepth);
      mpPeakList->mis_exta_run = dicvol_count > mMaxDicVolDepth;

      mpPeakList->mvcount.assign(mMaxDicVolDepth + 1, -1);
      // In the triclinic case, first try assigning a* and b* from the first reflections
      if(false) //mlattice==TRICLINIC)
         for(float minv=mVolumeMin;minv<mVolumeMax;minv+=vstep)
         {
            float maxv=minv+vstep;
            if(maxv>mVolumeMax)maxv=mVolumeMax;
            cout<<"Starting: V="<<minv<<"->"<<maxv<<endl;
            const float minr=1/(mLengthMax*mLengthMax);
            const float maxr=1/(mLengthMin*mLengthMin);
            const float stepr=(maxr-minr)/24.999;
            float p1,p2;
            for(unsigned int i=0;i<=5;i++)
            {
               switch(i)
               {// Try to find a and b from the first observed reflections
                  case 0: p1=mpPeakList->GetPeakList()[0].d2obs  ;p2=mpPeakList->GetPeakList()[1].d2obs  ; break;
                  case 1: p1=mpPeakList->GetPeakList()[0].d2obs  ;p2=mpPeakList->GetPeakList()[2].d2obs  ; break;
                  case 2: p1=mpPeakList->GetPeakList()[1].d2obs/2;p2=mpPeakList->GetPeakList()[0].d2obs  ; break;
                  case 3: p1=mpPeakList->GetPeakList()[1].d2obs/2;p2=mpPeakList->GetPeakList()[2].d2obs  ; break;
                  case 4: p1=mpPeakList->GetPeakList()[2].d2obs/2;p2=mpPeakList->GetPeakList()[0].d2obs  ; break;
                  case 5: p1=mpPeakList->GetPeakList()[2].d2obs/2;p2=mpPeakList->GetPeakList()[1].d2obs  ; break;
               }
               //if(i>0) exit(0);
               if(p1>p2) continue;
               cout<<"Trying #"<<i<<": a*="<<p1<<", b*="<<p2<<endl;
               float min3r=p2,
                     max3r=maxr;//:TODO: use larger value to take angles into account ?
               const float step3r=(max3r-min3r)/(ceil((max3r-min3r)/stepr)-.001);
               vector<unsigned int> vdepth(mnpar-1);
               for(vector<unsigned int>::iterator pos=vdepth.begin();pos!=vdepth.end();) *pos++=0;
               vdepth[0]=3;
               vdepth[1]=3;
               for(float p3=min3r;p3<max3r;p3+=step3r)
               {
                  //cout<<"    p3="<<p3<<endl;
                  const float max4r=p3+step3r;
                  const float step4r=max4r/(ceil(max4r/stepr)-.001);
                  for(float p4=0;p4<max4r;p4+=step4r)
                  {
                     //cout<<"      p4="<<p4<<endl;
                     float max5r=(p2+stepr);
                     const float step5r=max5r/(ceil(max5r/stepr)-.001);
                     for(float p5=0;p5<max5r;p5+=step5r)
                     {
                        float max6r=(p1+stepr);
                        const float step6r=max6r/(ceil(max6r/stepr)-.001);
                        for(float p6=-max6r;p6<max6r;p6+=step6r)
                        {
                           //cout<<"          p6="<<p6<<"/"<<p1<<"/"<<p3<<endl;
                           dpar.par[1]=stepr*pow(float(0.51),int(vdepth[0]));
                           dpar.par[2]=stepr*pow(float(0.51),int(vdepth[1]));
                           dpar.par[3]=step3r*0.51;
                           dpar.par[4]=step4r*0.51;
                           dpar.par[5]=step5r*0.51;
                           dpar.par[6]=step6r*0.51;

                           par0.par[0]=0;
                           par0.par[1]=p1;
                           par0.par[2]=p2;
                           par0.par[3]=p3+step3r/2;
                           par0.par[4]=p4+step4r/2;
                           par0.par[5]=p5+step5r/2;
                           par0.par[6]=p6+step6r/2;
                           //for(int i=0;i<=6;++i)cout<<par0.par[i]<<",";
                           //cout<<endl;
                           //for(int i=0;i<=6;++i)cout<<dpar.par[i]<<",";
                           //cout<<endl;
                           RDicVol(par0,dpar,0,nbCalc,minv,maxv,vdepth);
                        }
                     }
                  }
               }
               cout<<"Finished trying: a*="<<p1<<" A, b*="<<p2<<" A, "<<nbCalc
                  <<" unit cells tested, "<<nbCalc/chrono.seconds()<<" tests/s,   Elapsed time="
                  <<chrono.seconds()<<"s, Best score="<<mBestScore<<", "<<stopOnScore<<", "<<breakDepth<<endl;
               breakDepth=false;
               if(stopOnDepth>0)
                  for(unsigned int i=stopOnDepth; i<mvNbSolutionDepth.size();++i)
                     if(mvNbSolutionDepth[i]>1) {breakDepth=true;break;}
               if((mBestScore>stopOnScore)&&(breakDepth)) break;
            }//cases
            cout<<"Finished triclinic QUICK tests for: V="<<minv<<"->"<<maxv<<" A^3, "<<nbCalc
               <<" unit cells tested, "<<nbCalc/chrono.seconds()<<" tests/s,   Elapsed time="
               <<chrono.seconds()<<"s, Best score="<<mBestScore<<endl;
            if((mBestScore>stopOnScore)&&(breakDepth)) break;
         }//volume
      if((mBestScore<stopOnScore)||(!breakDepth))
      for(float minv=mVolumeMin;minv<mVolumeMax;minv+=vstep)
      {
         float maxv=minv+vstep;
         if(maxv>mVolumeMax)maxv=mVolumeMax;
         cout<<"Starting: V="<<minv<<"->"<<maxv<<endl;
         switch(mlattice)
         {
            case TRICLINIC:
            {
               const unsigned int nbstep=25;
               vector<float> v1=linspace(mLengthMin,mLengthMax,nbstep);
               const float lstep=v1[1]-v1[0];
               for(unsigned int i1=0;i1<(nbstep-1);++i1)
               {
                  const float p1 =(1/(v1[i1]*v1[i1])+1/(v1[i1+1]*v1[i1+1]))/2;
                  const float dp1=(1/(v1[i1]*v1[i1])-1/(v1[i1+1]*v1[i1+1]))/2;
                  //cout<<"p1="<<p1<<endl;
                  //cout<<(v1[i1+1]-mLengthMin)/lstep+2.001<<endl;
                  const unsigned int nb2=int((v1[i1+1]-mLengthMin)/lstep+2.001);
                  vector<float> v2=linspace(mLengthMin,v1[i1+1],nb2);
                  //for(unsigned int i2=0;i2<nb2;++i2) cout<<1/v2[i2]/v2[i2]<<"   ";
                  //cout<<endl;
                  for(unsigned int i2=0;i2<(nb2-1);++i2)
                  {
                     const float p2 =(1/(v2[i2]*v2[i2])+1/(v2[i2+1]*v2[i2+1]))/2;
                     const float dp2=(1/(v2[i2]*v2[i2])-1/(v2[i2+1]*v2[i2+1]))/2;
                     //cout<<"  p2="<<p2<<endl;
                     const unsigned int nb3=int((v2[i2+1]-mLengthMin)/lstep+2.001);
                     vector<float> v3=linspace(mLengthMin,v2[i2+1],nb3);
                     for(unsigned int i3=0;i3<(nb3-1);++i3)
                     {
                        const float p3 =(1/(v3[i3]*v3[i3])+1/(v3[i3+1]*v3[i3+1]))/2;
                        const float dp3=(1/(v3[i3]*v3[i3])-1/(v3[i3+1]*v3[i3+1]))/2;

                        //const float vmax3=v1[i1+1]*v2[i2+1]*v3[i3+1];
                        //const float vmin3=v1[i1]*v2[i2]*v3[i3]/5;
                        const float vmin3=1/sqrt((p1+dp1)*(p2+dp2)*(p3+dp3));
                        const float vmax3=1/sqrt((p1-dp1)*(p2-dp2)*(p3-dp3))*2; // *2 - sufficient margin to take into account angles
                        if(vmax3<minv) continue;
                        if(vmin3>maxv) continue;

                        //cout<<"  p3="<<p3<<endl;
                        //char buf[200];
                        //sprintf(buf,"a1=%6.4f +/-%6.4f (%6.3f-%6.3f) a2=%6.4f +/-%6.4f (%6.3f-%6.3f) a3=%6.4f +/-%6.4f (%6.3f-%6.3f), V=%6.2f-%6.2f",
                        //        p1,dp1,1/sqrt(p1+dp1),1/sqrt(p1-dp1),p2,dp2,1/sqrt(p2+dp2),1/sqrt(p2-dp2),
                        //        p3,dp3,1/sqrt(p3+dp3),1/sqrt(p3-dp3),vmin3,vmax3);
                        //cout<<buf<<endl;
                        const unsigned int nb4=int((p1+dp1)/(4*dp1)+2.001);
                        vector<float> v4=linspace(0,p1+dp1,nb4);
                        for(unsigned int i4=0;i4<(nb4-1);++i4)
                        {
                           const float p4 =(v4[i4+1]+v4[i4])/2;
                           const float dp4=(v4[i4+1]-v4[i4])/2;
                           //cout<<"      p4="<<p4<<endl;
                           const unsigned int nb5=int((p2+dp2)/(4*dp2)+2.001);
                           vector<float> v5=linspace(0,p2+dp2,nb5);
                           for(unsigned int i5=0;i5<(nb5-1);++i5)
                           {
                              const float p5 =(v5[i5+1]+v5[i5])/2;
                              const float dp5=(v5[i5+1]-v5[i5])/2;
                              //cout<<"        p5="<<p5<<endl;

                              float vmax6=(p1+dp1)+(p2+dp2)-(p4-dp4)-(p5-dp5);
                              if(vmax6>(p1+dp1))  vmax6=p1+dp1;
                              if(vmax6<0) continue;
                              const unsigned int nb6=int((2*vmax6)/(4*dp1)+2.001);
                              vector<float> v6=linspace(-vmax6,vmax6,nb6);
                              //cout<<"        p6="<<nb6<<","<<vmax6<<endl;
                              for(unsigned int i6=0;i6<(nb6-1);++i6)
                              {
                                 const float p6 =(v6[i6+1]+v6[i6])/2;
                                 const float dp6=(v6[i6+1]-v6[i6])/2;
                                 //cout<<"          p6="<<p6<<endl;

                                 //char buf[200];
                                 //sprintf(buf,"%6.4f-%6.4f  %6.4f-%6.4f  %6.4f-%6.4f  %6.4f-%6.4f  %6.4f-%6.4f  %6.4f-%6.4f  ",
                                 //       p1-dp1,p1+dp1,p2-dp2,p2+dp2,p3-dp3,p3+dp3,p4-dp4,p4+dp4,p5-dp5,p5+dp5,p6-dp6,p6+dp6);
                                 //cout<<"Testing: "<<buf<<endl;
                                 //cout<<i1<<","<<i2<<","<<i3<<","<<i4<<","<<i5<<","<<i6<<endl;
                                 dpar.par[1]=dp1;
                                 dpar.par[2]=dp2;
                                 dpar.par[3]=dp3;
                                 dpar.par[4]=dp4;
                                 dpar.par[5]=dp5;
                                 dpar.par[6]=dp6;

                                 par0.par[0]=0;
                                 par0.par[1]=p1;
                                 par0.par[2]=p2;
                                 par0.par[3]=p3;
                                 par0.par[4]=p4;
                                 par0.par[5]=p5;
                                 par0.par[6]=p6;

                                 RDicVol(par0,dpar,0,nbCalc,minv,maxv);
                              }
                           }
                        }
                     }
                  }
               }
               break;
            }
            case MONOCLINIC:
            {
               RecUnitCell parlarge,//parm: smallest reciprocal, largest direct cell
                           parsmall;//parp: largest reciprocal, smallest direct cell
               vector<float> parlarged,parsmalld;
               latstep=(mLengthMax-mLengthMin)/24.999;
               for(float x4=0;x4<mCosAngMax+cosangstep;x4+=cosangstep)
               {
                  const float sinbeta=sqrt(abs(1-x4*x4));
                  float x1=mLengthMin;
                  for(;x1<mLengthMax;x1+=latstep)
                  {
                     float x2=mLengthMin;
                     for(;x2<mLengthMax;x2+=latstep)
                     {
                        float x3=x1;
                        const float x3step=(mLengthMax-x1)/(ceil((mLengthMax-x1)/latstep)-0.001);
                        for(;x3<mLengthMax;x3+=x3step) //x3+=(latstep+x3*sin4)
                        {
                           if((x3*x4)>x1) break;// | c * cos(beta) | <a
                           dpar.par[1]=(1/(x1)-1/(x1+latstep))*0.5/sinbeta;
                           dpar.par[2]=(1/(x2)-1/(x2+latstep))*0.5/sinbeta;
                           dpar.par[3]=(1/(x3)-1/(x3+x3step ))*0.5/sinbeta;
                           dpar.par[4]=cosangstep*0.5;

                           par0.par[0]=0;
                           par0.par[1]=(1/(x1)+1/(x1+latstep))*0.5/sinbeta;
                           par0.par[2]=(1/(x2)+1/(x2+latstep))*0.5/sinbeta;
                           par0.par[3]=(1/(x3)+1/(x3+x3step ))*0.5/sinbeta;
                           par0.par[4]=x4+cosangstep*.5;

                           const float smallv=x1*x2*x3*sinbeta;
                           if(smallv>maxv) break;
                           const float largev=(x1+latstep)*(x2+latstep)*(x3+latstep)*(sinbeta+cosangstep);
                           if(largev<minv) continue;
                           /*
                           char buf[200];
                           sprintf(buf,"a=%5.2f-%5.2f b=%5.2f-%5.2f c=%5.2f-%5.2f alpha=%5.2f-%5.2f beta=%5.2f-%5.2f gamma=%5.2f-%5.2f V=%5.2f-%5.2f",
                                    parsmalld[0],parlarged[0],parsmalld[1],parlarged[1],parsmalld[2],parlarged[2],parsmalld[3]*RAD2DEG,parlarged[3]*RAD2DEG,
                                    parsmalld[4]*RAD2DEG,parlarged[4]*RAD2DEG,parsmalld[5]*RAD2DEG,parlarged[5]*RAD2DEG,parsmalld[6],parlarged[6]);
                           cout<<buf<<"   VM="<<maxv<<", x3="<<x3<<endl;
                           */
                           RDicVol(par0,dpar,0,nbCalc,minv,maxv);
                        }//x3
                        //if(((parsmalld[6]>maxv)&&(x3==x1))||(parlarged[1]>mLengthMax)) break;
                     }//x2
                  }//x1
                  // Test if we have one solution before going to the next angle range
                  for(list<pair<RecUnitCell,double> >::iterator pos=mvSolution.begin();pos!=mvSolution.end();++pos)
                  {
                     const double score=pos->second;//Score(*mpPeakList,pos->first,mNbSpurious);
                     if(score>bestscore) {bestscore=score;bestpos=pos;}
                  }
                  bool breakDepth=false;
                  if(stopOnDepth>0)
                     for(unsigned int i=stopOnDepth; i<mvNbSolutionDepth.size();++i)
                        if(mvNbSolutionDepth[i]>1) {breakDepth=true;break;}
                  if((bestscore>stopOnScore)&&(breakDepth)) break;
               }//x4
               break;
            }
            case ORTHOROMBIC:
            {
               if(false)
               {
                  // Test 7.677350  5.803770  10.313160   V=480
                  //const float a=7.75,b=5.75,c=10.25;
                  // Test 6.062000 16.779400 16.8881 v=1750
                  const float a=6.25,b=16.75,c=16.75;
                  dpar.par[1]=(1/(a-.25)-1/(a+.25))*0.5;
                  dpar.par[2]=(1/(b-.25)-1/(b+.25))*0.5;
                  dpar.par[3]=(1/(c-.25)-1/(c+.25))*0.5;
                  par0.par[0]=0;
                  par0.par[1]=1/a;
                  par0.par[2]=1/b;
                  par0.par[3]=1/c;
                  RDicVol(par0,dpar,0,nbCalc,minv,maxv);
                  break;
               }
               latstep=(mLengthMax-mLengthMin)/24.999;
               for(float x1=mLengthMin;x1<mLengthMax;x1+=latstep)
               {
                  for(float x2=x1;x2<mLengthMax;x2+=latstep)
                  {
                     for(float x3=x2;x3<mLengthMax;x3+=latstep)
                     {
                        dpar.par[1]=(1/(x1)-1/(x1+latstep))*0.5;
                        dpar.par[2]=(1/(x2)-1/(x2+latstep))*0.5;
                        dpar.par[3]=(1/(x3)-1/(x3+latstep))*0.5;

                        par0.par[0]=0;
                        par0.par[1]=(1/(x1)+1/(x1+latstep))*0.5;
                        par0.par[2]=(1/(x2)+1/(x2+latstep))*0.5;
                        par0.par[3]=(1/(x3)+1/(x3+latstep))*0.5;

                        const float vmin=x1*x2*x3,vmax=(x1+latstep)*(x2+latstep)*(x3+latstep);
                        if(vmin>maxv) break;
                        if(vmax>=minv) RDicVol(par0,dpar,0,nbCalc,minv,maxv);
                     }
                     if((x1*x2*x2)>maxv) break;
                  }
                  if((x1*x1*x1)>maxv) break;
               }
               break;
            }
            case HEXAGONAL:
            {
               vector<float> parlarged,parsmalld;// Small & large UC in direct space
               latstep=(mLengthMax-mLengthMin)/24.999;
               for(float x1=mLengthMin;;x1+=latstep)
               {
                  for(float x2=mLengthMin;x2<(mLengthMax+latstep);x2+=latstep)
                  {
                     dpar.par[1]=(1/(x1)-1/(x1+latstep))*0.5;
                     dpar.par[2]=(1/(x2)-1/(x2+latstep))*0.5;

                     par0.par[0]=0;
                     par0.par[1]=(1/(x1)+1/(x1+latstep))*0.5;
                     par0.par[2]=(1/(x2)+1/(x2+latstep))*0.5;

                     RecUnitCell parlarge=par0,parsmall=par0;
                     for(unsigned int i=0;i<6;++i) {parlarge.par[i]-=dpar.par[i];parsmall.par[i]+=dpar.par[i];}
                     parlarged=parlarge.DirectUnitCell();
                     parsmalld=parsmall.DirectUnitCell();
                     /*
                     char buf[200];
                     sprintf(buf,"a=%5.2f-%5.2f b=%5.2f-%5.2f c=%5.2f-%5.2f alpha=%5.2f-%5.2f beta=%5.2f-%5.2f gamma=%5.2f-%5.2f V=%5.2f-%5.2f",
                              parsmalld[0],parlarged[0],parsmalld[1],parlarged[1],parsmalld[2],parlarged[2],parsmalld[3]*RAD2DEG,parlarged[3]*RAD2DEG,
                              parsmalld[4]*RAD2DEG,parlarged[4]*RAD2DEG,parsmalld[5]*RAD2DEG,parlarged[5]*RAD2DEG,parsmalld[6],parlarged[6]);
                     */
                     if((parsmalld[6]<maxv)&&(parlarged[6]>minv))
                     {
                        //cout<<buf<<endl;
                        RDicVol(par0,dpar,0,nbCalc,minv,maxv);
                     }
                     //else cout<<buf<<" BREAK"<<endl;
                  }
                  if(parlarged[0]>mLengthMax) break;
               }
               break;
            }
            case RHOMBOHEDRAL:   //:TODO:
            {
               latstep=(mLengthMax-mLengthMin)/24.999;
               for(float x1=mLengthMin;x1<(mLengthMax+latstep);x1+=latstep)
               {
                  for(float x2=0;x2<mCosAngMax+cosangstep;x2+=cosangstep)
                  {
                     dpar.par[1]=latstep/2*1.1;
                     dpar.par[2]=cosangstep/2*1.1;

                     par0.par[0]=0;
                     par0.par[1]=x1-latstep/2*1.1;
                     par0.par[2]=x2-cosangstep/2*1.1;
                     vector<float> par=par0.DirectUnitCell();
                     if((par[6]<maxv)&&(par[6]>minv))
                     {
                        RDicVol(par0,dpar,0,nbCalc,minv,maxv);
                     }
                  }
               }
               break;
            }
            case TETRAGONAL:
            {
               vector<float> parlarged,parsmalld;// Small & large UC in direct space
               latstep=(mLengthMax-mLengthMin)/24.999;
               for(float x1=mLengthMin;x1<mLengthMax;x1+=latstep)
               {
                  for(float x2=mLengthMin;x2<mLengthMax;x2+=latstep)
                  {
                     dpar.par[1]=(1/(x1)-1/(x1+latstep))*0.5;
                     dpar.par[2]=(1/(x2)-1/(x2+latstep))*0.5;

                     par0.par[0]=0;
                     par0.par[1]=(1/(x1)+1/(x1+latstep))*0.5;
                     par0.par[2]=(1/(x2)+1/(x2+latstep))*0.5;

                     RecUnitCell parlarge=par0,parsmall=par0;
                     for(unsigned int i=0;i<6;++i) {parlarge.par[i]-=dpar.par[i];parsmall.par[i]+=dpar.par[i];}
                     parlarged=parlarge.DirectUnitCell();
                     parsmalld=parsmall.DirectUnitCell();
                     /*
                     char buf[200];
                     sprintf(buf,"a=%5.2f-%5.2f b=%5.2f-%5.2f c=%5.2f-%5.2f alpha=%5.2f-%5.2f beta=%5.2f-%5.2f gamma=%5.2f-%5.2f V=%5.2f-%5.2f",
                              parsmalld[0],parlarged[0],parsmalld[1],parlarged[1],parsmalld[2],parlarged[2],parsmalld[3]*RAD2DEG,parlarged[3]*RAD2DEG,
                              parsmalld[4]*RAD2DEG,parlarged[4]*RAD2DEG,parsmalld[5]*RAD2DEG,parlarged[5]*RAD2DEG,parsmalld[6],parlarged[6]);
                     */
                     if((parsmalld[6]<maxv)&&(parlarged[6]>minv))
                     {
                        RDicVol(par0,dpar,0,nbCalc,minv,maxv);
                     }
                     if(parsmalld[6]>maxv) break;
                  }
                  if((x1*mLengthMin*mLengthMin)>maxv) break;
               }
               break;
            }
            case CUBIC:
            {
               latstep=(mLengthMax-mLengthMin)/24.999;
               cout<<mLengthMax<<","<<mLengthMin<<","<<latstep<<endl;
               for(float x1=mLengthMin;x1<(mLengthMax+latstep);x1+=latstep)
               {
                  dpar.par[1]=(1/(x1)-1/(x1+latstep))*0.5;

                  par0.par[0]=0;
                  par0.par[1]=(1/(x1)+1/(x1+latstep))*0.5;

                  const float vmin=x1*x1*x1,vmax=(x1+latstep)*(x1+latstep)*(x1+latstep);
                  if(vmin>maxv)break;
                  if(vmax>minv) RDicVol(par0,dpar,0,nbCalc,minv,maxv);
               }
               break;
            }
         }
         cout<<"Finished: V="<<minv<<"->"<<maxv<<" A^3, "<<nbCalc
            <<" unit cells tested, "<<nbCalc/chrono.seconds()<<" tests/s,   Elapsed time="
            <<chrono.seconds()<<"s"<<endl;
         //info-debug.
         cout << "global_p: " << global_p << endl;
         cout << "global_pmax: " << global_pmax << endl;
         cout << "global_pkmax: " << global_pkmax << endl;
         cout << "control depth global_pk_integralmax: " << global_pk_integralmax << endl;
         cout << "global_kmax: " << global_kmax << endl;
         double excess_ratio = (double)global_check_inside_count/global_is_inside_count;
         cout << "global_check_inside_count/global_is_inside_count=" << global_check_inside_count << "/" << global_is_inside_count << "=" << excess_ratio << endl;
         // cout << "global_min_ratio: " << global_min_ratio << endl;
         // cout << "global_max_ratio: " << global_max_ratio << endl;

         cout << "calls_per level:" << endl;
         for (int level = 0; level < 20; ++level)
         {
            if (calls_per_level[level] == 0)
               break;
            cout << level << ": " << calls_per_level[level] << "; ";
         }
         cout << endl;
         cout << "sum: " << accumulate(calls_per_level.begin(), calls_per_level.end(), 0) << endl;
         cout << "repeated_calls_per level:" << endl;
         for (int level = 0; level < 20; ++level)
         {
            if (calls_per_level[level] == 0)
               break;
            cout << level << ": " << repeated_calls_per_level[level] << "; ";
         }
         cout << endl;
         cout << "sum: " << accumulate(repeated_calls_per_level.begin(), repeated_calls_per_level.end(), 0) << endl;
         cout << "sum_passed_PK_per level/sum_PK_per level:" << endl;
         for (int level = 0; level < 20; ++level)
         {
            double ratio = sum_passed_PK_per_level[level]/sum_PK_per_level[level];
            if (sum_PK_per_level[level] == 0)
               break;
            cout << level << ": " << sum_passed_PK_per_level[level] << "/" << sum_PK_per_level[level] << "=" << ratio << "; ";
         }
         cout << endl;
         for(list<pair<RecUnitCell,double> >::iterator pos=mvSolution.begin();pos!=mvSolution.end();++pos)
         {
            const float score=pos->second;//Score(*mpPeakList,pos->first,mNbSpurious);
            if(score>bestscore) {bestscore=score;bestpos=pos;}
         }
         bool breakDepth=false;
         if(stopOnDepth>0)
            for(unsigned int i=stopOnDepth; i<mvNbSolutionDepth.size();++i)
               if(mvNbSolutionDepth[i]>1) {breakDepth=true;break;}
         if((bestscore>stopOnScore)&&(breakDepth)) break;
      }
      /*
      {// Tag spurious lines
         vector<int> vSpuriousScore;
         for(vector<PeakList::hkl>::const_iterator pos=mpPeakList->GetPeakList().begin();pos!=mpPeakList->GetPeakList().end();++pos)
            vSpuriousScore.push_back(pos->stats);
         sort(vSpuriousScore.begin(),vSpuriousScore.end());
         const int threshold=vSpuriousScore[vSpuriousScore.size()/2]*5;
         for(vector<PeakList::hkl>::iterator pos=mpPeakList->mvHKL.begin();pos!=mpPeakList->mvHKL.end();++pos)
            if(pos->stats > threshold) pos->isSpurious=true;
            else pos->isSpurious=false;
         mpPeakList->Print(cout);
      }
      */

     // analyzing run
      if (mpPeakList->active_depth == mMaxDicVolDepth)
      {
         vector<double> PK_integral_sorted = mpPeakList->mvPK_integral;
         sort(PK_integral_sorted.begin(), PK_integral_sorted.end(), greater<double>());
         int mvSolution_size = 1000 < PK_integral_sorted.size() ? 1000 : PK_integral_sorted.size();
         mpPeakList->mPK_solution_critical = PK_integral_sorted[mvSolution_size-1];
      }
      if (!mpPeakList->mis_exta_run)
      {
         this->AnalyzeLevel();
      }
   }
   this->ReduceSolutions(true);
   for(list<pair<RecUnitCell,double> >::iterator pos=mvSolution.begin();pos!=mvSolution.end();++pos)
   {
      vector<float> par=pos->first.DirectUnitCell();
      cout<<__FILE__<<":"<<__LINE__<<" Solution ? a="<<par[0]<<", b="<<par[1]<<", c="<<par[2]
          <<", alpha="<<par[3]*RAD2DEG<<", beta="<<par[4]*RAD2DEG<<", gamma="<<par[5]*RAD2DEG
          <<", V="<<par[6]<<", PK="<<pos->second <<endl;
   }
      cout<<nbCalc<<"unit cells tested, "<<nbCalc/chrono.seconds()<<" tests/s,   Elapsed time="
          <<chrono.seconds()<<"s"<<endl;

   // switch off
   // bestscore=0;bestpos=mvSolution.end();
   // for(list<pair<RecUnitCell,double> >::iterator pos=mvSolution.begin();pos!=mvSolution.end();++pos)
   // {
   //    const float score=Score(*mpPeakList,pos->first,mNbSpurious);
   //    if(score>bestscore) {bestpos=pos;bestscore=score;}
   //    vector<float> par=pos->first.DirectUnitCell();
   //    cout<<__FILE__<<":"<<__LINE__<<" Solution ? a="<<par[0]<<", b="<<par[1]<<", c="<<par[2]
   //        <<", alpha="<<par[3]*RAD2DEG<<", beta="<<par[4]*RAD2DEG<<", gamma="<<par[5]*RAD2DEG
   //        <<", V="<<par[6]<<", score="<<score<<endl;
   // }
   // if(bestpos!=mvSolution.end())
   // {
   //    vector<float> par=bestpos->first.DirectUnitCell();
   //    //my_debug
   //    // get reciprocal cell
   //    cout << endl << "!!!!!mine"  << endl;
   //    cout << "reciprocal unit cell parameters" << endl;
   //    for (int i = 0; i < 7; ++i)
   //       cout << i << ": " << (bestpos->first.par)[i] << endl;
   //    // bad position
   //    // cout <<"global_p: " << global_p << endl;
   //    // cout <<"global_pkmax: " << global_pkmax << endl;


   //    cout<<__FILE__<<":"<<__LINE__<<" BEST ? a="<<par[0]<<", b="<<par[1]<<", c="<<par[2]
   //          <<", alpha="<<par[3]*RAD2DEG<<", beta="<<par[4]*RAD2DEG<<", gamma="<<par[5]*RAD2DEG
   //          <<", V="<<par[6]<<", score="<<bestscore<<endl;

   //    cout<<nbCalc<<"unit cells tested, "<<nbCalc/chrono.seconds()<<" tests/s,   Elapsed time="
   //        <<chrono.seconds()<<"s"<<endl;
   // }
   // end switch off
}

bool SimilarRUC(const RecUnitCell &c0,const RecUnitCell &c1, const double delta=0.005)
{
   if(c0.mNbSpurious != c1.mNbSpurious) return false;
   vector<float> par0=c0.DirectUnitCell();
   vector<float> par1=c1.DirectUnitCell();
   float diff=0;
   for(unsigned int i=0;i<6;++i) diff += abs(par0[i]-par1[i]);
   return (diff/6)<delta;
}

bool compareRUCScore(std::pair<RecUnitCell,double> &p1, std::pair<RecUnitCell,double> &p2)
{
   return p1.second > p2.second;
}

void CellExplorer::ReduceSolutions(const bool updateReportThreshold)
{
   const bool verbose=false;
   std::list<std::pair<RecUnitCell,double> > vSolution2;
   // TODO: take into account number of spurious lines for cutoff value.
   // keep only solutions above mBestScore/5
   
   // switch off
   // for(list<pair<RecUnitCell,double> >::iterator pos=mvSolution.begin();pos!=mvSolution.end();)
   // {
   //    if(pos->second<(mBestScore/5)) pos=mvSolution.erase(pos);
   //    else ++pos;
   // }
   // if(updateReportThreshold&& ((mBestScore/5)>mMinScoreReport))
   // {
   //    cout<<"CellExplorer::ReduceSolutions(): update threshold for report from "
   //        <<mMinScoreReport<<" to "<<mBestScore/5<<endl;
   //    mMinScoreReport=mBestScore/5;
   // }
   // end switch off


   while(mvSolution.size()>0)
   {
      vSolution2.push_back(mvSolution.front());
      mvSolution.pop_front();
      vector<float> par=vSolution2.back().first.DirectUnitCell();
      if(verbose)
         cout<<__FILE__<<":"<<__LINE__<<" SOLUTION: a="<<par[0]<<", b="<<par[1]<<", c="<<par[2]
               <<", alpha="<<par[3]*RAD2DEG<<", beta="<<par[4]*RAD2DEG<<", gamma="<<par[5]*RAD2DEG
               <<", V="<<par[6]<<", score="<<vSolution2.back().second<<",   SIMILAR TO:"<<endl;
      for(list<pair<RecUnitCell,double> >::iterator pos=mvSolution.begin();pos!=mvSolution.end();)
      {
         if(SimilarRUC(pos->first,vSolution2.back().first))
         {
            par=pos->first.DirectUnitCell();
            if(verbose)
               cout<<__FILE__<<":"<<__LINE__<<"        1: a="<<par[0]<<", b="<<par[1]<<", c="<<par[2]
                     <<", alpha="<<par[3]*RAD2DEG<<", beta="<<par[4]*RAD2DEG<<", gamma="<<par[5]*RAD2DEG
                     <<", V="<<par[6]<<", score="<<pos->second<<"       ("<<mvSolution.size()<<")"<<endl;
            if(vSolution2.back().first.mlattice==pos->first.mlattice)
            {
               if(pos->second>vSolution2.back().second) vSolution2.back()=*pos;
            }
            else if(vSolution2.back().first.mlattice<pos->first.mlattice) vSolution2.back()=*pos;
            pos=mvSolution.erase(pos);
         }
         else
         {
            par=pos->first.DirectUnitCell();
            if(verbose)
               cout<<__FILE__<<":"<<__LINE__<<"        0: a="<<par[0]<<", b="<<par[1]<<", c="<<par[2]
                     <<", alpha="<<par[3]*RAD2DEG<<", beta="<<par[4]*RAD2DEG<<", gamma="<<par[5]*RAD2DEG
                     <<", V="<<par[6]<<", score="<<pos->second<<"       ("<<mvSolution.size()<<")"<<endl;
            ++pos;
         }
      }
   }
   mvSolution=vSolution2;
   mvSolution.sort(compareRUCScore);

   // keep at most 100 solutions, update mDicVolDepthReport and mMinScoreReport if necessary
   int n_solutions = 20;
   if(mvSolution.size()>n_solutions)
   {
      mvSolution.resize(n_solutions);
      // switch off
      // if(updateReportThreshold && (mvSolution.back().second>mMinScoreReport))
      // {
      //    cout<<"CellExplorer::ReduceSolutions(): update threshold for report from "
      //        <<mMinScoreReport<<" to "<<mvSolution.back().second<<endl;
      //    mMinScoreReport=mvSolution.back().second;
      // }
      // end switch off
   }
}

int get_children_count(CrystalSystem cs)
{
   switch(cs)
   {
   case TRICLINIC:     return 64;
   case MONOCLINIC:    return 16;
   case ORTHOROMBIC:   return 8;
   case HEXAGONAL:     return 4;
   case RHOMBOHEDRAL:  return 4;
   case TETRAGONAL:    return 4;
   case CUBIC:         return 2;
   }
}

// double accumulate_PK(vector<short> const & vcriterion, int count)
// {
//    double ans = 0;
//    for (int i = 0; i != count; ++i)
//    {
//       ans += exp(0.002*vcriterion[i]);
//    }
//    return ans;
// }

void optimize_dicvol(vector<vector<bool> > & vvis_good, const int children_count, const int control_depth)
{
   for (int depth = control_depth; depth > 0; --depth)
   {
      int i = 0;
      vector<bool> buffer;
      vector<bool> new_deep_is_good;
      for (int j = 0; j != vvis_good[depth-1].size(); ++j)
      {
         if (vvis_good[depth-1][j])
         {
            buffer.assign(vvis_good[depth].begin() + i, vvis_good[depth].begin() + i + children_count);
            i += children_count;
            if (find(buffer.begin(), buffer.end(), true) == buffer.end())
            {
               vvis_good[depth-1][j] = false;
            }
            else
            {
               new_deep_is_good.insert(new_deep_is_good.end(), buffer.begin(), buffer.end());
            }
         }
      }
      vvis_good[depth].swap(new_deep_is_good);
   }
}

void CellExplorer::AnalyzeLevel()
{
   // make it global
   int control_depth = 6;
   // make it CellExplorer member
   size_t max_level_size = mMaxLevelSize;
   // size_t max_level_size = 200000000;
   const int active_depth = mpPeakList->active_depth;
   size_t max_ancestor_size = max_level_size / get_children_count(mlattice);
   // vector<vector<short> >::const_iterator pvcriterion = mpPeakList->mvvcriterion.begin() + active_depth;
   vector<vector<double> >::const_iterator pvPK = mpPeakList->mvvPK.begin() + active_depth;
   // vector<short> vcriterion_sorted = *pvcriterion;
   vector<double> vPK_sorted = *pvPK;
   // sort(vcriterion_sorted.begin(), vcriterion_sorted.end(), greater<short>());
   sort(vPK_sorted.begin(), vPK_sorted.end(), greater<double>());


   double sum_PK_all = accumulate(vPK_sorted.begin(), vPK_sorted.end(), 0.0);
   double max_passability = accumulate(vPK_sorted.begin(), vPK_sorted.begin() + min(vPK_sorted.size(), max_ancestor_size), 0.0)/sum_PK_all;
   double passability = -1;
   size_t ancestor_size = 0;
   if (active_depth == 0 || max_passability < mpPeakList->critical_passability*1.01)
   {
      passability = max_passability;
      mpPeakList->critical_passability = passability;
      if (active_depth <= control_depth)
      {
         mpPeakList->critical_depth = active_depth;
      }
      ancestor_size = min(vPK_sorted.size(), max_ancestor_size);
   }
   else
   {
      size_t reduced_ancestor_size = (1.0/pow(3, active_depth - mpPeakList->critical_depth) + 0.05)*max_ancestor_size;
      double reduced_passability = accumulate(vPK_sorted.begin(), vPK_sorted.begin() + min(vPK_sorted.size(), reduced_ancestor_size), 0.0)/sum_PK_all;
      if (active_depth <= control_depth && reduced_passability < mpPeakList->critical_passability)
      {
         passability = mpPeakList->critical_passability;
         double sum_PK_needed = passability*sum_PK_all;
         double sum_PK = 0;
         for(int count = 0; count != vPK_sorted.size(); ++count)
         {
            sum_PK += vPK_sorted[count];
            if (sum_PK >= sum_PK_needed)
            {
               ancestor_size = count;
               break;
            }
         }
      }
      else
      {
         passability = reduced_passability;
         ancestor_size = reduced_ancestor_size;
      }
   }
   double min_PK = vPK_sorted[ancestor_size - 1];
   // do not calculate PK = 0
   min_PK = max(min_PK, 1.0e-30);
   mpPeakList->mvvis_good[active_depth].resize(pvPK->size());
   for (int i = 0; i != pvPK->size(); ++i)
   {
      mpPeakList->mvvis_good[active_depth][i] = (*pvPK)[i] > min_PK || ((*pvPK)[i] == min_PK && i < ancestor_size);
      sum_PK_per_level[active_depth] += mpPeakList->mvvPK[active_depth][i];
      if (mpPeakList->mvvis_good[active_depth][i])
      {
         sum_passed_PK_per_level[active_depth] += mpPeakList->mvvPK[active_depth][i];
      }
   }


   if (active_depth <= control_depth)
   {
      optimize_dicvol(mpPeakList->mvvis_good, get_children_count(mlattice), min(control_depth, active_depth));
   }

   // number of ancestors does not exactly match ancestor_size, but probably close to it
   cout << __LINE__ << ": Analyze level start" << endl;
   cout << "active_depth: " << active_depth << "\nancestor_size: " << ancestor_size << endl;
   cout << "PK_min: " << min_PK << endl;
   cout << "Analyze level end" << endl;

}


void CellExplorer::Init()
{
   // Prepare global optimisation
   //for(unsigned int i=0;i<mpPeakList->nb;++i)
   //   cout<<__FILE__<<":"<<__LINE__<<":d*="<<mpPeakList->mvdobs[i]<<", d*^2="<<mpPeakList->mvd2obs[i]<<endl;
   srand(time(NULL));
   vector<pair<RecUnitCell,float> >::iterator pos;
   const float min_latt=1./mLengthMax;
   const float max_latt=1./mLengthMin;
   const float amp_crossp=abs(cos(mAngleMax));
   //mMin[0]=-.002;mAmp[0]=.004;
   mMin[0]=.00;mAmp[0]=.00;
   switch(mlattice)
   {
      case TRICLINIC:
         mMin[1]=min_latt;mAmp[1]=max_latt-min_latt;
         mMin[2]=min_latt;mAmp[2]=max_latt-min_latt;
         mMin[3]=min_latt;mAmp[3]=max_latt-min_latt;
         mMin[4]=0;mAmp[4]=amp_crossp;
         mMin[5]=0;mAmp[5]=amp_crossp;
         mMin[6]=0;mAmp[6]=amp_crossp;
         mnpar=7;
	 break;
      case MONOCLINIC:
         mMin[1]=min_latt;mAmp[1]=max_latt-min_latt;
         mMin[2]=min_latt;mAmp[2]=max_latt-min_latt;
         mMin[3]=min_latt;mAmp[3]=max_latt-min_latt;
         mMin[4]=0;mAmp[4]=amp_crossp;
	 mnpar=5;
         break;
      case ORTHOROMBIC:
         mMin[1]=min_latt;mAmp[1]=max_latt-min_latt;
         mMin[2]=min_latt;mAmp[2]=max_latt-min_latt;
         mMin[3]=min_latt;mAmp[3]=max_latt-min_latt;
	 mnpar=4;
         break;
      case HEXAGONAL:
         mMin[1]=min_latt;mAmp[1]=max_latt-min_latt;
         mMin[2]=min_latt;mAmp[2]=max_latt-min_latt;
         mnpar=3;
         break;
      case RHOMBOHEDRAL:
         mMin[1]=min_latt;mAmp[1]=max_latt-min_latt;
         mMin[2]=-amp_crossp;mAmp[2]=2*amp_crossp;
         mnpar=3;
         break;
      case TETRAGONAL:
         mMin[1]=min_latt;mAmp[1]=max_latt-min_latt;
         mMin[2]=min_latt;mAmp[2]=max_latt-min_latt;
         mnpar=3;
         break;
      case CUBIC:
         mMin[1]=min_latt;mAmp[1]=max_latt-min_latt;
         mnpar=2;
         break;
   }
   //for(unsigned int k=0;k<mnpar;++k) cout<<"par["<<k<<"]: "<<mMin[k]<<"->"<<mMin[k]+mAmp[k]<<endl;

   unsigned int nb1=0,nb2=0;
   switch(mlattice)
   {
      case TRICLINIC:
      {
         nb1=3;
         nb2=3;
         break;
      }
      case MONOCLINIC:
      {
         nb1=3;
         nb2=1;
         break;
      }
      case ORTHOROMBIC:
      {
         nb1=3;
         break;
      }
      case HEXAGONAL:
      {
         nb1=2;
         break;
      }
      case RHOMBOHEDRAL:
      {
         nb1=2;
         break;
      }
      case TETRAGONAL:
      {
         nb1=2;
         break;
      }
      case CUBIC:
      {
         nb1=1;
         break;
      }
   }
   this->ResetParList();
   {
      RefinablePar tmp("Zero",mRecUnitCell.par+0,-0.01,0.01,gpRefParTypeObjCryst,REFPAR_DERIV_STEP_ABSOLUTE,
                       true,false,true,false);
      tmp.SetDerivStep(1e-4);
      this->AddPar(tmp);
   }
   char buf [50];
   string str="Reciprocal unit cell par #";
   for(unsigned int i=0;i<nb1;++i)
   {
      sprintf(buf,"%i",i);
      RefinablePar tmp(str+(string)buf,mRecUnitCell.par+i+1,
                       0.01,1.,gpRefParTypeObjCryst,REFPAR_DERIV_STEP_ABSOLUTE,
                       false,false,true,false);
      tmp.SetDerivStep(1e-4);
      this->AddPar(tmp);
   }
   for(unsigned int i=nb1;i<(nb1+nb2);++i)
   {
      sprintf(buf,"%i",i);
      RefinablePar tmp(str+(string)buf,mRecUnitCell.par+i+1,
                       0.0,0.5,gpRefParTypeObjCryst,REFPAR_DERIV_STEP_ABSOLUTE,
                       false,false,true,false);
      tmp.SetDerivStep(1e-4);
      this->AddPar(tmp);
   }
}

void CellExplorer::LoadPrior()
{// mlattice
   vector<string> par_group_names;
   string lattice_name;
   switch(mlattice)
   {
   case TRICLINIC:
      lattice_name = "TRICLINIC";
      par_group_names.push_back("a");
      par_group_names.push_back("cos_alpha");
      break;
   case MONOCLINIC:
      lattice_name = "MONOCLINIC";
      par_group_names.push_back("a");
      par_group_names.push_back("b");
      par_group_names.push_back("cos_beta");
      break;
   case ORTHOROMBIC:
      lattice_name = "ORTHOROMBIC";
      par_group_names.push_back("a");
      break;
   case HEXAGONAL:
      lattice_name = "HEXAGONAL";
      par_group_names.push_back("a");
      par_group_names.push_back("c");
      break;
   case RHOMBOHEDRAL:
      lattice_name = "RHOMBOHEDRAL";
      par_group_names.push_back("a");
      par_group_names.push_back("c");
      break;
   case TETRAGONAL:
      lattice_name = "TETRAGONAL";
      par_group_names.push_back("a");
      par_group_names.push_back("c");
      break;
   case CUBIC:
      lattice_name = "CUBIC";
      par_group_names.push_back("a");
      break;
   }

   cout << "Importing prior for " << lattice_name << endl; 
   string pdf_dirname = "./FOX_data/Prior/" + mdatabase + "/" + lattice_name + "/";
   for (int i = 0; i != par_group_names.size(); ++i)
   {
      string filename = pdf_dirname + par_group_names[i] + "-PDF.txt";
      ifstream is(filename.c_str());
      if(!is)
      {
         cout<<"Cannot find a priory density file:"<< filename <<endl;
         exit(0);
      }
      // crutch
      bool is_angle = (par_group_names[i].rfind("cos_", 0) == 0);


      string header;
      getline(is, header);
      // string header = "# a, density. dx=0.05, size=2984\n";
      int entry_start = header.find("dx=") + 3;
      int entry_end = header.find(",", entry_start);
      string substring = header.substr(entry_start, entry_end - entry_start).c_str();
      istringstream istr(substring);
      float dx = 0.0;
      istr >> dx;
      // not working with locale
      // float dx = atof(header.substr(entry_start, entry_end - entry_start).c_str());
      float scale = 1/dx;
      float param, pdf;
      std::vector<float> pdfs;
      float xmax = 0.0;
      while(true)
      {// :TODO: use readline to make sure when the end is reached
         is >> param;
         if(is.good()==false) break;
         is >> pdf;
         if(is.good()==false) break;
         pdfs.push_back(pdf);
         xmax = param;
         if(is.good()==false) break;
         if (!is_angle && (param > 3*mLengthMax)) break;
      }
      cout << "Imported " << pdfs.size() << " pdf values for " << par_group_names[i] << endl;
      mPrior[par_group_names[i]] = FunctionTable(pdfs, scale, xmax, true);
      // remember pdfs is empty
      }

   cout << "checking. for a pdf(10)= " << mPrior["a"][10.0] << endl;
   cout << "checking. for cos_beta pdf(0.3)= " << mPrior["cos_beta"][0.3] << endl;
   cout << "scale for a " << mPrior["a"].scale << endl;
}



}//namespace
