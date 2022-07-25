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

   modified by Daniil Demidov
*/

/*
*  source file for Indexing classes & functions
*
*/
#ifndef _OBJCRYST_INDEXING_H_
#define _OBJCRYST_INDEXING_H_

#include <math.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <list>
#include <map>
#include <string>
#include <time.h>
#include "ObjCryst/RefinableObj/RefinableObj.h"
#include "ObjCryst/RefinableObj/LSQNumObj.h"
#include "ObjCryst/ObjCryst/PowderPattern.h"
namespace ObjCryst
{
/** Different lattice types.
*
*/
enum CrystalSystem
{ TRICLINIC, MONOCLINIC, ORTHOROMBIC, HEXAGONAL, RHOMBOHEDRAL, TETRAGONAL, CUBIC};

enum CrystalCentering
{ LATTICE_P,LATTICE_I,LATTICE_A,LATTICE_B,LATTICE_C,LATTICE_F};

/** Estimate volume from number of peaks at a given dmin
* See J. Appl. Cryst. 20 (1987), 161
*
* \param dmin,dmax: 1/d limits between which the number of reflections has been observed
* \param nbrefl: number of observed reflections
* \param kappa: estimated percentage of reflections actually observed (default=1.0)
*/
float EstimateCellVolume(const float dmin, const float dmax, const float nbrefl,
                         const CrystalSystem system,const CrystalCentering centering,const float kappa=1);

/** Lightweight class describing the reciprocal unit cell, for the fast computation of d*_hkl^2.
*
*
*/
class RecUnitCell
{
   public:
      RecUnitCell(const float zero=0,const float par0=0,const float par1=0,const float par2=0,
                  const float par3=0,const float par4=0,const float par5=0,CrystalSystem lattice=CUBIC,
                  const CrystalCentering cent=LATTICE_P, const unsigned int nbspurious=0);
      RecUnitCell(const RecUnitCell &old);
      void operator=(const RecUnitCell &rhs);
      // access to ith parameter
      //float operator[](const unsigned int i);
      /// Compute d*^2 for hkl reflection
      /// if deriv != -1, compute derivate versus the corresponding parameter
      ///
      /// If derivhkl=1,2,3, compute derivative versus h,k or l.
      float hkl2d(const float h,const float k,const float l,REAL *derivpar=NULL,const unsigned int derivhkl=0) const;
      /// Compute d*^2 for one hkl reflection: this functions computes a d*^2 range (min,max)
      /// for a given range of unit cell parameter (given in the delta parameter) around the
      /// current parameters.
      ///
      /// Used for DicVol algorithm
      void hkl2d_delta(const float h,const float k,const float l,const RecUnitCell &delta, float & dmin, float &dmax) const;
      /** Compute real space unit cell from reciprocal one
      *
      *\param equiv: if true, return real unit cell \e equivalent to the one computed from the reciprocal one,
      * so that alpha, beta and gamma are larger or equal to pi/2, and minimum. This is done
      * by adding multiples of \b a to \b b and multiples of \b a and \b b to \b c.
      */
      std::vector<float> DirectUnitCell(const bool equiv=false)const;
      /** The 6 parameters defining 1/d_hkl^2 = d*_hkl^2, for different crystal classes, from:
      * comments are in error. consult hkl2d - mine
      * d*_hkl^2 = zero + a*^2 h^2 + b*^2 k^2 + c*^2 l^2 + 2 a*.b* hk + 2 b*.c* kl + 2 a*.c* hl
      *
      * for triclinic:
      *   d*_hkl^2 = par[0] + par[1] h^2 + par[1] k^2 + par[2]^2 l^2 + par[3] hk + par[4] kl + par[5] hl
      * for monoclinic:
      *   d*_hkl^2 = zero + par[0]^2 h^2 + par[1]^2 k^2 + par[2]^2 l^2 + par[0]*par[2]*par[3] hl
      * for orthorombic:
      *   d*_hkl^2 = zero + par[0]^2 h^2 + par[1]^2 k^2 + par[2]^2 l^2
      * for hexagonal:
      *   d*_hkl^2 = zero + par[0]^2 h^2 + par[0]^2 k^2 + par[2]^2 l^2 + sqrt(3)/2*par[0]^2 hk
      * for rhombohedral:
      *   d*_hkl^2 = zero + par[0]^2 h^2 + par[1]^2 k^2 + par[2]^2 l^2 + par[3] (hk + kl + hl)
      * for quadratic:
      *   d*_hkl^2 = zero + par[0]^2 h^2 + par[0]^2 k^2 + par[1]^2 l^2
      * for cubic
      *   d*_hkl^2 = zero + par[0]^2 (h^2 + k^2 + l^2)
      */
      REAL par[7];
      CrystalSystem mlattice;
      CrystalCentering mCentering;
      /// The number of spurious lines used to match this RecUnitCell
      unsigned int mNbSpurious;
      // a priory probability
      double P;
};

/** Class for getting a value of a function given float argument
*/
struct FunctionTable
{
    FunctionTable();
    // construct from data and scale. if rob is set, swaps instead of copying
    FunctionTable(std::vector<float> & data, float scale, float xmax, bool rob=false);
    // return data[scale*x]
    float operator[](float x) const;
    std::vector<float> data;
    float scale;
    float xmax;
};


/** Class to store positions of observed reflections.
*
*
*/
class PeakList
{
   public:
      PeakList();
      PeakList(const PeakList &old);
      void operator=(const PeakList &rhs);
      ~PeakList();
      void ImportDhklDSigmaIntensity(std::istream &is,float defaultsigma=.001);
      void ImportDhklIntensity(std::istream &is);
      void ImportDhkl(std::istream &is);
      // pure mine so no collisions. Currently only distance (sigma) of (fi/rho < phi) interval (K_interval)
      // sigma for 1/d is calculated from d K_interval sigma so a little bit different from 1/d K_interval
      void ImportDhklDSigma(std::istream &is);
      void Import2ThetaIntensity(std::istream &is, const float wavelength=1.5418);
      void ImportIntegrals(std::istream &is);
      /** Generate a list of simulated peak positions, from given lattice parameters
      *
      * \param zero: the zero, given for d*^2 (in Angstroems^-2)
      * \param a,b,c,alpha,beta,gamma: lattice parameters
      * \param deg: if true, angles are in degrees instead of radians
      * \param nb: the number of peak positions to generate
      * \param nbspurious: number of spurious lines to be included in the list
      * \param sigma: the maximum relative error for d* - the d* values will be within [d_calc*(1-sigma) ;d_calc*(1+sigma)]
      * \param percentMissing: percentage (between 0 and 1) of missing reflections - maximum allowed 0.9
      * \param verbose: print some info
      * \return: the volume of the simulated unit cell
      */
      float Simulate(float zero, float a, float b, float c,
                    float alpha, float beta, float gamma,
                    bool deg, unsigned int nb=20, unsigned int nbspurious=0,
                    float sigma=0, float percentMissing=0, const bool verbose=false);
      void ExportDhklDSigmaIntensity(std::ostream &out)const;
      /// Add one peak
      ///\param d: 1/d for this peak (Angstroem)
      void AddPeak(const float d, const float iobs=1.0,const float dobssigma=0.0,const float iobssigma=0.0,
                   const int h=0,const int k=0, const int l=0,const float d2calc=0);
      void RemovePeak(unsigned int i);
      void Print(std::ostream &os) const;
      // d2 (1/d**2, Rho, F)
      struct integral
      {
         integral(const float d2=0, const float Rho=0, const float F=0, const float Rho_ext=0);
         float d2, Rho, F, Rho_ext;
      };
      // reset local variables
      void clear_locals() const;
      //sorted linear array of integrals
      vector<integral> integrals;
      // integrals[coef*d2].d2 = d2
      float coef;
      /// maximum d2 for which integral is defined
      float d2max;
      /// One set of Miller indices, a possible indexation for a reflection
      struct hkl0
      {
         hkl0(const int h=0,const int k=0, const int l=0);
         /// Miller indices
         int h,k,l;
      };
      struct thkl
      {
         thkl(int h=0, int k=0, int l=0, int count=0, float d2min=1, float d2max=1, double K=1, double dRho=-1, double P=-1, bool is_disabled=false);

         int h, k, l;
         // number of intersecting experimental peaks
         mutable int count;
         mutable double K, dRho, P;
         mutable float d2min, d2max;
         mutable bool is_disabled;
      };

      // Full list of calculated HKL positions
      mutable vector<thkl> vthkl;

      // crystallographic axis
      struct Ax
      {
         // for a dummy axis P = 1
         explicit Ax(vector<thkl *> const & vthkl, double P=1, double N_with1=0,  double N_without1=0);

         // pointer to thkl from vthkl belonging to this axis
         vector<thkl *> vthkl;
         // maximal PKa for this axis on the path from this to root (depth=0). For local use only.
         mutable double P;
         mutable double N_with1;
         mutable double N_without1;

      };
      
      // vector of axes. 
      // mutable vector<vector<thkl *> > vax;
      mutable vector<Ax> vax;


      /** One observed diffraction line, to be indexed
      *
      */
      struct hkl
      {
         hkl(const float dobs=1.0,const float iobs=0.0,const float dobssigma=0.0,const float iobssigma=0.0,
             const int h=0,const int k=0, const int l=0,const float d2calc=0);
         /// Observed peak position 1/d
         float dobs;
         /// Error on peak position
         float dobssigma;
         /// Observed peak position 1/d^2
         float d2obs;
         /// Min value for observed peak position 1/(d+disgma/2)^2
         float d2obsmin;
         /// Min value for observed peak position 1/(d-disgma/2)^2
         float d2obsmax;
         /// Observed peak intensity
         float iobs;
         /// Error on observed peak intensity
         float iobssigma;
         /// Miller indices, after line is indexed
         mutable int h,k,l;
         // Possible thkl pointers
         // mutable vector<thkl *> vDicVolTHKL;
         // /// Possible Miller indices, stored during a dichotomy search.
         // mutable std::list<hkl0> vDicVolHKL;

         /// Is this line indexed ?
         mutable bool isIndexed;
         /// Is this an impurity line ?
         mutable bool isSpurious;
         /// Indexing statistics
         mutable unsigned long stats;
         /// Calculated position, 1/d^2
         mutable float d2calc;
         /// 1/d^2 difference, obs-calc
         mutable float d2diff;
      };
      /// Get peak list
      vector<hkl> & GetPeakList();
      /// Get peak list
      const vector<hkl> & GetPeakList()const ;
      /// Predict peak positions
      /// Best h,k,l for each observed peak (for least-squares refinement)
      /// This is stored by the Score function, optionnally.
      mutable vector<hkl>mvHKL;
      /// Full list of calculated HKL positions for a given solution, up to a given resolution
      /// After finding a candidate solution, use score with pPredictedHKL=&mvPredictedHKL
      mutable list<hkl> mvPredictedHKL;

      // current number of cells tested (DichoIndexed called), per depth.
      mutable vector<int> mvcount;
      // 500*ln(PK), per count per depth
      mutable vector<vector<short> > mvvcriterion;
      // parallel version
      mutable vector<vector<double> > mvvPK;
      // is the cell good, per count per depth
      mutable vector<vector<bool> > mvvis_good;
      // vector of PK_integral for control depth
      mutable vector<double> mvPK_integral;
      // is it the extra run for determination of candidates
      mutable bool mis_exta_run;
      // critical PK used to choose relatively large (1000) list of possible candidates
      mutable double mPK_solution_critical;
      mutable int active_depth;
      mutable int critical_depth;
      mutable double critical_passability;
};

// Information available from parent
struct ParentInfo
{
   // ndhkl - number of dhkls
   explicit ParentInfo(int ndhkl = 0, int nax = 0);

   // Infrormation about dhkl
   struct DhklInfo
   {
      // pointers in vthkl for thkl possibly matching this dhkl. vDicVolTHKL analogue
      vector<PeakList::thkl * > vpthkl;
   };

   vector<DhklInfo> vdhkl_info;

   struct AxInfo
   {
      explicit AxInfo(int n=0, double P=1, double N_with1=0, double N_without1=0);
      vector<double> vK;
      double P;
      double N_with1, N_without1;
   };

   vector<AxInfo> vax_info;

   // Pointer on the best (maximum K) cell from the control (depth >= contol_depth) block. With K
   // initialized with zero
   mutable std::pair<RecUnitCell, double> * pbest_cell;

   // be careful with interpretation for zero depth. However only the value on control depth is significant
   mutable double PK_integral;
   // usual one-cell PK
   double PK;
};


/// Compute score for a candidate RecUnitCell and a PeakList
float Score(const PeakList &dhkl, const RecUnitCell &ruc, const unsigned int nbSpurious=0,
            const bool verbose=false,const bool storehkl=false,
            const bool storePredictedHKL=false);

/** Algorithm class to find the correct indexing from observed peak positions.
*
*/
class CellExplorer:public RefinableObj
{
   public:
      CellExplorer(const PeakList &dhkl, const CrystalSystem lattice, const unsigned int nbSpurious);
      void Evolution(unsigned int ng,const bool randomize=true,const float f=0.7,const float cr=0.5,unsigned int np=100);
      void SetLengthMinMax(const float min,const float max);
      void SetAngleMinMax(const float min,const float max);
      void SetVolumeMinMax(const float min,const float max);
      void SetNbSpurious(const unsigned int nb);
      void SetMaxLevelSize(const unsigned int nb);
      /// Allowed error on 1/d (squared!), used for dicvol
      void SetD2Error(const float err);
      void SetMinMaxZeroShift(const float min,const float max);
      void SetCrystalSystem(const CrystalSystem system);
      void SetCrystalCentering(const CrystalCentering cent);
      virtual const string& GetClassName() const;
      virtual const string& GetName() const;
      virtual void Print() const;
      virtual unsigned int GetNbLSQFunction() const;
      virtual const CrystVector_REAL& GetLSQCalc(const unsigned int) const;
      virtual const CrystVector_REAL & GetLSQObs(const unsigned int) const;
      virtual const CrystVector_REAL & GetLSQWeight(const unsigned int) const;
      virtual const CrystVector_REAL & GetLSQDeriv(const unsigned int, RefinablePar &);
      virtual void BeginOptimization(const bool allowApproximations=false, const bool enableRestraints=false);
      void LSQRefine(int nbCycle=1, bool useLevenbergMarquardt=true, const bool silent=false);
      /// Run DicVOl algorithm, store only solutions with score >minScore or depth>=minDepth,
      /// stop at the end of one volume interval (~400 A^3) if best score>stopOnScore,
      /// or if one solution was found at depth>=stopOnDepth
      ///
      /// If stopOnDepth==0, do not stop for any depth
      void DicVol(const float minScore=10,const unsigned int minDepth=3,const float stopOnScore=50.0,const unsigned int stopOnDepth=6);
      /** Sort all solutions by score, remove duplicates
      *
      * \param updateReportThreshold: if true, when too many solutions are produced,
      * the threshold above which solutions are reported will be updated to get less solutions.
      */

      void AnalyzeLevel();
      void ReduceSolutions(const bool updateReportThreshold=false);
      double GetBestScore()const;
      const std::list<std::pair<RecUnitCell,double> >& GetSolutions()const;
      std::list<std::pair<RecUnitCell,double> >& GetSolutions();
   private:
      unsigned int RDicVol(RecUnitCell uc0, RecUnitCell uc1, unsigned int depth,unsigned long &nbCalc,const float minV,const float maxV,
                           vector<unsigned int> vdepth=vector<unsigned int>(), ParentInfo  const * pparent_info = 0);
      void Init();
      void LoadPrior();
      /* collection of tables for a priory densities for different parameters.
      possible keys: "a", "b", "c", "cos_alpha", "cos_beta" and (never used) "cos_gamma"
      */
      std::map<std::string, FunctionTable> mPrior;
      /// Max number of obs reflections to use
      std::list<std::pair<RecUnitCell, double> > mvSolution;
      unsigned int mnpar;
      const PeakList *mpPeakList;
      float mLengthMin,mLengthMax;
      float mAngleMin,mAngleMax;
      float mVolumeMin,mVolumeMax;
      float mZeroShiftMin,mZeroShiftMax;
      /// Min values for all parameters (7=unit cell +zero)
      float mMin[7];
      /// Max amplitude (max=min+amplitude) for all parameters
      /// All parameters are treated as periodic for DE (??)
      float mAmp[7];
      // database used for statistics on cells
      // options: "ICSD", hope "CSD" and "PDB" will be added
      std::string mdatabase;
      /// Lattice type for which we search
      CrystalSystem mlattice;
      /// Centering type
      CrystalCentering mCentering;
      unsigned int mNbSpurious;
      int mMaxLevelSize;
      float mD2Error;
      LSQNumObj mLSQObj;
      mutable CrystVector_REAL mObs;
      mutable CrystVector_REAL mCalc;
      mutable CrystVector_REAL mWeight;
      mutable CrystVector_REAL mDeriv;
      /// Reciprocal unit cell used for least squares refinement
      RecUnitCell mRecUnitCell;
      /// Current best score
      float mBestScore;
      /// Number of solutions found during dicvol search, at each depth.
      std::vector<unsigned int> mvNbSolutionDepth;
      float mMinScoreReport;
      unsigned int mMaxDicVolDepth,mDicVolDepthReport;
      /// Stored value of cos(max ang) for tricilinic search - we do
      /// not want to recompute the cos at every dicvol iteration
      mutable float mCosAngMax;
      /// Number of exceptions caught during LSQ, in a given search - above 20 LSQ is disabled
      unsigned int mNbLSQExcept;
};


}//namespace
#endif
