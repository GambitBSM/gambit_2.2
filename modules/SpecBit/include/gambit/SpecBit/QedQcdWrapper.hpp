//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  This class is used to wrap the QedQcd object used by SoftSUSY
///  and FlexibleSUSY in a Gambit Spectrum object. This is to enable
///  access to the parameters of the SM defined as a low-energy effective theory 
///  (as opposed to correspending information defined in a UV model). 
///  Parameters defined this way are often used as input to a physics calculator.
///
///  *********************************************
///
///  Authors: 
///  <!-- add name and date if you modify -->
///   
///  \author Ben Farmer
///          (benjamin.farmer@fysik.su.se)
///  \date 2015 Mar 
///
///  *********************************************

#ifndef __QedQcd_wrap_hpp__
#define __QedQcd_wrap_hpp__

#include "lowe.h" ///TODO: wrap using BOSS at some point, i.e. get this from FlexibleSUSY or SoftSUSY
#include "gambit/Utils/Spectrum.hpp"

namespace Gambit {

   typedef softsusy::QedQcd QedQcdModel;

   class QedQcdWrapper;

   //this contains scale and scheme dependent stuff
   class QedQcd_MSbarPars : public RunparDer<QedQcdModel>
   {
      REDO_TYPEDEFS(QedQcdModel)
      private:
         QedQcdWrapper& my_parent;

         // Limits for running
         double softup;
         double hardup; // Be careful of order in constructor!
         virtual double hard_upper() const {return hardup;}
         virtual double soft_upper() const {return softup;}
         virtual double soft_lower() const {return 0.;}
         virtual double hard_lower() const {return 0.;}
         virtual void RunToScaleOverride(double);

         /* MAPS */

         //mass
         static fmap_plain mass_map_extra;
         static fmap_plain fill_mass_map_extra();
         virtual const fmap_plain& get_mass_map_extra() const;  // MUST have this name for overload 
 
         //mass0
         static fmap_plain mass0_map_extra;
         static fmap_plain fill_mass0_map_extra();
         virtual const fmap_plain& get_mass0_map_extra() const;  // MUST have this name for overload 

   public:
         QedQcd_MSbarPars(QedQcdWrapper&, QedQcdModel&); 
         virtual ~QedQcd_MSbarPars();
         virtual double GetScale() const;
         virtual void SetScale(double scale);
           
   };
 
   /// Class for accessing physical spectrum
   class QedQcd_Phys : public PhysDer<QedQcdModel> 
   {
      REDO_TYPEDEFS(QedQcdModel)
      private:
         QedQcdWrapper &my_parent;

         /* MAPS */
         static fmap PoleMass_map;
         static fmap fill_PoleMass_map();
         virtual const fmap& get_PoleMass_map() const;  // MUST have this name for overload 

      public:
         QedQcd_Phys(QedQcdWrapper&, QedQcdModel&); 
         virtual ~QedQcd_Phys();
   };

    
   class QedQcdWrapper : public Spectrum 
   {
      friend class QedQcd_MSbarPars;
      friend class QedQcd_Phys;
      private:
         //int index_offset;
         //virtual int get_index_offset() const {return index_offset;}
         QedQcdModel qedqcd; // We will make an internal copy of this object

      public:
         /// Override of clone function
         DEFINE_CLONE(QedQcdWrapper)

         /// Internal instances of the derived "inner" classes
         QedQcd_Phys qedqcd_ph;
         QedQcd_MSbarPars qedqcd_msbar_pars;

         // Constructors/destructors
         QedQcdWrapper(bool switch_index_convention=false);
         QedQcdWrapper(const QedQcdModel&, bool switch_index_convention=false);
         QedQcdWrapper(const QedQcdWrapper&);
         virtual ~QedQcdWrapper();
      
         virtual int get_index_offset() const;   
         virtual int get_numbers_stable_particles() const; 
   };
 


} // end Gambit namespace

#endif
