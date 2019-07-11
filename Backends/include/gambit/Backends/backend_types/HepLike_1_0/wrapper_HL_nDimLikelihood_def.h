#ifndef __wrapper_HL_nDimLikelihood_def_HepLike_1_0_h__
#define __wrapper_HL_nDimLikelihood_def_HepLike_1_0_h__

#include <string>
#include <vector>

#include "identification.hpp"

namespace CAT_3(BACKENDNAME,_,SAFE_VERSION)
{
   
   
   // Member functions: 
   inline void HL_nDimLikelihood::Read()
   {
      get_BEptr()->Read();
   }
   
   inline double HL_nDimLikelihood::GetChi2(::std::vector<double, std::allocator<double> > theory)
   {
      return get_BEptr()->GetChi2(theory);
   }
   
   inline double HL_nDimLikelihood::GetLikelihood(::std::vector<double, std::allocator<double> > theory)
   {
      return get_BEptr()->GetLikelihood(theory);
   }
   
   inline double HL_nDimLikelihood::GetLogLikelihood(::std::vector<double, std::allocator<double> > theory)
   {
      return get_BEptr()->GetLogLikelihood(theory);
   }
   
   inline void HL_nDimLikelihood::Profile()
   {
      get_BEptr()->Profile();
   }
   
   inline double HL_nDimLikelihood::GetChi2_profile(double theory, ::std::basic_string<char, std::char_traits<char>, std::allocator<char> > arg_1)
   {
      return get_BEptr()->GetChi2_profile(theory, arg_1);
   }
   
   inline double HL_nDimLikelihood::GetLikelihood_profile(double theory, ::std::basic_string<char, std::char_traits<char>, std::allocator<char> > axis)
   {
      return get_BEptr()->GetLikelihood_profile(theory, axis);
   }
   
   inline double HL_nDimLikelihood::GetLogLikelihood_profile(double theory, ::std::basic_string<char, std::char_traits<char>, std::allocator<char> > X)
   {
      return get_BEptr()->GetLogLikelihood_profile(theory, X);
   }
   
   
   // Wrappers for original constructors: 
   inline HL_nDimLikelihood::HL_nDimLikelihood() :
      HL_Data(__factory0())
   {
      get_BEptr()->set_wptr(this);
      get_BEptr()->set_delete_wrapper(false);
   }
   
   inline HL_nDimLikelihood::HL_nDimLikelihood(::std::basic_string<char, std::char_traits<char>, std::allocator<char> > s) :
      HL_Data(__factory1(s))
   {
      get_BEptr()->set_wptr(this);
      get_BEptr()->set_delete_wrapper(false);
   }
   
   // Special pointer-based constructor: 
   inline HL_nDimLikelihood::HL_nDimLikelihood(Abstract_HL_nDimLikelihood* in) :
      HL_Data(in)
   {
      get_BEptr()->set_wptr(this);
      get_BEptr()->set_delete_wrapper(false);
   }
   
   // Copy constructor: 
   inline HL_nDimLikelihood::HL_nDimLikelihood(const HL_nDimLikelihood& in) :
      HL_Data(in.get_BEptr()->pointer_copy__BOSS())
   {
      get_BEptr()->set_wptr(this);
      get_BEptr()->set_delete_wrapper(false);
   }
   
   // Assignment operator: 
   inline HL_nDimLikelihood& HL_nDimLikelihood::operator=(const HL_nDimLikelihood& in)
   {
      if (this != &in)
      {
         get_BEptr()->pointer_assign__BOSS(in.get_BEptr());
      }
      return *this;
   }
   
   
   // Destructor: 
   inline HL_nDimLikelihood::~HL_nDimLikelihood()
   {
      if (get_BEptr() != 0)
      {
         get_BEptr()->set_delete_wrapper(false);
         if (can_delete_BEptr())
         {
            delete BEptr;
            BEptr = 0;
         }
      }
      set_delete_BEptr(false);
   }
   
   // Returns correctly casted pointer to Abstract class: 
   inline Abstract_HL_nDimLikelihood* HL_nDimLikelihood::get_BEptr() const
   {
      return dynamic_cast<Abstract_HL_nDimLikelihood*>(BEptr);
   }
   
}


#include "gambit/Backends/backend_undefs.hpp"

#endif /* __wrapper_HL_nDimLikelihood_def_HepLike_1_0_h__ */
