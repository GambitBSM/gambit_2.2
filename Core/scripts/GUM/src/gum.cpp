//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Class mathematica_variable, needed to overload
///  constructor and assignment operators to send
///  messages throught WSTP
///
///  ***********************************************
///
///  Authors
///  =======
///
///  (add name and date if you modify)
///
///  \author Tomas Gonzalo
///          (t.e.gonzalo@fys.uio.no)
///  \date 2017 May
///
///  ***********************************************

#include <iostream>
#include <sstream>
#include "../../../../cmake/include/gambit/cmake/cmake_variables.hpp"

#ifdef HAVE_MATHEMATICA
#include MATHEMATICA_WSTP_H

extern WSENV WSenv;

using namespace std;

// Setup the WS environment and start the Kernel
int MStart(void *&pHandle)
{
  WSenv = WSInitialize(0);
  if(WSenv == (WSENV)0)
    return 0;

  int WSerrno;
  std::stringstream WSTPflags;
  #ifdef __APPLE__
  WSTPflags << "-linkname " << MATHEMATICA_KERNEL << " -mathlink";
  #else
  WSTPflags << "-linkname math -mathlink";
  #endif

  pHandle = WSOpenString(WSenv, WSTPflags.str().c_str(), &WSerrno);
  if((WSLINK)pHandle == (WSLINK)0 || WSerrno != WSEOK)
    return 0;

  return 1;
}

// Send a function through WSTP
int MPutFunction(void *pHandle, std::string name, int args)
{
  return WSPutFunction((WSLINK)pHandle, name.c_str(), args);
}

// Overloaded functions to get data through WSTP
int MGetInteger(void *pHandle, int* val) { return WSGetInteger(WSLINK(pHandle), val); }
int MGetFloat(void *pHandle, float* val) { return WSGetReal32(WSLINK(pHandle), val); }
int MGetDouble(void *pHandle, double* val) { return WSGetReal64(WSLINK(pHandle), val); } 
int MGetBool(void *pHandle, bool* val) 
{ 
  const char *val2;
  int ret = WSGetString(WSLINK(pHandle), &val2); 
  *val = (std::string(val2) == "True");
  return ret;
}
int MGetChar(void *pHandle, char* val)
{ 
  const char *val2;
  int ret = WSGetString(WSLINK(pHandle), &val2);
  *val = val2[0];
  return ret;
}
int MGetString(void *pHandle, std::string* val) 
{ 
  const char *val2;
  int ret = WSGetString(WSLINK(pHandle), &val2);
  *val = std::string(val2);
  return ret;
} 
// FIXME: This won't work in python because of the template
//template <typename T> inline int WSGetVariable(WSLINK WSlink, std::vector<T>* val)
//{
//  long int dim;
//  if(!WSCheckFunction(WSlink, "List", &dim))
//    return 0;
//  for(int i=0; i<dim; i++) 
//  {
//    T value;
//    if(!WSGetVariable(WSlink, &value))
//      return 0;
//    val->push_back(value);
//  }
//  return 1;
//}  

// Overloaded functions to put data through WSTP
int MPutInteger(void *pHandle, int val) { return WSPutInteger32(WSLINK(pHandle), val); }
int MPutFloat(void *pHandle, float val) { return WSPutReal32(WSLINK(pHandle), val); }
int MPutDouble(void *pHandle, double val) { return WSPutReal64(WSLINK(pHandle), val); }
int MPutBool(void *pHandle, bool val) 
{ 
  if(val)
    return WSPutSymbol(WSLINK(pHandle), "True");
  else
    return WSPutSymbol(WSLINK(pHandle), "False");
}
int MPutChar(void *pHandle, char val)
{
  return WSPutString(WSLINK(pHandle), std::string(&val).c_str());
}
int MPutString(void *pHandle, std::string val) { return WSPutString(WSLINK(pHandle), val.c_str()); }
//FIXME: This won't work because in python because of the template
//template <typename T> inline int WSPutVariable(WSLINK WSlink, std::vector<T> val)
//{
//  if(!WSPutFunction(WSlink, "List", val.size()))
//    return 0; 
//  for(auto it = val.begin(); it != val.end(); it++)
//    if(!WSPutVariable(WSlink, *it))
//      return 0;
//  return 1;
//}

// Function to close the link
int MEnd(void *pHandle)
{

  WSPutFunction((WSLINK)pHandle, "Exit", 0);
  WSClose((WSLINK)pHandle);
  WSDeinitialize(WSenv);

  return 1;
}

// Rubbish main function
int main()
{
  return 1;
}

#endif /* HAVE_MATHEMATICA */
