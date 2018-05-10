#ifndef FEYNRULES_H
#define FEYNRULES_H

#include "wstp.h"
#include <math.h>
#include <stdio.h>
#include <sstream>
#include <iostream>
#include <vector>

#include "cmake_variables.hpp"
#include "options.hpp"

class FeynRules 
{
  private:
    
    void *pHandle;
    WSLINK link;
    std::string name;
    
  public:
  
    void set_name(std::string modelname)
    {
      name = modelname;
    }

    void create_wstp_link()
    {
      int WSerrno;
      WSENV WSenv = WSInitialize(0);
      
      if(WSenv == (WSENV)0)
      {
        std::cout << "Unable to initialize WSTP environment" << std::endl;
      }
      else
      {
        std::cout << "The environment is initialized successfully... ";
      }
      
      std::stringstream WSTPflags;
      
      // This opens a WSTP connection
      #ifdef __APPLE__
        WSTPflags << "-linkname " << MATHEMATICA_KERNEL << " -mathlink";
      #else
        WSTPflags << "-linkname math -mathlink";
      #endif
      
      pHandle = WSOpenString(WSenv, WSTPflags.str().c_str(), &WSerrno);
      
      link = (WSLINK)pHandle;
      
      if(link == (WSLINK)0 || WSerrno != WSEOK)
      {
        std::cout << "Unable to create link to the Kernel" << std::endl;
        WSNewPacket(link);
      }
      else
      {
        std::cout << "WSTP link started" << std::endl;
      }
      
    }
    
    void close_wstp_link()
    {
        std::string command = "Quit[];";
        send_to_math(command);
        WSClose(link);
        std::cout << "WSTP link closed successfully." << std::endl;
    }
    
    // Wait to receive a packet from the kernel
    void wait_for_packet()
    {
      int pkt;
      while( (pkt = WSNextPacket(link), pkt) && pkt != RETURNPKT)
      {
        WSNewPacket(link);
        if (WSError(link))
        {
          std::cout << "Error reading packet from WSTP" << std::endl;
        }
      }
    }
    
    // Send a string to be evaluated in Mathematica via WSTP
    void send_to_math(std::string &input)
    {
    
      WSNewPacket(link);
      WSPutFunction(link, "ToExpression", 1);
      WSPutString(link, input.c_str());
        
      wait_for_packet();
      input = "";
      
    }
    
    // Load package
    void load_feynrules();
    
    // Load model
    void load_model(std::string model);
    
    // Load restriction e.g. diagonal CKM, etc.
    void load_restriction(std::string model, std::string rst);
    
    // Hermiticity check
    void check_herm();
    
    // Get model name
    void get_modelname(std::string&);
    
    // Particle list
    void get_partlist(std::vector<Particle>&);
  
    // Parameters list
    void get_paramlist(std::vector<Parameter>&);
    
    // Set gauge: unitary, or feynman
    void set_gauge(std::string);
    
    // Outputs: CalcHEP, MadGraph, ... 
    void write_ch_output();
    void write_madgraph_output();
    
};

// Everything
void all_feynrules(Options, std::vector<Particle>&, std::vector<Parameter>&, Outputs&);
#endif
