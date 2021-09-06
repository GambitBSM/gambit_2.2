//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Utility functions for bibtex files
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Tomas Gonzalo
///          (gonzalo@physik.rwth-aachen.de)
///  \date 2021 Sep
///
///  *********************************************


#include <vector>
#include <fstream>
#include <regex>

#include "gambit/Utils/util_types.hpp"
#include "gambit/Utils/bibtex_functions.hpp"
#include "gambit/Utils/stream_overloads.hpp"

namespace Gambit
{

  // Constructor
  BibTeX::BibTeX(str bibtex_file)
  {
    _filename = bibtex_file;
    std::ifstream file(bibtex_file);
    str line;

    if(file.fail())
    {
      std::ostringstream errmsg;
      errmsg << "Error opening bibtex file " << bibtex_file
           << ". Please make sure that bibtex fiel exists." << std::endl;
      utils_error().raise(LOCAL_INFO,errmsg.str());
    }

    // Loop over lines
    str entry;
    while(getline(file, line))
    {
      // If it's a comment ignore it
      if(line[0] == '#') continue;

      if(line[0] == '@')
      {
        std::smatch matched_entry;
        std::regex_search(line, matched_entry, std::regex("@.*\\{(.*),"));
        entry = matched_entry[1].str();
        if(std::find(_bibtex_entries.begin(),_bibtex_entries.end(),entry) == _bibtex_entries.end())
          _bibtex_entries.push_back(entry);
        continue;
      }

      // If there's a closing bracket, reset entry
      if(line[0] == '}')
      {
        entry = "";
        continue;
      }

      // Otherwise, if entry has been initialized fill bibtex data
      if(entry != "")
      {
        str key = line.substr(0, line.find('='));
        str data = line.substr(line.find('=')+1, std::string::npos);
        _bibtex_data[entry].insert(sspair(key,data));
      }
    }
  }

  // Return input file name
  const str BibTeX::filename() const { return _filename; }

  // Get the list of bibtex entries
  const std::vector<str> BibTeX::getBibTeXEntries() const { return _bibtex_entries; }

  // Get the bibtex data
  const std::map<str,std::map<str,str>> BibTeX::getBibTeXData() const { return _bibtex_data; }

  // Drop a bibtex file with all stored data
  void BibTeX::dropBibTeXFile(str output_filename) const
  {
    std::ofstream out(output_filename);

    for(const auto& entry : _bibtex_data)
    {
      out << "@article{" << entry.first << "," << std::endl;

      for(const auto& data : entry.second)
      {
        out << data.first << " = " << data.second << std::endl;
      }
      out << "}" << std::endl << std::endl;

    }

  }

  // Drop a bibtex file only with given set of keys
  void BibTeX::dropBibTeXFile(std::vector<str> &citation_keys, str output_filename) const
  {
    std::ofstream out(output_filename);

    for(str key : citation_keys)
    {
      if(std::find(_bibtex_entries.begin(), _bibtex_entries.end(), key) == _bibtex_entries.end())
      {
        std::ostringstream errmsg;
        errmsg << "Error writing bibtex file. Could not find key " << key << "in BibTeX database." << std::endl;
        utils_error().raise(LOCAL_INFO,errmsg.str());
      }
      std::map<str,str> entry = _bibtex_data.at(key);
      out << "@article{" << key << "," << std::endl;

      for(const auto& data : entry)
      {
        out << data.first << " = " << data.second << std::endl;
      }
      out << "}" << std::endl << std::endl;

    }

  }
  // Drop a sample tex file citing a given set of keys
  void BibTeX::dropTeXFile(std::vector<str> &citation_keys, str tex_filename, str bibtex_filename) const
  {
    std::ofstream out(tex_filename);
 
    out << "\\documentclass{article}" << std::endl << std::endl;

    out << "\\title{Your Paper}" << std::endl;
    out << "\\author{You}" << std::endl << std::endl;

    out << "\\begin{document}" << std::endl;
    out << "\\maketitle" << std::endl << std::endl;

    out << "This is a sample TeX File citing all relevant references" << std::endl;
    out << "In this scan you have used the following tools~\\cite{";

    for(str key : citation_keys)
    {
      if(std::find(_bibtex_entries.begin(), _bibtex_entries.end(), key) == _bibtex_entries.end())
      {
        std::ostringstream errmsg;
        errmsg << "Error writing bibtex file. Could not find key " << key << "in BibTeX database." << std::endl;
        utils_error().raise(LOCAL_INFO,errmsg.str());
      }
      out << key << ",";
    }
    out << "}." << std::endl << std::endl;

    out << "\\bibliography{" << bibtex_filename.erase(bibtex_filename.find(".bib")) << "}" << std::endl << std::endl;

    out << "\\end{document}";

  }
}

