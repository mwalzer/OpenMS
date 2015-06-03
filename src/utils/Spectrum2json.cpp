// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer $
// $Author: Mathias Walzer $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/FORMAT/IndexedMzMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/FORMAT/CsvFile.h>

#include <QByteArray>
#include <QFile>
#include <QString>
#include <QFileInfo>

//~ #include <QIODevice>
#include <fstream>
#include <vector>
#include <map>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
    @page UTILS_Spectrum2json Spectrum2json

    @brief will export a spectrum from a mzML file into a json object

    <CENTER>
      <table>
        <tr>
        <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
        <td VALIGN="middle" ROWSPAN=3> \f$ \longrightarrow \f$ Spectrum2json \f$ \longrightarrow \f$</td>
        <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
        </tr>
        <tr>
        <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_FileFilter </td>
        </tr>
        <tr>
        <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref javascript </td>
        </tr>
      </table>
    </CENTER>

    TBA

    <B>The command line parameters of this tool are:</B>
    @verbinclude UTILS_Spectrum2json.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude UTILS_Spectrum2json.html

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES
class TOPPSpectrum2json :
  public TOPPBase
{
public:
  TOPPSpectrum2json() :
    TOPPBase("Spectrum2json", "Will extract several qp from several run/sets in a tabular format.", false)
  {
  }

protected:
  void registerOptionsAndFlags_()
  {
    registerInputFile_("in", "<file>", String(), "Input mzml file", true);
    setValidFormats_("in", ListUtils::create<String>("mzML"));
    registerIntOption_("spectrum", "<scannumber>", -1, "The scan number of the target spectrum.", false); //set to required breaks, since toppbase checks and complains (there is no value to indicate it is missing)
    registerStringOption_("sequence", "<aminoacidsequence>", String(), "The name of the target runs or sets to be exported from. If empty, from all will be exported.", false);
    registerOutputFile_("out", "<file>", "", "Output txt file with json of given spectrum.", false);
    setValidFormats_("out", ListUtils::create<String>("txt"));
    registerFlag_("force_regular_read", "will read teh mzml file in regular fashion (slow), even if it is indexed", true);
  }

  ExitCodes main_(int, const char**)
  {
    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------
    String in = getStringOption_("in");
    int spec = getIntOption_("spectrum");
    String seq = getStringOption_("sequence");
    bool slow = getFlag_("force_regular_read");
    String out = getStringOption_("out");

    if (spec < 0) return ILLEGAL_PARAMETERS;
    stringstream ss;

    //-------------------------------------------------------------
    // reading input
    //------------------------------------------------------------
    IndexedMzMLFile ifile;
    ifile.setSkipXMLChecks(true);
    ifile.openFile(in);
    if (ifile.getParsingSuccess() && !slow)
    {
      LOG_DEBUG << "Found valid index for the mzML file " << in << std::endl;
      OpenMS::Interfaces::SpectrumPtr p = ifile.getSpectrumById(spec);
      const std::vector<double>& mz_array = p->getMZArray()->data;
      const std::vector<double>& int_array = p->getIntensityArray()->data;

      // PRECONDITION
      assert(mz_array.size() == int_array.size());
      assert(mz_array.size() > 2);

      if (seq.empty())
      {
        ss << "{ \"sequence\": null," << std::endl;
      }
      else
      {
        ss << "{ \"sequence\": \"" << seq << "\"," << std::endl;
      }
      ss << "\t\"scanNum\":" << spec+1 << "," << std::endl;
      ss << "\t\"fileName\":" << in << "," << std::endl;
      ss << "\t\"peaks\":[" << std::endl;
      for (size_t i = 0; i < int_array.size(); ++i)
      {
          ss << "\t\t[" << mz_array[i] << "," << int_array[i] << "]," << std::endl;
      }
      ss << "\t]" << std::endl;
      ss << "}" << std::endl;
    }
    else
    {
      LOG_DEBUG << "Could not detect a valid index for the mzML file " << in << "\nEither the index is not present or is not correct." << std::endl;
      MzMLFile file;
      file.setLogType(log_type_);
      //speedup stuff
      file.getOptions().setSkipXMLChecks(true);

      MSExperiment<> exp;
      file.load(in, exp);
      if (!exp.empty() && !exp[spec].empty())
      {
          if (seq.empty())
          {
            ss << "{ \"sequence\": null," << std::endl;
          }
          else
          {
            ss << "{ \"sequence\": \"" << seq << "\"," << std::endl;
          }
          ss << "\t\"scanNum\":" << spec+1 << "," << std::endl;
          ss << "\t\"fileName\":" << in << "," << std::endl;
          if (exp[spec].getMSLevel() > 1)
          {
            ss << "\t\"precursorMz\":" << exp[spec].getPrecursors().front().getMZ() << "," << std::endl;
            ss << "\t\"charge\":" << exp[spec].getPrecursors().front().getCharge() << "," << std::endl;
          }
          ss << "\t\"peaks\":[" << std::endl;
          for (size_t i = 0; i < exp[spec].size(); ++i)
          {
              ss << "\t\t[" << exp[spec][i].getMZ() << "," << exp[spec][i].getIntensity() << "]," << std::endl;
          }
          ss << "\t]" << std::endl;
          ss << "}" << std::endl;
      }
      else
      {
        LOG_DEBUG << "Could not extract spectrum " << spec << " from file " << in << std::endl;
        return ILLEGAL_PARAMETERS;
      }
    }    
    if (out.empty())
    {
      std::cout << ss.str() << std::endl;
    }
    else
    {
      TextFile txt;
      txt.addLine(ss.str());
      txt.store(out);
    }
    return EXECUTION_OK;
  }

};
int main(int argc, const char** argv)
{
  TOPPSpectrum2json tool;
  return tool.main(argc, argv);
}

/// @endcond
