// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Author: Mathias Walzer, Sven Nahnsen $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CHEMISTRY/EnzymaticDigestion.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/QcMLFile.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>
#include <OpenMS/FORMAT/ControlledVocabulary.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/FORMAT/MzIdentMLFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/KERNEL/SpectrumHelper.h>

#include <QFileInfo>
//~ #include <boost/regex.hpp>

#include <vector>
#include <map>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
    @page UTILS_QCCalculator QCCalculator

    @brief Calculates basic quality parameters from MS experiments and compiles data for subsequent QC into a qcML file.

    <CENTER>
      <table>
        <tr>
        <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
        <td VALIGN="middle" ROWSPAN=3> \f$ \longrightarrow \f$ QCCalculator \f$ \longrightarrow \f$</td>
        <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
        </tr>
        <tr>
        <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_FeatureFinderCentroided </td>
        <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref UTILS_QCMerger </td>
        </tr>
        <tr>
        <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_XTandemAdapter </td>
        <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref UTILS_QCExporter </td>
        </tr>
      </table>
    </CENTER>

    The calculated quality parameters or data compiled as attachments for easy plotting input include file origin, spectra distribution, aquisition details, ion current stability ( & TIC ), id accuracy statistics and feature statistics.
    The MS experiments base name is used as name to the qcML element that is comprising all quality parameter values for the given run (including the given downstream analysis data).

    - @p id produces quality parameter values for the identification file; this file should contain either only the final psm to each spectrum (1 PeptideHit per identified spectrum) or have the PeptideHits sorted to 'best' first, where 'best' depends on the use case.
    - @p feature produces quality parameter values for the feature file; this file can be either mapped or unmapped, the latter reulting in less metrics available.
    - @p consensus produces quality parameter values for the consensus file;
    some quality parameter calculation are only available if both feature and ids are given.
    - @p remove_duplicate_features only needed when you work with a set of merged features. Then considers duplicate features only once.

    Output is in qcML format (see parameter @p out) which can be viewed directly in a modern browser (chromium, firefox, safari).

    @note For mzid in-/out- put, due to legacy reason issues you are temporarily asked to use IDFileConverter as a wrapper.
    <B>The command line parameters of this tool are:</B>
    @verbinclude UTILS_QCCalculator.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude UTILS_QCCalculator.html

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wshadow"

class TOPPQCCalculator :
  public TOPPBase
{
public:
  TOPPQCCalculator() :
    TOPPBase("QCCalculator",
      "Calculates basic quality parameters from MS experiments and subsequent analysis data as identification or feature detection.",
      false,
      {{ "Walzer M, Pernas LE, Nasso S, Bittremieux W, Nahnsen S, Kelchtermans P,  Martens, L",
         "qcML: An Exchange Format for Quality Control Metrics from Mass Spectrometry Experiments",
         "Molecular & Cellular Proteomics 2014; 13(8)" , "10.1074/mcp.M113.035907"
      }})
  {
  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "raw data input file (this is relevant if you want to look at MS1, MS2 and precursor peak information)");
    setValidFormats_("in", ListUtils::create<String>("mzML,mgf"));
    registerOutputFile_("out", "<file>", "", "Your qcML file.");
    setValidFormats_("out", ListUtils::create<String>("qcML"));
    registerInputFile_("id", "<file>", "", "Input idXML file containing the identifications. Your identifications will be exported in an easy-to-read format", false);
    setValidFormats_("id", ListUtils::create<String>("idXML,mzid"));
    registerInputFile_("feature", "<file>", "", "feature input file (this is relevant for most QC issues)", false);
    setValidFormats_("feature", ListUtils::create<String>("featureXML"));
    registerInputFile_("consensus", "<file>", "", "consensus input file (this is only used for charge state deconvoluted output. Use the consensusXML output form the DeCharger)", false);
    setValidFormats_("consensus", ListUtils::create<String>("consensusXML"));
    registerFlag_("remove_duplicate_features", "This flag should be set, if you work with a set of merged features.");
    //~ registerFlag_("MS1", "This flag should be set, if you want to work with MS1 stats.");
    //~ registerFlag_("MS2", "This flag should be set, if you want to work with MS2 stats.");
  }

//  void calculateSNident (PeptideHit& hit, MSSpectrum& spec)
//  {
    // TODO
//  }

  float calculateSNmedian (MSSpectrum& spec, bool norm = true)
  {
    if (spec.size() == 0) return 0;
    float median = 0;
    float maxi = 0;
    spec.sortByIntensity();

    if (spec.size() % 2 == 0)
    {
      median = (spec[spec.size() / 2 - 1].getIntensity() + spec[spec.size() / 2].getIntensity()) / 2;
    }
    else
    {
      median = spec[spec.size() / 2].getIntensity();
    }
    maxi = spec.back().getIntensity();
    if (!norm)
    {
      float sn_by_max2median = maxi / median;
      return sn_by_max2median;
    }

    float sign_int= 0;
    float nois_int = 0;
    size_t sign_cnt= 0;
    size_t nois_cnt = 0;
    for (MSSpectrum::const_iterator pt = spec.begin(); pt != spec.end(); ++pt)
    {
      if (pt->getIntensity() <= median)
      {
        ++nois_cnt;
        nois_int += pt->getIntensity();
      }
      else
      {
        ++sign_cnt;
        sign_int += pt->getIntensity();
      }
    }
    float sn_by_max2median_norm = (sign_int / sign_cnt) / (nois_int / nois_cnt);

    return sn_by_max2median_norm;
  }

  QcMLFile::QualityParameter fillQualityParameter (string id, string name, string ref, string acc)
  {
    QcMLFile::QualityParameter qp;
    qp.id = id;
    qp.name = name;
    qp.cvRef = ref;
    qp.cvAcc = acc;

    return qp;
  }

  QcMLFile::Attachment fillAttachment (string id, string name, string qp, string ref, string acc, vector<OpenMS::String> col_names)
  {
    QcMLFile::Attachment at;
    at.id = id;
    at.name = name;
    at.qualityRef = qp;
    at.cvRef = ref;
    at.cvAcc = acc;
    at.colTypes = col_names;

    return at;
  }

  string fetchCVTermNameOrDefault (ControlledVocabulary& cv, string cvAcc, string def_arg)
  {
    try
    {
      const ControlledVocabulary::CVTerm& term = cv.getTerm(cvAcc);
      return term.name;
    }
    catch (...)
    {
      return def_arg;
    }
  }

  ExitCodes main_(int, const char**) override
  {
    vector<ProteinIdentification> prot_ids;
    vector<PeptideIdentification> pep_ids;

    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------
    String inputfile_id = getStringOption_("id");
    String inputfile_feature = getStringOption_("feature");
    String inputfile_consensus = getStringOption_("consensus");
    String inputfile_raw = getStringOption_("in");
    String outputfile_name = getStringOption_("out");

    //~ bool Ms1(getFlag_("MS1"));
    //~ bool Ms2(getFlag_("MS2"));
    bool remove_duplicate_features(getFlag_("remove_duplicate_features"));

    //-------------------------------------------------------------
    // fetch vocabularies
    //------------------------------------------------------------
    ControlledVocabulary cv;
    cv.loadFromOBO("PSI-MS", File::find("/CV/psi-ms.obo"));
    cv.loadFromOBO("QC", File::find("/CV/qc-cv.obo"));

    QcMLFile qcmlfile;

    //-------------------------------------------------------------
    // MS acquisition
    //------------------------------------------------------------
    String base_name = QFileInfo(QString::fromStdString(inputfile_raw)).baseName();
    MSExperiment exp;
    cout << "Reading spectrum file ..." << endl;

    FileHandler fh;
    FileTypes::Type in_type = fh.getType(inputfile_raw);
    if (!fh.loadExperiment(inputfile_raw, exp, in_type))
    {
      writeLog_("Unsupported or corrupt input file. Aborting!");
      printUsage_();
      return ILLEGAL_PARAMETERS;
    }

    //---prep input
    exp.sortSpectra();
    UInt min_mz = std::numeric_limits<UInt>::max();
    UInt max_mz = 0;
    std::map<Size, UInt> mslevelcounts;

    qcmlfile.registerRun(base_name,base_name); //TODO use UIDs

    QcMLFile::QualityMetric qm = QcMLFile::QualityMetric();
    qm.value = exp.getInstrument().getName();
    qm.id = base_name + "_instrument_name";
    qm.name = "instrument model";
    qm.cvAcc = "MS:1000031";
    qm.cvRef = "MS";
    qcmlfile.addRunQualityMetric(base_name, qm);

//    qp = fillQualityParameter(base_name + "_date", "completion time", "MS", "MS:1000747");
//    qp.value = exp.getDateTime().getDate();
//    qcmlfile.addRunQualityParameter(base_name, qp);

    //---precursors and SN
    qm = QcMLFile::QualityMetric();
        qm.id = base_name + "_precursors";
        qm.name = "precursors";
        qm.cvAcc = "QC:0000044";
        qm.cvRef = "QC";
        qm.content_value = String(exp.size());
        qm.content_id = base_name + "_precursors_table";
        qm.content_name = "precursortable";
        qm.content_cvAcc = "linkedlist";
        qm.content_cvRef = "QC";
    for (Size i = 0; i < exp.size(); ++i)
    {
      mslevelcounts[exp[i].getMSLevel()]++;
      if (exp[i].getMSLevel() == 2)
      {
        if (exp[i].getPrecursors().front().getMZ() < min_mz)
        {
          min_mz = exp[i].getPrecursors().front().getMZ();
        }
        if (exp[i].getPrecursors().front().getMZ() > max_mz)
        {
          max_mz = exp[i].getPrecursors().front().getMZ();
        }
        qm.content["MS:1000894_[sec]"].push_back(exp[i].getRT());
        qm.content["MS:1000040"].push_back(exp[i].getPrecursors().front().getMZ());
        qm.content["MS:1000041"].push_back(exp[i].getPrecursors().front().getCharge());
        qm.content["S/N"].push_back(calculateSNmedian(exp[i]));
        qm.content["peak count"].push_back(exp[i].size());
      }
    }
    qcmlfile.addRunQualityMetric(base_name, qm);

//    //---aquisition results qp
//    qp_acc = "QC:0000006";
//    qp = fillQualityParameter(base_name + "_ms1aquisition",
//                              fetchCVTermNameOrDefault(cv, qp_acc, "number of ms1 spectra"),
//                              "QC",
//                              qp_acc);
//    qp.value = String(mslevelcounts[1]);
//    qcmlfile.addRunQualityParameter(base_name, qp);

//    qp_acc = "QC:0000007";
//    qp = fillQualityParameter(base_name + "_ms2aquisition",
//                              fetchCVTermNameOrDefault(cv, qp_acc, "number of ms2 spectra"),
//                              "QC",
//                              qp_acc);
//    qp.value = String(mslevelcounts[2]);
//    qcmlfile.addRunQualityParameter(base_name, qp);

//    qp_acc = "QC:0000008";
//    qp = fillQualityParameter(base_name + "_chromaquisition",
//                              fetchCVTermNameOrDefault(cv, qp_acc, "number of chromatograms"),
//                              "QC",
//                              qp_acc);
//    qp.value = String(exp.getChromatograms().size());
//    qcmlfile.addRunQualityParameter(base_name, qp);


//    at_acc = "QC:0000009";
//    at = fillAttachment(base_name + "_mzrange",
//                        fetchCVTermNameOrDefault(cv, at_acc, "MS MZ aquisition ranges"),
//                        msaq_ref,
//                        "QC",
//                        at_acc,
//                        {"QC:0000010","QC:0000011"});
//    std::vector<String> rowmz;
//    rowmz.push_back(String(min_mz));
//    rowmz.push_back(String(max_mz));
//    at.tableRows.push_back(rowmz);
//    qcmlfile.addRunAttachment(base_name, at);

//    at_acc = "QC:0000012";
//    at = fillAttachment(base_name + "_rtrange",
//                        fetchCVTermNameOrDefault(cv, at_acc, "MS RT aquisition ranges"),
//                        msaq_ref,
//                        "QC",
//                        at_acc,
//                        {"QC:0000013","QC:0000014"});
//    std::vector<String> rowrt;
//    rowrt.push_back(String(exp.begin()->getRT()));
//    rowrt.push_back(String(exp.getSpectra().back().getRT()));
//    at.tableRows.push_back(rowrt);
//    qcmlfile.addRunAttachment(base_name, at);


//    //---ion current stability ( & tic ) qp
//    at_acc = "QC:0000022";
//    at = fillAttachment(base_name + "_tics",
//                        fetchCVTermNameOrDefault(cv, at_acc, "MS TICs"),
//                        msaq_ref,
//                        "QC",
//                        at_acc,
//                        {"MS:1000894_[sec]","MS:1000285"});
//    Size below_10k = 0;
//    std::vector<OpenMS::Chromatogram> chroms = exp.getChromatograms();
//    if (!chroms.empty()) //real TIC from the mzML
//    {
//      for (Size t = 0; t < chroms.size(); ++t)
//      {
//        if (chroms[t].getChromatogramType() == ChromatogramSettings::TOTAL_ION_CURRENT_CHROMATOGRAM)
//        {
//          for (Size i = 0; i < chroms[t].size(); ++i)
//          {
//            double sum = chroms[t][i].getIntensity();
//            if (sum < 10000)
//            {
//              ++below_10k;
//            }
//            std::vector<String> row;
//            row.push_back(chroms[t][i].getRT() * 60);
//            row.push_back(sum);
//            at.tableRows.push_back(row);
//          }
//          break;  // what if there are more than one? should generally not be though ...
//        }
//      }
//      qcmlfile.addRunAttachment(base_name, at);

//      qp = QcMLFile::QualityParameter();
//      qp.id = base_name + "_ticslump"; ///< Identifier
//      qp.cvRef = "QC"; ///< cv reference
//      qp.cvAcc = "QC:0000023";
//      qp.value = String((100 / exp.size()) * below_10k);
//      try
//      {
//        const ControlledVocabulary::CVTerm& term = cv.getTerm(qp.cvAcc);
//        qp.name = term.name; ///< Name
//      }
//      catch (...)
//      {
//        qp.name = "percentage of tic slumps"; ///< Name
//      }
//      qcmlfile.addRunQualityParameter(base_name, qp);
//    }

//    // -- reconstructed TIC or RIC from the MS1 intensities
//    at_acc = "QC:0000056";
//    at = fillAttachment(base_name + "_rics",
//                        fetchCVTermNameOrDefault(cv, at_acc, "MS RICs"),
//                        msaq_ref,
//                        "QC",
//                        at_acc,
//                        {"MS:1000894_[sec]","MS:1000285","S/N","peak count"});
//    Size prev = 0;
//    below_10k = 0;
//    Size jumps = 0;
//    Size drops = 0;
//    Size fact = 10;
//    for (Size i = 0; i < exp.size(); ++i)
//    {
//      if (exp[i].getMSLevel() == 1)
//      {
//        UInt sum = 0;
//        for (Size j = 0; j < exp[i].size(); ++j)
//        {
//          sum += exp[i][j].getIntensity();
//        }
//        if (prev > 0 && sum > fact * prev)  // no jumps after complete drops (or [re]starts)
//        {
//          ++jumps;
//        }
//        else if (sum < fact*prev)
//        {
//          ++drops;
//        }
//        if (sum < 10000)
//        {
//          ++below_10k;
//        }
//        prev = sum;
//        std::vector<String> row;
//        row.push_back(exp[i].getRT());
//        row.push_back(sum);
//        row.push_back(calculateSNmedian(exp[i]));
//        row.push_back(exp[i].size());
//        at.tableRows.push_back(row);
//      }
//    }
//    qcmlfile.addRunAttachment(base_name, at);

//    qp_acc = "QC:0000057";
//    qp = fillQualityParameter(base_name + "_ricslump", fetchCVTermNameOrDefault(cv, qp_acc, "percentage of ric slumps"), "QC", qp_acc);
//    qp.value = String((100 / exp.size()) * below_10k);
//    qcmlfile.addRunQualityParameter(base_name, qp);

//    qp_acc = "QC:0000059";
//    qp = fillQualityParameter(base_name + "_ricjump",
//                              fetchCVTermNameOrDefault(cv, qp_acc, "IS-1A"),
//                              "QC",
//                              qp_acc);
//    qp.value = String(jumps);
//    qcmlfile.addRunQualityParameter(base_name, qp);

//    qp_acc = "QC:0000060";
//    qp = fillQualityParameter(base_name + "_ricdump",
//                              fetchCVTermNameOrDefault(cv, qp_acc, "IS-1B"),
//                              "QC",
//                              qp_acc);
//    qp.value = String(drops);
//    qcmlfile.addRunQualityParameter(base_name, qp);

//    //---injection times MSn
//    at_acc = "QC:0000018";
//    at = fillAttachment(base_name + "_ms2inj",
//                        fetchCVTermNameOrDefault(cv, at_acc, "MS2 injection time"),
//                        msaq_ref,
//                        "QC",
//                        at_acc,
//                        {"MS:1000894_[sec]","MS:1000927"});
//    for (Size i = 0; i < exp.size(); ++i)
//    {
//      if (exp[i].getMSLevel() > 1 && !exp[i].getAcquisitionInfo().empty())
//      {
//        for (Size j = 0; j < exp[i].getAcquisitionInfo().size(); ++j)
//        {
//          if (exp[i].getAcquisitionInfo()[j].metaValueExists("MS:1000927"))
//          {
//            std::vector<String> row;
//            row.push_back(String(exp[i].getRT()));
//            row.push_back(exp[i].getAcquisitionInfo()[j].getMetaValue("MS:1000927"));
//            at.tableRows.push_back(row);
//          }
//        }
//      }
//    }
//    if (!at.tableRows.empty())
//    {
//      qcmlfile.addRunAttachment(base_name, at);
//    }

//    //-------------------------------------------------------------
//    // MS  id
//    //------------------------------------------------------------
//    if (inputfile_id != "")
//    {
//      FileTypes::Type in_type = FileHandler::getTypeByFileName(inputfile_id);

//      if (in_type == FileTypes::MZIDENTML)
//      {
//        MzIdentMLFile().load(inputfile_id, prot_ids, pep_ids);   //, document_id);
//      }
//      else
//      {
//        IdXMLFile().load(inputfile_id, prot_ids, pep_ids);
//      }
//      cerr << "idXML read ended. Found " << pep_ids.size() << " peptide identifications." << endl;

//      ProteinIdentification::SearchParameters params = prot_ids[0].getSearchParameters();
//      vector<String> var_mods = params.variable_modifications;
//      //~ boost::regex re("(?<=[KR])(?=[^P])");

//      string msid_ref = base_name + "_msid";
//      qp_acc = "QC:0000025";
//      qp = fillQualityParameter(msid_ref, fetchCVTermNameOrDefault(cv, qp_acc, "percentage of ric slumps"), "QC", qp_acc);
//      qcmlfile.addRunQualityParameter(base_name, qp);

//      at_acc = "QC:0000026";
//      at = fillAttachment(base_name + "_idsetting",
//                          fetchCVTermNameOrDefault(cv, at_acc, "MS id settings"),
//                          msid_ref,
//                          "QC",
//                          at_acc,
//                          {"MS:1001013","MS:1001016","MS:1001020"});
//      std::vector<String> row;
//      row.push_back(String(prot_ids.front().getSearchParameters().db));
//      row.push_back(String(prot_ids.front().getSearchParameters().db_version));
//      row.push_back(String(prot_ids.front().getSearchParameters().taxonomy));
//      at.tableRows.push_back(row);
//      qcmlfile.addRunAttachment(base_name, at);

//      //---count missed cleavage numbers
//      UInt spectrum_count = 0;
//      Size peptide_hit_count = 0;
//      UInt runs_count = 0;
//      Size protein_hit_count = 0;
//      set<String> peptides;
//      set<String> proteins;
//      Size missedcleavages = 0;
//      for (Size i = 0; i < pep_ids.size(); ++i)
//      {
//        if (!pep_ids[i].empty())
//        {
//          ++spectrum_count;
//          peptide_hit_count += pep_ids[i].getHits().size();
//          const vector<PeptideHit>& temp_hits = pep_ids[i].getHits();
//          for (Size j = 0; j < temp_hits.size(); ++j)
//          {
//            peptides.insert(temp_hits[j].getSequence().toString());
//          }
//        }
//      }
//      for (set<String>::iterator it = peptides.begin(); it != peptides.end(); ++it)
//      {
//        for (String::const_iterator st = it->begin(); st != it->end() - 1; ++st)
//        {
//          if (*st == 'K' || *st == 'R')
//          {
//            ++missedcleavages;
//          }
//        }
//      }
//      qp_acc = "QC:0000037";
//      qp = fillQualityParameter(base_name + "_misscleave", fetchCVTermNameOrDefault(cv, qp_acc, "total number of missed cleavages"), "QC", qp_acc);
//      qp.value = missedcleavages;
//      qcmlfile.addRunQualityParameter(base_name, qp);

//      for (Size i = 0; i < prot_ids.size(); ++i)
//      {
//        ++runs_count;
//        protein_hit_count += prot_ids[i].getHits().size();
//        const vector<ProteinHit>& temp_hits = prot_ids[i].getHits();
//        for (Size j = 0; j < temp_hits.size(); ++j)
//        {
//          proteins.insert(temp_hits[j].getAccession());
//        }
//      }
//      qp_acc = "QC:0000032";
//      qp = fillQualityParameter(base_name + "_totprot", fetchCVTermNameOrDefault(cv, qp_acc, "total number of identified proteins"), "QC", qp_acc);
//      qp.value = protein_hit_count;
//      qcmlfile.addRunQualityParameter(base_name, qp);

//      qp_acc = "QC:0000033";
//      qp = fillQualityParameter(base_name + "_totuniprot", fetchCVTermNameOrDefault(cv, qp_acc, "total number of uniquely identified proteins"), "QC", qp_acc);
//      qp.value = String(proteins.size());
//      qcmlfile.addRunQualityParameter(base_name, qp);

//      qp_acc = "QC:0000029";
//      qp = fillQualityParameter(base_name + "_psms", fetchCVTermNameOrDefault(cv, qp_acc, "total number of PSM"), "QC", qp_acc);
//      qp.value = String(spectrum_count);
//      qcmlfile.addRunQualityParameter(base_name, qp);

//      qp_acc = "QC:0000030";
//      qp = fillQualityParameter(base_name + "_totpeps", fetchCVTermNameOrDefault(cv, qp_acc, "total number of identified petides"), "QC", qp_acc);
//      qp.value = String(peptide_hit_count);
//      qcmlfile.addRunQualityParameter(base_name, qp);

//      qp_acc = "QC:0000030";
//      qp = fillQualityParameter(base_name + "_totunipeps", fetchCVTermNameOrDefault(cv, qp_acc, "total number of uniquely identified petides"), "QC", qp_acc);
//      qp.value = String(peptides.size());
//      qcmlfile.addRunQualityParameter(base_name, qp);

//      at_acc = "QC:0000038";
//      at = fillAttachment(base_name + "_massacc",
//                          fetchCVTermNameOrDefault(cv, at_acc, "delta ppm tables"),
//                          msid_ref,
//                          "QC",
//                          at_acc,
//      {"RT","MZ","Score","PeptideSequence","Charge","TheoreticalWeight","delta_ppm","S/N"});
//      for (UInt w = 0; w < var_mods.size(); ++w)
//      {
//        at.colTypes.push_back(String(var_mods[w]).substitute(' ', '_'));
//      }

//      std::vector<double> deltas;
//      //~ prot_ids[0].getSearchParameters();
//      for (vector<PeptideIdentification>::iterator it = pep_ids.begin(); it != pep_ids.end(); ++it)
//      {
//        if (!it->getHits().empty())
//        {
//          std::vector<String> row;
//          row.push_back(it->getRT());
//          row.push_back(it->getMZ());
//          PeptideHit tmp = it->getHits().front();  //N.B.: depends on score & sort
//          vector<UInt> pep_mods;
//          for (UInt w = 0; w < var_mods.size(); ++w)
//          {
//            pep_mods.push_back(0);
//          }
//          for (AASequence::ConstIterator z =  tmp.getSequence().begin(); z != tmp.getSequence().end(); ++z)
//          {
//            Residue res = *z;
//            String temp;
//            if (res.isModified() && res.getModificationName() != "Carbamidomethyl")
//            {
//              temp = res.getModificationName() + " (" + res.getOneLetterCode()  + ")";
//              //cout<<res.getModification()<<endl;
//              for (UInt w = 0; w < var_mods.size(); ++w)
//              {
//                if (temp == var_mods[w])
//                {
//                  //cout<<temp;
//                  pep_mods[w] += 1;
//                }
//              }
//            }
//          }

//          row.push_back(tmp.getScore());
//          row.push_back(tmp.getSequence().toString().removeWhitespaces());
//          row.push_back(tmp.getCharge());
//          row.push_back(String((tmp.getSequence().getMonoWeight() + tmp.getCharge() * Constants::PROTON_MASS_U) / tmp.getCharge()));
//          double dppm = /* std::abs */ (OpenMS::getDeltaPpm(((tmp.getSequence().getMonoWeight() + tmp.getCharge() * Constants::PROTON_MASS_U) / tmp.getCharge()), it->getMZ()));
//          row.push_back(String(dppm));
////          row.push_back(String(calculateSNident(tmp)));
//          deltas.push_back(dppm);
//          for (UInt w = 0; w < var_mods.size(); ++w)
//          {
//            row.push_back(pep_mods[w]);
//          }
//          at.tableRows.push_back(row);
//        }
//      }
//      qcmlfile.addRunAttachment(base_name, at);

//      qp_acc = "QC:0000040";
//      qp = fillQualityParameter(base_name + "_mean_delta", fetchCVTermNameOrDefault(cv, qp_acc, "mean delta ppm"), "QC", qp_acc);
//      qp.value = String(OpenMS::Math::mean(deltas.begin(), deltas.end()));
//      qcmlfile.addRunQualityParameter(base_name, qp);

//      qp_acc = "QC:0000041";
//      qp = fillQualityParameter(base_name + "_median_delta", fetchCVTermNameOrDefault(cv, qp_acc, "median delta ppm"), "QC", qp_acc);
//      qp.value = String(OpenMS::Math::median(deltas.begin(), deltas.end(), false));
//      qcmlfile.addRunQualityParameter(base_name, qp);

//      qp_acc = "QC:0000035";
//      qp = fillQualityParameter(base_name + "_ratio_id", fetchCVTermNameOrDefault(cv, qp_acc, "id ratio"), "QC", qp_acc);
//      qp.value = String(double(pep_ids.size()) / double(mslevelcounts[2]));
//      qcmlfile.addRunQualityParameter(base_name, qp);
//    }

//    //-------------------------------------------------------------
//    // MS quantitation
//    //------------------------------------------------------------
//    FeatureMap map;
//    String msqu_ref = base_name + "_msqu";
//    if (inputfile_feature != "")
//    {
//      FeatureXMLFile f;
//      f.load(inputfile_feature, map);

//      cout << "Read featureXML file..." << endl;

//      //~ UInt fiter = 0;
//      map.sortByRT();
//      map.updateRanges();

//      qp_acc = "QC:0000045";
//      qp = fillQualityParameter(msqu_ref, fetchCVTermNameOrDefault(cv, qp_acc, "MS quantification result details"), "QC", qp_acc);
//      qcmlfile.addRunQualityParameter(base_name, qp);

//      qp_acc = "QC:0000046";
//      qp = fillQualityParameter(base_name + "_feature_count", fetchCVTermNameOrDefault(cv, qp_acc, "MS quantification result details"), "QC", qp_acc);
//      qp.value = String(map.size());
//      qcmlfile.addRunQualityParameter(base_name, qp);
//    }

//    if (inputfile_feature != "" && !remove_duplicate_features)
//    {

//      at_acc = "QC:0000047";
//      at = fillAttachment(base_name + "_features",
//                          fetchCVTermNameOrDefault(cv, at_acc, "features"),
//                          msqu_ref,
//                          "QC",
//                          at_acc,
//                          {"RT","MZ","Intensity","Charge","Quality","FWHM","IDs"});
//      UInt fiter = 0;
//      UInt ided = 0;
//      map.sortByRT();
//      //ofstream out(outputfile_name.c_str());
//      while (fiter < map.size())
//      {
//        std::vector<String> row;
//        row.push_back(map[fiter].getMZ());
//        row.push_back(map[fiter].getRT());
//        row.push_back(map[fiter].getIntensity());
//        row.push_back(map[fiter].getCharge());
//        row.push_back(map[fiter].getOverallQuality());
//        row.push_back(map[fiter].getWidth());
//        row.push_back(map[fiter].getPeptideIdentifications().size());
//        if (map[fiter].getPeptideIdentifications().size() > 0)
//        {
//          ++ided;
//        }
//        fiter++;
//        at.tableRows.push_back(row);
//      }
//      qcmlfile.addRunAttachment(base_name, at);

//      qp_acc = "QC:0000058";
//      qp = fillQualityParameter(base_name + "_idfeature_count", fetchCVTermNameOrDefault(cv, qp_acc, "number of identified features"), "QC", qp_acc);
//      qp.value = ided;
//      qcmlfile.addRunQualityParameter(base_name, qp);

//    }
//    else if (inputfile_feature != "" && remove_duplicate_features)
//    {
//      at_acc = "QC:0000047";
//      at = fillAttachment(base_name + "_features",
//                          fetchCVTermNameOrDefault(cv, at_acc, "features"),
//                          msqu_ref,
//                          "QC",
//                          at_acc,
//                          {"RT","MZ","Intensity","Charge"});
//      FeatureMap map, map_out;
//      FeatureXMLFile f;
//      f.load(inputfile_feature, map);
//      UInt fiter = 0;
//      map.sortByRT();
//      while (fiter < map.size())
//      {
//        FeatureMap map_tmp;
//        for (UInt k = fiter; k <= map.size(); ++k)
//        {
//          if (abs(map[fiter].getRT() - map[k].getRT()) < 0.1)
//          {
//            //~ cout << fiter << endl;
//            map_tmp.push_back(map[k]);
//          }
//          else
//          {
//            fiter = k;
//            break;
//          }
//        }
//        map_tmp.sortByMZ();
//        UInt retif = 1;
//        map_out.push_back(map_tmp[0]);
//        while (retif < map_tmp.size())
//        {
//          if (abs(map_tmp[retif].getMZ() - map_tmp[retif - 1].getMZ()) > 0.01)
//          {
//            cout << "equal RT, but mass different" << endl;
//            map_out.push_back(map_tmp[retif]);
//          }
//          retif++;
//        }
//      }
//      qcmlfile.addRunAttachment(base_name, at);
//    }
//    if (inputfile_consensus != "")
//    {
//      cout << "Reading consensusXML file..." << endl;
//      ConsensusXMLFile f;
//      ConsensusMap map;
//      f.load(inputfile_consensus, map);
//      //~ String CONSENSUS_NAME = "_consensus.tsv";
//      //~ String combined_out = outputfile_name + CONSENSUS_NAME;
//      //~ ofstream out(combined_out.c_str());

//      at_acc = "QC:xxxxxxxx";
//      at = fillAttachment(base_name + "_consensuses",
//                          fetchCVTermNameOrDefault(cv, at_acc, "consensuspoints"),
//                          msqu_ref,
//                          "QC",
//                          at_acc,
//      {"Native_spectrum_ID","DECON_RT_(sec)","DECON_MZ_(Th)","DECON_Intensity","Feature_RT_(sec)","Feature_MZ_(Th)","Feature_Intensity","Feature_Charge"});
//      for (ConsensusMap::const_iterator cmit = map.begin(); cmit != map.end(); ++cmit)
//      {
//        const ConsensusFeature& CF = *cmit;
//        for (ConsensusFeature::const_iterator cfit = CF.begin(); cfit != CF.end(); ++cfit)
//        {
//          std::vector<String> row;
//          FeatureHandle FH = *cfit;
//          row.push_back(CF.getMetaValue("spectrum_native_id"));
//          row.push_back(CF.getRT()); row.push_back(CF.getMZ());
//          row.push_back(CF.getIntensity());
//          row.push_back(FH.getRT());
//          row.push_back(FH.getMZ());
//          row.push_back(FH.getCharge());
//          at.tableRows.push_back(row);
//        }
//      }
//      qcmlfile.addRunAttachment(base_name, at);
//    }


    //-------------------------------------------------------------
    // finalize
    //------------------------------------------------------------
    qcmlfile.store(outputfile_name);
    return EXECUTION_OK;
  }

};

#pragma clang diagnostic pop

int main(int argc, const char** argv)
{
  TOPPQCCalculator tool;
  return tool.main(argc, argv);
}

/// @endcond
