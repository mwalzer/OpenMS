// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
// $Maintainer: Mathias Walzer$
// $Authors: Mathias Walzer$
// --------------------------------------------------------------------------

#ifndef OPENMS_METADATA_MSQUANTIFICATIONS_H
#define OPENMS_METADATA_MSQUANTIFICATIONS_H

#include <OpenMS/METADATA/ExperimentalSettings.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
//~ #include <OpenMS/KERNEL/ConsensusFeature.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/METADATA/DataProcessing.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

#include <vector>
#include <map>

namespace OpenMS
{
  class OPENMS_DLLAPI MSQuantifications :
    public ExperimentalSettings
  {
public:
    /// @name Base type definitions
    //@{
    /// typedef docu
    typedef CVTermList ParamGroupList;         // userparams are exclusively inside the CVTermList's MetaInfoInterface

    enum QUANT_TYPES {MS1LABEL = 0, MS2LABEL, LABELFREE, SIZE_OF_QUANT_TYPES};       // derived from processing applied
    static const std::string NamesOfQuantTypes[SIZE_OF_QUANT_TYPES];
    //@}
    
    struct AnalysisSummary
    {
      AnalysisSummary()
      {
      }

      AnalysisSummary(const AnalysisSummary & rhs) :
        cv_params_(rhs.cv_params_)
      {
        user_params_ = rhs.user_params_;
        quant_type_ = rhs.quant_type_;
      }

      virtual ~AnalysisSummary()
      {
      }

      AnalysisSummary & operator=(const AnalysisSummary & rhs)
      {
        if (&rhs != this)
        {
          cv_params_ = rhs.cv_params_;
          user_params_ = rhs.user_params_;
          quant_type_ = rhs.quant_type_;
        }
        return *this;
      }

      MetaInfo user_params_;
      CVTermList cv_params_;
      QUANT_TYPES quant_type_;
    };

    struct Assay
    {
      Assay()
      {
      }

      Assay(const Assay & rhs)
      {
        uid_ = rhs.uid_;
        rfg_ref_ = rhs.rfg_ref_;
        mods_ = rhs.mods_;
      }

      virtual ~Assay()
      {
      }

      Assay & operator=(const Assay & rhs)
      {
        if (&rhs != this)
        {
          uid_ = rhs.uid_;
          rfg_ref_ = rhs.rfg_ref_;
          mods_ = rhs.mods_;
        }
        return *this;
      }

      UInt64 uid_;
      UInt64 rfg_ref_;
      std::vector<std::pair<String, DoubleReal> > mods_;
    };

    //~ InputFiles: //~ searchdb tracking version,releasedate,#entries,dbname via paramgrouplist
    //~ struct ParamGroupList
    //~ {
    //~ ParamGroupList()
    //~ {
    //~ }

    //~ ParamGroupList(const ParamGroupList& rhs)
    //~ :	cv_params(rhs.cv_params)
    //~ {
    //~ }

    //~ ~ParamGroupList()
    //~ {
    //~ }

    //~ ParamGroupList& operator = (const ParamGroupList& rhs)
    //~ {
    //~ if (&rhs != this)
    //~ {
    //~ cv_params = rhs.cv_params;
    //~ user_params = rhs.user_params;
    //~ }
    //~ return *this;
    //~ }

    //~ MetaInfoInterface user_params;
    //~ CVTermList cv_params;
    //~ };
    
    // TODO handle referencing from consensusmaps to featuremaps/rawfiles
    // TODO add ContactPerson or something to (Consensus)FeatureMap or DataProcessing (see below)
    // TODO rewrite OpenMS::DataProcessing - data not yet linked in openms core formats - below should go in analysissummary of MSQuantifications - input/output not possible to be carried along
    //~ if(DataProcessing::NamesOfProcessingAction[*it] == String("Quantitation"))
    //~ {
    //~ if (processing.getSoftware().getName()==String("SILACAnalyzer"))
    //~ {
    //~ experiment_type = MS1LABEL;
    //~ }
    //~ else if (processing.getSoftware().getName()==String("ITRAQAnalyzer"))
    //~ {
    //~ experiment_type = MS2LABEL;
    //~ }
    //~ else
    //~ {
    //~ experiment_type = LABELFREE;
    //~ }
    //~ }
    //~ QUANT_TYPES experiment_type = MS1LABEL;

    /// Constructor
    MSQuantifications();

    /// Destructor
    ~MSQuantifications();

    /// Copy constructor
    MSQuantifications(const MSQuantifications & source);

    /// Assignment operator
    MSQuantifications & operator=(const MSQuantifications & source);

    /// Equality operator
    bool operator==(const MSQuantifications & rhs) const;

    /// Equality operator
    bool operator!=(const MSQuantifications & rhs) const;

    // getter & setter
    const std::vector<DataProcessing> getDataProcessingList() const;
    const std::map< UInt64, MSQuantifications::Assay > & getAssays() const;
    std::map< UInt64, MSQuantifications::Assay> & getAssays();
    const std::map< UInt64, ConsensusMap > & getConsensusMaps() const;
    std::map< UInt64, ConsensusMap > & getConsensusMaps();
    //~ void setConsensusMaps(const std::map< UInt64, ConsensusMap > & );
    const std::map< UInt64, FeatureMap<> > & getFeatureMaps() const;
    std::map< UInt64, FeatureMap<> > & getFeatureMaps();
    const AnalysisSummary & getAnalysisSummary() const;
    AnalysisSummary & getAnalysisSummary();
    const std::map< UInt64, std::set<ExperimentalSettings> > & getRawFiles() const;
    const std::vector< UInt64 > & getSourceFiles() const;
    void setDataProcessingList(std::vector<DataProcessing> & dpl);
    void setAnalysisSummaryQuantType(QUANT_TYPES r);
    
    std::set<UInt64> getDataProcessingInRefs(UInt64 dp_ref) const;
    UInt64 getDataProcessingOutRefs(UInt64 dp_ref) const;
    const UInt64 fromWhichInput(const UInt64 & feat) const;
    
    /// registerers
    void addConsensusMap(ConsensusMap & m, std::vector<UInt64> file_uids);
    void registerFeatureMap(FeatureMap<> & m, UInt64 rawfile_uid);        
    const std::pair< std::vector<UInt64>,UInt64 > registerExperimentMap(MSExperiment<Peak1D> & exp, std::vector<std::vector<std::pair<String, DoubleReal> > > labels = (std::vector<std::vector<std::pair<String, DoubleReal> > >()));
    
    const std::pair< std::vector<UInt64>,UInt64 > registerExperimentMap(ExperimentalSettings & es, std::vector<DataProcessing>& dps, std::vector<std::vector<std::pair<String, DoubleReal> > > labels = (std::vector<std::vector<std::pair<String, DoubleReal> > >()));
    
    const UInt64 addExperiment( std::vector<UInt64> & assay_uids, MSExperiment<Peak1D> & exp);
    
    const UInt64 addExperiment( std::vector<UInt64> & assay_uids ,ExperimentalSettings & es, std::vector<DataProcessing>& dps);
    
    ///for mzquantml consumption
    void stubFeatureMap(FeatureMap<> & m);        
    void stubExperimentMap(FeatureMap<> & m);        

private:
    AnalysisSummary analysis_summary_;
    std::map< UInt64, std::set<ExperimentalSettings> > raw_files_group_; //implicite: raw->mzml ! so:ExperimentalSettings.getUniqueId()->ms-raw xsd:ID
    std::vector< UInt64 > source_files_; // entry is a featuremap/consensusmap UID, also present in either feature or consensus map containing the file info in DocumentIdentifier Interface 
    // TODO @mths : add missing searchdatabase files, identification files and method files

    std::vector<DataProcessing> data_processings_;
    std::multimap< UInt64, UInt64 > in_data_processings_; 
    std::map< UInt64, UInt64 > out_data_processings_; 

    std::map< UInt64, Assay > assays_;

    std::map< UInt64, FeatureMap<> > feature_maps_;
    std::map< UInt64, UInt64 > feature_to_raw_;
    std::map< UInt64, ConsensusMap > consensus_maps_;
    std::map< UInt64, std::vector<UInt64> > consensus_to_features_;

    //~ std::vector<MetaInfo> bibliographic_reference_;
    //~ std::map<String,ConsensusFeature::Ratio > ratio_calculations_;

  };

} // namespace OpenMS

#endif // OPENMS_METADATA_MSQUANTIFICATIONS_H
