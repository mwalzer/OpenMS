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

#include <OpenMS/METADATA/MSQuantifications.h>
#include <set>
#include <iostream>

using namespace std;

namespace OpenMS
{
  const std::string MSQuantifications::NamesOfQuantTypes[] = {"MS1LABEL", "MS2LABEL", "LABELFREE"};

  /// Constructor
  MSQuantifications::MSQuantifications() :
    ExperimentalSettings()
  {
  }  
  
  /// Copy constructor
  MSQuantifications::MSQuantifications(const MSQuantifications & source) :
    ExperimentalSettings(source)
  {
  }

  MSQuantifications::~MSQuantifications()
  {
  }

  /// Assignment operator
  MSQuantifications & MSQuantifications::operator=(const MSQuantifications & source)
  {
    if (&source == this)
      return *this;

    //~ reassign members
    ExperimentalSettings::operator=(source);
    //~ PersistentObject::operator=(source);

    return *this;
  }

  /// Equality operator
  bool MSQuantifications::operator==(const MSQuantifications & rhs) const
  {
    return ExperimentalSettings::operator==(rhs);
  }

  /// Equality operator
  bool MSQuantifications::operator!=(const MSQuantifications & rhs) const
  {
    return !(operator==(rhs));
  }

  /// getter & setter
  void MSQuantifications::setDataProcessingList(std::vector<DataProcessing> & dpl)
  {
    data_processings_ = dpl;
  }

  const std::vector<DataProcessing> MSQuantifications::getDataProcessingList() const
  {
    return data_processings_;
  }

  const std::map< UInt64, MSQuantifications::Assay> & MSQuantifications::getAssays() const
  {
    return assays_;
  }

  std::map< UInt64, MSQuantifications::Assay> & MSQuantifications::getAssays()
  {
    return assays_;
  }

  //~ std::map<String,ConsensusFeature::Ratio>& MSQuantifications::getRatioCalculations()
  //~ {
  //~ return ratio_calculations_;
  //~ }

  const std::map< UInt64, FeatureMap<> > & MSQuantifications::getFeatureMaps() const
  {
    return feature_maps_;
  }

  const std::map< UInt64, ConsensusMap > & MSQuantifications::getConsensusMaps() const
  {
    return consensus_maps_;
  }

  //~ void MSQuantifications::setConsensusMaps(const std::vector<ConsensusMap> & consensus_maps)
  //~ {
      //~ consensus_maps_ = consensus_maps;
  //~ }

  std::map< UInt64, ConsensusMap > & MSQuantifications::getConsensusMaps()
  {
    return consensus_maps_;
  }

  const MSQuantifications::AnalysisSummary & MSQuantifications::getAnalysisSummary() const
  {
    return analysis_summary_;
  }

  MSQuantifications::AnalysisSummary & MSQuantifications::getAnalysisSummary()
  {
    return analysis_summary_;
  }
  
  const std::map< UInt64, std::set<ExperimentalSettings> > & MSQuantifications::getRawFiles() const
  {
    return raw_files_group_;
  }
  
  const std::vector< UInt64 > & MSQuantifications::getSourceFiles() const
  {
    return source_files_;
  }

  void MSQuantifications::setAnalysisSummaryQuantType(MSQuantifications::QUANT_TYPES r)
  {
    analysis_summary_.quant_type_ = r;
  }
  
  std::set<UInt64> MSQuantifications::getDataProcessingInRefs(UInt64 dp_ref) const
  {
    std::pair<std::multimap<UInt64, UInt64>::const_iterator,std::multimap<UInt64, UInt64>::const_iterator> er = in_data_processings_.equal_range(dp_ref);
    std::set<UInt64> ret;
    for (std::multimap<UInt64, UInt64>::const_iterator it=er.first; it!=er.second; ++it)
    {
      ret.insert(it->second);
    }
    return ret;
  }
  
  UInt64 MSQuantifications::getDataProcessingOutRefs(UInt64 dp_ref) const
  {
    std::map<UInt64, UInt64>::const_iterator ret = out_data_processings_.find(dp_ref);
    if (ret != out_data_processings_.end())
    {
      return ret->second;
    }
    return 0; // TODO @mths : this could be nicer
  }
 
  const std::vector<UInt64> MSQuantifications::getFeatureMapUIDs() const
  {
    std::vector<UInt64> fmuids;
    fmuids.reserve(feature_maps_.size());
    for (std::map<UInt64, FeatureMap<> >::const_iterator it = feature_maps_.begin(); it != feature_maps_.end(); ++it) 
    {
      fmuids.push_back(it->first);
    }
    return fmuids;
  } 
  
  std::vector<FeatureMap <> > MSQuantifications::getFeatureMapVector() const
  {
    std::vector<FeatureMap<> > fms;
    fms.reserve(feature_maps_.size());
    for (std::map<UInt64, FeatureMap<> >::const_iterator it = feature_maps_.begin(); it != feature_maps_.end(); ++it) 
    {
      fms.push_back(it->second);
    }
    return fms;
  }
  
  /// registerers
  const std::pair< std::vector<UInt64>, UInt64 > MSQuantifications::registerExperimentMap(MSExperiment<Peak1D> & exp, std::vector<std::vector<std::pair<String, DoubleReal> > > labels)
  {
    //TODO check if exp.getSpectra().front() not empty, or better get and merge for each MS1-n
    return this->registerExperimentMap(exp.getExperimentalSettings(),exp.getSpectra().front().getDataProcessing(),labels);
  }
  
  const std::pair< std::vector<UInt64>, UInt64 > MSQuantifications::registerExperimentMap(ExperimentalSettings & es, std::vector<DataProcessing>& dps, std::vector<std::vector<std::pair<String, DoubleReal> > > labels)
  {
    UInt64 rfg_uid = UniqueIdGenerator::getUniqueId();
    std::vector<UInt64> auids;
    for (std::vector<std::vector<std::pair<String, DoubleReal> > >::const_iterator lit = labels.begin(); lit != labels.end(); ++lit)
    {
      Assay a;
      a.uid_ = UniqueIdGenerator::getUniqueId();
      a.rfg_ref_ = rfg_uid;
      auids.push_back(a.uid_);
      a.mods_ = (*lit);
      assays_.insert(std::pair<UInt64, Assay >(a.uid_, a));
    }
    
    if (labels.empty())
    {
      Assay a;
      a.uid_ = UniqueIdGenerator::getUniqueId();
      a.rfg_ref_ = rfg_uid;
      auids.push_back(a.uid_);
      assays_.insert(std::pair<UInt64, Assay >(a.uid_, a));
    }
  
    std::set<ExperimentalSettings> es_set;
    es_set.insert(es);
    raw_files_group_.insert( std::pair<UInt64, std::set<ExperimentalSettings> > (rfg_uid, es_set ) );
    
    data_processings_.insert(data_processings_.end(), dps.begin(), dps.end()); 
    for (std::vector<DataProcessing>::const_iterator dpit = dps.begin(); dpit != dps.end(); ++dpit)
    {
      in_data_processings_.insert( std::pair<UInt64, UInt64>(dpit->getUniqueId(), es.getUniqueId()) ); 
      out_data_processings_.insert( std::pair<UInt64, UInt64>(dpit->getUniqueId(), es.getUniqueId()) ); 
    }
    
    raw_to_assays_.insert(std::make_pair(es.getUniqueId(),auids));

    return std::make_pair(auids,es.getUniqueId());
  }

  const UInt64 MSQuantifications::addExperiment( std::vector<UInt64> & assay_uids, MSExperiment<Peak1D> & exp)
  {
    //TODO check if exp.getSpectra().front() not empty, or better get and merge for each MS1-n
    return this->addExperiment(assay_uids, exp.getExperimentalSettings(), exp.getSpectra().front().getDataProcessing());   
  }

  const UInt64 MSQuantifications::addExperiment( std::vector<UInt64> & assay_uids, ExperimentalSettings & es, std::vector<DataProcessing>& dps)
  {
    //~ UInt64 rfg_uid = UniqueIdGenerator::getUniqueId();
    for (std::vector<UInt64>::const_iterator ait = assay_uids.begin(); ait != assay_uids.end(); ++ait)
    {
      std::map< UInt64, Assay >::const_iterator assay_it = assays_.find(*ait);
      if (assay_it != assays_.end())
      {  
        std::map< UInt64, std::set<ExperimentalSettings> >::iterator rfg_it = raw_files_group_.find(assay_it->second.rfg_ref_);
        if (rfg_it != raw_files_group_.end())
        { 
          rfg_it->second.insert(es); //if not exisits
        }
        else
        {
          //TODO @mths : warn
        }
      }
      else
      {
        //TODO @mths : warn
      }
    }
    
    data_processings_.insert(data_processings_.end(), dps.begin(), dps.end());
    for (std::vector<DataProcessing>::const_iterator dpit = dps.begin(); dpit != dps.end(); ++dpit)
    {
      in_data_processings_.insert( std::pair<UInt64, UInt64>(dpit->getUniqueId(), es.getUniqueId()) ); 
      out_data_processings_.insert( std::pair<UInt64, UInt64>(dpit->getUniqueId(), es.getUniqueId()) ); 
    }
    
    raw_to_assays_.insert(std::make_pair(es.getUniqueId(),assay_uids));
    
    return es.getUniqueId();
  }
  
  void MSQuantifications::registerFeatureMap(FeatureMap<> & m, UInt64 rawfile_uid)
  {
    feature_maps_.insert(std::pair<UInt64, FeatureMap<> >(m.getUniqueId(), m)); //TODO @mths : whatif experiment from featuremap is not registered
    featuremap_to_raw_.insert(std::pair<UInt64, UInt64>(m.getUniqueId(), rawfile_uid)); 
    source_files_.push_back(m.getUniqueId()); //TODO @tue this needs filename as well 
    data_processings_.insert(data_processings_.end(), m.getDataProcessing().begin(), m.getDataProcessing().end()); 
    for (std::vector<DataProcessing>::const_iterator dpit = m.getDataProcessing().begin(); dpit != m.getDataProcessing().end(); ++dpit)
    {
      in_data_processings_.insert( std::pair<UInt64, UInt64>(dpit->getUniqueId(), rawfile_uid) ); 
      out_data_processings_.insert( std::pair<UInt64, UInt64>(dpit->getUniqueId(), m.getUniqueId()) ); 
    }
  }

  void MSQuantifications::addConsensusMap(ConsensusMap & m, std::vector<UInt64> file_uids)
  {
    consensus_maps_.insert(std::pair<UInt64, ConsensusMap >(m.getUniqueId(), m));
    consensus_to_features_.insert(std::pair<UInt64, std::vector<UInt64> >(m.getUniqueId(), file_uids));
  }  
  
  void MSQuantifications::registerInRefs(const UInt64 & dp_ref, const UInt64 & object_ref)
  {
    in_data_processings_.insert(std::make_pair(dp_ref,object_ref));
  }
  
  void MSQuantifications::registerOutRefs(const UInt64 & dp_ref, const UInt64 & object_ref)
  {
    out_data_processings_.insert(std::make_pair(dp_ref,object_ref));
  }  
  
  const UInt64 MSQuantifications::featureMapWhichRaw(const UInt64 & what) const
  {
    UInt64 ret(0);
    std::map< UInt64, UInt64 >::const_iterator it = featuremap_to_raw_.find(what);
    if (it != featuremap_to_raw_.end())
    {
      return it->second;
    }
    return ret;
    // TODO @mths : add consensus_to_features_
  }  
  
  const std::vector<UInt64> MSQuantifications::rawWhichAssays(const UInt64 & what) const
  {
    std::vector<UInt64> ret;
    std::map< UInt64, std::vector<UInt64> >::const_iterator it = raw_to_assays_.find(what);
    if (it != raw_to_assays_.end())
      return it->second;
    return ret;
    // TODO @mths : add consensus_to_features_
  }
  
  void MSQuantifications::consume_raw_file_groups(std::map<UInt64, std::set<ExperimentalSettings> > & rfgs)
  {
    raw_files_group_.swap(rfgs);
  }     
  
  void MSQuantifications::consumeDataProcessingList(std::vector<DataProcessing> & dps)
  {
    data_processings_.swap(dps);
  }
  
  void MSQuantifications::consumeAssays(std::map<UInt64, MSQuantifications::Assay> & asys)
  {
    assays_.swap(asys);
  }
  
  void MSQuantifications::consumeFeatureMap(FeatureMap<> & fm, UInt64 & rfref)
  {
    UInt64 fid = fm.getUniqueId();
    feature_maps_.insert(std::make_pair(fid, fm)); // TODO @mths : create empty and swap to avoid copy!
    featuremap_to_raw_.insert(std::make_pair(fid, rfref));
  }
  
  void MSQuantifications::simpleMerge(MSQuantifications& rhs)
  {
    // analysis_summary_ check if equal?!
    raw_files_group_.insert(rhs.raw_files_group_.begin(), rhs.raw_files_group_.end());
    source_files_.insert(source_files_.end(), rhs.source_files_.begin(), rhs.source_files_.end()); 
    data_processings_.insert(data_processings_.end(), rhs.data_processings_.begin(), rhs.data_processings_.end()); // TODO merge software?
    assays_.insert(rhs.assays_.begin(), rhs.assays_.end());

    feature_maps_.insert(rhs.feature_maps_.begin(), rhs.feature_maps_.end());
    consensus_maps_.insert(rhs.consensus_maps_.begin(), rhs.consensus_maps_.end());

    featuremap_to_raw_.insert(rhs.featuremap_to_raw_.begin(), rhs.featuremap_to_raw_.end());
    raw_to_assays_.insert(rhs.raw_to_assays_.begin(), rhs.raw_to_assays_.end());
    consensus_to_features_.insert(rhs.consensus_to_features_.begin(), rhs.consensus_to_features_.end());
    
    in_data_processings_.insert(rhs.in_data_processings_.begin(), rhs.in_data_processings_.end());
    out_data_processings_.insert(rhs.out_data_processings_.begin(), rhs.out_data_processings_.end()); 
  }



} //namespace OpenMS
