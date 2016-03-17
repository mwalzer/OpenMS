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
// $Maintainer: Nico Pfeifer $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <boost/math/special_functions/fpclassify.hpp>

#include <sstream>
#include <iostream>

using namespace std;

namespace OpenMS
{

  SpectrumIdentification::SpectrumIdentification() :
    MetaInfoInterface(),
    id_(),
    hits_(),
    significance_threshold_(0.0),
    score_type_(),
    higher_score_better_(true),
    base_name_(),
    mz_(std::numeric_limits<double>::quiet_NaN()),
    rt_(std::numeric_limits<double>::quiet_NaN())
  {
  }

  SpectrumIdentification::SpectrumIdentification(const SpectrumIdentification& rhs) :
    MetaInfoInterface(rhs),
    id_(rhs.id_),
    hits_(rhs.hits_),
    significance_threshold_(rhs.significance_threshold_),
    score_type_(rhs.score_type_),
    higher_score_better_(rhs.higher_score_better_),
    base_name_(rhs.base_name_),
    mz_(rhs.mz_),
    rt_(rhs.rt_)
  {
    setExperimentLabel( rhs.getExperimentLabel() );
  }

  SpectrumIdentification::~SpectrumIdentification()
  {
  }

  SpectrumIdentification& SpectrumIdentification::operator=(const SpectrumIdentification& rhs)
  {
    if (this == &rhs)
    {
      return *this;
    }

    MetaInfoInterface::operator=(rhs);
    id_ = rhs.id_;
    hits_ = rhs.hits_;
    significance_threshold_ = rhs.significance_threshold_;
    score_type_ = rhs.score_type_;
    higher_score_better_ = rhs.higher_score_better_;
    setExperimentLabel( rhs.getExperimentLabel() );
    base_name_ = rhs.base_name_;
    mz_ = rhs.mz_;
    rt_ = rhs.rt_;

    return *this;
  }

  // Equality operator
  bool SpectrumIdentification::operator==(const SpectrumIdentification& rhs) const
  {
    return MetaInfoInterface::operator==(rhs)
           && id_ == rhs.id_
           && hits_ == rhs.getHits()
           && significance_threshold_ == rhs.getSignificanceThreshold()
           && score_type_ == rhs.score_type_
           && higher_score_better_ == rhs.higher_score_better_
           && getExperimentLabel() == rhs.getExperimentLabel()
           && base_name_ == rhs.base_name_
           && (mz_ == rhs.mz_ || (!this->hasMZ() && !rhs.hasMZ())) // might be NaN, so comparing == will always be false
           && (rt_ == rhs.rt_ || (!this->hasRT() && !rhs.hasRT()));// might be NaN, so comparing == will always be false
  }

  // Inequality operator
  bool SpectrumIdentification::operator!=(const SpectrumIdentification& rhs) const
  {
    return !(*this == rhs);
  }

  double SpectrumIdentification::getRT() const
  {
    return rt_;
  }

  void SpectrumIdentification::setRT(double rt)
  {
    rt_ = rt;
  }

  bool SpectrumIdentification::hasRT() const
  {
    return !boost::math::isnan(rt_);
  }

  double SpectrumIdentification::getMZ() const
  {
    return mz_;
  }

  void SpectrumIdentification::setMZ(double mz)
  {
    mz_ = mz;
  }

  bool SpectrumIdentification::hasMZ() const
  {
    return !boost::math::isnan(mz_);
  }

  const std::vector<SpectrumMatch>& SpectrumIdentification::getHits() const
  {
    return hits_;
  }

  std::vector<SpectrumMatch>& SpectrumIdentification::getHits()
  {
    return hits_;
  }

  void SpectrumIdentification::insertHit(const SpectrumMatch& hit)
  {
    hits_.push_back(hit);
  }

  void SpectrumIdentification::setHits(const std::vector<SpectrumMatch>& hits)
  {
    hits_ = hits;
  }

  double SpectrumIdentification::getSignificanceThreshold() const
  {
    return significance_threshold_;
  }

  void SpectrumIdentification::setSignificanceThreshold(double value)
  {
    significance_threshold_ = value;
  }

  const String& SpectrumIdentification::getScoreType() const
  {
    return score_type_;
  }

  void SpectrumIdentification::setScoreType(const String& type)
  {
    score_type_ = type;
  }

  bool SpectrumIdentification::isHigherScoreBetter() const
  {
    return higher_score_better_;
  }

  void SpectrumIdentification::setHigherScoreBetter(bool value)
  {
    higher_score_better_ = value;
  }

  const String& SpectrumIdentification::getIdentifier() const
  {
    return id_;
  }

  void SpectrumIdentification::setIdentifier(const String& id)
  {
    id_ = id;
  }

  const String& SpectrumIdentification::getBaseName() const
  {
    return base_name_;
  }

  void SpectrumIdentification::setBaseName(const String& base_name)
  {
    base_name_ = base_name;
  }

  const String SpectrumIdentification::getExperimentLabel() const
  {
    // implement as meta value in order to reduce bloat of PeptideIdentification object
    //  -> this is mostly used for pepxml at the moment which allows each peptide id to belong to a different experiment
    if (metaValueExists("experiment_label"))
    {
      return getMetaValue("experiment_label").toString();
    }
    else
    {
      return "";
    }
  }

  void SpectrumIdentification::setExperimentLabel(const String& label)
  {
    // do not store empty label (default value)
    if (!label.empty())
    {
      setMetaValue("experiment_label", label);
    }
  }

  void SpectrumIdentification::assignRanks()
  {
    if (hits_.empty())
    {
      return;
    }
    UInt rank = 1;
    sort();
    vector<SpectrumMatch>::iterator lit = hits_.begin();
    double last_score = lit->getScore();
    while (lit != hits_.end())
    {
      if ((double)lit->getScore() != last_score)
      {
        ++rank;
        last_score = lit->getScore();
      }
      lit->setRank(rank);
      ++lit;
    }
  }

  void SpectrumIdentification::sort()
  {
    if (higher_score_better_)
    {
      std::stable_sort(hits_.begin(), hits_.end(), SpectrumMatch::ScoreMore());
    }
    else
    {
      std::stable_sort(hits_.begin(), hits_.end(), SpectrumMatch::ScoreLess());
    }
  }

  void SpectrumIdentification::sortByRank()
  {
    std::sort(hits_.begin(), hits_.end(), SpectrumMatch::RankLess());
  }

  bool SpectrumIdentification::empty() const
  {
    return id_ == ""
           && hits_.empty()
           && significance_threshold_ == 0.0
           && score_type_ == ""
           && higher_score_better_ == true
           && base_name_ == "";
  }

  std::vector<SpectrumMatch> SpectrumIdentification::getReferencingHits(const std::vector<SpectrumMatch>& hits, const std::set<String>& accession)
  {
    std::vector<SpectrumMatch> filtered;
    for (std::vector<SpectrumMatch>::const_iterator h_it = hits.begin(); h_it != hits.end(); ++h_it)
    {
      set<String> hit_accessions = h_it->extractProteinAccessions();
      set<String> intersect;
      set_intersection(hit_accessions.begin(), hit_accessions.end(), accession.begin(), accession.end(), std::inserter(intersect, intersect.begin()));
      if (!intersect.empty())
      {
        filtered.push_back(*h_it);
      }
    }
    return filtered;
  }

  /// re-implemented from MetaValueInterface as a precaution against deprecated usage of "RT" and "MZ" values
  const DataValue& SpectrumIdentification::getMetaValue(const String& name) const
  {
    if (name == "RT" || name == "MZ")
    { // this line should never the triggered. Set a breakpoint, find out who called getMetaValue() and replace with PeptideIdentification.getRT()/.getMZ() !!!!
      std::cerr << "\n\nUnsupported use of MetavalueInferface for 'RT' detected in " << __FILE__ << ":" << __LINE__ << ". Please notify the developers, so they can remove outdated code!\n\n";
      exit(1);
    }
    return MetaInfoInterface::getMetaValue(name);
  }

  /// re-implemented from MetaValueInterface as a precaution against deprecated usage of "RT" and "MZ" values
  void SpectrumIdentification::setMetaValue(const String& name, const DataValue& value)
  {
    if (name == "RT" || name == "MZ")
    { // this line should never the triggered. Set a breakpoint, find out who called getMetaValue() and replace with PeptideIdentification.getRT()/.getMZ() !!!!
      std::cerr << "\n\nUnsupported use of MetavalueInferface for 'RT' detected in " << __FILE__ << ":" << __LINE__ << ". Please notify the developers, so they can remove outdated code!\n\n";
      exit(1);
    }
    MetaInfoInterface::setMetaValue(name, value);
  }

} // namespace OpenMS
