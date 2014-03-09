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
// $Maintainer: Mathias Walzer $
// $Authors: Mathias Walzer$
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/HANDLERS/MzQuantMLHandler.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/METADATA/DataProcessing.h>
#include <OpenMS/DATASTRUCTURES/DateTime.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/Feature.h>
#include <set>
#include <vector>
#include <map>
#include <iostream>
#include <algorithm>

using namespace std;

namespace OpenMS
{
  namespace Internal
  {
    MzQuantMLHandler::MzQuantMLHandler(const MSQuantifications& msq, const String& filename, const String& version, const ProgressLogger& logger) :
      XMLHandler(filename, version),
      logger_(logger),
      msq_(0),
      cmsq_(&msq)
    {
      cv_.loadFromOBO("MS", File::find("/CV/psi-ms.obo")); //TODO unimod -> then automatise CVList writing
    }

    MzQuantMLHandler::MzQuantMLHandler(MSQuantifications& msq, /* FeatureMap& feature_map, */ const String& filename, const String& version, const ProgressLogger& logger) :
      XMLHandler(filename, version),
      logger_(logger),
      msq_(&msq),
      cmsq_(0)
    {
      cv_.loadFromOBO("MS", File::find("/CV/psi-ms.obo"));
    }

    MzQuantMLHandler::~MzQuantMLHandler()
    {
    }

    void MzQuantMLHandler::startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes)
    {
      tag_ = sm_.convert(qname);
      open_tags_.push_back(tag_);

      static set<String> to_ignore;
      if (to_ignore.empty())
      {
        to_ignore.insert("CvList"); // for now static set of obos.
        to_ignore.insert("Cv"); // for now static set of obos.
        to_ignore.insert("ProteinGroupList"); // for now no proteins or groups
        to_ignore.insert("ProteinList"); // .
        to_ignore.insert("Protein"); // .
        to_ignore.insert("StudyVariableList"); // We can't deal with these right now, but that is coming
        to_ignore.insert("StudyVariable"); // .
        to_ignore.insert("Assay_refs"); // .

        to_ignore.insert("FeatureList"); // we only need to see the features and datamatrices rows
        to_ignore.insert("AssayList"); // we only need to see the assays
        to_ignore.insert("DataProcessingList"); // we only need to see the DataProcessings
        to_ignore.insert("SoftwareList"); // we only need to see the Softwares
        to_ignore.insert("InputFiles"); // we only need to see the Files
        to_ignore.insert("Label"); // we only need to see the Modifications
        to_ignore.insert("DataType"); // we only need to see the Modifications
        to_ignore.insert("ColumnIndex"); // we only need to see the inside characters
        to_ignore.insert("DataMatrix"); // we only need to see the inside characters
      }

      if (to_ignore.find(tag_) != to_ignore.end())
      {
        return;
      }

      //determine parent tag
      String parent_tag;
      if (open_tags_.size() > 1)
      {
        parent_tag = *(open_tags_.end() - 2);
      }
      String parent_parent_tag;
      if (open_tags_.size() > 2)
      {
        parent_parent_tag = *(open_tags_.end() - 3);
      }

      static const XMLCh* s_value = xercesc::XMLString::transcode("value");
      static const XMLCh* s_type = xercesc::XMLString::transcode("type");
      static const XMLCh* s_name = xercesc::XMLString::transcode("name");
      static const XMLCh* s_unit_accession = xercesc::XMLString::transcode("unitAccession");
      static const XMLCh* s_cv_ref = xercesc::XMLString::transcode("cvRef");
      static const XMLCh* s_accession = xercesc::XMLString::transcode("accession");

      if (tag_ == "cvParam")
      {
        String value, unit_accession, cv_ref;
        optionalAttributeAsString_(value, attributes, s_value);
        optionalAttributeAsString_(unit_accession, attributes, s_unit_accession);
        optionalAttributeAsString_(cv_ref, attributes, s_cv_ref); //TODO
        handleCVParam_(parent_parent_tag, parent_tag, attributeAsString_(attributes, s_accession), attributeAsString_(attributes, s_name), value, attributes, cv_ref, unit_accession);
      }
      else if (tag_ == "MzQuantML")
      {
        // handle version and experiment type
      }
      else if (tag_ == "AnalysisSummary")
      {
        // handle version and experiment type
        // consists of cvParam only -> handleCvParam
      }
      else if (tag_ == "RawFilesGroup")
      {
        current_id_ = attributeAsString_(attributes, "id");
        std::set<ExperimentalSettings> exp_set;
        current_files_.insert(std::make_pair(current_id_, exp_set));
      }
      else if (tag_ == "RawFile")
      {
        ExperimentalSettings es;
        es.setLoadedFilePath(attributeAsString_(attributes, "location")+attributeAsString_(attributes, "name"));
        es.setUniqueId(parseUID_(attributeAsString_(attributes, "id")));
        //here would be the place to start looking for additional experimentalsettings readin
        current_files_[current_id_].insert(es);
      }

/*      else if (tag_ == "DataProcessing")
      {
        int order = asInt_(attributeAsString_(attributes, "order"));
        current_dp_ = std::make_pair(order, DataProcessing());
        current_pas_.clear();
        //~ order
        DataValue sw_ref(attributeAsString_(attributes, "software_ref"));
        current_dp_.second.setMetaValue("software_ref", sw_ref);
      }
      else if (tag_ == "ProcessingMethod")
      {
        //order gets implicity imposed by set<ProcessingAction> - so nothing to do here
      }
      else if (tag_ == "Software")
      {
        current_id_ = attributeAsString_(attributes, "id");
        current_sws_.insert(std::make_pair(current_id_, Software()));
        String vers = attributeAsString_(attributes, "version");
        current_sws_[current_id_].setVersion(vers);
      }
      else if (tag_ == "userParam")
      {
        String type = "";
        optionalAttributeAsString_(type, attributes, s_type);
        String value = "";
        optionalAttributeAsString_(value, attributes, s_value);
        handleUserParam_(parent_parent_tag, parent_tag, attributeAsString_(attributes, s_name), type, value);
      }
      else if (tag_ == "Assay")
      {
        current_assay_ = MSQuantifications::Assay();
        current_assay_.uid_ = attributeAsString_(attributes, "id");
        if (current_assay_.uid_.hasPrefix("a_"))
        {
          current_assay_.uid_ =  current_assay_.uid_.substr(2, current_assay_.uid_.size());
        }
        current_id_ = attributeAsString_(attributes, "rawFilesGroup_ref");
        current_assay_.raw_files_ = current_files_[current_id_];
        //TODO CVhandling
      }
      else if (tag_ == "Modification")
      {
        if (parent_tag == "Label")
        {
          String massdelta_string;
          optionalAttributeAsString_(massdelta_string, attributes, "massDelta");
          String residue;
          optionalAttributeAsString_(residue, attributes, "residues");
          if (massdelta_string != "145")
          {
            current_assay_.mods_.push_back(std::make_pair(residue, massdelta_string.toDouble()));
          }
          //TODO CVhandling
        }
        else
        {
          error(LOAD, "MzQuantMLHandler::startElement: Unhandable element found: '" + tag_ + "' in tag '" + parent_tag + "', ignoring.");
        }
      }
      else if (tag_ == "Ratio")
      {
        current_id_ = attributeAsString_(attributes, "id");
        String num = attributeAsString_(attributes, "numerator_ref");
        if (num.hasPrefix("a_"))
        {
          num = num.substr(2, num.size());
        }
        String den = attributeAsString_(attributes, "denominator_ref");
        if (den.hasPrefix("a_"))
        {
          den = den.substr(2, den.size());
        }
        ConsensusFeature::Ratio r;
        r.denominator_ref_ = den;
        r.numerator_ref_ = num;
        r_rtemp_.insert(std::make_pair(current_id_, r));
      }
      else if (tag_ == "PeptideConsensusList")
      {
        current_id_ = attributeAsString_(attributes, "id"); //needed in all PeptideConsensus elements
      }
      else if (tag_ == "PeptideConsensus")
      {
        ConsensusFeature current_cf;
        current_cf_id_ = attributeAsString_(attributes, "id");
        int c = attributeAsInt_(attributes, "charge");
        current_cf.setCharge(c);
        //TODO read searchDatabase map from inputfiles
        String searchDatabase_ref;
        if (optionalAttributeAsString_(searchDatabase_ref, attributes, "SearchDatabase_ref"))
        {
          current_cf.setMetaValue("SearchDatabase_ref", DataValue(searchDatabase_ref));
        }
        cm_cf_ids_.insert(std::make_pair(current_id_, current_cf_id_));
        cf_cf_obj_.insert(std::make_pair(current_cf_id_, current_cf));
      }
      else if (tag_ == "EvidenceRef")
      {
        //~ String searchDatabase_ref;
        //~ if (optionalAttributeAsString_(searchDatabase_ref,attributes,"SearchDatabase_ref") )
        //~ {
        //~ current_cf.setMetaValue("SearchDatabase_ref",DataValue(searchDatabase_ref));
        //~ }
        //~ String identificationFile_ref;
        //~ if (optionalAttributeAsString_(identificationFile_ref,attributes,"identificationFile_ref") )
        //~ {
        //~ current_cf.setMetaValue("identificationFile_ref",DataValue(identificationFile_ref)); //TODO add identificationFile_ref to PeptideIdentification
        //~ }
        //~ StringList id_refs;
        //~ if (optionalAttributeAsStringList_(id_refs,attributes,"id_refs") )
        //~ {
        //~ for (StringList::const_iterator it; it != id_refs.end(); ++it) //complete loop wont work! TODO add id_refs to PeptideIdentification
        //~ {
        //~ current_cf.setMetaValue("identificationFile_ref",DataValue(*it));
        //~ }
        //~ }

        String f_ref = attributeAsString_(attributes, "feature_ref"); // models which features will be included in this consensus feature - idependent from id(is optional)
        f_cf_ids_.insert(std::make_pair(f_ref, current_cf_id_));

        //~ StringList a_refs = attributeAsStringList_(attributes,"assay_refs"); // what to do with these??
        //~ for (StringList::const_iterator it = a_refs.begin(); it != a_refs.end(); ++it)
        //~ {
        //~ }
      }
      else if (tag_ == "Feature")
      {
        current_id_ = attributeAsString_(attributes, "id");
        DoubleReal rt = attributeAsDouble_(attributes, "rt");
        DoubleReal mz = attributeAsDouble_(attributes, "mz");
        FeatureHandle fh;
        fh.setRT(rt);
        fh.setMZ(mz);
        int c;
        if (optionalAttributeAsInt_(c, attributes, "charge"))
        {
          fh.setCharge(c);
        }
        f_f_obj_.insert(std::make_pair(current_id_, fh)); // map_index was lost!! TODO artificial ones produced in ConsensusMap assembly
      }
      else if (tag_ == "FeatureQuantLayer" || tag_ == "RatioQuantLayer" || tag_ == "MS2AssayQuantLayer")
      {
        //TODO Column attribute index!!!
        current_col_types_.clear();
      }
      else if (tag_ == "Column")
      {
        current_count_ = Size(attributeAsInt_(attributes, "index"));
      }
      else if (tag_ == "Row")
      {
        current_id_ = attributeAsString_(attributes, "object_ref");
        current_row_.clear();
      }
      else
      {
        error(LOAD, "MzQuantMLHandler::startElement: Unkown element found: '" + tag_ + "' in tag '" + parent_tag + "', ignoring.");
      } */        
    }

    void MzQuantMLHandler::characters(const XMLCh* const chars, const XMLSize_t /*length*/)
    {
/*       //if there is data between the element tags - !attention if element is derived from a xsd:list type, each list entry is a charecters call :(

      if (tag_ == "PeptideSequence")
      {
        AASequence p(sm_.convert(chars));
        PeptideHit ph = PeptideHit(0, 0, cf_cf_obj_[current_cf_id_].getCharge(), p);
        cf_cf_obj_[current_cf_id_].getPeptideIdentifications().back().insertHit(ph); // just moments before added
        return;
      }
      else if (tag_ == "Row")
      {
        String r = sm_.convert(chars);
        r.trim();
        if (!r.empty()) // always two notifications for a row, only the first one contains chars - dunno why
        {
          std::vector<String> splits;
          r.split(" ", splits);
          for (std::vector<String>::iterator it = splits.begin(); it != splits.end(); ++it)
          {
            current_row_.push_back(it->toDouble());
          }
        }
      }
      else if (tag_ == "ColumnIndex")
      {
        //overwrites current_col_types_ with the ratio_refs or the assay_refs
        String r = sm_.convert(chars);
        //clear must have happened earlyer in QuantLayer tag
        r.trim();
        if (!r.empty()) // always two notifications for a row, only the first one contains chars - dunno why
        {
          r.split(" ", current_col_types_);
        }
      }
      else
      {
        String transcoded_chars2 = sm_.convert(chars);
        transcoded_chars2.trim();
        if (transcoded_chars2 != "")
          warning(LOAD, "MzQuantMLHandler::characters: Unkown character section found: '" + tag_ + "', ignoring: " + transcoded_chars2);
      }
 */
    }
    void MzQuantMLHandler::endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname)
    {
/*       static set<String> to_ignore;
      if (to_ignore.empty())
      {
        to_ignore.insert("Cv");
      }

      tag_ = sm_.convert(qname);

      //determine parent tag
      String parent_tag;
      if (open_tags_.size() > 1)
      {
        parent_tag = *(open_tags_.end() - 2);
      }
      String parent_parent_tag;
      if (open_tags_.size() > 2)
      {
        parent_parent_tag = *(open_tags_.end() - 3);
      }

      //close current tag
      open_tags_.pop_back();

      if (to_ignore.find(tag_) != to_ignore.end())
      {
        return;
      }

      // no ProcessingMethod endElement action so each userParam under Dataprocessing will be one processingaction - no other way for core-lib compability yet
      if (tag_ == "DataProcessing")
      {
        current_dp_.second.setProcessingActions(current_pas_);
        current_orderedps_.insert(current_dp_);
        return;
      }
      else if (tag_ == "DataProcessingList")
      {
        std::vector<DataProcessing> dps;
        for (std::map<int, DataProcessing>::const_iterator it = current_orderedps_.begin(); it != current_orderedps_.end(); ++it)
        {
          dps.push_back(it->second);
        }
        for (std::vector<DataProcessing>::iterator it = dps.begin(); it != dps.end(); ++it)
        {
          if (it->metaValueExists("software_ref") && current_sws_.find(it->getMetaValue("software_ref")) != current_sws_.end())
          {
            it->setSoftware(current_sws_[it->getMetaValue("software_ref")]);
          }
        }
        msq_->setDataProcessingList(dps);
      }
      else if (tag_ == "Assay")
      {
        msq_->getAssays().push_back(current_assay_);
      }
      else if (tag_ == "ColumnDefinition")
      {
        //TODO check all current_col_types_[] are not empty
      }
      else if (tag_ == "Row")
      {
        if (current_col_types_.size() != current_row_.size())
        {
          warning(LOAD, String("Unknown/unmatching row content in Row element of '") + parent_tag + "'.");
          return;
        }

        if (parent_parent_tag == "RatioQuantLayer")
        {
          for (Size i = 0; i < current_row_.size(); ++i)
          {
            ConsensusFeature::Ratio r(r_rtemp_[current_col_types_[i]]);
            r.ratio_value_ = current_row_[i];
            cf_cf_obj_[current_id_].addRatio(r);
          }
        }

        if (parent_parent_tag == "MS2AssayQuantLayer")
        {
          ConsensusFeature ms2cf;
          ms2cf.setMZ(f_f_obj_[current_id_].getMZ());
          ms2cf.setRT(f_f_obj_[current_id_].getRT());
          cf_cf_obj_.insert(std::make_pair(current_id_, ms2cf));
          for (Size i = 0; i < current_row_.size(); ++i)
          {
            FeatureHandle fh;
            fh.setRT(f_f_obj_[current_id_].getMZ());
            fh.setMZ(f_f_obj_[current_id_].getRT());
            fh.setIntensity(current_row_[i]);
            fh.setUniqueId(i);
            cf_cf_obj_[current_id_].insert(fh);
          }
        }

        if (parent_parent_tag == "FeatureQuantLayer")
        {
          for (Size i = 0; i < current_row_.size(); ++i)
          {
            if (current_col_types_[i] == "MS:1001141")
            {
              f_f_obj_[current_id_].setIntensity(current_row_[i]);
            }
            else if (current_col_types_[i] == "width")
            {
              //TODO featurehandle have no width
            }
          }
        }
      }
      //~ assemble consensusmap MS2LABEL
      else if (tag_ == "MS2AssayQuantLayer")
      {
        ConsensusMap cm;
        for (std::map<String, ConsensusFeature>::iterator it = cf_cf_obj_.begin(); it != cf_cf_obj_.end(); ++it)
        {
          cm.push_back(it->second);
        }
        f_cf_ids_.clear();
        cm_cf_ids_.clear();
        msq_->addConsensusMap(cm);
      }
      //~ assemble consensusmap MS1LABEL
      else if (tag_ == "FeatureList") // TODO what if there are more than one FeatureQuantLayer?
      {
        //~ assemble consensus features
        for (std::map<String, String>::iterator it = f_cf_ids_.begin(); it != f_cf_ids_.end(); ++it)
        {
          cf_cf_obj_[it->second].insert(f_f_obj_[it->first]);
        }

        //~ assemble consensus maps
        ConsensusMap cm;
        std::multimap<String, String>::const_iterator last = cm_cf_ids_.begin();
        for (std::multimap<String, String>::const_iterator it = cm_cf_ids_.begin(); it != cm_cf_ids_.end(); ++it)
        {
          if (it->first != last->first)
          {
            msq_->addConsensusMap(cm);
            cm = ConsensusMap();
            last = it;
          }
          cm.push_back(cf_cf_obj_[it->second]);
        }
        if (!f_cf_ids_.empty()) //in case of MS2QuantLayer we do not need that and so after datamatrix f_cf_ids_ get cleared so we know here.
        {
          msq_->addConsensusMap(cm);
        }
      }
      else
        warning(LOAD, String("MzQuantMLHandler::endElement: Unkown element found: '" + tag_ + "', ignoring."));
 */
    }
    
    void MzQuantMLHandler::handleCVParam_(const String& parent_parent_tag, const String& parent_tag, const String& accession, const String& name, const String& value, const xercesc::Attributes& /* attributes */, const String& /* cv_ref */, const String& /* unit_accession */)
    {
       //Abort on unknown terms
      if (!cv_.exists(accession)) // TODO @mths: more permissive cv handling - do not directly abort
      {
        warning(LOAD, String("Unknown cvParam '") + accession + "' in tag '" + parent_tag + "'.");
        return;
      }
      else
      {
        const ControlledVocabulary::CVTerm& term = cv_.getTerm(accession);
        //obsolete CV terms
        if (term.obsolete)
        {
          warning(LOAD, String("Obsolete CV term '") + accession + " - " + cv_.getTerm(accession).name + "' used in tag '" + parent_tag + "'.");
        }
        //check if term name and parsed name match
        String parsed_name = name;
        parsed_name.trim();
        String correct_name = term.name;
        correct_name.trim();
        if (parsed_name != correct_name)
        {
          warning(LOAD, String("Name of CV term not correct: '") + term.id + " - " + parsed_name + "' should be '" + correct_name + "'");
        }
        if (term.obsolete)
        {
          warning(LOAD, String("Obsolete CV term '") + accession + " - " + cv_.getTerm(accession).name + "' used in tag '" + parent_tag + "'.");
        }
        //values used in wrong places and wrong value types
        if (value != "")
        {
          if (term.xref_type == ControlledVocabulary::CVTerm::NONE)
          {
            //Quality CV does not state value type :(
            if (!accession.hasPrefix("PATO:"))
            {
              warning(LOAD, String("The CV term '") + accession + " - " + cv_.getTerm(accession).name + "' used in tag '" + parent_tag + "' must not have a value. The value is '" + value + "'.");
            }
          }
          else
          {
            switch (term.xref_type)
            {
            //string value can be anything
            case ControlledVocabulary::CVTerm::XSD_STRING:
              break;

            //int value => try casting
            case ControlledVocabulary::CVTerm::XSD_INTEGER:
            case ControlledVocabulary::CVTerm::XSD_NEGATIVE_INTEGER:
            case ControlledVocabulary::CVTerm::XSD_POSITIVE_INTEGER:
            case ControlledVocabulary::CVTerm::XSD_NON_NEGATIVE_INTEGER:
            case ControlledVocabulary::CVTerm::XSD_NON_POSITIVE_INTEGER:
              try
              {
                value.toInt();
              }
              catch (Exception::ConversionError&)
              {
                warning(LOAD, String("The CV term '") + accession + " - " + cv_.getTerm(accession).name + "' used in tag '" + parent_tag + "' must have an integer value. The value is '" + value + "'.");
                return;
              }
              break;

            //double value => try casting
            case ControlledVocabulary::CVTerm::XSD_DECIMAL:
              try
              {
                value.toDouble();
              }
              catch (Exception::ConversionError&)
              {
                warning(LOAD, String("The CV term '") + accession + " - " + cv_.getTerm(accession).name + "' used in tag '" + parent_tag + "' must have a floating-point value. The value is '" + value + "'.");
                return;
              }
              break;

            //date string => try conversion
            case ControlledVocabulary::CVTerm::XSD_DATE:
              try
              {
                DateTime tmp;
                tmp.set(value);
              }
              catch (Exception::ParseError&)
              {
                warning(LOAD, String("The CV term '") + accession + " - " + cv_.getTerm(accession).name + "' used in tag '" + parent_tag + "' must be a valid date. The value is '" + value + "'.");
                return;
              }
              break;

            default:
              warning(LOAD, String("The CV term '") + accession + " - " + cv_.getTerm(accession).name + "' used in tag '" + parent_tag + "' has the unknown value type '" + ControlledVocabulary::CVTerm::getXRefTypeName(term.xref_type) + "'.");
              break;
            }
          }
        }
        //no value, although there should be a numerical value
        else if (term.xref_type != ControlledVocabulary::CVTerm::NONE && term.xref_type != ControlledVocabulary::CVTerm::XSD_STRING)
        {
          warning(LOAD, String("The CV term '") + accession + " - " + cv_.getTerm(accession).name + "' used in tag '" + parent_tag + "' should have a numerical value. The value is '" + value + "'.");
          return;
        }
      }

      if (parent_tag == "AnalysisSummary")
      {
        if (accession == "MS:1001834" or accession == "MS:1002019"  or accession == "MS:1002020")
        {
          MSQuantifications::QUANT_TYPES quant_type = MSQuantifications::LABELFREE;
          msq_->setAnalysisSummaryQuantType(quant_type); 
        }
        // TODO @anyone : register more elaborate analysis summary with cv param stack and extension of analysissummary 
        //~ cvp_stack_.push_back(cv_.getTerm(accession)); //set value!
      }    
/*      if (parent_tag == "DataType" && parent_parent_tag == "Column")
      {
        if (current_count_ >= current_col_types_.size())
        {
          current_col_types_.resize(current_count_ + 1, "");
        }
        current_col_types_[current_count_] = accession; //TODO real cv handling here (i.e. translate name into decision string for the "row-loop")
      }
      else if (parent_parent_tag == "Label")
      {
        //TODO
        if (accession == "MOD:01522")
          current_assay_.mods_.push_back(std::make_pair<String, DoubleReal>("114", DoubleReal(114)));
        else if (accession == "MOD:01523")
          current_assay_.mods_.push_back(std::make_pair<String, DoubleReal>("115", DoubleReal(115)));
        else if (accession == "MOD:01524")
          current_assay_.mods_.push_back(std::make_pair<String, DoubleReal>("116", DoubleReal(116)));
        else if (accession == "MOD:01525")
          current_assay_.mods_.push_back(std::make_pair<String, DoubleReal>("117", DoubleReal(117)));

      } */
      else
        warning(LOAD, String("Unhandled cvParam '") + name + "' in tag '" + parent_tag + "'.");

    }
    
    void MzQuantMLHandler::handleUserParam_(const String& parent_parent_tag, const String& parent_tag, const String& name, const String& type, const String& value)
    {
/*       //create a DataValue that contains the data in the right type
      DataValue data_value;
      //float type
      if (type == "xsd:double" || type == "xsd:float")
      {
        data_value = DataValue(value.toDouble());
      }
      //integer type
      else if (type == "xsd:byte" || type == "xsd:decimal" || type == "xsd:int" || type == "xsd:integer" || type == "xsd:long" || type == "xsd:negativeInteger" || type == "xsd:nonNegativeInteger" || type == "xsd:nonPositiveInteger" || type == "xsd:positiveInteger" || type == "xsd:short" || type == "xsd:unsignedByte" || type == "xsd:unsignedInt" || type == "xsd:unsignedLong" || type == "xsd:unsignedShort")
      {
        data_value = DataValue(value.toInt());
      }
      //everything else is treated as a string
      else
      {
        data_value = DataValue(value);
      }

      if (parent_parent_tag == "")
      {
        //~ TODO: dummy
        warning(LOAD, String("The user param '") + name + "' used in tag '" + parent_tag + "' has no valid grand parent.'");
      }
      //find the right MetaInfoInterface
      if (parent_tag == "ProcessingMethod")
      {
        //~ value is softwarename - will get handled elsewhere
        int x = std::distance(DataProcessing::NamesOfProcessingAction, std::find(DataProcessing::NamesOfProcessingAction, DataProcessing::NamesOfProcessingAction + DataProcessing::SIZE_OF_PROCESSINGACTION, name));
        DataProcessing::ProcessingAction a = static_cast<DataProcessing::ProcessingAction>(x); // ugly and depends on NamesOfProcessingAction^=ProcessingAction-definitions - see TODO rewrite DataProcessing!
        current_pas_.insert(a);
      }
      else if (parent_tag == "Software")
      {
        if (value == "")
        {
          current_sws_[current_id_].setName(name);
        }
        else
        {
          current_sws_[current_id_].setMetaValue(name, data_value);
        }
      }
      else if (parent_tag == "AnalysisSummary")
      {
        if (name == "QuantType")
        {
          const std::string* match = std::find(MSQuantifications::NamesOfQuantTypes, MSQuantifications::NamesOfQuantTypes + MSQuantifications::SIZE_OF_QUANT_TYPES, value);
          int i = (MSQuantifications::NamesOfQuantTypes + MSQuantifications::SIZE_OF_QUANT_TYPES == match) ? -1 : std::distance(MSQuantifications::NamesOfQuantTypes, match);
          MSQuantifications::QUANT_TYPES quant_type = static_cast<MSQuantifications::QUANT_TYPES>(i); //this is so not nice and soooo unsafe why enum in the first place?!
          msq_->setAnalysisSummaryQuantType(quant_type);
        }
        else
        {
          msq_->getAnalysisSummary().user_params_.setValue(name, data_value);
        }
      }
      else if (parent_tag == "RatioCalculation")
      {
        r_rtemp_[current_id_].description_.push_back(name);
      }
      else if (parent_tag == "Feature")
      {
        if (name == "feature_index")
        {
          f_f_obj_[current_id_].setUniqueId(UInt64(value.toInt()));
        }
        else if (name == "map_index")
        {
          f_f_obj_[current_id_].setMapIndex(UInt64(value.toInt()));
        }
      }
      else
        warning(LOAD, String("Unhandled userParam '") + name + "' in tag '" + parent_tag + "'.");
 */
    }
    
    void MzQuantMLHandler::writeTo(std::ostream& os)
    {
      //~ TODO logger_.startProgress(0,exp.size(),"storing mzQuantML file");
      String line; //everyone walk the line!!!

      //---header---
      os << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n";
      os << "<MzQuantML xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\"http://psidev.info/psi/pi/mzQuantML/1.0.0-rc3 ../../schema/mzQuantML_1_0_0-rc3.xsd\" xmlns=\"http://psidev.info/psi/pi/mzQuantML/1.0.0-rc3\" id=\"OpenMS_" << String(UniqueIdGenerator::getUniqueId()) << "\" version=\"1.0.0\"" << " creationDate=\"" << DateTime::now().get() << "\">\n";

      //---CVList---
      os << "\t<CvList>\n";
      os << " \t\t<Cv id=\"PSI-MS\" fullName=\"Proteomics Standards Initiative Mass Spectrometry Vocabularies\"  uri=\"http://psidev.cvs.sourceforge.net/viewvc/*checkout*/psidev/psi/psi-ms/mzML/controlledVocabulary/psi-ms.obo\" version=\"3.41.0\"/>\n";
      os << "\t\t<Cv id=\"PSI-MOD\" fullName=\"Proteomics Standards Initiative Protein Modifications Vocabularies\" uri=\"http://psidev.cvs.sourceforge.net/psidev/psi/mod/data/PSI-MOD.obo\" version=\"1.2\"/>\n";
      os << "\t\t<Cv id=\"UO\" fullName=\"Unit Ontology\" uri=\"http://obo.cvs.sourceforge.net/*checkout*/obo/obo/ontology/phenotype/unit.obo\"/>\n";
      os << "\t</CvList>\n";

      //---AnalysisSummary---
      os << "\t<AnalysisSummary>\n";
      cmsq_->getAnalysisSummary().quant_type_;
      //~ os << "\t\t<userParam name=\"QuantType\" value=\"";
      //~ os << String(MSQuantifications::NamesOfQuantTypes[cmsq_->getAnalysisSummary().quant_type_]);
      switch (cmsq_->getAnalysisSummary().quant_type_)
      {
      case 0:
        os << "\t\t<cvParam cvRef=\"PSI-MS\" accession=\"MS:1002018\" name=\"MS1 label-based analysis\"/>\n";
        os << "\t\t<cvParam cvRef=\"PSI-MS\" accession=\"MS:1001837\" name=\"SILAC quantitation analysis\"/>\n";
        os << "\t\t<cvParam cvRef=\"PSI-MS\" accession=\"MS:1002001\" name=\"MS1 label-based raw feature quantitation\" value=\"true\"/>\n";
        os << "\t\t<cvParam cvRef=\"PSI-MS\" accession=\"MS:1002002\" name=\"MS1 label-based peptide level quantitation\" value=\"true\"/>\n";
        os << "\t\t<cvParam cvRef=\"PSI-MS\" accession=\"MS:1002003\" name=\"MS1 label-based protein level quantitation\" value=\"false\"/>\n";
        os << "\t\t<cvParam cvRef=\"PSI-MS\" accession=\"MS:1002004\" name=\"MS1 label-based proteingroup level quantitation\" value=\"false\"/>\n";
        break;

      case 1:
        os << "\t\t<cvParam cvRef=\"PSI-MS\" accession=\"MS:1002023\" name=\"MS2 tag-based analysis\"/>\n";
        os << "\t\t<cvParam cvRef=\"PSI-MS\" accession=\"MS:1002024\" name=\"MS2 tag-based analysis feature level quantitation\" value=\"true\"/>\n";
        os << "\t\t<cvParam cvRef=\"PSI-MS\" accession=\"MS:1002025\" name=\"MS2 tag-based peptide level quantitation\" value=\"true\"/>\n";
        os << "\t\t<cvParam cvRef=\"PSI-MS\" accession=\"MS:1002026\" name=\"MS2 tag-based analysis protein level quantitation\" value=\"false\"/>\n";
        os << "\t\t<cvParam cvRef=\"PSI-MS\" accession=\"MS:1002027\" name=\"MS2 tag-based analysis protein group level quantitation\" value=\"false\"/>\n";
        break;

      case 2:
        os << "\t\t<cvParam accession=\"MS:1001834\" cvRef=\"PSI-MS\" name=\"LC-MS label-free quantitation analysis\"/>\n";
        os << "\t\t<cvParam accession=\"MS:1002019\" cvRef=\"PSI-MS\" value=\"false\" name=\"label-free raw feature quantitation\"/>\n";
        os << "\t\t<cvParam accession=\"MS:1002020\" cvRef=\"PSI-MS\" value=\"true\" name=\"label-free peptide level quantitation\"/>\n";
        //~ TODO
        //~ os << "\t\t<cvParam accession=\"MS:1002021\" cvRef=\"PSI-MS\" value=\"true\" name=\"label-free protein level quantitation\"/>\n";
        //~ os << "\t\t<cvParam accession=\"MS:1002022\" cvRef=\"PSI-MS\" value=\"false\" name=\"label-free proteingroup level quantitation\"/>\n";
        break;

      case 3:
        break; //why SIZE_OF_QUANT_TYPES anyway?
      }
      //~ writeUserParam_(dataprocessinglist_tag, cmsq_->getAnalysisSummary().getUserParams(), UInt(2));
      //~ writeCVParams_(dataprocessinglist_tag, (cmsq_->getAnalysisSummary().getCVTerms(), UInt(2));
      os << "\t</AnalysisSummary>\n";
      //---AnalysisSummary---
      
      //---InputFiles---
      os << "\t<InputFiles>\n"; 
      
      // TODO @mths : add raw files to mzML 
      for (std::map< UInt64, std::set<ExperimentalSettings> >::const_iterator rfgit = cmsq_->getRawFiles().begin(); rfgit != cmsq_->getRawFiles().end(); ++rfgit)
      {
        os << "\t\t<RawFilesGroup id=\"" << "f_" << String(rfgit->first) << "\">\n";
        for (std::set<ExperimentalSettings>::const_iterator sit=rfgit->second.begin(); sit!=rfgit->second.end(); ++sit)
        {
          os << "\t\t\t<RawFile location=\"" << sit->getSourceFiles().front().getPathToFile() << "\" name=\"" << sit->getSourceFiles().front().getNameOfFile() << "\" id=\"" << "f_" << String(sit->getUniqueId()) <<  "\"/>\n";
        }
        os << "\t\t</RawFilesGroup>\n";
      }
      // TODO @mths : add <IdentificationFiles> & <SearchDatabase>
      os << "\t</InputFiles>\n"; 
      //---InputFiles---

      
      //---Software & DataProcessing---
      String softwarelist_tag;
      softwarelist_tag += "\t<SoftwareList>\n";

      String dataprocessinglist_tag;
      dataprocessinglist_tag += "\t<DataProcessingList>\n";
      // TODO Software DefaultTag for each file: OpenMS
      Size order_d = 0;

      String idfile_tag, idfile_ref, searchdb_ref;

      std::vector<DataProcessing> pl = cmsq_->getDataProcessingList();
      for (std::vector<DataProcessing>::const_iterator dit = pl.begin(); dit != pl.end(); ++dit)
      {
        String sw_ref;
        sw_ref = "sw_" + String(UniqueIdGenerator::getUniqueId());
        softwarelist_tag += "\t\t<Software id=\"" +  sw_ref + "\" version=\"" + String(dit->getSoftware().getVersion()) + "\">\n";
        writeCVParams_(softwarelist_tag, dit->getSoftware().getCVTerms(), UInt(3)); 
        
        // TODO @mths : refactor this block to a private function
        if (dit->getSoftware().getCVTerms().empty())
        {
          softwarelist_tag += "\t\t\t<userParam name=\"" + String(dit->getSoftware().getName()) + "\"/>\n";
        }
        if (dit->getSoftware().getName() == "SILACAnalyzer")
        {
          softwarelist_tag += "\t\t\t<cvParam cvRef=\"PSI-MS\" accession=\"MS:1001831\" name=\"SILACAnalyzer\"/>\n";
        }
        if (dit->getSoftware().getName() == "ITRAQAnalyzer")
        {
          softwarelist_tag += "\t\t\t<cvParam cvRef=\"PSI-MS\" accession=\"MS:1001831\" name=\"ITRAQAnalyzer\"/>\n";
        }
        softwarelist_tag += "\t\t</Software>\n";
        ++order_d;
        
        
        dataprocessinglist_tag += "\t\t<DataProcessing id=\"dp_" + String(UniqueIdGenerator::getUniqueId()) + "\" software_ref=\"" + sw_ref + "\" order=\"" + String(order_d) + "\">\n";
        Size order_c = 0;
        
        std::vector<UInt64> dp_in = cmsq_->getDataProcessingInRefs(dit->getUniqueId());
        if (!dp_in.empty())
        {
          String ior;
          for (std::vector<UInt64>::iterator init = dp_in.begin(); init != dp_in.end(); ++init)
          {
            ior += "f_";
            ior += String(*init);
            dataprocessinglist_tag += " ";
          }
          dataprocessinglist_tag += "\t\t\t<InputObjects_ref>" + ior.trim() + "\t\t\t</InputObjects_ref>\n";
        }
        UInt64 dp_out = cmsq_->getDataProcessingOutRefs(dit->getUniqueId());
        if (dp_out != 0)
        {
          dataprocessinglist_tag += "\t\t\t<OutputObjects_ref>f_" + String(dp_out) + "\t\t\t</OutputObjects_ref>\n";
        }
        
        for (std::set<DataProcessing::ProcessingAction>::const_iterator pit = dit->getProcessingActions().begin(); pit != dit->getProcessingActions().end(); ++pit)
        {
          ++order_c;
          dataprocessinglist_tag += "\t\t\t<ProcessingMethod order=\"" + String(order_c) + "\">\n";
          dataprocessinglist_tag += "\t\t\t\t<userParam name=\"" + String(DataProcessing::NamesOfProcessingAction[*pit]) + "\" value=\"" + String(dit->getSoftware().getName()) + "\" />\n";
          dataprocessinglist_tag += "\t\t\t</ProcessingMethod>\n";
        }
        dataprocessinglist_tag += "\t\t</DataProcessing>\n";
      }
      dataprocessinglist_tag += "\t</DataProcessingList>\n";

      softwarelist_tag += "\t</SoftwareList>\n";
      os << softwarelist_tag << dataprocessinglist_tag;
      //---Software & DataProcessing---

      // TODO @mths: add BibliographicReference
      
      // ---Assay --- 
      os << "\t<AssayList>\n";
      
      std::map< UInt64, MSQuantifications::Assay > as_rfg = cmsq_->getAssays();
      for (std::map< UInt64, MSQuantifications::Assay >::iterator agit = as_rfg.begin(); agit!=as_rfg.end(); ++agit)
      {
        os << "\t\t<Assay rawFilesGroup_ref=\"" << "f_" << String(agit->second.rfg_ref_) << "\" id=\"" << "f_" << String(agit->first) << "\">\n";
        if (!agit->second.mods_.empty())
        {
          // TODO @mths : add the labels CVs
        }
        else
        {
          os << "\t\t\t<Label>\n";
          os << "\t\t\t\t<Modification>\n";
          os << "\t\t\t\t\t<cvParam accession=\"MS:1002038\" cvRef=\"PSI-MS\" name=\"unlabeled sample\"/>\n";
          os << "\t\t\t\t</Modification>\n";
          os << "\t\t\t</Label>\n";
        }
        os << "\t</Assay>\n"; 
      }
      os << "\t</AssayList>\n";
      // ---Assay ---
      
      // ---StudyVariables---
      
      // TODO @mths: add Tech. reps. : <cvParam accession="MS:1001808" cvRef="PSI-MS" name="technical replicate" value=""/>
      os << "\t<StudyVariableList>\n";
      for (std::map< UInt64, MSQuantifications::Assay >::iterator agit = as_rfg.begin(); agit!=as_rfg.end(); ++agit)
      {
        os << "\t\t<StudyVariable name=\"" << "sv_" << String(agit->first) << "\" id=\"" << "sv_" << String(agit->first) << "\">\n";
        os << "\t\t\t<Assay_refs>" << "f_" << String(agit->first) << "</Assay_refs>\n";
        os << "\t\t</StudyVariable>\n";
      }
      os << "\t</StudyVariableList>\n";
      // ---StudyVariables---
      
      // ---PeptideConsensus---
      
      // TODO @mths: add 
      
      // ---PeptideConsensus---
      
      // ---Ratios tag---
      String ratio_xml;
      switch (cmsq_->getAnalysisSummary().quant_type_)
      {
      case 0:
/*         //~ register ratio elements in numden_r_ids_ and r_r_obj_
        for (std::vector<ConsensusMap>::const_iterator mit = cmsq_->getConsensusMaps().begin(); mit != cmsq_->getConsensusMaps().end(); ++mit)
        {
          //~ std::vector< std::vector<UInt64> > cmid;
          for (ConsensusMap::const_iterator cit = mit->begin(); cit != mit->end(); ++cit)
          {
            std::vector<ConsensusFeature::Ratio> rv = cit->getRatios();
            //~ for (std::vector<ConsensusFeature::Ratio>::const_iterator rit = cit->getRatios().begin(); rit != cit->getRatios().end(); ++rit)
            for (Size i = 0; i < rv.size(); ++i)
            {
              ConsensusFeature::Ratio robj(rv[i]);
              //~ String rd = rit->numerator_ref_ + rit->denominator_ref_; // add ratiocalculation params too?
              String rd = robj.numerator_ref_ + robj.denominator_ref_; // add ratiocalculation params too?
              String tid = String(UniqueIdGenerator::getUniqueId());
              numden_r_ids_.insert(std::make_pair(rd, tid));

              //~ ConsensusFeature::Ratio robj(*rit); this segfaults!!! why???? oh, why?!?!?!?!
              r_r_obj_.insert(std::make_pair(tid, robj));
            }
          }
        }

        ratio_xml += "\t<RatioList>\n";
        for (std::map<String, String>::const_iterator rit = numden_r_ids_.begin(); rit != numden_r_ids_.end(); ++rit)
        {
          ratio_xml += "\t\t<Ratio id=\"r_" + String(rit->second) + "\" numerator_ref=\"a_" + String(r_r_obj_[rit->second].numerator_ref_) + "\" denominator_ref=\"a_" + String(r_r_obj_[rit->second].denominator_ref_) + "\" >\n";
          ratio_xml += "\t\t\t<RatioCalculation>\n";
          for (std::vector<String>::const_iterator dit = r_r_obj_[rit->second].description_.begin(); dit != r_r_obj_[rit->second].description_.end(); ++dit)
          {
            ratio_xml += "\t\t\t\t<userParam name=\"" + String(*dit) + "\"/>\n";
          }
          ratio_xml += "\t\t\t\t<cvParam cvRef=\"PSI-MS\" accession=\"MS:1001848\" name=\"simple ratio of two values\"/>\n";
          ratio_xml += "\t\t\t</RatioCalculation>\n";
          ratio_xml += "\t\t\t<NumeratorDataType>\n\t\t\t\t<cvParam accession=\"MS:1001847\" cvRef=\"PSI-MS\" name=\"reporter ion intensity\"/>\n\t\t\t</NumeratorDataType>\n";
          ratio_xml += "\t\t\t<DenominatorDataType>\n\t\t\t\t<cvParam accession=\"MS:1001847\" cvRef=\"PSI-MS\" name=\"reporter ion intensity\"/>\n\t\t\t</DenominatorDataType>\n";
          ratio_xml += "\t\t</Ratio>\n";
        }
        ratio_xml += "\t</RatioList>\n";
 */        break;

      case 1:
        break; //TODO for SILACAnalyzer to produce some ratios

      case 2:
        break; //no tool yet

      case 3:
        break; //why SIZE_OF_QUANT_TYPES anyway?
      }
      // ---Ratios tag---

      //--------------------------------------------------------------------------------------------
      // Proteins and Proteingroups
      //--------------------------------------------------------------------------------------------
      // TODO - omitted as there are no ids yet
      

      // ---Features and QuantLayers---
      String feature_tag;
      for (std::map< UInt64, FeatureMap<> >::const_iterator fmit = cmsq_->getFeatureMaps().begin(); fmit != cmsq_->getFeatureMaps().end(); ++fmit)
      {
        writeFeatureMap_(feature_tag,fmit->second,fmit->first,1);
      }
      os << feature_tag; 
      os << "</MzQuantML>\n";
    }

    void MzQuantMLHandler::writeCVParams_(String& s, const Map<String, std::vector<CVTerm> >& cvl, UInt indent)
    {
      String inden((size_t)indent, '\t');
      for (std::map<String, std::vector<CVTerm> >::const_iterator jt = cvl.begin(); jt != cvl.end(); ++jt)
      {
        for (std::vector<CVTerm>::const_iterator kt =  (*jt).second.begin(); kt !=  (*jt).second.end(); ++kt)
        {
          s += inden;
          s += "<cvParam cvRef=\"" + kt->getCVIdentifierRef() + "\" accession=\"" + (*jt).first + "\" name=\"" + kt->getName();
          if (kt->hasValue())
          {
            s += "\" value=\"" + kt->getValue().toString() + "\"/>\n"; // value is OpenMS::DataValue
          }
          else
          {
            s +=     "\"/>\n";
          }
        }
      }
    }

    void MzQuantMLHandler::writeUserParams_(std::ostream& os, const MetaInfoInterface& meta, UInt indent)
    {
      String h;
      writeUserParams_(h, meta, indent);
      os << h;
    }

    void MzQuantMLHandler::writeUserParams_(String& s, const MetaInfoInterface& meta, UInt indent)
    {
      if (meta.isMetaEmpty())
      {
        return;
      }
      std::vector<String> keys;
      meta.getKeys(keys);

      for (Size i = 0; i != keys.size(); ++i)
      {
        s += String(indent, '\t') + "<userParam name=\"" + keys[i] + "\" unitName=\"";

        DataValue d = meta.getMetaValue(keys[i]);
        //determine type
        if (d.valueType() == DataValue::INT_VALUE)
        {
          s += "xsd:integer";
        }
        else if (d.valueType() == DataValue::DOUBLE_VALUE)
        {
          s += "xsd:double";
        }
        else //string or lists are converted to string
        {
          s += "xsd:string";
        }
        s += "\" value=\"" + (String)(d) + "\"/>" + "\n";
      }
    }

    UInt64 MzQuantMLHandler::parseUID_(String xml_id_string);
    {
      std::vector< String > substrings;
      xml_id_string.split("_", substrings);
      try
      {
        substrings[1].toInt();
      }
      catch (Exception::ConversionError&)
      {
        warning(LOAD, String("xml uid could not be parsed");
        return;
      }
    }
    
    void MzQuantMLHandler::writeFeatureMap_(String& feature_tag, const FeatureMap<> & fm, const UInt64 & ar, UInt indentation_level)
    {
      feature_tag = String(indentation_level, '\t') + "<FeatureList id=\"fl_" + String(ar) + "\" rawFilesGroup_ref=\"f_" + String(ar) + " \">\n";
      String fql = String(indentation_level, '\t') + "\t<FeatureQuantLayer id=\"fql_" + String(UniqueIdGenerator::getUniqueId()) + "\">\n";
      fql += String(indentation_level, '\t') + "\t\t<ColumnDefinition>\n";
      fql += String(indentation_level, '\t') + "\t\t\t<Column index=\"0\">\n";
      fql += String(indentation_level, '\t') + "\t\t\t\t<DataType>\n";
      fql += String(indentation_level, '\t') + "\t\t\t\t\t<cvParam cvRef=\"PSI-MS\" accession=\"MS:1001141\" name=\"intensity of precursor ion\"/>\n";
      fql += String(indentation_level, '\t') + "\t\t\t\t</DataType>\n";
      fql += String(indentation_level, '\t') + "\t\t\t</Column>	\n";
      fql += String(indentation_level, '\t') + "\t\t\t<Column index=\"1\">\n";
      fql += String(indentation_level, '\t') + "\t\t\t\t<DataType>\n";
      fql += String(indentation_level, '\t') + "\t\t\t\t\t<cvParam cvRef=\"PSI-MS\" accession=\"PATO:0000001\" name=\"quality\"/>\n";
      fql += String(indentation_level, '\t') + "\t\t\t\t</DataType>\n";
      fql += String(indentation_level, '\t') + "\t\t\t</Column>	\n";
      fql += String(indentation_level, '\t') + "\t\t\t<Column index=\"2\">\n";
      fql += String(indentation_level, '\t') + "\t\t\t\t<DataType>\n";
      fql += String(indentation_level, '\t') + "\t\t\t\t\t<cvParam cvRef=\"PSI-MS\" accession=\"MS:1000086\" name=\"full width at half-maximum\"/>\n";
      fql += String(indentation_level, '\t') + "\t\t\t\t</DataType>\n";
      fql += String(indentation_level, '\t') + "\t\t\t</Column>	\n";
      fql += String(indentation_level, '\t') + "\t\t</ColumnDefinition>\n";
      fql += String(indentation_level, '\t') + "\t\t<DataMatrix>\n";

      for (FeatureMap<>::const_iterator fit = fm.begin(); fit != fm.end(); ++fit)
      {
        feature_tag += "\t<Feature charge=\"" + String(fit->getCharge()) + "\" mz=\"" + String(fit->getMZ()) + "\" rt=\"" + String(fit->getRT()) + "\" id=\"ft_" + String(fit->getUniqueId()) + "\">\n";
        feature_tag += "\t\t<MassTrace>" + fit->getConvexHull().getBoundingBox().toString() + "</MassTrace>\n"; // alt: ANALYSIS/TARGETED/OfflinePrecursorIonSelection.h:444: Size lower_rt = features[feat].getConvexHull().getBoundingBox().minX();
        feature_tag += "\t</Feature>\n";
        
        fql += "\t\t<Row object_ref=\"ft_" + String(fit->getUniqueId()) + "\">" + String(fit->getIntensity()) + " " + String(fit->getOverallQuality()) + " " + String(fit->getWidth()) + "</Row>\n";
      }
      
      feature_tag += "\t<cvParam accession=\"MS:1001826\" cvRef=\"PSI-MS\" name=\"mass trace reporting: rectangles\"/>\n";
      
      fql += String(indentation_level, '\t') + "\t\t</DataMatrix>\n";
      fql += String(indentation_level, '\t') + "\t</FeatureQuantLayer>\n";
      
      feature_tag += fql;
      feature_tag += String(indentation_level, '\t') + "</FeatureList>\n";
    }

    void MzQuantMLHandler::writeConsensusMap_(String& consensus_tag, const ConsensusMap & cm, UInt indentation_level)
    {
      consensus_tag = String(indentation_level, '\t') + "<PeptideConsensusList id=\"cl_" + String(cm.getUniqueId()) + "\" finalResult=\"false\">\n";
      String fql = String(indentation_level, '\t') + "\t<StudyVariableQuantLayer id=\"fql_" + String(UniqueIdGenerator::getUniqueId()) + "\">\n";
      fql += String(indentation_level, '\t') + "\t\t<ColumnDefinition>\n";
      fql += String(indentation_level, '\t') + "\t\t\t<Column index=\"0\">\n";
      fql += String(indentation_level, '\t') + "\t\t\t\t<DataType>\n";
      fql += String(indentation_level, '\t') + "\t\t\t\t\t<cvParam cvRef=\"PSI-MS\" accession=\"MS:1000040\" name=\"m/z\"/>\n";
      fql += String(indentation_level, '\t') + "\t\t\t\t</DataType>\n";
      fql += String(indentation_level, '\t') + "\t\t\t</Column>	\n";
      fql += String(indentation_level, '\t') + "\t\t\t<Column index=\"1\">\n";
      fql += String(indentation_level, '\t') + "\t\t\t\t<DataType>\n";
      fql += String(indentation_level, '\t') + "\t\t\t\t\t<cvParam cvRef=\"PSI-MS\" accession=\"MS:1000894\" name=\"retention time\"/>\n";
      fql += String(indentation_level, '\t') + "\t\t\t\t</DataType>\n";
      fql += String(indentation_level, '\t') + "\t\t\t</Column>	\n";
      fql += String(indentation_level, '\t') + "\t\t\t<Column index=\"2\">\n";
      fql += String(indentation_level, '\t') + "\t\t\t\t<DataType>\n";
      fql += String(indentation_level, '\t') + "\t\t\t\t\t<cvParam cvRef=\"PSI-MS\" accession=\"MS:1001141\" name=\"intensity of precursor ion\"/>\n";
      fql += String(indentation_level, '\t') + "\t\t\t\t</DataType>\n";
      fql += String(indentation_level, '\t') + "\t\t\t</Column>	\n";
      fql += String(indentation_level, '\t') + "\t\t\t<Column index=\"3\">\n";
      fql += String(indentation_level, '\t') + "\t\t\t\t<DataType>\n";
      fql += String(indentation_level, '\t') + "\t\t\t\t\t<cvParam cvRef=\"PSI-MS\" accession=\"PATO:0000001\" name=\"quality\"/>\n";
      fql += String(indentation_level, '\t') + "\t\t\t\t</DataType>\n";
      fql += String(indentation_level, '\t') + "\t\t\t</Column>	\n";
      fql += String(indentation_level, '\t') + "\t\t\t<Column index=\"4\">\n";
      fql += String(indentation_level, '\t') + "\t\t\t\t<DataType>\n";
      fql += String(indentation_level, '\t') + "\t\t\t\t\t<cvParam cvRef=\"PSI-MS\" accession=\"MS:1000086\" name=\"full width at half-maximum\"/>\n";
      fql += String(indentation_level, '\t') + "\t\t\t\t</DataType>\n";
      fql += String(indentation_level, '\t') + "\t\t\t</Column>	\n";
      fql += String(indentation_level, '\t') + "\t\t</ColumnDefinition>\n";
      fql += String(indentation_level, '\t') + "\t\t<DataMatrix>\n";

      for (ConsensusMap::const_iterator cit = cm.begin(); cit != cm.end(); ++cit)
      {
        consensus_tag += "\t<PeptideConsensus charge=\"" + String(cit->getCharge()) + "\" id=\"ct_" + String(cit->getUniqueId()) + "\">\n";
        for (std::set<FeatureHandle, FeatureHandle::IndexLess>::const_iterator fhit = cit->getFeatures().begin(); fhit != cit->getFeatures().end(); ++fhit)
        {
          consensus_tag += "\t\t<EvidenceRef feature_ref=\"f_" + String(fhit->getUniqueId()) + "\" assay_refs=\"" + String("none") + "\">\n";           
        }
        consensus_tag += "\t</PeptideConsensus>\n";
        
        fql += "\t\t<Row object_ref=\"ft_" + String(cit->getUniqueId()) + "\">" + String(cit->getMZ()) + " " + String(cit->getRT()) + " " + String(cit->getIntensity()) + " " + String(cit->getQuality()) + " " + String(cit->getWidth()) + "</Row>\n";
      }
      
      consensus_tag += "\t<cvParam accession=\"MS:1001826\" cvRef=\"PSI-MS\" name=\"mass trace reporting: rectangles\"/>\n";
      
      fql += String(indentation_level, '\t') + "\t\t</DataMatrix>\n";
      fql += String(indentation_level, '\t') + "\t</PeptideConsensusQuantLayer>\n";
      
      consensus_tag += fql;
      consensus_tag += String(indentation_level, '\t') + "</PeptideConsensusList>\n";
    }

  } //namespace Internal
} // namespace OpenMS
