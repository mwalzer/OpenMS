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
// $Maintainer: Sandro Andreotti $
// $Authors: Sandro Andreotti, Andreas Bertsch$
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_DENOVO_DENOVOIDENTIFICATION_H
#define OPENMS_ANALYSIS_DENOVO_DENOVOIDENTIFICATION_H

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <vector>

namespace OpenMS
{
  /**
    @brief Base class for de novo identification



        @ingroup Analysis_DeNovo
  */
  class OPENMS_DLLAPI DeNovoIdentification :
    public DefaultParamHandler
  {
public:

    /** @name Constructors and destructors
    */
    //@{
    /// default constructor
    DeNovoIdentification();

    /// destructor
    virtual ~DeNovoIdentification();

    /// copy constructor
    DeNovoIdentification(const DeNovoIdentification & rhs);
    //@}

    /// assignment operator
    DeNovoIdentification & operator=(const DeNovoIdentification & rhs);

    /// performs an ProteinIdentification run on a RichPeakMap
    virtual void getIdentifications(std::vector<SpectrumIdentification> & ids, const RichPeakMap & exp) = 0;

    /// performs an ProteinIdentification run on a PeakSpectrum
    virtual void getIdentification(SpectrumIdentification & id, const RichPeakSpectrum & spectrum) = 0;

  };

} // namespace OpenMS

#endif // OPENMS_ANALYSIS_DENOVO_DENOVOIDENTIFICATION_H
