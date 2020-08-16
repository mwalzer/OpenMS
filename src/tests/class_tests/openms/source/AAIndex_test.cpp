// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/AAIndex.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

AASequence seq1 = AASequence::fromString("ALEGDEK");
AASequence seq2 = AASequence::fromString("GTVVTGR");
AASequence seq3 = AASequence::fromString("EHVLLAR");


START_TEST(AASequenceIndeces, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

//sequence spec_id protein_id mass GB500 arginin_count KHAG800101 VASM830103 NADH010106 NADH010107 WILM950102 ROBB760107 OOBM850104 FAUJ880111 FINA770101 ARGP820102 M F H Q Y target_log
//ALEGDEK 15 0587  761.368 1337.53 0  129.3 1.145   31  565  1.5200000 -6.60000e+00  -3.240000 1  7.18  5.23 0 0 0 0 0 2.08623342
//GTVVTGR 15 0587  689.394 1442.70 1  383.2 1.042  241  403  7.1800000 -3.00000e-01 -16.010000 1  5.55  5.02 0 0 0 0 0 1.35346120
//EHVLLAR 15 0587  837.494 1442.70 1  318.5 1.259  171  190 18.1300000  3.00000e-01  -9.970000 2  7.73  9.34 0 0 1 0 0 5.22075034

TOLERANCE_ABSOLUTE(0.01)

START_SECTION(static double calculateGB(const AASequence& seq, double T=500.0) )
  TEST_REAL_SIMILAR(AAIndex::calculateGB(seq1), 1337.53)
  TEST_REAL_SIMILAR(AAIndex::calculateGB(seq2), 1442.70)
  TEST_REAL_SIMILAR(AAIndex::calculateGB(seq3), 1442.70)

  TEST_NOT_EQUAL(AAIndex::calculateGB(seq1,100.0), 1337.53)
  TEST_NOT_EQUAL(AAIndex::calculateGB(seq2,100.0), 1442.70)
  TEST_NOT_EQUAL(AAIndex::calculateGB(seq3,100.0), 1442.70)
END_SECTION

START_SECTION(static double aliphatic(char aa))
  TEST_REAL_SIMILAR(AAIndex::aliphatic('A'),1.0)
  TEST_REAL_SIMILAR(AAIndex::aliphatic('B'),0.0)
END_SECTION

START_SECTION(static double acidic(char aa))
  TEST_REAL_SIMILAR(AAIndex::acidic('D'),1.0)
  TEST_REAL_SIMILAR(AAIndex::acidic('A'),0.0)
END_SECTION

START_SECTION(static double basic(char aa))
  TEST_REAL_SIMILAR(AAIndex::basic('K'),1.0)
  TEST_REAL_SIMILAR(AAIndex::basic('A'),0.0)
END_SECTION

START_SECTION(static double polar(char aa))
  TEST_REAL_SIMILAR(AAIndex::polar('S'),1.0)
  TEST_REAL_SIMILAR(AAIndex::polar('A'),0.0)
END_SECTION

START_SECTION(static double getKHAG800101(char aa))
 TEST_REAL_SIMILAR(AAIndex::getKHAG800101('A'),49.1)
END_SECTION

START_SECTION(static double getVASM830103(char aa))

  TEST_REAL_SIMILAR(AAIndex::getVASM830103('A'),0.159)
  TEST_REAL_SIMILAR(AAIndex::getVASM830103('R'),0.194)
  TEST_REAL_SIMILAR(AAIndex::getVASM830103('N'),0.385)
  TEST_REAL_SIMILAR(AAIndex::getVASM830103('D'),0.283)
  TEST_REAL_SIMILAR(AAIndex::getVASM830103('C'),0.187)
  TEST_REAL_SIMILAR(AAIndex::getVASM830103('Q'),0.236)
  TEST_REAL_SIMILAR(AAIndex::getVASM830103('E'),0.206)
  TEST_REAL_SIMILAR(AAIndex::getVASM830103('G'),0.049)
  TEST_REAL_SIMILAR(AAIndex::getVASM830103('H'),0.233)
  TEST_REAL_SIMILAR(AAIndex::getVASM830103('I'),0.581)
  TEST_REAL_SIMILAR(AAIndex::getVASM830103('L'),0.083)
  TEST_REAL_SIMILAR(AAIndex::getVASM830103('K'),0.159)
  TEST_REAL_SIMILAR(AAIndex::getVASM830103('M'),0.198)
  TEST_REAL_SIMILAR(AAIndex::getVASM830103('F'),0.682)
  TEST_REAL_SIMILAR(AAIndex::getVASM830103('P'),0.366)
  TEST_REAL_SIMILAR(AAIndex::getVASM830103('S'),0.150)
  TEST_REAL_SIMILAR(AAIndex::getVASM830103('T'),0.074)
  TEST_REAL_SIMILAR(AAIndex::getVASM830103('W'),0.463)
  TEST_REAL_SIMILAR(AAIndex::getVASM830103('Y'),0.737)
  TEST_REAL_SIMILAR(AAIndex::getVASM830103('V'),0.301)

END_SECTION

START_SECTION(static double getNADH010106(char aa))
  TEST_REAL_SIMILAR(AAIndex::getNADH010106('A'),5.0)
  TEST_REAL_SIMILAR(AAIndex::getNADH010106('R'),-57.0)
  TEST_REAL_SIMILAR(AAIndex::getNADH010106('N'),-77.0)
  TEST_REAL_SIMILAR(AAIndex::getNADH010106('D'),45.0)
  TEST_REAL_SIMILAR(AAIndex::getNADH010106('C'),224.0)
  TEST_REAL_SIMILAR(AAIndex::getNADH010106('Q'),-67.0)
  TEST_REAL_SIMILAR(AAIndex::getNADH010106('E'),-8.0)
  TEST_REAL_SIMILAR(AAIndex::getNADH010106('G'),-47.0)
  TEST_REAL_SIMILAR(AAIndex::getNADH010106('H'),-50.0)
  TEST_REAL_SIMILAR(AAIndex::getNADH010106('I'),83.0)
  TEST_REAL_SIMILAR(AAIndex::getNADH010106('L'),82.0)
  TEST_REAL_SIMILAR(AAIndex::getNADH010106('K'),-38.0)
  TEST_REAL_SIMILAR(AAIndex::getNADH010106('M'),83.0)
  TEST_REAL_SIMILAR(AAIndex::getNADH010106('F'),117.0)
  TEST_REAL_SIMILAR(AAIndex::getNADH010106('P'),-103.0)
  TEST_REAL_SIMILAR(AAIndex::getNADH010106('S'),-41.0)
  TEST_REAL_SIMILAR(AAIndex::getNADH010106('T'),79.0)
  TEST_REAL_SIMILAR(AAIndex::getNADH010106('W'),130.0)
  TEST_REAL_SIMILAR(AAIndex::getNADH010106('Y'),27.0)
  TEST_REAL_SIMILAR(AAIndex::getNADH010106('V'),117.0)
END_SECTION

START_SECTION(static double getNADH010107(char aa))
  TEST_REAL_SIMILAR(AAIndex::getNADH010107('A'),-2.0)
  TEST_REAL_SIMILAR(AAIndex::getNADH010107('R'),-41.0)
  TEST_REAL_SIMILAR(AAIndex::getNADH010107('N'),-97.0)
  TEST_REAL_SIMILAR(AAIndex::getNADH010107('D'),248.0)
  TEST_REAL_SIMILAR(AAIndex::getNADH010107('C'),329.0)
  TEST_REAL_SIMILAR(AAIndex::getNADH010107('Q'),-37.0)
  TEST_REAL_SIMILAR(AAIndex::getNADH010107('E'),117.0)
  TEST_REAL_SIMILAR(AAIndex::getNADH010107('G'),-66.0)
  TEST_REAL_SIMILAR(AAIndex::getNADH010107('H'),-70.0)
  TEST_REAL_SIMILAR(AAIndex::getNADH010107('I'),28.0)
  TEST_REAL_SIMILAR(AAIndex::getNADH010107('L'),36.0)
  TEST_REAL_SIMILAR(AAIndex::getNADH010107('K'),115.0)
  TEST_REAL_SIMILAR(AAIndex::getNADH010107('M'),62.0)
  TEST_REAL_SIMILAR(AAIndex::getNADH010107('F'),120.0)
  TEST_REAL_SIMILAR(AAIndex::getNADH010107('P'),-132.0)
  TEST_REAL_SIMILAR(AAIndex::getNADH010107('S'),-52.0)
  TEST_REAL_SIMILAR(AAIndex::getNADH010107('T'),174.0)
  TEST_REAL_SIMILAR(AAIndex::getNADH010107('W'),179.0)
  TEST_REAL_SIMILAR(AAIndex::getNADH010107('Y'),-7.0)
  TEST_REAL_SIMILAR(AAIndex::getNADH010107('V'),114.0)
END_SECTION

START_SECTION(static double getWILM950102(char aa))
  TEST_REAL_SIMILAR(AAIndex::getWILM950102('A'),2.62)
  TEST_REAL_SIMILAR(AAIndex::getWILM950102('R'),1.26)
  TEST_REAL_SIMILAR(AAIndex::getWILM950102('N'),-1.27)
  TEST_REAL_SIMILAR(AAIndex::getWILM950102('D'),-2.84)
  TEST_REAL_SIMILAR(AAIndex::getWILM950102('C'),0.73)
  TEST_REAL_SIMILAR(AAIndex::getWILM950102('Q'),-1.69)
  TEST_REAL_SIMILAR(AAIndex::getWILM950102('E'),-0.45)
  TEST_REAL_SIMILAR(AAIndex::getWILM950102('G'),-1.15)
  TEST_REAL_SIMILAR(AAIndex::getWILM950102('H'),-0.74)
  TEST_REAL_SIMILAR(AAIndex::getWILM950102('I'),4.38)
  TEST_REAL_SIMILAR(AAIndex::getWILM950102('L'),6.57)
  TEST_REAL_SIMILAR(AAIndex::getWILM950102('K'),-2.78)
  TEST_REAL_SIMILAR(AAIndex::getWILM950102('M'),-3.12)
  TEST_REAL_SIMILAR(AAIndex::getWILM950102('F'),9.14)
  TEST_REAL_SIMILAR(AAIndex::getWILM950102('P'),-0.12)
  TEST_REAL_SIMILAR(AAIndex::getWILM950102('S'),-1.39)
  TEST_REAL_SIMILAR(AAIndex::getWILM950102('T'),1.81)
  TEST_REAL_SIMILAR(AAIndex::getWILM950102('W'),5.91)
  TEST_REAL_SIMILAR(AAIndex::getWILM950102('Y'),1.39)
  TEST_REAL_SIMILAR(AAIndex::getWILM950102('V'),2.30)
END_SECTION

START_SECTION(static double getROBB760107(char aa))
  TEST_REAL_SIMILAR(AAIndex::getROBB760107('A'),0.0)
  TEST_REAL_SIMILAR(AAIndex::getROBB760107('R'),1.1)
  TEST_REAL_SIMILAR(AAIndex::getROBB760107('N'),-2.0)
  TEST_REAL_SIMILAR(AAIndex::getROBB760107('D'),-2.6)
  TEST_REAL_SIMILAR(AAIndex::getROBB760107('C'),5.4)
  TEST_REAL_SIMILAR(AAIndex::getROBB760107('Q'),2.4)
  TEST_REAL_SIMILAR(AAIndex::getROBB760107('E'),3.1)
  TEST_REAL_SIMILAR(AAIndex::getROBB760107('G'),-3.4)
  TEST_REAL_SIMILAR(AAIndex::getROBB760107('H'),0.8)
  TEST_REAL_SIMILAR(AAIndex::getROBB760107('I'),-0.1)
  TEST_REAL_SIMILAR(AAIndex::getROBB760107('L'),-3.7)
  TEST_REAL_SIMILAR(AAIndex::getROBB760107('K'),-3.1)
  TEST_REAL_SIMILAR(AAIndex::getROBB760107('M'),-2.1)
  TEST_REAL_SIMILAR(AAIndex::getROBB760107('F'),0.7)
  TEST_REAL_SIMILAR(AAIndex::getROBB760107('P'),7.4)
  TEST_REAL_SIMILAR(AAIndex::getROBB760107('S'),1.3)
  TEST_REAL_SIMILAR(AAIndex::getROBB760107('T'),0.0)
  TEST_REAL_SIMILAR(AAIndex::getROBB760107('W'),-3.4)
  TEST_REAL_SIMILAR(AAIndex::getROBB760107('Y'),4.8)
  TEST_REAL_SIMILAR(AAIndex::getROBB760107('V'),2.7)
END_SECTION

START_SECTION(static double getOOBM850104(char aa))
  TEST_REAL_SIMILAR(AAIndex::getOOBM850104('A'),-2.49)
  TEST_REAL_SIMILAR(AAIndex::getOOBM850104('R'),2.55)
  TEST_REAL_SIMILAR(AAIndex::getOOBM850104('N'),2.27)
  TEST_REAL_SIMILAR(AAIndex::getOOBM850104('D'),8.86)
  TEST_REAL_SIMILAR(AAIndex::getOOBM850104('C'),-3.13)
  TEST_REAL_SIMILAR(AAIndex::getOOBM850104('Q'),1.79)
  TEST_REAL_SIMILAR(AAIndex::getOOBM850104('E'),4.04)
  TEST_REAL_SIMILAR(AAIndex::getOOBM850104('G'),-0.56)
  TEST_REAL_SIMILAR(AAIndex::getOOBM850104('H'),4.22)
  TEST_REAL_SIMILAR(AAIndex::getOOBM850104('I'),-10.87)
  TEST_REAL_SIMILAR(AAIndex::getOOBM850104('L'),-7.16)
  TEST_REAL_SIMILAR(AAIndex::getOOBM850104('K'),-9.97)
  TEST_REAL_SIMILAR(AAIndex::getOOBM850104('M'),-4.96)
  TEST_REAL_SIMILAR(AAIndex::getOOBM850104('F'),-6.64)
  TEST_REAL_SIMILAR(AAIndex::getOOBM850104('P'),5.19)
  TEST_REAL_SIMILAR(AAIndex::getOOBM850104('S'),-1.60)
  TEST_REAL_SIMILAR(AAIndex::getOOBM850104('T'),-4.75)
  TEST_REAL_SIMILAR(AAIndex::getOOBM850104('W'),-17.84)
  TEST_REAL_SIMILAR(AAIndex::getOOBM850104('Y'),9.25)
  TEST_REAL_SIMILAR(AAIndex::getOOBM850104('V'),-3.97)
END_SECTION

START_SECTION(static double getFAUJ880111(char aa))
  TEST_REAL_SIMILAR(AAIndex::getFAUJ880111('A'),0.)
  TEST_REAL_SIMILAR(AAIndex::getFAUJ880111('R'),1.)
  TEST_REAL_SIMILAR(AAIndex::getFAUJ880111('N'),0.)
  TEST_REAL_SIMILAR(AAIndex::getFAUJ880111('D'),0.)
  TEST_REAL_SIMILAR(AAIndex::getFAUJ880111('C'),0.)
  TEST_REAL_SIMILAR(AAIndex::getFAUJ880111('Q'),0.)
  TEST_REAL_SIMILAR(AAIndex::getFAUJ880111('E'),0.)
  TEST_REAL_SIMILAR(AAIndex::getFAUJ880111('G'),0.)
  TEST_REAL_SIMILAR(AAIndex::getFAUJ880111('H'),1.)
  TEST_REAL_SIMILAR(AAIndex::getFAUJ880111('I'),0.)
  TEST_REAL_SIMILAR(AAIndex::getFAUJ880111('L'),0.)
  TEST_REAL_SIMILAR(AAIndex::getFAUJ880111('K'),1.)
  TEST_REAL_SIMILAR(AAIndex::getFAUJ880111('M'),0.)
  TEST_REAL_SIMILAR(AAIndex::getFAUJ880111('F'),0.)
  TEST_REAL_SIMILAR(AAIndex::getFAUJ880111('P'),0.)
  TEST_REAL_SIMILAR(AAIndex::getFAUJ880111('S'),0.)
  TEST_REAL_SIMILAR(AAIndex::getFAUJ880111('T'),0.)
  TEST_REAL_SIMILAR(AAIndex::getFAUJ880111('W'),0.)
  TEST_REAL_SIMILAR(AAIndex::getFAUJ880111('Y'),0.)
  TEST_REAL_SIMILAR(AAIndex::getFAUJ880111('V'),0.)
END_SECTION

START_SECTION(static double getFINA770101(char aa))
  TEST_REAL_SIMILAR(AAIndex::getFINA770101('A'),1.08)
  TEST_REAL_SIMILAR(AAIndex::getFINA770101('R'),1.05)
  TEST_REAL_SIMILAR(AAIndex::getFINA770101('N'),0.85)
  TEST_REAL_SIMILAR(AAIndex::getFINA770101('D'),0.85)
  TEST_REAL_SIMILAR(AAIndex::getFINA770101('C'),0.95)
  TEST_REAL_SIMILAR(AAIndex::getFINA770101('Q'),0.95)
  TEST_REAL_SIMILAR(AAIndex::getFINA770101('E'),1.15)
  TEST_REAL_SIMILAR(AAIndex::getFINA770101('G'),0.55)
  TEST_REAL_SIMILAR(AAIndex::getFINA770101('H'),1.00)
  TEST_REAL_SIMILAR(AAIndex::getFINA770101('I'),1.05)
  TEST_REAL_SIMILAR(AAIndex::getFINA770101('L'),1.25)
  TEST_REAL_SIMILAR(AAIndex::getFINA770101('K'),1.15)
  TEST_REAL_SIMILAR(AAIndex::getFINA770101('M'),1.15)
  TEST_REAL_SIMILAR(AAIndex::getFINA770101('F'),1.10)
  TEST_REAL_SIMILAR(AAIndex::getFINA770101('P'),0.71)
  TEST_REAL_SIMILAR(AAIndex::getFINA770101('S'),0.75)
  TEST_REAL_SIMILAR(AAIndex::getFINA770101('T'),0.75)
  TEST_REAL_SIMILAR(AAIndex::getFINA770101('W'),1.10)
  TEST_REAL_SIMILAR(AAIndex::getFINA770101('Y'),1.10)
  TEST_REAL_SIMILAR(AAIndex::getFINA770101('V'),0.95)
END_SECTION

START_SECTION(static double getARGP820102(char aa))
  TEST_REAL_SIMILAR(AAIndex::getARGP820102('A'),1.18)
  TEST_REAL_SIMILAR(AAIndex::getARGP820102('R'),0.20)
  TEST_REAL_SIMILAR(AAIndex::getARGP820102('N'),0.23)
  TEST_REAL_SIMILAR(AAIndex::getARGP820102('D'),0.05)
  TEST_REAL_SIMILAR(AAIndex::getARGP820102('C'),1.89)
  TEST_REAL_SIMILAR(AAIndex::getARGP820102('Q'),0.72)
  TEST_REAL_SIMILAR(AAIndex::getARGP820102('E'),0.11)
  TEST_REAL_SIMILAR(AAIndex::getARGP820102('G'),0.49)
  TEST_REAL_SIMILAR(AAIndex::getARGP820102('H'),0.31)
  TEST_REAL_SIMILAR(AAIndex::getARGP820102('I'),1.45)
  TEST_REAL_SIMILAR(AAIndex::getARGP820102('L'),3.23)
  TEST_REAL_SIMILAR(AAIndex::getARGP820102('K'),0.06)
  TEST_REAL_SIMILAR(AAIndex::getARGP820102('M'),2.67)
  TEST_REAL_SIMILAR(AAIndex::getARGP820102('F'),1.96)
  TEST_REAL_SIMILAR(AAIndex::getARGP820102('P'),0.76)
  TEST_REAL_SIMILAR(AAIndex::getARGP820102('S'),0.97)
  TEST_REAL_SIMILAR(AAIndex::getARGP820102('T'),0.84)
  TEST_REAL_SIMILAR(AAIndex::getARGP820102('W'),0.77)
  TEST_REAL_SIMILAR(AAIndex::getARGP820102('Y'),0.39)
  TEST_REAL_SIMILAR(AAIndex::getARGP820102('V'),1.08)
END_SECTION


END_TEST
