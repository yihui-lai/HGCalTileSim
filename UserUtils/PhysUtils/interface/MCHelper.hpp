/**
 * @file    MCHelper.hpp
 * @brief   Helper functions and variables for MC truth information crawling.
 * @author  [Yi-Mu "Enoch" Chen](https://github.com/yimuchen)
 */
#ifndef USERUTILS_PHYSUTILS_MCHELPER_HPP
#define USERUTILS_PHYSUTILS_MCHELPER_HPP

#include "DataFormats/Candidate/interface/Candidate.h"

namespace usr {

/**
 * @defgroup MCHelper MC Truth Helpers
 * @brief    Collections of helpers to get parse the truth information
 * @ingroup  PhysUtils
 * @{
 */

/**
 * @brief
 * Standard model particle ids in enum format (add negative sign for
 * anti-particle)
 */
enum PID
{
  DOWN_QUARK_ID    = 1,
  UP_QUARK_ID      = 2,
  STRANGE_QUARK_ID = 3,
  CHARM_QUARK_ID   = 4,
  BOTTOM_QUARK_ID  = 5,
  TOP_QUARK_ID     = 6,

  ELEC_ID    = 11,
  ELECNU_ID  = 12,
  MUON_ID    = 13,
  MUONNU_ID  = 14,
  TAU_ID     = 15,
  TAUNU_ID   = 16,

  GLUON_ID       = 21,
  GAMMA_ID       = 22,
  Z_BOSON_ID     = 23,
  W_BOSON_ID     = 24,
  HIGGS_BOSON_ID = 25
};

/*------------------------------------------------------------------------------
 *  MC truth topology crawling function.
   ---------------------------------------------------------------------------*/
const reco::Candidate* GetDirectMother( const reco::Candidate*, int );
const reco::Candidate* GetDaughter( const reco::Candidate*, int );

/** @} */

}/* usr */


#endif/* end of include guard: USERUTILS_PHYSUTILS_MCHELPER_HPP */
