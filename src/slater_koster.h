/***************************************************************************
* Slater-Koster method
*
* Copyright (C) 2023 Hidemaro Suwa <suwamaro@phys.s.u-tokyo.ac.jp>
*
*****************************************************************************/

#ifndef _SLATER_KOSTER_H_
#define _SLATER_KOSTER_H_

#include "tbc.h"

namespace TBC
{
  class SlaterKoster {
  public:
    SlaterKoster();
    ~SlaterKoster(){}
    mat_type t_SK() const { return t_SK_; }

  private:
    double ts = 1.5;
    double tp = - 1.0;
    double td = 0.25;
    mat_type t_SK_;
  };

  SlaterKoster::SlaterKoster(){
    vec_type v_SK({td, tp, ts, tp, td});
    t_SK_ = arma::diagmat(v_SK);
  }
}

#endif // _SLATER_KOSTER_H_
