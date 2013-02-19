/**
 * @file   AcceptanceRatio.hxx
 * @author Suvayu Ali <Suvayu.Ali@cern.ch>
 * @date   Thu Nov 30 12:47:46 2012
 *
 * @brief  This class implements acceptance ratio function.
 *
 *         This code is based on a skeleton code generated by
 *         RooClassFactory from RooFit. The functional form and the
 *         integral was coded in by hand.
 *
 */


#ifndef __ACCEPTANCERATIO_HXX
#define __ACCEPTANCERATIO_HXX

#include <RooAbsReal.h>
#include <RooRealProxy.h>
#include <RooCategoryProxy.h>
#include <RooAbsReal.h>
#include <RooAbsCategory.h>
#include <RooConstVar.h>


class AcceptanceRatio : public RooAbsReal {
public:

  AcceptanceRatio();
  AcceptanceRatio(const char *name, const char *title,
		  RooAbsReal& time, RooAbsReal& turnon,
		  RooAbsReal& offset, RooAbsReal& beta);
  AcceptanceRatio(const AcceptanceRatio& other, const char* name=0);
  virtual ~AcceptanceRatio();
  virtual TObject* clone(const char* newname) const;
  AcceptanceRatio& operator=(const AcceptanceRatio& other);

protected:

  Double_t evaluate() const;

  RooRealProxy _time;
  RooRealProxy _turnon;
  RooRealProxy _offset;
  RooRealProxy _beta;

private:

  ClassDef(AcceptanceRatio, 1); // Implements acceptance ratio function.
};

#endif	// __ACCEPTANCERATIO_HXX