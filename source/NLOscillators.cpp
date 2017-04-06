/*
 *  X75Oscillators.cpp
 *
 *  Some ODE based oscillators, Oswald Berthold
 *  Derived from Dan Stowells ChaosUGens
 *  Partly based on Lance Putnam's chaos UGen work.
 *  Available under the terms of the GNU Public License (GPL).
 *
 */

#include "SC_PlugIn.h"
#define TWOPI 6.283185307179586
#define PI 3.141592653589793
#define RECPI 0.3183098861837907
#define RECTWOPI 0.1591549430918953
#define ONESIXTH 0.1666666666666667

// InterfaceTable contains pointers to functions in the host (server).
static InterfaceTable *ft;

struct NonLinear : public Unit {
  double x0, y0, xn, yn, xnm1, ynm1;
  float counter;
  //bool stable;
};

// declare struct to hold unit generator state
struct VanDerPol : public NonLinear {
  double frac;
};
struct Duffing : public NonLinear {
  double frac;
};

// declare unit generator functions 
extern "C"
{
  void load(InterfaceTable *inTable);

  void VanDerPol_next(VanDerPol *unit, int inNumSamples);
  void VanDerPol_Ctor(VanDerPol *unit);
  void Duffing_next(Duffing *unit, int inNumSamples);
  void Duffing_Ctor(Duffing *unit);
};

//////////////////////////////////////////////////////////////////

void VanDerPol_next(VanDerPol *unit, int inNumSamples)
{
  float *xout = ZOUT(0);
  float *yout = ZOUT(1);
  float *drv = ZIN(0);
  float freq = ZIN0(1);
  double e = ZIN0(2);
  double a = ZIN0(3);
  double h = ZIN0(4);
  double x0 = ZIN0(5);
  double y0 = ZIN0(6);
	
  double xn = unit->xn;
  double yn = unit->yn;
  float counter = unit->counter;
  double xnm1 = unit->xnm1;
  double ynm1 = unit->ynm1;
  double frac = unit->frac;
	
  float samplesPerCycle;
  double slope;
  if(freq < unit->mRate->mSampleRate){
    samplesPerCycle = unit->mRate->mSampleRate / sc_max(freq, 0.001f);
    slope = 1.f / samplesPerCycle;
  }
  else {
    samplesPerCycle = 1.f;
    slope = 1.f;
  }

  if((unit->x0 != x0) || (unit->y0 != y0)){
    xnm1 = xn;
    ynm1 = yn;
    unit->x0 = xn = x0;
    unit->y0 = yn = y0;
  }
	
  double dx = xn - xnm1;
  double dy = yn - ynm1;
	
  for (int i=0; i<inNumSamples; ++i) {
    if(counter >= samplesPerCycle){
      counter -= samplesPerCycle;
      frac = 0.f;

      xnm1 = xn;
      ynm1 = yn;
			
      double k1x, k2x, k3x, k4x, 
	k1y, k2y, k3y, k4y, 
	kxHalf, kyHalf;
			
      // 4th order Runge-Kutta approximate solution for differential equations
      k1x = h * ynm1;
      k1y = h * (e*ynm1*(1 - pow(xnm1, 2)) - xnm1 + (a * ZXP(drv)));
      kxHalf = k1x * 0.5;
      kyHalf = k1y * 0.5;
			
      k2x = h * (ynm1 + kyHalf);
      k2y = h * (e*(ynm1 + kyHalf)*(1 - pow(xnm1-kxHalf, 2)) - (xnm1+kxHalf) + (a * ZXP(drv)));
      kxHalf = k2x * 0.5;
      kyHalf = k2y * 0.5;
			
      k3x = h * (ynm1 + kyHalf);
      k3y = h * (e*(ynm1 + kyHalf)*(1 - pow(xnm1-kxHalf, 2)) - (xnm1+kxHalf) + (a * ZXP(drv)));

      k4x = h * (ynm1 + k3y);
      k4y = h * (e*(ynm1 + k3y)*(1 - pow(xnm1-k3x, 2)) - (xnm1+k3x) + (a * ZXP(drv)));

      xn = xn + (k1x + 2.0*(k2x + k3x) + k4x) * ONESIXTH;
      yn = yn + (k1y + 2.0*(k2y + k3y) + k4y) * ONESIXTH;
			
      dx = xn - xnm1;
      dy = yn - ynm1;
    }
    counter++;
    ZXP(xout) = (xnm1 + dx * frac) * 0.5f;
    ZXP(yout) = (ynm1 + dy * frac) * 0.5f;
    frac += slope;
  }
	
  unit->xn = xn;
  unit->yn = yn;
  unit->counter = counter;
  unit->xnm1 = xnm1;
  unit->ynm1 = ynm1;
  unit->frac = frac;
}

void VanDerPol_Ctor(VanDerPol* unit){
  SETCALC(VanDerPol_next);
	
  unit->x0 = unit->xn = unit->xnm1 = ZIN0(5);
  unit->y0 = unit->yn = unit->ynm1 = ZIN0(6);
  unit->counter = 0.f;
  unit->frac = 0.f;
	
  VanDerPol_next(unit, 1);
}

//////////////////////////////////////////////////////////////////

void Duffing_next(Duffing *unit, int inNumSamples)
{
  float *xout = ZOUT(0);
  float *yout = ZOUT(1);
  float *drv = ZIN(0);
  float freq = ZIN0(1);
  double alpha = ZIN0(2);
  double beta = ZIN0(3);
  double gamma = ZIN0(4);
  double delta = ZIN0(5);
  double h = ZIN0(6);
  double x0 = ZIN0(7);
  double y0 = ZIN0(8);
	
  double xn = unit->xn;
  double yn = unit->yn;
  float counter = unit->counter;
  double xnm1 = unit->xnm1;
  double ynm1 = unit->ynm1;
  double frac = unit->frac;
	
  float samplesPerCycle;
  double slope;
  if(freq < unit->mRate->mSampleRate){
    samplesPerCycle = unit->mRate->mSampleRate / sc_max(freq, 0.001f);
    slope = 1.f / samplesPerCycle;
  }
  else {
    samplesPerCycle = 1.f;
    slope = 1.f;
  }

  if((unit->x0 != x0) || (unit->y0 != y0)){
    xnm1 = xn;
    ynm1 = yn;
    unit->x0 = xn = x0;
    unit->y0 = yn = y0;
  }
	
  double dx = xn - xnm1;
  double dy = yn - ynm1;
	
  for (int i=0; i<inNumSamples; ++i) {
    if(counter >= samplesPerCycle){
      counter -= samplesPerCycle;
      frac = 0.f;

      xnm1 = xn;
      ynm1 = yn;

      double k1x, k2x, k3x, k4x, 
	k1y, k2y, k3y, k4y, 
	kxHalf, kyHalf;

      // 4th order Runge-Kutta approximate solution for differential equations
      k1x = h * ynm1;
      k1y = h * (-delta*ynm1 - alpha*xnm1 - beta*pow(xnm1, 3) + (gamma * ZXP(drv)));
      kxHalf = k1x * 0.5;
      kyHalf = k1y * 0.5;
			
      k2x = h * (ynm1 + kyHalf);
      k2y = h * (-delta*(ynm1 + kyHalf) - alpha*(xnm1+kxHalf) - beta*pow(xnm1+kxHalf, 3) + (gamma * ZXP(drv)));
      kxHalf = k2x * 0.5;
      kyHalf = k2y * 0.5;
			
      k3x = h * (ynm1 + kyHalf);
      k3y = h * (-delta*(ynm1 + kyHalf) - alpha*(xnm1+kxHalf) - beta*pow(xnm1+kxHalf, 3) + (gamma * ZXP(drv)));

      k4x = h * (ynm1 + k3y);
      k4y = h * (-delta*(ynm1 + k3y) - alpha*(xnm1+k3x) - beta * pow(xnm1+k3x, 3) + (gamma * ZXP(drv)));

      xn = xn + (k1x + 2.0*(k2x + k3x) + k4x) * ONESIXTH;
      yn = yn + (k1y + 2.0*(k2y + k3y) + k4y) * ONESIXTH;
			
      dx = xn - xnm1;
      dy = yn - ynm1;
    }
    counter++;
    ZXP(xout) = (xnm1 + dx * frac) * 0.5f;
    ZXP(yout) = (ynm1 + dy * frac) * 0.5f;
    frac += slope;
  }
	
  unit->xn = xn;
  unit->yn = yn;
  unit->counter = counter;
  unit->xnm1 = xnm1;
  unit->ynm1 = ynm1;
  unit->frac = frac;
}

void Duffing_Ctor(Duffing* unit){
  SETCALC(Duffing_next);
	
  unit->x0 = unit->xn = unit->xnm1 = ZIN0(5);
  unit->y0 = unit->yn = unit->ynm1 = ZIN0(6);
  unit->counter = 0.f;
  unit->frac = 0.f;
	
  Duffing_next(unit, 1);
}


//////////////////////////////////////////////////////////////////


// the load function is called by the host when the plug-in is loaded
//void load(InterfaceTable *inTable)
PluginLoad(NLOscillators)
{
  ft = inTable;
  DefineSimpleUnit(VanDerPol);
  DefineSimpleUnit(Duffing);
}

////////////////////////////////////////////////////////////////////
