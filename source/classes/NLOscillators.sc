/*

These UGens were written by Oswald Berthold, re-using code
from Dan Stowell's chaos UGens, re-using code from Lance
Putnam's chaos UGens. Released under the GNU Public License.

*/

// van der Pol Oscillator
VanDerPol : MultiOutUGen {
	const <equation="x' = y\ny' = e(1-x^2)x' - x + a in(0)"; // driven model, with special case a=0
	*ar { arg drv=0, freq=22050, e=0.2, a=0.0, h=0.05, xi=0.1, yi=0, mul=1.0, add=0.0;
		^this.multiNew('audio', drv, freq, e, a, h, xi, yi).madd(mul, add)
	}	
	init { arg ... theInputs;
		inputs = theInputs;
		^this.initOutputs(2, rate);
	}
	*categories {	^ #["UGens>Generators>Chaotic"]	}
}

// Duffing Oscillator
Duffing : MultiOutUGen {
	const <equation="x' = y\ny' = -delta y - alpha x - beta x^3 + gamma in(0)";
	*ar { arg drv=0, freq=22050, alpha=0.0, beta=0.0, gamma=0.0, delta=0.1, h=0.05, xi=0.1, yi=0, mul=1.0, add=0.0;
		^this.multiNew('audio', drv, freq, alpha, beta, gamma, delta, h, xi, yi).madd(mul, add)
	}	
	init { arg ... theInputs;
		inputs = theInputs;
		^this.initOutputs(2, rate);
	}
	*categories {	^ #["UGens>Generators>Chaotic"]	}
}
