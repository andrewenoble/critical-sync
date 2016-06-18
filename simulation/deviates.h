// Ranq2 is somewhat faster than Ran but has a shorter period 
// (see NR 3rd ed.).
// The NR function Normaldev has been modified here to generate  
// deviates from a standard normal distribution.
struct SNormaldev : Ranq2 { 

   	SNormaldev(Ullong i) : Ranq2(i) {} 
   	
   	Doub dev() {
  		Doub u,v,x,y,q;
  		do {
 			u = doub();
 			v = 1.7156*(doub()-0.5);
 			x = u - 0.449871;
 			y = abs(v) + 0.386595;
 			q = SQR(x) + y*(0.19600*y-0.25472*x);
  		} while (q > 0.27597
 			&& (q > 0.27846 || SQR(v) > -4.*log(u)*SQR(u)));
  		return v/u;
   	}
};
