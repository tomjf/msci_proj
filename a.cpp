#include <iostream> 
#include <cmath>
#include <fstream>
#include <iomanip>

using namespace std;

void Euler (double, double, double, double, double);
void Leapfrog (double, double, double, double, double);
void RK4 (double, double, double, double, double);
void RK4doublependulum (double, double, double, double, double, double, double, double);

double calcdv (double, double, double);
double calcdv2pend( double, double, double, double, double, double);
double calcdomega2pend (double, double, double, double, double);
double calcEnergydoublependulum (double, double, double, double, double);
double energyanalytic(double, double, double);

ofstream dataeuler("dataEuler.txt");
ofstream dataleap ("dataLeapfrog.txt");
ofstream dataRK4 ("dataRK4.txt");
ofstream data2pend ("datadoublependulum.txt");
ofstream accuracyRK4 ("accuracyofRK4.txt");
ofstream dataeulerdamped ("dataEulerdamped.txt");
ofstream accuracyleapfrog ("accuracyofleapfrog.txt");
ofstream accuracyeuler ("accuracyofEuler.txt");
ofstream dataleapdamp ("leapfrogdamped.txt");
ofstream dataRK4damp ("RK4damped.txt");
ofstream accthreshleap ("accuracythresholdofleapfrog.txt"); 
ofstream accthreshRK4 ("accuracythresholdofRK4.txt");
ofstream accthresheuler ("accuracythersholdeuler.txt");

int main(void) 
{

  // define intial conditions for single pendulum and damping
  double u = 0.1, v=0.0, gamma = 0, h = 0.01;
  // define time in dimensionless units we would like to find the solution up to
  double imax = 100;
 
// Testing the Accuracy - COMMENT OUT TO COMPILE REST OF CODE
// This loop was used to vary the time step to study accuracy to find the order of h and to find the 1% accuracy threshold
//for (double dt=0.001;dt<=1;dt+=0.001)
//{
  //  double istop = imax/dt;
	//Euler (gamma,u,v,dt,istop);
	//Leapfrog(gamma, u, v, dt, istop);
	//RK4 (gamma, u, v, dt, istop);
	//cout << dt << endl;
//}

 // calculate the number of steps required at given step size and final value of dimensionless time
 double istop = imax/h;
 // calculate solution for single pendulum using Euler, Leapfrog and RK4 methods
 Euler(gamma, u,v,h,istop);
 Leapfrog(gamma,u,v,h,istop);
 RK4(gamma, u, v, h, istop);
 // define variables for double pendulum. Note that this is where G and R are varied
 double G = 1, R = 100, theta = 0.1, phi = 0, omega = 0, vdoublepend = 0;
 // calculate solution to double pendulum using RK4 method
 RK4doublependulum(h, G, R, theta, phi, omega, vdoublepend, istop);
 return 0; 
} 

// Euler method
void Euler(double gamma, double u, double v, double h, double istop)
{
	double vn, un, energy;
	// calculate energy and output initial values
	energy = v*v +u*u;
	dataeuler << "t" << "\t" << 0 << "\t" << "theta" <<  "\t" << u << "\t" << "w" <<  "\t" << v << "\t" << "E" <<  "\t" << energy << endl;
	//loop that calculates each step using Euler method
	for(int i=1;i<=istop;i++)
	{
		vn = v + ((calcdv(gamma, v, u))*h);
		un = u + (v*h);
		// reset new values ready for the next step
		u = un;
		v = vn;
		//record new values in the ouptut file
		energy = v*v + u*u;
		dataeuler << "t" << "\t" << i*h << "\t" << "theta" <<  "\t" << un << "\t" << "w" <<  "\t" << vn << "\t" << "E" <<  "\t" << energy << endl;
	}
	// accuracy calculations commented out to make program run faster put back in to verify accuracy
	accuracyeuler << setprecision(15) << h << "\t" << energy - 0.01 << endl;
	dataeulerdamped << setprecision(20) << h  << "\t" << energy << "\t" << energy-(((0.01/(1-((gamma*gamma)*0.25))*exp(-(gamma*(istop*h)))*(1 + (gamma/2)*sin(sqrt(4-(gamma*gamma))*(istop*h) - asin(gamma*0.5)))))) << endl;
	// acccuracy threshold
	double percentagedeviation = ((energy-energyanalytic(gamma,istop,h))/energyanalytic(gamma,istop,h));
	if(percentagedeviation <= 0.01) accthresheuler << setprecision(20) << h << "\t" << energy << "\t" << energyanalytic(gamma,istop,h)  << "\t" << percentagedeviation << endl;
}


// Leapfrog method
void Leapfrog(double gamma, double unminus1, double vnminus1, double h, double istop)
{
	 // calculate intial energy and output initial values of variables
	 double energy = vnminus1*vnminus1 + unminus1*unminus1;
	 dataleap << "t" << "\t" << 0 << "\t" << "theta" <<  "\t" << unminus1 << "\t" << "w" <<  "\t" << vnminus1 << "\t" << "E" <<  "\t" << energy << endl;
	 // calculate the first step which must be an Euler step
	 double vn0 = vnminus1 + ((calcdv(gamma,vnminus1,unminus1))*h);
	 double un0 = unminus1 + (vnminus1*h);
	 // output values after first step
	 energy = vn0*vn0 + un0*un0;
	 dataleap << "t" << "\t" << h << "\t" << "theta" <<  "\t" << un0 << "\t" << "w" <<  "\t" << vn0 << "\t" << "E" <<  "\t" << energy << endl;
	 // Loop that Calculates new values at each step from the Leapfrog method
	 for(int i=3;i<=istop;i++)
	 {
		 double vnplus1 = vnminus1 + (2*(calcdv(gamma, vn0, un0))*h);
		 double unplus1 = unminus1 + (2*vn0*h);
		 // reset new values ready for the next step
		 vnminus1 = vn0;
		 vn0 = vnplus1;
		 unminus1 = un0;
		 un0 = unplus1; 
		 energy = unplus1*unplus1 + vnplus1*vnplus1;
		 //record both values in the ouptut file
		 dataleap << "t" << "\t" << i*h << "\t" << "theta" <<  "\t" << un0 << "\t" << "w" <<  "\t" << vn0 << "\t" << "E" <<  "\t" << energy << endl;
	 }
	 // calculating accuracy RK4
	 dataleapdamp << setprecision(20) << h  << "\t" << energy << "\t" << energy-(((0.01/(1-((gamma*gamma)*0.25))*exp(-(gamma*(istop*h)))*(1 + (gamma/2)*sin(sqrt(4-(gamma*gamma))*(istop*h) - asin(gamma*0.5)))))) << endl;
	 double accuracy = energy-0.01;
	 accuracyleapfrog << setprecision(20) << h << "\t" << accuracy << endl;
	 double percentagedeviation = ((energy-energyanalytic(gamma,istop,h))/energyanalytic(gamma,istop,h));
	 if(percentagedeviation <= 0.01) accthreshleap << setprecision(20) << h << "\t" << energy << "\t" << energyanalytic(gamma,istop,h)  << "\t" << percentagedeviation << endl;
}

//RK4 method
void RK4(double gamma, double un, double vn, double h, double istop)
{
	double uk1, uk2, uk3, uk4, vk1, vk2, vk3, vk4, unplus1, vnplus1, energy;
	// output initial conditions to file
	energy = vn*vn + un*un;
	dataRK4 << "t" << "\t" << 0 << "\t" << "theta" <<  "\t" << un << "\t" << "w" <<  "\t" << vn << "\t" << "E" <<  "\t" << energy << endl;
	// loop that integrates in time via the RK4 method
	for(int i=1;i<=istop;i++)
	{
		//work out k1s
		vk1 = -h*((gamma*vn) + (un));
		uk1 = h*vn;
		//work out k2s
		vk2 = -h*((gamma*(vn+(0.5*vk1))) + (un+(0.5*uk1)));
		uk2 = h*(vn + (0.5*vk1));
		//work out k3s
		vk3 = -h*((gamma*(vn+(0.5*vk2))) + (un+(0.5*uk2)));
		uk3 = h*(vn + (0.5*vk2));
		//work out k4s
		vk4 = -h*((gamma*(vn+vk3)) + (un+uk3));
		uk4 = h*(vn + vk3);
		//work out n + 1 values
		vnplus1 = vn + ((1.0/6.0)*(vk1 + (2*vk2) + (2*vk3) + vk4));
		unplus1 = un + ((1.0/6.0)*(uk1 + (2*uk2) + (2*uk3) + uk4));
		//work out energy and output all values to file
		energy = (vnplus1*vnplus1 + unplus1*unplus1);
	    dataRK4 << "t" << "\t" << i*h << "\t" << "theta" <<  "\t" << unplus1 << "\t" << "w" <<  "\t" << vnplus1 << "\t" << "E" <<  "\t" << energy << endl;
		// reset new values ready for the next step
		vn = vnplus1;
		un = unplus1;
	}
	//calculating accuracy of RK4
	double accuracy = 0.01 - energy;
	accuracyRK4 << setprecision(20) << h << "\t" << accuracy << endl;
	dataRK4damp << setprecision(20) << h  << "\t" << energy- (energyanalytic(gamma, istop, h)-energy) << endl;
	double percentagedeviation = (sqrt((energyanalytic(gamma, istop, h)-energy)*(energyanalytic(gamma, istop, h)-energy))/energyanalytic(gamma,istop,h));
	if(percentagedeviation <= 0.0100) accthreshRK4 << setprecision(20) << h << "\t" << energyanalytic(gamma,istop,h)  << "\t" << percentagedeviation << endl;
}

void RK4doublependulum (double h, double G, double R, double theta, double phi, double omega, double v, double istop)
{
       // define all variables
	   double vk1, vk2, vk3, vk4, omegak1, omegak2, omegak3, omegak4;
       double thetak1, thetak2, thetak3, thetak4, phik1, phik2, phik3, phik4;
       double vn, omegan, thetan, phin, energy;
	   // outputting initial values to file
	   energy = calcEnergydoublependulum(R, theta, phi, omega, v);
	   data2pend << "t" << "\t" << 0 << "\t" << "v" << "\t" << v << "\t" << "phi" << "\t" << phi << "\t" << "w" << "\t" <<omega << "\t" << "theta" << "\t" << theta << "\t" << "energy" << "\t" << energy << endl;
       // loop that integrates the 4 coupled linear ODEs in time using the RK4 method
       for(int i=1;i<=istop;i++)
       {
              // finding all the k1's
              vk1 = (h*calcdv2pend(G,R,theta,phi,omega,v));
              phik1 = (h*v);
              omegak1 = (h*calcdomega2pend(G,R,theta,phi,omega));
              thetak1 = (h*omega);
              // finding all the k2's
              vk2 = (h*calcdv2pend(G,R,(theta+((0.5)*thetak1)),(phi+((0.5)*phik1)),(omega+((0.5)*omegak1)),(v+((0.5)*vk1))));
              phik2 = (h*(v + ((0.5)*vk1)));
              omegak2 = (h*calcdomega2pend(G,R,(theta+((0.5)*thetak1)),(phi+((0.5)*phik1)),(omega+((0.5)*omegak1))));
              thetak2 = (h*(omega + ((0.5)*omegak1)));
              // finding all the k3's
              vk3 = (h*calcdv2pend(G,R,(theta+((0.5)*thetak2)),(phi+((0.5)*phik2)),(omega+((0.5)*omegak2)),(v+((0.5)*vk2))));
              phik3 = (h*(v + ((0.5)*vk2)));
              omegak3 = (h*calcdomega2pend(G,R,(theta+((0.5)*thetak2)),(phi+((0.5)*phik2)),(omega+((0.5)*omegak2))));
              thetak3 = (h*(omega + ((0.5)*omegak2)));
              // finding all the k4's
              vk4 = (h*calcdv2pend(G,R,(theta+thetak3),(phi+phik3),(omega+omegak3),(v+vk3)));
              phik4 = (h*(v + vk3));
              omegak4 = (h*calcdomega2pend(G,R,(theta+thetak3),(phi+phik3),(omega+omegak3)));
              thetak4 = (h*(omega + omegak3));
              // finding the new values of all the variables
              vn = v + ((1.0/6.0)*(vk1 + (2*vk2) + (2*vk3) + vk4));
              phin = phi + ((1.0/6.0)*(phik1 + (2*phik2) + (2*phik3) + phik4));
              omegan = omega + ((1.0/6.0)*(omegak1 + (2*omegak2) + (2*omegak3) + omegak4));
              thetan = theta + ((1.0/6.0)*(thetak1 + (2*thetak2) + (2*thetak3) + thetak4));
			  // calculate the energy of the system
			  energy = calcEnergydoublependulum(R, thetan, phin, omegan, vn);
			  // output calculated values for this new step to file
			  data2pend << setprecision(15) << "t" << "\t" << i*h << "\t" << "v" << "\t" << vn << "\t" << "phi" << "\t" << phin << "\t" << "w" << "\t" <<omegan << "\t" << "theta" << "\t" << thetan << "\t" << "energy" << "\t" << energy << endl;
              // reset values ready for the calculation of the next step
              v = vn;
              phi = phin;
              omega = omegan;
              theta = thetan;
             
       }
}
 
// Calculate dv/dt for the double pendulum
double calcdv2pend (double G, double R, double theta, double phi, double omega, double v)
{
       double dv = ((R+1)*theta) + (-(R+1)*phi) + (G*(1-(1/R))*omega) - ((G/R)*v);
       return double (dv);
}
 
// Calculate domega/dt for the double pendulum
double calcdomega2pend (double G, double R, double theta, double phi, double omega)
{
       double domega = ((-(R+1)*theta) + (R*phi) - (G*omega));
       return double (domega);
}

//Fuction to calculate dv/dt
double calcdv(double gamma, double v, double u)
{
	double dv = -(gamma*v) - u;
	return double (dv);
}

// Function to calculate the Energy of the solution for the double pendulum by finding sum of kinetic and potential energies
double calcEnergydoublependulum (double R, double theta, double phi, double omega, double v)
{
	double Pm = (-cos(theta));
	double PM = (-R*(cos(theta) + cos(phi))); 
	double KEm = ((0.5*omega*omega));
	double KEM = ((0.5)*R*((omega*omega) + (v*v) +(2*omega*v*cos(theta-phi))));
	double E = Pm + PM + KEm + KEM;
	return double (E);
}

// Function to calculate analytical energy for the damped single pendulum
double energyanalytic (double gamma, double istop, double h)
{
	double energyanalytic = (((0.01/(1-((gamma*gamma)*0.25))*exp(-(gamma*(istop*h)))*(1 + (gamma/2)*sin(sqrt(4-(gamma*gamma))*(istop*h) - asin(gamma*0.5))))));
	return double (energyanalytic);
}
