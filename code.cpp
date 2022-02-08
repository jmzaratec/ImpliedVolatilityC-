// INCLUDE LIBRARIES
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <fstream>
using namespace std;

double GENERAL=0;

/// ESTRUC///
struct Payoff 
{
	virtual double operator () (double y)=0;
};

double f(double y, double X, double sigma, double T, double r)//f(double xi)
{
	
	return(X*(exp(-r*sigma*T)))*(exp(T*y));
	//return(xi*xi);
}

double integrate(int n, double a, double b, double X, double sigma, double T, double r)//integrate(int n, double a, double b) 
{
		

		////////////////////////////////////////////////////

		///////////////////////////////////////////////	


	
	double h=((b-a))/double(n);
	double sum1=0, sum2=0;
	
	for (int i=1;i<=(n/2)-1;i++){
		double yi=a +(h*2*i);
		sum1= sum1 + f(yi,X,sigma,T,r);
	}

	
	for (int i=1;i<=(n/2);i++){
		double yi=a +(h*((2*i)-1));
		sum2= sum2 + f(yi,X,sigma,T,r);
	}

	
	double result = (h/3)*(f(a,X,sigma,T,r)+f(b,X,sigma,T,r)+(2*sum1)+(4*sum2));
	
	cout<<result<<"\n";

		return(result);

}
		
//// NEW FUNCTION!!!
class payoffcalltest: public Payoff{
	double X;
public: 
	
	payoffcalltest(double X_):X(X_){}
	double operator()(double y){
	double a=(+X*(exp(y))-X);
	double b=0.;
	return max(a,b);
	}
};

class payoffputtest: public Payoff{
	double X;
public: 
	
	payoffputtest(double X_):X(X_){}
	double operator()(double y){
	double a=(-X*(exp(y))+X);
	double b=0.;
	return max(a,b);
	}
};

double payoffcall(double y, double X){

	double a=(X*(exp(y))-X);
	double b=0.;
	return max(a,b);
}

double A(double x,double r,double sigma,double T,double k){
	
	double pi=atan(1.0)*4;

	double a1,a2,a3;
	a1= 1/(sqrt(2*sigma*sigma*pi*T));
	a2= (-0.5*k*x)-(0.125*sigma*sigma*k*k*T)-(r*T);
	a3= a1*(exp(a2));
	return a3;
}

double B(double x,double y,double sigma,double T,double k){
	
	
	double a1,a2,a3;
	a1= -pow((x-y),2)/(2*sigma*sigma*T);
	a2= 0.5*k*y;
	a3= exp(a1+a2);
	return a3;
}

double getx(double S0,double X){
	double x= log(S0/X);
	return x;
}

double gety(double ST,double X){
	double y=log(ST/X);// any number between a and b
	return y;
}

double getk(double r,double sigma){
	double k= ((2*r)/(sigma*sigma))-1;
	return k;
}

/// NEW FUNCTION FXY
double fxy_objettest(double x,double y,double sigma,double T,double k, double X, Payoff& off){//adding object off 
	
	//payoffcalltest o1(X);
	
	double a3= B(x,y,sigma,T,k);
	double c=off(y);
	return a3*c;

}


double fxy(double x,double y,double sigma,double T,double k, double X){

/*//B
	double a1,a2,a3;
	a1= -pow((x-y),2)/(2*sigma*sigma*T);
	a2= 0.5*k*y;
	a3= exp(a1+a2);
//payoff	
	double a=(X*exp(y))-X;
	double b=0.;
	double c=max(a,b);
	
	*/
	double a3= B(x,y,sigma,T,k);
	double c=payoffcall(y,X);

	return a3*c;


}


double integralCall(double S0,double X, double sigma, double T, double r)//integrate(int n, double a, double b) 
{
	int n= 10000;
	double a=-10, b=10;
	double x, y,k;
	
	x= getx(S0,X);
	y= rand()% int(a) &&rand()%int(b);//gety (ST,X);
	k= getk(r,sigma);

	double h=((b-a))/double(n);
	double sum1=0, sum2=0;
	
	for (int i=1;i<=(n/2)-1;i++){
		double yi=a +(h*2*i);
		sum1= sum1 + fxy(x,yi,sigma,T,k,X);
		
	}

	
	for (int i=1;i<=(n/2);i++){
		double yi=a +(h*((2*i)-1));
		sum2= sum2 + fxy(x,yi,sigma,T,k,X);
	}

	
	double result = (h/3)*((fxy(x,a,sigma,T,k,X))+fxy(x,b,sigma,T,k,X)+(2*sum1)+(4*sum2));
	
	//cout<<result<<"\n";

		return(result);

}

// NEW FUNCTION INTEGRATE!!!
double integralCall_objetctest(double S0,double X, double sigma, double T, double r,Payoff& off)// add objet off
{
	int n= 177147;
	double a=-10, b=10;
	double x, y,k;
	

	x= getx(S0,X);
	y= rand()% int(a) &&rand()%int(b);//gety (ST,X);
	k= getk(r,sigma);

//	payoffcalltest o1(X);
	off(y);
	
	double h=((b-a))/double(n);
	double sum1=0, sum2=0;
	
	for (int i=1;i<=(n/2)-1;i++){
		double yi=a +(h*2*i);
		//sum1= sum1 + fxy(x,yi,sigma,T,k,X);
		sum1= sum1 + fxy_objettest(x,yi,sigma,T,k,X,off);//OFF
	}

	
	for (int i=1;i<=(n/2);i++){
		double yi=a +(h*((2*i)-1));
		//sum2= sum2 + fxy(x,yi,sigma,T,k,X);
		sum2= sum2 + fxy_objettest(x,yi,sigma,T,k,X,off);//OFF
	}

	
	//double result = (h/3)*((fxy(x,a,sigma,T,k,X))+fxy(x,b,sigma,T,k,X)+(2*sum1)+(4*sum2));
	double result = (h/3)*((fxy_objettest(x,a,sigma,T,k,X,off))+fxy_objettest(x,b,sigma,T,k,X,off)+(2*sum1)+(4*sum2));//OFF
	//cout<<result<<"\n";

		return(result);
}

double valueCallOption(double S0,double X, double sigma, double T, double r)//integrate(int n, double a, double b) 
{
	double x, y,k;
	double call=0;
	x= getx(S0,X);
	//y= gety (ST,X);
	k= getk(r,sigma);
	
	call = A(x,r,sigma,T,k)*integralCall(S0,X,sigma,T,r);

	cout<<"The value of call Option is: "<<call<<"\n";
	return call;

}

// NEW FUNCTION TEST
double valueCallOption_TEST(double S0,double X, double sigma, double T, double r,Payoff& off)//integrate(int n, double a, double b) 
{
	double x, y,k;
	double call=0;
	x= getx(S0,X);
	//y= gety (ST,X);
	k= getk(r,sigma);
	
	call = A(x,r,sigma,T,k)*integralCall_objetctest(S0,X,sigma,T,r,off);//ADD OFF

	//cout<<"The value of  Option TEST is: "<<call<<"\n";
	return call;

}




///////////////////////////////////////////////////////////////////////
//VECTOR!


class MVector
{	
  // storage for the new vector class 
  vector<double> v;
  public:
    // constructor
    explicit MVector(){}
    explicit MVector(int n):v(n){}
    explicit MVector(int n,double x):v(n,x){}
    // equate vectors;
    MVector& operator=(const MVector& X)
    {if(&X==this)return *this;v=X.v;return *this;}
    // access data in vector 
    double& operator[](int index){return v[index];}
    // access data in vector (const)
    double operator[](int index) const {return v[index];}
    // size of vector
    int size() const {return v.size();}

	void push_back(double x){
		v.push_back(x);
	}
	// Average of a vector
	MVector average (const MVector&x){
		double suma=0;
		
		for (int i=0;i<x.size();i++){
		suma+=x[i];
		}
		
		cout<<"Average: (_"<<suma/x.size()<<"_)\n";
				
	return x;
	}

	// Moving average of a vector

	MVector movingAverage (const MVector&x,double m){
		
		MVector mA;
	

		for (int i=0;i<=x.size()-m;i++){
		double suma=0;
		double averagem=0;
		for (int n=i;n<i+m;n++){
		suma+=x[n];
		//cout<<"cycle 1: suma"<<suma;
		}
		averagem=suma/m;
		mA.push_back(averagem);
		
		
		//cout<<"avergare"<<mA[i];
		}
		cout<<"Moving Average"<< " ("<<m<<"): \n";
	
		mA.show(mA);
		
	return mA;
	}

	//Show vector
	MVector show (const MVector&x){
		cout<<"V: (_";
	for (int i=0;i<x.size();i++)
	{
		if (i<1) cout<<x[i]; 
		else cout<<" , "<<x[i];
			
	}
	cout<<"_) \n";
				
	return x;
	}

	MVector expmovAverage (const MVector&x,double m,double alfa){
		
		MVector mA;
		MVector expmA;
		double exponential=0;

		for (int i=0;i<=x.size()-m;i++){
		double suma=0;
		double averagem=0;
		for (int n=i;n<i+m;n++){
		suma+=x[n];
		//cout<<"cycle 1: suma"<<suma;
		}
		averagem=suma/m;
		mA.push_back(averagem);
		}
		
		exponential = mA[0]*(alfa)+(1.-alfa)*mA[0];
		
		for (int i=0;i<mA.size();i++){
			exponential = mA[i]*(alfa)+(1.-alfa)*exponential;
						
			expmA.push_back(exponential);
			
		}
		
		cout<<"Exponential Moving Average"<< " ("<<m<<")"<<"alfa: ("<<alfa<<"): \n";
		expmA.show(expmA);
		
	return expmA;
	}


	MVector variance (const MVector&x){
		double suma=0;
		double sumasquares=0;
		double variance=0;
		double m;
		MVector VAR;
		m= x.size();


		for (int i=0;i<x.size();i++){
		suma+=x[i];
		sumasquares+= pow(x[i],2);
		}
		
		suma=pow(suma,2);

		variance = (sumasquares/(x.size()-1)) -	(suma/(pow(m,2)-m));
		
		VAR.push_back(variance);
		cout<<"Variance: (_"<<VAR[0]<<"_)\n";
				
	return VAR;
	}

	MVector convertToLog (const MVector&x){
		
		MVector Clog;
		double ln=0;
		

		for (int i=0;i<x.size();i++){
		
		ln= log(x[i]);
		Clog.push_back(ln);
		}
		
		cout<<"Convert vector into log(x[i])\n";
		Clog.show(Clog);
				
	return Clog;
	}
	
	MVector dlog(const MVector&x){
		
		MVector dlog;
		double dln=0;
		

		for (int i=1;i<x.size();i++){
		
		dln= (log(x[i])-log(x[i-1]));
		dlog.push_back(dln);
		}
		
		cout<<" dlog vector {log(x[i]) -log (x[i-1])}\n";
		dlog.show(dlog);
				
	return dlog;
	}

	MVector stddev (const MVector&x){	
	MVector std;
	double dev=0;

		dev=sqrt((std.variance(x)[0]));
		cout<<"Std dev (sigma): (_"<<dev<<"_)\n";
	return(std);
	}


}; // end class MVector

//////////////////////////////////////////////////////////
//SECANT!!
class Secant// public payoffcalltest
{

public:
		
	//double a,b,c,d,e;
	double S0;
	double X;
	double sigma;
	double T; 
	double r;
	double callmarket;
	
	
	double secantstart(double x0,double xminus1,Payoff& off){
	
	double x0old, max=20,tol=1.e-4; 
	int n;
	//cout<<" 0 "<<x0<<"-"; Show iteractions
	//cout<<xminus1<<"\n";Show iteractions

	////////////////////////
	/*
	std::ofstream output;
	 output.open("C:/Users/jm_zarate_c/Documents/Visual Studio 2010/Projects/ScientificComputing/5_Implied_volatilities/allcode/pivot.csv");

	if(!output.is_open())
 {
     std::cout << " File not opened \n"; // To validate if the file was opened.
    
	 // stop the program here
    throw;
  }
  */
	/////////////////////////////////////////




	for (n=0; n<max; n++)
	{
	x0old = x0;
		//(double S0,double X, double sigma, double T, double r,Payoff& off,double callmarket)
	//x0=ite(x0,xminus1,function(a,b,c,d,e,x0),function(a,b,c,d,e,xminus1));
	x0=ite(x0,xminus1,function(S0,X,x0,T,r,off,callmarket),function(S0,X,xminus1,T,r,off,callmarket));
	xminus1=x0old;
	//cout<<n<<" "<<x0<<"-"; //Show iteractions
	//cout<<xminus1<<"\n"; //Show iteractions	
	
	/////////////////////////////
	 
	//GENERAL=GENERAL+n;
	//cout<<GENERAL<<"\n";

	 // output <<n<<" , "<<x0<<"," <<xminus1<<","<<GENERAL<< std::endl;
	 
	 // std::cout << " File write successful \n";
  
  // CLOSE file
  

	



	


	///////////////////////////////////


	if (abs((x0 -xminus1))<= tol) break;
	}
	//////
	

	///////
	if (n>(max-1.))
		{
		cout<<"Validate data input \n";
		x0=0.;
		}
	else
		{
		x0; //cout<<"x="<<x0<<"\n";
		}
	return x0;
	//output.close();
	
	}


	
	double ite(double x0, double xminus1, double fx0, double fxminus1)
{
	double a1=0,a2=0,a3=0;

	a3= ((fx0-fxminus1) / (x0-xminus1));
	a2= fx0 / a3;
	a1= x0-a2;

	return (a1);
}

	double function(double S0,double X, double sigma, double T, double r,Payoff& off,double callmarket)//(double a2,double b2,double c2,double d2,double e2, double x)
{
	//	cout<<"variable"<<x<<" : "<<(((a2*x+b2)*x+c2)*x+d2)*x + e2<<"\n";
	//cout<<"The solution using  METHOD is:\n";
	
	//return(((a2*x+b2)*x+c2)*x+d2)*x + e2;
	//return((a2*(x*x*x*x))+(b2*(x*x*x))+(c2*(x*x))+(d2*x)+ (e2));
	double a;
	a = callmarket-valueCallOption_TEST(S0,X,sigma,T,r,off);

	return (a);
}

};


//int main(int argc, char** argv) { for VECTOR   & pull info!!!
//int main(int argc, char** argv){ for SECANT!!!!!!
int main(int argc, char* argv[]) // for options!!!
{
	
	cout<<"WELCOME\n"<<"1. Please insert Options prices in *.txt file\n";
	// system("pause");

//VECTORS MAIN!!
	
			


	//std::ifstream input("C:/Users/jm_zarate_c/Documents/Visual Studio 2010/Projects/ScientificComputing/4_Execises_class/4_newtonClasses/2_Vectors/datasource.txt");
	std::ifstream input("C:/Users/jm_zarate_c/Documents/Visual Studio 2010/Projects/ScientificComputing/5_Implied_volatilities/allcode/datasource.txt");
	Secant o5;

	MVector callMarketV;
	MVector S0V;
	MVector XV;
	MVector TV;
	MVector rV;
	MVector impsigmaV;
	MVector typeOptV;


	double a1=0;
	double a2=0;
	double a3=0;
	double a4=0;
	double a5=0;
	double a6=0;
	double initial;
	double initial2;
	

		//cout<<input.eof();


	do{
		input>>a1>>a2>>a3>>a4>>a5>>a6;
		callMarketV.push_back(a1);
		S0V.push_back(a2);
		XV.push_back(a3);
		TV.push_back(a4);
		rV.push_back(a5);
		typeOptV.push_back(a6);
	
	} while (input.eof()==0);
	
	
	//typeOptV.show(typeOptV);
 
  input.close();
  
	  //system("pause");
  
  /*
  
	double S0=30;
	double X=40;
	double sigma=0.3;
	double T=2; 
	double r=0.05;

	*/
	
	//valueCallOption(S0,X,sigma,T,r);
		
	//valueCallOption_TEST(S0V[0],XV[0],sigma,TV[0],rV[0],Ecall);
	//valueCallOption_TEST(S0,X,sigma,T,r,Ecall);
	//valueCallOption_TEST(S0,X,sigma,T,r,Eput);

	//	system("pause");
				// MAIN SECANT!!!
		//S0V.show(S0V);

initial = .1;//.3;//(-b/(2.*a))-1;
initial2 = 0.9;//.5;//(-b/(2.*a))+1;
	
	
	//payoffputtest Eput(XV[0]);

for (int i=0;i<S0V.size();i++){
	//cout<<S0V[i]<<" "<<XV[i]<<" "<<TV[i]<<" "<<"\n";
	double a;
	//cout<<impsigmaV[i];
	o5.callmarket = double(callMarketV[i]);
	o5.S0 =double(S0V[i]);
	o5.X = double(XV[i]);
	o5.T = double(TV[i]);
	o5.r = double(rV[i]);
	 
	//cout<<S0V[i]<<" "<<XV[i]<<" "<<TV[i]<<" "<<"\n";
	//payoffcalltest Eput(XV[i]);	
	//a= o5.secantstart(initial,initial2,Eput);

	if (double(typeOptV[i])==1){
		payoffcalltest Ecall(XV[i]);	
		a= o5.secantstart(initial,initial2,Ecall);
		
	}
	else if (double(typeOptV[i])==2){
		payoffputtest Eput(XV[i]);
		a= o5.secantstart(initial,initial2,Eput);
	}
	else {
		cout<<"Validate options type";
		}
	impsigmaV.push_back(a);
	//cout<<impsigmaV[i];
}

	// system("pause");
	 std::ofstream output;
	 output.open("C:/Users/jm_zarate_c/Documents/Visual Studio 2010/Projects/ScientificComputing/5_Implied_volatilities/allcode/impliedVolatilities.csv");
	 if(!output.is_open())
 {
     std::cout << " File not opened \n"; // To validate if the file was opened.
    
	 // stop the program here
    throw;
  }

	 output << "Value of Option" << " , " << "Spot Price" << " , " << "Strike Price" << " , " << "Maturity"<< " , " << " Interest rate "<<" , "<<"Type Call = 1-Put = 2"<<" , "<<"Implied Volatility"  <<std::endl;
	   for (int i=0;i<S0V.size();i++){
   	  output << callMarketV[i] << " , " << S0V[i] << " , " << XV[i] << " , " << TV[i]<< " , " << rV[i]<<" , "<<typeOptV[i]<<" , "<<impsigmaV[i]<<std::endl;
	 }
	 std::cout << " File write successful \n";
  
  // CLOSE file
  output.close();
  
  system("pause");
}