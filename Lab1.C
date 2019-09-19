#include <math.h>


float func(float x,int n){
	if (n <= 3){return pow(x,n);}
	else if (n == 4){return pow(sin(x),2);}
	else if (n == 5){return 300*exp(-pow((x-2),2)/0.00001)+3*exp(x-2);}
	else {return 0;}
}

string Method(int method){
	switch (method){
		case 0:
			return "Midpt";
			break;
		case 1:
			return "Trap";
			break;
		case 2:
			return "MCI";
			break;
		default:
			return "No method";
			break;
	}
}



void Lab1(){

	// ============================================================================================
	// ==============================  Easy Function Integration  =================================
	// ============================================================================================



	TCanvas * c1 = new TCanvas("c1","three easy functions",1000,700);
	c1->Divide(3,3);
	
	int method;
	vector<double> vec[15];
	vector<double> bins[15];
	vector<double> rmsvec[5];
	vector<double> uncertainty[5];
	float rms;
	float average;
	float piece, nbins, dx, deltax;
	float precision = 0.0001;
	float xmin = 0.0;
	float xmax = 1.0;
	int powers[5] = {0,1,3,4,5};
	int index;
	float integral;

	for (int funcnum=0;funcnum<5;funcnum++){
		if (funcnum == 3){
			xmin = 0;
			xmax = M_PI/2;
		}
		if (funcnum == 4){
			xmin = 0;
			xmax = 3;
		}

		// ==========================================================================================
		// ============================  Midpoint for each function  ===============================
		// ==========================================================================================
		
		method = 0;
		index = 3*funcnum + method;
		std::cout << " function " << funcnum << ", method " << Method(method) << std::endl;
		nbins = 2;
		deltax = 0.0;	
		do {
			integral = 0.0;
			dx = (xmax-xmin)/nbins;
			for (int n=0;n<nbins;n++){
		                piece = func(xmin+(0.5+n)*dx,powers[funcnum])*dx;
               			integral += piece;
        		}
			vec[index].push_back(integral);
			bins[index].push_back(nbins);
			deltax = abs(float(vec[index][vec[index].size()-1])-float(vec[index][vec[index].size()-2]));
			nbins += 20;
		} while(deltax > precision && nbins < pow(2,10));
		
		std::cout << "" << std::endl;


		// ==============================================================================================
		// ================================  Trap for each function  ====================================
		// ==============================================================================================

		method = 1;
		index = 3*funcnum + method;
		std::cout << " function " << funcnum << ", method " << Method(method) << std::endl;
		nbins = 1;
		deltax = 0.0;
		do {
                        dx = (xmax-xmin)/nbins;
			integral =0.0;
                        for (int n=0;n<nbins;n++){
                                piece = 0.5*(func(xmin+(n*dx),powers[funcnum])+func(xmin+((n+1)*dx),powers[funcnum]))*dx;
                                integral += piece;
                        }
                        vec[index].push_back(integral);
			bins[index].push_back(nbins);
                        deltax = abs(float(vec[index][vec[index].size()-1])-float(vec[index][vec[index].size()-2]));
                        nbins += 20;
                } while(deltax > precision && nbins < pow(2,10));

		std::cout << "" << std::endl;

		
		
		
		
		
		// ==============================================================================================
		// ====================================  MCI for each function  =================================
		// ==============================================================================================
		
		method = 2;
		index = 3*funcnum + method;
		std::cout << " function " << funcnum << ", method " << Method(method) << std::endl;
		nbins=1;
		int step=0;
		float intsum=0.0;
		float sum;
		integral =0.0;
		TRandom r;
		deltax=0.0;
		do {
			step += 1;
			sum = 0.0;
                        for (int n=0;n<nbins;n++){
				piece = func(r.Uniform(xmin,xmax),powers[funcnum]);
				sum += piece;
				rms += pow(piece,2);	
			}
			
			integral = sum/nbins*(xmax-xmin);
			vec[index].push_back(integral);
			bins[index].push_back(nbins);
	
			// calculating and filling rms vector
			rmsvec[funcnum].push_back(sqrt(rms/nbins));
			
			// calculating and filling uncertainty vector
			intsum += integral;
			average = intsum/step;
			uncertainty[funcnum].push_back(sqrt((rms-pow(sum,2))/nbins)/average);
			
			// calculate Delta X for precision purposes
			deltax = abs(vec[index][vec[index].size()-1]-vec[index][vec[index].size()-2]);
				
                        nbins += 20;
                } while(deltax > precision && nbins < pow(2,10));
	}


	std::cout << "" << std::endl;

	// =====================================================================================================
	// =============================== change all vectors into arrays ======================================
	// =====================================================================================================

	double variable[15][100];
	double variable2[15][100];
	double rmsarray[5][100];
	double uncertarray[5][100];
	
	for (int num=0; num <15; num++){
		for (int size=0;size<vec[num].size();size++){
			variable[num][size] = vec[num][size];
			variable2[num][size] = bins[num][size];
		}
	}
        for (int num=0; num <5; num++){
                for (int size=0;size<rmsvec[num].size();size++){
                        rmsarray[num][size] = rmsvec[num][size];
                	uncertarray[num][size] = uncertainty[num][size];
		}

        }




	// =====================================================================================================
	// ==================================== Draw all graphs ================================================
	// =====================================================================================================

	c1->cd(1);
	TGraph *g0 = new TGraph(int(vec[0].size()),variable2[0],variable[0]);
	g0->SetTitle("f1(x)=1 midpoint");
	g0->Draw();
	c1->cd(2);
	TGraph *g1 = new TGraph(int(vec[1].size()),variable2[1],variable[1]);
	g1->SetTitle("f1(x)=1 trap");
	g1->Draw();
	c1->cd(3);
	TGraph *g2 = new TGraph(int(vec[2].size()),variable2[2],variable[2]);
	g2->SetTitle("f1(x)=1 MonteCarlo");
	g2->Draw();
	c1->cd(4);
	TGraph *g3 = new TGraph(int(vec[3].size()),variable2[3],variable[3]);
	g3->SetTitle("f2(x)=x midpoint");
	g3->Draw();
	c1->cd(5);
	TGraph *g4 = new TGraph(int(vec[4].size()),variable2[4],variable[4]);
	g4->SetTitle("f2(x)=x trap");
	g4->Draw();
	c1->cd(6);
	TGraph *g5 = new TGraph(int(vec[5].size()),variable2[5],variable[5]);
	g5->SetTitle("f2(x)=x MonteCarlo");
	g5->Draw();
	c1->cd(7);
	TGraph *g6 = new TGraph(int(vec[6].size()),variable2[6],variable[6]);
	g6->SetTitle("f3(x)=x^3 midpoint");
	g6->Draw();
	c1->cd(8);
	TGraph *g7 = new TGraph(int(vec[7].size()),variable2[7],variable[7]);
	g7->SetTitle("f3(x)=x^3 trap");
	g7->Draw();
	c1->cd(9);
	TGraph *g8 = new TGraph(int(vec[8].size()),variable2[8],variable[8]);
	g8->SetTitle("f3(x)=x^3 MonteCarlo");
	g8->Draw();



	// Draw 3 basic integration methods for 1st complicated function on c2

	TCanvas *c2 = new TCanvas("c2","f(x) = sin^2(x) for x=[0,pi/2]",1200,300);
	c2->Divide(3,1);

	c2->cd(1);
	TGraph *g9 = new TGraph(int(vec[9].size()),variable2[9],variable[9]);
	g9->GetYaxis()->SetRangeUser(0.0,1.0);
	g9->SetTitle("midpoint");
	g9->Draw();
	c2->cd(2);
	TGraph *g10 = new TGraph(int(vec[10].size()),variable2[10],variable[10]);
	g10->GetYaxis()->SetRangeUser(0.0,1.0);
	g10->SetTitle("Trapezium");
	g10->Draw();
	c2->cd(3);
	TGraph *g11 = new TGraph(int(vec[11].size()),variable2[11],variable[11]);
	g11->GetYaxis()->SetRangeUser(0.0,1.0);
	g11->SetTitle("Monte Carlo");
	g11->Draw();



	// Draw 3 basic integration methods for 2nd complicated function on c3 
	
	TCanvas *c3 = new TCanvas("c3","f(x) = 300e^{-(x-2)^2/0.00001}+3e^{x-2}",1200,300);
	c3->Divide(3,1);

	c3->cd(1);
	TGraph *g12 = new TGraph(int(vec[12].size()),variable2[12],variable[12]);
	g12->Draw();
	c3->cd(2);
	TGraph *g13 = new TGraph(int(vec[13].size()),variable2[13],variable[13]);
	g13->Draw();
	c3->cd(3);
	TGraph *g14 = new TGraph(int(vec[14].size()),variable2[14],variable[14]);
	g14->Draw();



	// Draw graphs of rms on canvas c4

	TCanvas *c4 = new TCanvas("c4","rms",1200,800);
	c4-> Divide(3,2);

	c4->cd(1);
	TGraph *grms1 = new TGraph(int(rmsvec[0].size()),variable2[2],rmsarray[0]);
	grms1->SetTitle("f(x)=1");
	grms1->Draw();
	c4->cd(2);
        TGraph *grms2 = new TGraph(int(rmsvec[1].size()),variable2[5],rmsarray[1]);
        grms2->SetTitle("f(x)=x");
	grms2->Draw();
	c4->cd(3);
        TGraph *grms3 = new TGraph(int(rmsvec[2].size()),variable2[8],rmsarray[2]);
        grms3->SetTitle("f(x)=x^3");
	grms3->Draw();
	c4->cd(4);
        TGraph *grms4 = new TGraph(int(rmsvec[3].size()),variable2[11],rmsarray[3]);
        grms4->SetTitle("f(x)=sin^2 x");
	grms4->Draw();
	c4->cd(5);
        TGraph *grms5 = new TGraph(int(rmsvec[4].size()),variable2[14],rmsarray[4]);
        grms5->SetTitle("f(x)=300e^{-(x-2)^2/0.00001}+3e^{x-2}");
	grms5->Draw();




	// Draw graphs of relative uncertainty on canvas c5

	TCanvas *c5 = new TCanvas("c5","uncertainty",1200,800);
        c5-> Divide(3,2);

	c5->cd(1);
        TGraph *gunc1 = new TGraph(int(uncertainty[0].size()),variable2[2],uncertarray[0]);
        gunc1->SetTitle("f(x)=1");
        gunc1->Draw();
        c5->cd(2);
        TGraph *gunc2 = new TGraph(int(uncertainty[1].size()),variable2[5],uncertarray[1]);
        gunc2->SetTitle("f(x)=x");
        gunc2->Draw();
        c5->cd(3);
        TGraph *gunc3 = new TGraph(int(uncertainty[2].size()),variable2[8],uncertarray[2]);
        gunc3->SetTitle("f(x)=x^3");
        gunc3->Draw();
        c5->cd(4);
        TGraph *gunc4 = new TGraph(int(uncertainty[3].size()),variable2[11],uncertarray[3]);
        gunc4->SetTitle("f(x)=sin^2 x");
        gunc4->Draw();
        c5->cd(5);
        TGraph *gunc5 = new TGraph(int(uncertainty[4].size()),variable2[14],uncertarray[4]);
        gunc5->SetTitle("f(x)=300e^{-(x-2)^2/0.00001}+3e^{x-2}");
        gunc5->Draw();

}

