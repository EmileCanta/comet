#include "TFile.h"
#include "TH1.h"
#include "fstream"

void RootToAscii()
{
	TFile *input = new TFile("histogrammes/calibration_ge.root", "READ");

	TH1F *hist = (TH1F*)input->Get("hist");
	
	fstream output;
	output.open("histogrammes/calibration_ge.dat", ios::out);
	
	for(int i =0; i< 32000; i++)
	{
		output << i << " " << hist->GetBinContent(i) << endl;
	}

	output.close();
}
