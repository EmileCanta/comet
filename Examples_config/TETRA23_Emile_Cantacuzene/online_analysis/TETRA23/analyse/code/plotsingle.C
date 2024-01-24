void plotsingle()
{	
	UChar_t det;
	UInt_t energy;
	Double_t time;

	TFile *file = new TFile("/home/verney/online_analysis/TETRA23/TETRA23_RUN26.root", "READ");
	
	TTree *tree = (TTree*)file->Get("Narval_tree");

	TH1I *hist = new TH1I("hist", "hist", 32000, 0, 32000);

	tree->SetBranchAddress("Energy", &energy);
	tree->SetBranchAddress("Det_nbr", &det);
	tree->SetBranchAddress("Time", &time);

	for(int i = 0; i< tree->GetEntries(); i++)
	{
		tree->GetEntry(i);

		if(det = 15)
		{
			hist->Fill(0.32756506*energy+50.49);
		}
	}
	
	hist->Draw();

}
