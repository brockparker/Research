void readHeaderVarsFromRootFile(){
	const char inFileName[]="example_image.root";
	TFile f(inFileName);
	
	TTree *headerTree_0 = (TTree*) f.Get("headerTree_0");
	TTree *headerTree_1 = (TTree*) f.Get("headerTree_1");
	
	const int kMaxStringSize = 4096;
	
	char tempStr[kMaxStringSize];
	headerTree_0->SetBranchStatus("*",0);
	headerTree_0->SetBranchStatus("TEMPER",1);
	headerTree_0->SetBranchAddress("TEMPER",&tempStr);
	headerTree_0->GetEntry(0);
	
	char ncolStr[kMaxStringSize];
	headerTree_1->SetBranchStatus("*",0);
	headerTree_1->SetBranchStatus("CCDNCOL",1);
	headerTree_1->SetBranchAddress("CCDNCOL",&ncolStr);
	headerTree_1->GetEntry(0);

	float temp = atof(tempStr);
	float ncol = atof(ncolStr);
	cout << "Temperature: " << temp << endl;
	cout << "N cols:      " << ncol << endl;

	f.Close();
	exit();
}