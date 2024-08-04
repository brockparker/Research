#include <string.h>
#include <stdio.h>
#include "fitsio.h"

#include <iostream>
#include <sstream>
#include <cstdlib>
#include <cerrno>
#include <memory>

#include <sys/time.h>
#include <time.h>
#include <inttypes.h>
#include <fstream>
#include <unistd.h>
#include <getopt.h>    /* for getopt_long; standard getopt is in unistd.h */
#include <vector>
#include <algorithm>
#include <ctime>
#include <climits>
#include <cmath>
#include <iomanip>
#include <numeric>

#include "TFile.h"
#include "TTree.h"

#include "TCanvas.h"

#include "TH1F.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TSpline.h"

#include "globalConstants.h"

using namespace std;

struct options_t {
    bool saveSamples;
    bool saveFitsSamples;
    bool useZeroThr;
    bool autoZeroThr;
    bool saveFits;
    bool rawMode;
    bool useRunningBL;
    int blWRadius;
    bool useWholeImageAsOS;
    bool doNotSaveRoot;
    bool newRunningBL;
    bool fitCols;
    bool highCharge;
    bool dropFirst;
    bool dropNoise;
    double zeroCut;
    double gain;
    double trimQuantile;
    options_t():
        saveSamples(false), 
        saveFitsSamples(false), 
        useZeroThr(false), 
        autoZeroThr(false), 
        saveFits(false), 
        rawMode(false), 
        useRunningBL(false), 
        blWRadius(10),
        useWholeImageAsOS(false), 
        doNotSaveRoot(false),
        newRunningBL(false),
        fitCols(false),
        highCharge(false),
        dropFirst(false),
        dropNoise(false),
        zeroCut(100000),
        gain(1.0),
        trimQuantile(0.05) {};
};

// https://stackoverflow.com/questions/57151983/counting-elements-greater-than-a-number-in-vector
struct my_less_equal
{
    explicit my_less_equal(double _threshold) : threshold(_threshold){}
    bool operator()(double _other) const
    {
        return _other <= threshold;
    }
    double threshold;
};

const double kZeroThrSigmaCut = 4.0;

int deleteFile(const char *fileName){
    cout << yellow;
    cout << "Will overwrite: " << fileName << endl << endl;
    cout << normal;
    return unlink(fileName);
}

bool fileExist(const char *fileName){
    ifstream in(fileName,ios::in);

    if(in.fail()){
        //cout <<"\nError reading file: " << fileName <<"\nThe file doesn't exist!\n\n";
        in.close();
        return false;
    }

    in.close();
    return true;
}

/*========================================================
  ASCII progress bar
  ==========================================================*/
void showProgress(unsigned int currEvent, unsigned int nEvent) {

    const int nProgWidth=50;

    if ( currEvent != 0 ) {
        for ( int i=0;i<nProgWidth+8;i++)
            cout << "\b";
    }

    double percent = (double) currEvent/ (double) nEvent;
    int nBars = (int) ( percent*nProgWidth );

    cout << " |";
    for ( int i=0;i<nBars-1;i++)
        cout << "=";
    if ( nBars>0 )
        cout << ">";
    for ( int i=nBars;i<nProgWidth;i++)
        cout << " ";
    cout << "| " << setw(3) << (int) (percent*100.) << "%";
    cout << flush;

}

void printHelp(const char *exeName, bool printFullHelp=false){

    if(printFullHelp){
        cout << bold;
        cout << endl;
        cout << "This program process the raw Skipper CCD data. It computes overscan\n"
            << "mean for each sample and subtracts it line by line.\n"
            << "The output file will be a ROOT image containing a TTree with the pixels\n"
            << "values averaged over all the samples (after subtraction of the corresponding\n"
            << "overscan value). It's also possible to save an additional TTree that will\n"
            << "contain the individual values of all the samples.\n"
            << "If no <output filename> is provided a default name is used.\n";
        cout << normal;
    }
    cout << "==========================================================================\n";
    cout << red;
    cout << "\nUsage:\n";
    cout << "  "   << exeName << " <input file> -o <output filename> \n\n";
    cout << "\nOptions:\n";
    cout << "  -s for saving the individual values of all the samples.\n";
    cout << "  -i for saving a fits image file with the averaged pixels. It will be named \"proc_<input filename>.fits\".\n";
    cout << "  -S for saving a fits image file with all the samples (-i must be used). It will be named \"smpl_<input filename>.fits\".\n";
    cout << "  -n for NOT creating an output ROOT file. User must select \'-i\' or \'-S\' option.\n";
    cout << "  -w for ignoring the OS position and using the whole image as effective OS. ONLY works if image is mostly empty.\n";
    cout << "  -r raw mode. Raw image pixel values, no subtraction of the OS mean.\n";
    cout << "  -d for overwriting the output file if it exist.\n\n";
    cout << "  -z <zero threshold cut in ADC> for using pixels with skPix vale smaller than zeroCut in the OS mean.\n";
    cout << "  -a <auto zero threshold cut in OS SIGMAs> for using pixels with skPix vale smaller than OS_SD*zeroCut in the OS mean.\n";
    cout << "  -b for computing a running baseline using the empty pixels. This option mus be used together with \'-z\' or \'-a\' option. \n";
    cout << "  -R <window radius> sets the radius of the running baseline window. Only meaningful when used with the -b option. \n\n";
    cout << "  -B for computing a row-by-row baseline by fitting the 0-electron peak (new method). This option overrides most of the other options. \n";
    cout << "  -c for also fitting the column-by-column shifts when using the -B option. \n";
    cout << "  -g <gain> sets the gain assumed (ADU/e-/ssamp, default 1.0) when fitting the 0-electron peak when using the -B option.\n";
    cout << "  -Q <quantile> sets the fraction of pixels discarded when fitting the 0-electron peak when using the -B option. Use a large value if many pixels have crosstalk. \n\n";
    cout << "  -H if the first sample of a high-charge pixel is much bigger than the rest, use the first sample. \n";
    cout << "  -f ignore the first sample. \n";
    cout << "  -N ignore rows with a large number of noisy pixels (too many pixels in the tails of the 0e Gaussian). \n\n";
    cout << "  -q quiet\n";
    cout << "  -v verbose\n";
    cout << "  -h this help message\n";
    cout << normal;
    cout << blue;
    cout << "For any problems or bugs contact Javier Tiffenberg <javiert@fnal.gov>\n\n";
    cout << normal;
    cout << "==========================================================================\n\n";
}


std::string trim(const std::string& str, const std::string& whitespace = " \t\'"){ // removes leading and trailing spaces
    const auto strBegin = str.find_first_not_of(whitespace);
    if (strBegin == std::string::npos) return ""; // no content
    const auto strEnd = str.find_last_not_of(whitespace);
    const auto strRange = strEnd - strBegin + 1;
    return str.substr(strBegin, strRange);
}

inline bool isInteger(const char *c, int &n){
    if(strlen(c)==0) return false;
    string s = trim(c);
    if(((!isdigit(s[0])) && (s[0] != '-') && (s[0] != '+'))) return false ;
    std::istringstream iss(s);
    float f;
    iss >> noskipws >> f; // noskipws considers leading whitespace invalid
    if( (iss.eof() && !iss.fail()) == false ) return false; // it's not a number
    if(f != int(f) ) return false; // it's not integer

    n = int(f);
    return true;
}

class basicImageInfo_t
{
    public:
        int nsmp;
        int ssmp;
        int imcol;
        int ncol;
        int nrow;
        int npre;
        int trst;
        int runid;
        basicImageInfo_t(int nsmp=-1, int ssmp=-1, int imcol=-1, int ncol=-1, int nrow=-1, int npre=-1, int trst=-1, int runid=-1):nsmp(nsmp), ssmp(ssmp), imcol(imcol), ncol(ncol), nrow(nrow), npre(npre), trst(trst), runid(runid){};
        ~basicImageInfo_t(){};
        bool hasMinimalBasicInfo();
};

bool basicImageInfo_t::hasMinimalBasicInfo(){
    if(nsmp*ssmp*ncol*npre*trst > 0) return true;
    return false;
}

std::ostream& operator<< (std::ostream& o, const basicImageInfo_t& mt) {
    o << "nsmp  " << mt.nsmp  << endl;
    o << "ssmp  " << mt.ssmp  << endl;
    o << "imcol " << mt.imcol  << endl;
    o << "ncol  " << mt.ncol  << endl;
    o << "nrow  " << mt.nrow  << endl;
    o << "npre  " << mt.npre  << endl;
    o << "trst  " << mt.trst  << endl;
    o << "runid " << mt.runid << endl;
    return o;
}

int readIntFlagFromHdr(fitsfile *fptr, const char* keyName, const int hdu, int &intVal){
    int status = 0;   /*  CFITSIO status value MUST be initialized to zero!  */

    int hdutype;
    char keyValue[1024] = "";
    char comment[1024]  = "";

    // Get general data from ext=0 (hdu=1) header
    fits_movabs_hdu(fptr, 1, &hdutype, &status);
    // if hdu 1 does not contain image data (e.g. if this is a compressed FITS file), check the next HDU
    int naxis = 0;
    fits_get_img_dim(fptr, &naxis, &status);
    if (hdutype != IMAGE_HDU || naxis == 0) fits_movabs_hdu(fptr, 2, &hdutype, &status);
    if(status!=0) return -2;

    fits_read_keyword(fptr, keyName, keyValue, comment, &status);
    if(status==0){ // key exist
        if(isInteger(keyValue,intVal)==false) return -1;
    }
    else return -2;

    return 0;
}

int tryToReadBasicInfoFromHdr(fitsfile *fptr, const int hdu, basicImageInfo_t &bii){
    int status = 0;   /*  CFITSIO status value MUST be initialized to zero!  */

    int hdutype;
    char keyValue[1024] = "";
    char comment[1024]  = "";

    // =========================================== 
    // Try first to get data from hdu header and 
    // then try to get the remaining possibly missing 
    // info from ext=0 (hdu=1).
    // ===========================================
    std::vector<int> hdusToRead;
    if(hdu != 1) hdusToRead.push_back(hdu);
    hdusToRead.push_back(1);

    for (unsigned int i = 0; i < hdusToRead.size(); ++i)
    {
        fits_movabs_hdu(fptr, hdusToRead[i], &hdutype, &status);
        if(status!=0) return -2;

        if(bii.nsmp<0){
            fits_read_keyword(fptr, "NSAMP", keyValue, comment, &status);
            if(status==0){ // key exist
                if(isInteger(keyValue,bii.nsmp)==false) return -1;
            }
            status = 0;
        }

        if(bii.ssmp<0){
            fits_read_keyword(fptr, "SSAMP", keyValue, comment, &status);
            if(status==0){ // key exist
                if(isInteger(keyValue,bii.ssmp)==false) return -1;
            }
            status = 0;
        }

        if(bii.imcol<0){
            fits_read_keyword(fptr, "NCOL", keyValue, comment, &status);
            if(status==0){ // key exist
                if(isInteger(keyValue,bii.imcol)==false) return -1; 
            }
            status = 0;
        }

        if(bii.ncol<0){
            fits_read_keyword(fptr, "CCDNCOL", keyValue, comment, &status);
            if(status==0){ // key exist
                if(isInteger(keyValue,bii.ncol)==false) return -1; 
            }
            status = 0;
        }

        if(bii.nrow<0){
            fits_read_keyword(fptr, "CCDNROW", keyValue, comment, &status);
            if(status==0){ // key exist
                if(isInteger(keyValue,bii.nrow)==false) return -1;
            }
            status = 0;
        }

        if(bii.npre<0){
            fits_read_keyword(fptr, "CCDNPRES", keyValue, comment, &status);
            if(status==0){ // key exist
                if(isInteger(keyValue,bii.npre)==false) return -1;
            }
            status = 0;
        }

        // No-critical variables for processing
        if(bii.trst<0){
            fits_read_keyword(fptr, "TRACESTP", keyValue, comment, &status);
            if(status==0){ // key exist
                if(isInteger(keyValue,bii.trst)==false) return -1;
            }
            status = 0;
        }

        if(bii.runid<0){
            fits_read_keyword(fptr, "RUNID", keyValue, comment, &status);
            isInteger(keyValue,bii.runid); //no check..
            status = 0;
        }

    }
    fits_movabs_hdu(fptr, hdu, &hdutype, &status);

    return 0;
}

int fitsHeaderToTree(fitsfile *fptr, TFile *outRootFile){
    int status = 0;   /*  CFITSIO status value MUST be initialized to zero!  */
    int nhdu   = 0;
    outRootFile->cd();
    int hdutype;
    const int maxStringSize = 256;
    char keyName[maxStringSize];
    char keyValue[maxStringSize];
    char comment[maxStringSize];

    fits_get_num_hdus(fptr, &nhdu, &status); // get the number of HDUs
    int headerIndex = 0;
    for(int eI=1; eI<=nhdu; ++eI){  /* Main loop through each extension */
        long totpix = 0;
        int hdutype, bitpix, naxis = 0;
        long naxes[9] = {1, 1, 1, 1, 1, 1, 1, 1, 1};
        fits_movabs_hdu(fptr, eI, &hdutype, &status);
        fits_get_img_param(fptr, 9, &bitpix, &naxis, naxes, &status); /* get image dimensions and total number of pixels in image */
        totpix = naxes[0] * naxes[1] * naxes[2] * naxes[3] * naxes[4] * naxes[5] * naxes[6] * naxes[7] * naxes[8];

        if (hdutype == IMAGE_HDU && naxis > 0 && totpix > 0){
            int nKeys = 0;
            fits_get_hdrspace(fptr, &nKeys, 0, &status);
            std::vector<string> vKeyName;
            std::vector<string> vKeyValue;
            std::vector<string> vComment;
            for (int i = 0; i < nKeys; ++i){
                fits_read_keyn(fptr, i, keyName, keyValue, comment, &status);
                vKeyName.push_back(keyName);
                vKeyValue.push_back(keyValue);
                if(vKeyValue.back()[0] == '\'') vKeyValue.back().erase(0,1);
                if(vKeyValue.back()[vKeyValue.back().size()-1] == '\'') vKeyValue.back().erase(vKeyValue.back().size()-1,1);
                vComment.push_back(comment);
            }
            ostringstream headerTreeName;
            headerTreeName << "headerTree_" << headerIndex; 
            TTree headerTree(headerTreeName.str().c_str(),headerTreeName.str().c_str());
            for (int i = 0; i < nKeys; ++i){
                headerTree.Branch((vKeyName[i]).c_str(),(void*)(vKeyValue[i].c_str()),"string/C",maxStringSize);
            }
            headerTree.Fill();
            headerTree.Write();  

            headerIndex++;
        }
    }
    return 0;
}


int procSkipperImage(const char *inFile, const string outFileBaseName, basicImageInfo_t bii, const options_t opts, const string outDir=""){

    fitsfile *infptr;

    int status = 0;  
    int nhdu = 0;
    long totpix = 0;

    fits_open_file(&infptr, inFile, READONLY, &status); /* Open the input file */
    if (status != 0) return(status);

    fits_get_num_hdus(infptr, &nhdu, &status); // get the number of HDUs
    unsigned int nImHdu = nhdu;
    for(int eI=1; eI<=nhdu; ++eI){  /* Loop each extension to see how many have images */

        int hdutype, bitpix, naxis = 0;
        long naxes[9] = {1, 1, 1, 1, 1, 1, 1, 1, 1};
        fits_movabs_hdu(infptr, eI, &hdutype, &status);
        fits_get_img_param(infptr, 9, &bitpix, &naxis, naxes, &status); /* get image dimensions and total number of pixels in image */
        totpix = naxes[0] * naxes[1] * naxes[2] * naxes[3] * naxes[4] * naxes[5] * naxes[6] * naxes[7] * naxes[8];

        if (hdutype != IMAGE_HDU || naxis == 0 || totpix == 0){
            /* not an image hdu */
            nImHdu--;
        }
    }

    TCanvas *c1;
    if (gVerbosity>1 && opts.newRunningBL) {
        c1 = new TCanvas("c1", "c1", 800,600);
        c1->Print((outDir+outFileBaseName+".pdf[").c_str());
    }

    TFile *outRootFile;
    if(!opts.doNotSaveRoot) outRootFile = new TFile((outDir+outFileBaseName+".root").c_str(),"RECREATE");

    Int_t x;
    Int_t y;
    Int_t ohdu;
    Int_t runID;
    Double_t pix;

    TTree skPixTree("skPixTree", "skPixTree");
    skPixTree.Branch("x",     &x,     "x/I");
    skPixTree.Branch("y",     &y,     "y/I");
    skPixTree.Branch("ohdu",  &ohdu,  "ohdu/I");
    skPixTree.Branch("runID", &runID, "runID/I");
    skPixTree.Branch("pix",   &pix,   "pix/D");

    Double_t *pixExt = new Double_t[nImHdu];
    ostringstream pixExtArg;
    pixExtArg << "pix[" << nImHdu << "]/D";
    TTree skTablePixTree("skTablePixTree",       "skTablePixTree");
    skTablePixTree.Branch("x",     &x,           "x/I");
    skTablePixTree.Branch("y",     &y,           "y/I");
    skTablePixTree.Branch("runID", &runID,       "runID/I");
    skTablePixTree.Branch("pix",   &(pixExt[0]), pixExtArg.str().c_str());

    int kMaxNSpl = -1;
    Int_t nSpl;
    Double_t *splPix;

    TTree splPixTree("splPixTree", "splPixTree");
    splPixTree.Branch("x",      &x,     "x/I");
    splPixTree.Branch("y",      &y,     "y/I");
    splPixTree.Branch("ohdu",   &ohdu,  "ohdu/I");
    splPixTree.Branch("runID",  &runID, "runID/I");
    splPixTree.Branch("nSpl",   &nSpl,  "nSpl/I");
    splPixTree.Branch("skPix",  &pix,   "skPix/D");
    // this branch cannot be initialized yet, because we don't know nsamp: splPixTree.Branch("pix",    &(splPix[0]), "pix[nSpl]/D");

    Int_t    spl;
    Double_t osm;
    TTree osMeanTree("osMeanTree", "osMeanTree");
    osMeanTree.Branch("osm",   &osm,   "osm/D");
    osMeanTree.Branch("y",     &y,     "y/I");
    osMeanTree.Branch("spl",   &spl,   "spl/I");
    osMeanTree.Branch("ohdu",  &ohdu,  "ohdu/I");
    osMeanTree.Branch("runID", &runID, "runID/I");

    if(nhdu>0) readIntFlagFromHdr(infptr, "RUNID", 1, runID); // Read runID from header

    fitsfile *outMeanfptr;  
    fitsfile *outSmplfptr;  
    if(opts.saveFits){
        std::string outMeanFitsFile = outDir+"proc_"+outFileBaseName+".fits";
        fits_create_file(&outMeanfptr, outMeanFitsFile.c_str(), &status);/* Create the output file */
        if (status != 0) return(status);

        if(opts.saveFitsSamples){
            std::string outSmplFitsFile = outDir+"smpl_"+outFileBaseName+".fits";
            fits_create_file(&outSmplfptr, outSmplFitsFile.c_str(), &status);/* Create the output file */
            if (status != 0) return(status);
        } 
    } 

    std::vector<double*> vProcessedImg;
    std::vector<long>    vProcessedImgNPix;
    std::vector<int>     vProcessedImgNCols;

    ostringstream autoThrSummaryOSS;
    autoThrSummaryOSS << setw(5) <<  "HDU" << setw(12) << "Peak" << setw(12) << "SD" << setw(12) << "Cut" << setw(17) << "\%of zero pix" << endl;

    ohdu = 0; //first image HDU is 1
    for(int eI=1; eI<=nhdu; ++eI){  /* Main loop through each extension */
        //FITS variables
        long fpixel[2] = {1,1};
        double nulval = 0.;
        int anynul = 0;
        int hdutype, bitpix, naxis = 0;
        long naxes[9] = {1, 1, 1, 1, 1, 1, 1, 1, 1};

        fits_movabs_hdu(infptr, eI, &hdutype, &status);

        fits_get_img_param(infptr, 9, &bitpix, &naxis, naxes, &status); /* get image dimensions and total number of pixels in image */
        totpix = naxes[0] * naxes[1] * naxes[2] * naxes[3] * naxes[4] * naxes[5] * naxes[6] * naxes[7] * naxes[8];

        if (hdutype != IMAGE_HDU || naxis == 0 || totpix == 0){ /* Don't process if it's not an image HDU */
            /* ignore */
            /*
               if(opts.saveFits){
               fits_copy_hdu(infptr, outMeanfptr, 0, &status);
               if (status != 0) return(status);
               if(opts.saveFitsSamples){
               fits_copy_hdu(infptr, outSmplfptr, 0, &status);
               if (status != 0) return(status); 
               }
               }
               */
        }
        else{
            ohdu++;
            int headerStatus = tryToReadBasicInfoFromHdr(infptr, eI, bii);
            if(headerStatus != 0){
                cerr << red << "\n\nERROR:";
                cerr << "Could not retrieve basic information from image header\nWill not continue.\n" << normal << endl;
                return -100;
            }

            /* create output image */

            const int nRows = naxes[1];
            int nSamples;
            int imCols;
            if (bii.imcol==naxes[0]) { //this must be an already-processed FITS file - ignore the NSAMP value in the header
                imCols = naxes[0];
                nSamples = 1;
            } else {
                nSamples  = bii.nsmp;
                imCols    = naxes[0]/nSamples;
            }

            if (kMaxNSpl<0) { //if we haven't initialized splPix, do it now
                kMaxNSpl = nSamples;
                splPix = new Double_t[kMaxNSpl];
                splPixTree.Branch("pix",    &(splPix[0]), "pix[nSpl]/D");
            }

            if(nSamples > kMaxNSpl){ //this would happen if nSamples varies between HDUs for some reason
                cerr << red << "\nERROR: nSamples > kMaxNSpl !!!\n\nWill not continue.\n\n;" <<normal << endl;
                break;
            }
            long totpixSkp = totpix/nSamples; //number of pixels in the output image

            //define and initialize output data structures
            Double_t* fullOutArray; //baseline-subtracted value for each input pixel
            double*   outMeanArray      = new double[totpixSkp]; //the output value (baseline-subtracted and averaged) for each pixel
            double*   rblMeanArray; // Output mean array corrected using running baseline 

            fill(outMeanArray, outMeanArray+totpixSkp, 0);
            if(opts.useRunningBL){
                rblMeanArray = new double[totpixSkp];
                fill(rblMeanArray, rblMeanArray+totpixSkp, 0);
            }  

            double *osMeanV     = new double[nRows*nSamples]; //OS mean for each sample of each row

            if (opts.newRunningBL || opts.rawMode) { //new baseline subtraction algorithm, or no subtraction
                printf("%d samples, %d cols, %d rows\n", nSamples, imCols, nRows);

                const int chunkPixels = 10000000;
                const int chunkRows = max(1, chunkPixels/(nSamples*imCols)); //number of rows to read at a time

                //we don't use osMeanV, so fill with dummy values
                fill(osMeanV, osMeanV+(nRows*nSamples), 0);

                //for now: if we need to write samples, just copy the raw samples
                if (opts.saveFitsSamples || opts.saveSamples) {
                    double* inArray = new double[totpix]; //all pixels in the input image
                    //read the whole HDU into inArray
                    fits_read_pix(infptr, TDOUBLE, fpixel, totpix, &nulval, inArray, &anynul, &status);
                    fullOutArray = new Double_t[totpix]; 
                    memcpy(fullOutArray, inArray, totpix*sizeof(double));
                    delete[] inArray;
                }

                const int skipRows = 10; // skip the first rows, we expect a transient at start of readout
                int skipCols = 10; // skip the first columns, we expect a transient after the vertical shift
                if (imCols<15) skipCols = 5;
                //TODO: check against image dimensions

                float trimQuantile = opts.trimQuantile; //how much of the low-end tail to discard - this must be larger than the expected crosstalk fraction, smaller than the expected size of the 0 e- peak
                float gain = bii.ssmp*opts.gain;
                float fitRange = 0.7*gain; //how many ADUs to use in the initial fit range - this must be smaller than the expected gain, large enough to capture some of the shape of the Gaussian
                float extraRange = 2*gain; //how many ADUs to use in the histogram range, beyond the initial fit range - this can be large

                /*
                   double* sampleMeanArray = new double[nSamples](); //average value for each sample: this measures the transient after each horizontal shift
                   int nPix = 0;
                   for (int r = skipRows; r < nRows; ++r){
                   for (int c = skipCols; c < imCols; ++c){
                   for (int s = 0; s < nSamples; ++s){ // loop on samples
                   sampleMeanArray[s] += inArray[(r*imCols+c)*nSamples+s];
                   }
                   nPix++;
                   }
                   }
                   for (int s = 0; s < nSamples; ++s){ // loop on samples
                   sampleMeanArray[s] /= nPix;
                   }
                   */

                double* rawMeanArray = new double[totpixSkp]; //the average value for each pixel
                fill(rawMeanArray, rawMeanArray+totpixSkp, 0);

                {//calculate pixel averages
                    int nHighCharge = 0; //number of pixels where the high-charge override was applied
                    double* inChunkArray = new double[chunkRows*imCols*nSamples]; //all pixels in the input image
                    int firstRowThisChunk = 0;
                    while (firstRowThisChunk < nRows) {
                        int rowsRemaining = nRows - firstRowThisChunk;

                        int rowsThisChunk = chunkRows;
                        if (rowsRemaining<rowsThisChunk) rowsThisChunk = rowsRemaining;

                        long pixelsThisChunk = rowsThisChunk*imCols*nSamples;
                        fpixel[1] = firstRowThisChunk+1;
                        //read the rows
                        fits_read_pix(infptr, TDOUBLE, fpixel, pixelsThisChunk, &nulval, inChunkArray, &anynul, &status);
                        for (int iRowThisChunk = 0; iRowThisChunk < rowsThisChunk; ++iRowThisChunk){
                            int r = firstRowThisChunk + iRowThisChunk;

                            for (int c = 0; c < imCols; ++c){
                                if (opts.dropFirst) {
                                    rawMeanArray[r*imCols+c] = accumulate(inChunkArray+(iRowThisChunk*imCols+c)*nSamples+1, inChunkArray+(iRowThisChunk*imCols+c+1)*nSamples, 0.0);
                                    rawMeanArray[r*imCols+c] /= (nSamples-1);
                                } else {
                                    rawMeanArray[r*imCols+c] = accumulate(inChunkArray+(iRowThisChunk*imCols+c)*nSamples, inChunkArray+(iRowThisChunk*imCols+c+1)*nSamples, 0.0);
                                    rawMeanArray[r*imCols+c] /= nSamples;
                                }

                                //compare the first sample to the pixel mean
                                //if the first sample is much bigger than the mean, use the first sample
                                if (opts.highCharge) {
                                    if (inChunkArray[(iRowThisChunk*imCols+c)*nSamples]-rawMeanArray[r*imCols+c] > 100*gain) {
                                        //printf("%f\n",inChunkArray[(iRowThisChunk*imCols+c)*nSamples]-rawMeanArray[r*imCols+c]);
                                        rawMeanArray[r*imCols+c] = inChunkArray[(iRowThisChunk*imCols+c)*nSamples];
                                        nHighCharge++;
                                    }
                                }
                            }
                        }
                        firstRowThisChunk += rowsThisChunk;
                    }
                    if (nHighCharge!=0) {
                        cout << "nHighCharge " << nHighCharge << endl;
                    }
                    fpixel[1] = 1;//reset to the start of the image
                    delete[] inChunkArray;
                }

                if (!opts.rawMode) {

                    //calculate image mean
                    double imageMean = 0.0;
                    for (int r = 0; r < nRows; ++r){
                        for (int c = 0; c < imCols; ++c){
                            imageMean += rawMeanArray[r*imCols+c];
                        }
                    }
                    imageMean /= (nRows*imCols);


                    Double_t *fittedRows = new Double_t[nRows];
                    //Double_t *fittedRows = (Double_t *)malloc(nRows*sizeof(Double_t));
                    Double_t *rowBaselines = new Double_t[nRows];
                    Double_t *rowErrs = new Double_t[nRows];
                    Double_t *rowBlErrs = new Double_t[nRows];

                    char title[100];

                    const int minPixToFit = 400;
                    const int minPixFor0eFit = 10;
                    int nRowsToFit = 1 + minPixToFit/(imCols - skipCols);
                    int nPixToFit = nRowsToFit*(imCols - skipCols);
                    int nGoodRows = 0;
                    int nBadRows = 0;

                    for (int r = 0; r < nRows-nRowsToFit+1; ++r){
                        Double_t *rowVals = new Double_t[nPixToFit];
                        for (int dr = 0; dr < nRowsToFit; ++dr){
                            for (int c = skipCols; c < imCols; ++c){
                                rowVals[dr*(imCols - skipCols) + c-skipCols] = rawMeanArray[(r+dr)*imCols+c];
                            }
                        }
                        sort(rowVals, rowVals + nPixToFit);
                        int lowEnd = nPixToFit * trimQuantile;
                        //printf("row %d, %d-th value: %f\n", r, lowEnd, rowVals[lowEnd]);

                        sprintf(title, "HDU %d, row %d", ohdu, r);
                        TH1 * hrow = new TH1F("hrow", title, 500, rowVals[lowEnd] - extraRange, rowVals[lowEnd] + fitRange + extraRange);
                        for (int c = 0; c < nPixToFit; ++c){
                            hrow->Fill(rowVals[c]);
                        }
                        double rangeMin, rangeMax;
                        rangeMin = rowVals[lowEnd]-0.05*gain;
                        rangeMax = max(rangeMin+fitRange, rowVals[lowEnd+minPixFor0eFit])+0.05*gain;
                        if (gVerbosity>1) printf("row %d, mean: %f, range: %f %f\n", r, hrow->GetMean(), rangeMin, rangeMax);

                        TFitResultPtr result = hrow->Fit("gaus", "QSL", "", rangeMin, rangeMax);
                        if (gVerbosity>1) printf("first %f %f %f\n", result->Parameter(0), result->Parameter(1), result->Parameter(2));
                        int fitStatus = result;
                        double mean, sigma;
                        if (fitStatus==0) {
                            mean = result->Parameter(1);
                            sigma = result->Parameter(2);

                            result = hrow->Fit("gaus", "QSL", "", mean-2*sigma, mean+2*sigma);
                            fitStatus = result;
                            if (gVerbosity>1 && fitStatus==0) printf("old %f %f new %f %f\n", mean, sigma, result->Parameter(1), result->Parameter(2));
                        } else {
                            if (gVerbosity>1) printf("failed initial fit\n");
                        }
                        if (fitStatus==0) {
                            mean = result->Parameter(1);
                            sigma = result->Parameter(2);

                            if (opts.dropNoise) {
                                //count number of "noisy" pixels (too far away from the mean)
                                // the cut is very ad-hoc, but requiring <5% of pixels to be outside of +/- gain/3 is equiv to a Gaussian noise of gain/6.
                                Double_t *part1 = partition_point(rowVals,rowVals + nPixToFit, my_less_equal(mean - gain/3));
                                Double_t *part2 = partition_point(rowVals,rowVals + nPixToFit, my_less_equal(mean + gain/3));
                                Double_t *part3 = partition_point(rowVals,rowVals + nPixToFit, my_less_equal(mean + gain/2));
                                int tooLow = part1-rowVals;
                                int tooHigh = part3-part2;
                                if (tooLow+tooHigh > 0.05*nPixToFit) {
                                    //printf("%d: %d %d, %f \n", r, tooLow, tooHigh, ((double)(tooLow+tooHigh))/nPixToFit);
                                    fitStatus = -1;
                                    if (gVerbosity>1) printf("failed noise check\n");
                                }
                            }

                            if (mean<rowVals[0] || mean>rowVals[lowEnd]+fitRange) {
                                fitStatus=-1;
                                if (gVerbosity>1) printf("failed sanity check\n");
                            }
                        } else {
                            if (gVerbosity>1) printf("failed refit\n");
                        }
                        //}
                        //TODO: more checks - integral of the Gaussian?

                        if (fitStatus==0 && result->Parameter(1)>mean-sigma && result->Parameter(1)<mean+sigma) {//check that the new fitted mean is close to the first one
                            //printf("row %d, fit status %d, params: amplitude=%f, mean=%f, sigma=%f\n", r, fitStatus, result->Parameter(0), result->Parameter(1), result->Parameter(2));
                            fittedRows[nGoodRows] = r+0.5*(nRowsToFit-1);
                            rowBaselines[nGoodRows] = result->Parameter(1);
                            rowErrs[nGoodRows] = 0;
                            rowBlErrs[nGoodRows] = result->ParError(1);
                            nGoodRows++;
                        } else {
                            nBadRows++;
                            if (gVerbosity>1) printf("could not fit row %d\n", r);
                        }
                        if (gVerbosity>1 && (r<5 || fitStatus!=0)) {
                            hrow->Draw();
                            c1->Print((outDir+outFileBaseName+".pdf").c_str());
                        }
                        delete hrow;
                        delete[] rowVals;
                }
                printf("fitted %d rows, failed to fit %d rows\n", nGoodRows, nBadRows);

                TGraphErrors *rowGraph;
                TSpline3 *rowSpline;

                bool useRowFit = (nGoodRows>=4);
                if (useRowFit) {
                    rowGraph = new TGraphErrors(nGoodRows, fittedRows, rowBaselines, rowErrs, rowBlErrs);
                    if (gVerbosity>1) {
                        rowGraph->Draw("ALP");
                        c1->Print((outDir+outFileBaseName+".pdf").c_str());
                    }
                    rowSpline = new TSpline3("rowBaseline", rowGraph, "b2e2");
                    if (gVerbosity>1) {
                        rowSpline->Draw("CP");
                        c1->Print((outDir+outFileBaseName+".pdf").c_str());
                    }
                } else {
                    printf("was not able to fit enough rows; subtracting the image mean instead of doing a row-by-row baseline subtraction\n");
                }

                for (int r = 0; r < nRows; ++r){
                    double rowMean = imageMean;
                    for (int c = 0; c < imCols; ++c){
                        if (useRowFit) {
                            //rowMean = rowSpline->Eval(r);
                            rowMean = rowSpline->Eval(r + (c-0.5*imCols)/imCols); //interpolate, assuming that the row-by-row fit is correct for the center of each row
                        }
                        outMeanArray[r*imCols+c] = rawMeanArray[r*imCols+c] - rowMean;
                        if (opts.saveFitsSamples || opts.saveSamples) for (int sm = 0; sm < nSamples; ++sm) fullOutArray[(r*imCols+c)*nSamples + sm] -= rowMean; //For now just subtract the mean of the row regardless of the sample number
                    }
                }
                if (useRowFit) {
                    delete rowSpline;
                    delete rowGraph;
                }

                if (opts.fitCols && useRowFit && nRows>=200) {
                    double*   rowSubtractedArray      = new double[totpixSkp];
                    std::copy(outMeanArray, outMeanArray+totpixSkp, rowSubtractedArray);
                    int nGoodCols = 0;
                    int nBadCols = 0;
                    Double_t *fittedCols = new Double_t[imCols];
                    Double_t *colBaselines = new Double_t[imCols];
                    Double_t *colErrs = new Double_t[imCols];
                    Double_t *colBlErrs = new Double_t[imCols];
                    for (int c = 0; c < imCols; ++c){
                        Double_t *colVals = new Double_t[nRows-skipRows];

                        for (int r = skipRows; r < nRows; ++r){
                            colVals[r-skipRows] = rowSubtractedArray[r*imCols+c];
                        }

                        sort(colVals, colVals + nRows-skipRows);
                        int lowEnd = (nRows-skipRows) * trimQuantile;
                        //printf("col %d, %d-th value: %f\n", c, lowEnd, colVals[lowEnd]);

                        TH1 * hcol = new TH1F("hcol", "hcol", 500, colVals[lowEnd] - extraRange, colVals[lowEnd] + fitRange + extraRange);
                        for (int r = 0; r < nRows-skipRows; ++r){
                            hcol->Fill(colVals[r]);
                        }
                        TFitResultPtr result = hcol->Fit("gaus", "QSL", "", colVals[lowEnd], colVals[lowEnd]+fitRange);
                        int fitStatus = result;
                        if (fitStatus==0) {
                            double mean = result->Parameter(1);
                            double sigma = result->Parameter(2);
                            result = hcol->Fit("gaus", "QSL", "", mean-2*sigma, mean+2*sigma);
                            fitStatus = result;
                        }
                        //printf("fit status %d, params: amplitude=%f, mean=%f, sigma=%f\n", fitStatus, result->Parameter(0), result->Parameter(1), result->Parameter(2));
                        //if (c<skipCols*2) {
                        //hcol->Draw();
                        //c1->Print((outDir+outFileBaseName+".pdf").c_str());
                        //}
                        //TODO: more checks - integral of the Gaussian?
                        if (fitStatus==0) {
                            fittedCols[nGoodCols] = c;
                            colBaselines[nGoodCols] = result->Parameter(1);
                            colErrs[nGoodCols] = 0;
                            colBlErrs[nGoodCols] = result->ParError(1);
                            nGoodCols++;
                        } else {
                            nBadCols++;
                            if (gVerbosity>1) printf("could not fit col %d\n", c);
                        }
                        delete hcol;
                        delete[] colVals;
                    }
                    printf("fitted %d cols, failed to fit %d cols\n", nGoodCols, nBadCols);

                    TGraphErrors *colGraph = new TGraphErrors(nGoodCols, fittedCols, colBaselines, colErrs, colBlErrs);
                    if (gVerbosity>1) {
                        colGraph->Draw("ALP");
                        colGraph->GetXaxis()->SetRangeUser(-1,20);
                        c1->Print((outDir+outFileBaseName+".pdf").c_str());
                    }
                    TSpline3 *colSpline = new TSpline3("colBaseline", colGraph);
                    if (gVerbosity>1) {
                        colSpline->Draw("CP");
                        c1->Print((outDir+outFileBaseName+".pdf").c_str());
                    }
                    TSpline5 *colSpline5 = new TSpline5("colBaseline5", colGraph);
                    if (gVerbosity>1) {
                        colSpline5->Draw("CP");
                        c1->Print((outDir+outFileBaseName+".pdf").c_str());
                    }

                    //double*   colSubtractedArray      = new double[totpixSkp];
                    for (int r = 0; r < nRows; ++r){
                        for (int c = 0; c < imCols; ++c){
                            outMeanArray[r*imCols+c] = rowSubtractedArray[r*imCols+c] - colSpline->Eval(c);
                            //colSubtractedArray[r*imCols+c] = rawMeanArray[r*imCols+c] - colSpline->Eval(c);
                            if (opts.saveFitsSamples || opts.saveSamples) for (int sm = 0; sm < nSamples; ++sm) fullOutArray[(r*imCols+c)*nSamples + sm] -= colSpline->Eval(c); //For now just subtract the mean of the row regardless of the sample number
                        }
                    }
                    delete colSpline;
                    delete colGraph;

                    delete[] fittedCols;
                    delete[] colBaselines;

                    delete[] rowSubtractedArray;
                    //delete[] colSubtractedArray;
                }

                //double maxRow = *(max_element(fittedRows, fittedRows+nGoodRows));
                //double minRow = *(min_element(fittedRows, fittedRows+nGoodRows));

                //TH1 * hrunning = new TH1F("hcol", "hcol", 200, binVal-gain, binVal+gain);
                delete[] fittedRows;
                delete[] rowBaselines;
            }



            //for (int r = 0; r < nRows; ++r){
            //    for (int c = 0; c < imCols; ++c){
            //        for (int s = 0; s < nSamples; ++s){
            //            fullOutArray[nSamples*(r*imCols+c)+s] = inArray[nSamples*(r*imCols+c)+s] - sampleMeanArray[s];
            //        }
            //    }
            //}

            if (opts.rawMode) {
                for (int r = 0; r < nRows; ++r){
                    for (int c = 0; c < imCols; ++c){
                        outMeanArray[r*imCols+c] = rawMeanArray[r*imCols+c];
                    }
                }
            }

            delete[] rawMeanArray;
            //delete[] sampleMeanArray;

        } else { //old baseline subtraction algorithms
            const int trst      = bii.trst; //tracer pixels??
            const int nTrcrPix  = (bii.trst>0)? int(ceil(bii.ncol/bii.trst+0.5)) + (int(ceil(bii.ncol/bii.trst+0.5))&1) : 0; // round up to the nearest even number

            const int nPreScan  = bii.npre;
            const int nOS       = (opts.useWholeImageAsOS)? imCols - nPreScan : imCols - (bii.ncol + nTrcrPix)/2 - nPreScan -2; //The trailing '-2' is to ignore the first OS col (it may have charge due CTI), and to account for the fact that NPRESCAN is typically 1 less than the number of prescan pixels
            const int lOsStart  = imCols - nOS;
            const int nMeanTrim = (opts.useWholeImageAsOS)? nOS/10 : 4; // discard the 10% tails when computing stable mean : hardcoded 4 pix when using OS

            double zeroThr = opts.zeroCut; //this is overwritten if autoZeroThr is used

            if(nOS<0){
                cerr << red << "\nERROR: no overscan !!!\n\nOld baseline subtraction algorithm will not work. Try -B.\nWill not continue.\n\n;" <<normal << endl;
                break;
            }

            double* inArray = new double[totpix]; //all pixels in the input image
            //read the whole HDU into inArray
            fits_read_pix(infptr, TDOUBLE, fpixel, totpix, &nulval, inArray, &anynul, &status);

            fullOutArray = new Double_t[totpix]; 

            //first pass: use the trimmed mean of nontracer OS pixels in this row
            double osAuxV[nOS];
            for (int r = 0; r < nRows; ++r){ // loop on rows
                for (int s = 0; s < nSamples; ++s){ // loop on samples
                    /* compute stable mean for the OS pixels */
                    int nTrcOS = 0;
                    for (int o = 0; o < nOS; ++o){
                        if( ((o+lOsStart)%trst) == (trst-1) ){
                            osAuxV[o] = -1000; // negative value for the tracer pixels so they are ignored in the robust mean 
                            ++nTrcOS;
                        }
                        else osAuxV[o] = inArray[nSamples*(o+r*imCols+lOsStart)+s];
                    } 

                    sort(osAuxV, osAuxV+nOS);
                    osMeanV[r*nSamples+s] = (nOS>2*nMeanTrim)? accumulate(osAuxV+nMeanTrim+nTrcOS, osAuxV+(nOS-nMeanTrim), 0.0) : accumulate(osAuxV, osAuxV, 0.0);
                    //if we can, trim off the tracer pixels and the nMeanTrim largest and smallest values; otherwise, ???
                    osMeanV[r*nSamples+s] /= (nOS>4)? (nOS-nMeanTrim*2-nTrcOS) : (nOS-nTrcOS);
                }
            }

            if(opts.useZeroThr || opts.autoZeroThr){
                vector< vector<double> > lSkOsFirstPass(nRows, std::vector<double>(nOS,0) );

                for (int r = 0; r < nRows; ++r){
                    // subtract OS mean for each sample in the OS
                    for (int o = 0; o < nOS; ++o){
                        for (int s = 0; s < nSamples; ++s){
                            lSkOsFirstPass[r][o] += inArray[nSamples*(o+r*imCols+lOsStart)+s] - osMeanV[r*nSamples+s];
                        }
                        lSkOsFirstPass[r][o] /= nSamples;
                    }
                }

                if(opts.autoZeroThr){
                    vector<double> osPix(nOS*nRows);
                    auto p = 0;
                    int nTrcOS = 0;
                    for (int r = 0; r < nRows; ++r){
                        for (long o = 0; o < nOS; ++o){
                            if( ((o+lOsStart)%trst) == (trst-1) ){
                                ++nTrcOS;
                            }
                            else{
                                osPix[p] = lSkOsFirstPass[r][o];
                                ++p;
                            }
                        }
                    }

                    // Find the position of the zero e- peak
                    sort(osPix.begin(), osPix.end() - nTrcOS);
                    auto kPeakFindWind = (osPix.size()-nTrcOS)/2; // Assume that at least 1/2 of the pixels in the OS are empty
                    //we find the peak by measuring the distance between the oi'th value and the (oi+window)'th value
                    //when that distance is minimized, the window is straddling the peak
                    auto minDist = 100000.0;
                    unsigned int minOi = 0;
                    for (unsigned int oi = 0; oi < osPix.size()-nTrcOS-kPeakFindWind; ++oi)
                    {
                        double dist = osPix[oi+kPeakFindWind] - osPix[oi];
                        if(dist < minDist){
                            minDist = dist;
                            minOi = oi;
                        }
                    }
                    //now we take the peak position as the mean of values inside the window
                    double zeroPeakPos = accumulate(osPix.begin() + minOi, osPix.begin() + minOi +kPeakFindWind, 0.0);
                    zeroPeakPos /= kPeakFindWind;

                    // Compute sigma from the lower half of the gaussian of the zero e- peak
                    double sqDiff = 0;
                    int nHalfGaus = 0;
                    for (unsigned int oi = max(int(osPix.size()-nTrcOS)/50, 5); oi < osPix.size()-nTrcOS; ++oi) // discard the first 2% (or 5 entries) of the dist to make SD more stable
                    {
                        if(osPix[oi]>zeroPeakPos) break;
                        sqDiff += pow(osPix[oi]-zeroPeakPos, 2);
                        ++nHalfGaus;
                    }
                    const double zeroPeakSD = sqrt(sqDiff/nHalfGaus);
                    zeroThr = zeroPeakPos + opts.zeroCut*zeroPeakSD;
                    auto nLessThanZeroCut = std::lower_bound(osPix.begin(), osPix.end(), zeroThr) - osPix.begin();
                    autoThrSummaryOSS << setw(5) <<  ohdu << setw(12) << zeroPeakPos << setw(12) << zeroPeakSD << setw(12) << zeroThr << setw(17) << nLessThanZeroCut*1.0/(osPix.size()-nTrcOS) << endl;
                    /* It would be good to make a gaussianity test. Goodnes of a fit? */
                }
                /* End of zeroThr processing */

                //recalculate mean as the mean of zero, nontracer OS pixels in this row
                for (int r = 0; r < nRows; ++r){ // loop on rows
                    for (int s = 0; s < nSamples; ++s){ // loop on samples
                        osMeanV[r*nSamples+s] = 0;

                        int lNZero = 0; //number of nonzero pixels
                        for (long o = 0; o < nOS; ++o){
                            if(((o+lOsStart)%trst) == (trst-1)) continue; // ignore tracer pixels
                            if(lSkOsFirstPass[r][o]<zeroThr){
                                osMeanV[r*nSamples+s] += inArray[nSamples*(o+r*imCols+lOsStart)+s];
                                ++lNZero;
                            }
                        }

                        osMeanV[r*nSamples+s] /= lNZero;
                    }
                }
            }

            /* Subtract OS mean for each sample */
            for (int r = 0; r < nRows; ++r){ // loop on rows
                for (int s = 0; s < nSamples; ++s){ // loop on samples
                    // subtract OS mean for each pixel on the row for the current sample
                    double m = 0;
                    double b = osMeanV[r*nSamples+s];
                    if(r>0){ // it's not the first row, so interpolate between the current and previous row OS means
                        const double prevOsMean = osMeanV[(r-1)*nSamples+s];
                        m = (osMeanV[r*nSamples+s] - prevOsMean)/(imCols-1);
                        b = prevOsMean;
                    }
                    for (int c = 0; c < imCols; ++c){ // subtract the calculated mean from the pixels of this row
                        double corrMean = m*c+b;
                        fullOutArray[nSamples*(r*imCols+c)+s] = inArray[nSamples*(r*imCols+c)+s] - corrMean;
                        //if you don't want interpolation: use osMeanV[r*nSamples+s] instead of corrMean
                        outMeanArray[r*imCols+c] += fullOutArray[nSamples*(r*imCols+c)+s];
                    }
                }
                for (int c = 0; c < imCols; ++c){
                    outMeanArray[r*imCols+c] /= nSamples;
                }

                if(gVerbosity) showProgress(r + (ohdu-1)*nRows*(1 + opts.useRunningBL), nhdu*nRows*(1 + opts.useRunningBL));
            }


            /* START of Improved Baseline subtraction */
            if(opts.useRunningBL){
                for (int r = 0; r < nRows; ++r){ // loop on rows
                    for (int s = 0; s < nSamples; ++s){ // loop on samples
                        double rowSum[imCols]; //running sum of nonzero pixels
                        double rowNZeroSum[imCols]; //running count of nonzero pixels
                        if(outMeanArray[imCols*r]<zeroThr){
                            rowSum[0] = fullOutArray[nSamples*(r*imCols)+s];
                            rowNZeroSum[0] = 1;
                        }
                        else{
                            rowSum[0] = 0;
                            rowNZeroSum[0] = 0;
                        }
                        for (int c = 1; c < imCols; ++c){ //fill the running sum and count
                            if(outMeanArray[imCols*r + c]<zeroThr){
                                rowSum[c]      = rowSum[c-1] + fullOutArray[nSamples*(r*imCols+c)+s];
                                rowNZeroSum[c] = rowNZeroSum[c-1] + 1;
                            }
                            else{
                                rowSum[c]      = rowSum[c-1];
                                rowNZeroSum[c] = rowNZeroSum[c-1];
                            }
                        }

                        for (int c = 1; c < imCols; ++c){ //use the running sum and count to calculate the mean in a window around each pixel

                            int wEnd = c + opts.blWRadius;
                            int wBeg = c - opts.blWRadius;
                            if(wEnd>=imCols){
                                wEnd = imCols-1;
                                wBeg = wEnd - opts.blWRadius*2;
                            }
                            if(wBeg<0){
                                wBeg = 0;
                                wEnd = opts.blWRadius*2;
                            }
                            const double bl = (rowSum[wEnd]-rowSum[wBeg])/(rowNZeroSum[wEnd]-rowNZeroSum[wBeg]);
                            // cout << r << "\t" << bl << "\t" << wBeg << "\t" << wEnd << "\t" << (rowNZeroSum[wEnd]-rowNZeroSum[wBeg]) << "\t" << (rowSum[wEnd]-rowSum[wBeg]) << "\t" << rowNZeroSum[wEnd] << endl;
                            fullOutArray[nSamples*(r*imCols+c)+s] -= bl; //subtract the window mean
                            //TODO: this means col 0 is not subtracted?
                            rblMeanArray[imCols*r + c] += fullOutArray[nSamples*(r*imCols+c)+s];
                        }
                        rblMeanArray[imCols*r] += fullOutArray[nSamples*(r*imCols)+s];//since col 0 is skipped
                    }
                    for (int c = 0; c < imCols; ++c) rblMeanArray[imCols*r + c] /= nSamples;

                    if(gVerbosity) showProgress(r + ((ohdu-1)*2+1)*nRows, nhdu*nRows*2);
                }
            }
            /* END of Improved Baseline subtraction   */
            delete[] inArray;
        }

        if(!opts.doNotSaveRoot) {
            if(opts.saveSamples){
                nSpl = nSamples;
                for (long j = 0; j < totpixSkp; ++j){
                    x = j%imCols;
                    y = j/imCols;
                    pix = 0;
                    for (int s = 0; s < nSamples; ++s){
                        splPix[s] = fullOutArray[j*nSamples + s];
                        pix += splPix[s];
                    }
                    pix /= nSamples;
                    splPixTree.Fill();
                }
                splPixTree.Write();
            }

            for (long j = 0; j < totpixSkp; ++j){
                x = j%imCols;
                y = j/imCols;
                if(opts.useRunningBL) pix = rblMeanArray[j];
                else pix = outMeanArray[j];
                skPixTree.Fill();
            }

            for (int r = 0; r < nRows; ++r){ // loop on rows
                for (int s = 0; s < nSamples; ++s){ // loop on samples
                    osm = osMeanV[r*nSamples+s];
                    y   = r;
                    spl = s;
                    osMeanTree.Fill(); 
                }
            }
        }

        if(opts.saveFits){ // Create fits image with the averaged pixels
            long naxesOut[9] = {1, 1, 1, 1, 1, 1, 1, 1, 1};
            std::copy(naxes, naxes+9, naxesOut);
            naxesOut[0] /= nSamples;
            fits_create_img(outMeanfptr, -32, naxis, naxesOut, &status);
            if (status != 0) return(status);

            if(opts.useRunningBL) fits_write_pix(outMeanfptr, TDOUBLE, fpixel, totpixSkp, rblMeanArray, &status);
            else fits_write_pix(outMeanfptr, TDOUBLE, fpixel, totpixSkp, outMeanArray, &status);

            /* Copy header info */
            int nkeys=0;
            fits_get_hdrspace(infptr, &nkeys, NULL, &status);
            char card[FLEN_CARD];
            for (int k = 1; k <= nkeys; ++k) {
                fits_read_record(infptr, k, card, &status);
                if (fits_get_keyclass(card) == TYP_USER_KEY) fits_write_record(outMeanfptr, card, &status);
            }

            if (status != 0) return(status);
            if(opts.saveFitsSamples){ // Create fits image with all the samples
                double*   singleSampleArray = new double[totpixSkp]; //for each pixel of the output image, just the current sample
                for (int si = 0; si < nSamples; ++si){
                    for (long j = 0; j < totpixSkp; ++j) singleSampleArray[j] = fullOutArray[j*nSamples+si];
                    fits_create_img(outSmplfptr, -32, naxis, naxesOut, &status);
                    if (status != 0) return(status);
                    fits_write_pix(outSmplfptr, TDOUBLE, fpixel, totpixSkp, singleSampleArray, &status);
                    if (status != 0) return(status);
                }
                delete[] singleSampleArray;
            }
        }

        if(opts.useRunningBL) vProcessedImg.push_back(rblMeanArray);
        else             vProcessedImg.push_back(outMeanArray);
        vProcessedImgNPix.push_back(totpixSkp);
        vProcessedImgNCols.push_back(imCols);

        if (!opts.newRunningBL || opts.saveFitsSamples) {//the new algo does not allocate this array unless needed
            delete[] fullOutArray;
        }
        delete[] osMeanV;
    }
}

if(vProcessedImgNPix.size()==0){
    cerr << "\nERROR: no images found on the FITS file.\n\n";
    return -1;
}

for (unsigned int amp = 1; amp < vProcessedImgNPix.size(); ++amp){ // Check that all the extensions (amplifiers) have the same number of pixels
    if(vProcessedImgNPix[amp]!=vProcessedImgNPix[0]){
        cerr << "\nERROR: extensions do not have the same number of pixels.\n\n";
        return -1;
    }
}

if( vProcessedImg.size() != (nImHdu) ){
    cerr << "\nERROR: the number of extensions is not \"nImHdu\".\n\n";
    cout << vProcessedImg.size() <<  " " << nImHdu << endl;
    return -1;
}
const int imCols = vProcessedImgNCols[0]; // TO DO: should check that all extensions have the same number of cols
const int nAmps  = nImHdu;
for (long j = 0; j < vProcessedImgNPix[0]; ++j){
    x = j%imCols;
    y = j/imCols;
    for (int amp = 0; amp < nAmps; ++amp){ // Loop on amplifiers (extensions)
        pixExt[amp] = vProcessedImg[amp][j];
    }
    skTablePixTree.Fill();
}

for (unsigned int amp = 0; amp < vProcessedImg.size(); ++amp) delete[] (vProcessedImg[amp]);

if (gVerbosity>1 && opts.newRunningBL) {
    c1->Print((outDir+outFileBaseName+".pdf]").c_str());
    delete c1;
}

if(!opts.doNotSaveRoot){
    skPixTree.Write();
    osMeanTree.Write();
    skTablePixTree.Write();
    fitsHeaderToTree(infptr, outRootFile);
    outRootFile->Close();
    delete outRootFile;
}

fits_close_file(infptr,   &status);
if (status != 0) return(status);  
if(opts.saveFits){  
    fits_close_file(outMeanfptr,  &status);
    if (status != 0) return(status);
    if(opts.saveFitsSamples){
        fits_close_file(outSmplfptr,  &status);
        if (status != 0) return(status);
    }  
}

delete[] pixExt;

if(gVerbosity){
    showProgress(1, 1);
}
if(opts.autoZeroThr) cout << "\n\nAuto zero threshold summary:\n" << autoThrSummaryOSS.str() << endl;

return status;
}


int checkForExistingFileAndHandle(const char* fileName, const bool owOutFileFlag, const std::string fileType = "\b"){

    if(fileExist(fileName)){ // output file
        cout << yellow << "\nThe output "<< fileType << " file exist.\n" << normal;
        if(owOutFileFlag){
            cout << yellow << "Will overwrite the output file.\n\n" << normal;
            deleteFile(fileName);
        }
        else{
            cout << yellow << "Please provide a different name or use the \"-d\" option.\n\n" << normal;
            cout << red    << "Will NOT continue.\n" << normal;
            return 5;
        }
    }
    return 0;
}

int processCommandLineArgs(const int argc, char *argv[], string &inFile, string &outFileBaseName, options_t &opts, string &outDir, basicImageInfo_t &bii){

    if(argc == 1) return 1;

    bool outFileFlag   = false;
    bool owOutFileFlag = false;
    int opt=0;
    string outFile;
    while ( (opt = getopt(argc, argv, "a:o:z:R:g:Q:windsSrbBcqHfNvh?")) != -1) {
        switch (opt) {
            case 'o':
                if(!outFileFlag){
                    outFile = optarg;
                    outFileFlag = true;
                }
                else{
                    cerr << red << "\nError, can not set more than one output file!\n\n" << normal;
                    return 2;
                }
                break;
            case 'd': /* overwrite output fie is exist */
                owOutFileFlag = true;
                break;
            case 's':
                opts.saveSamples = true;
                break;
            case 'i':
                opts.saveFits = true;
                break;
            case 'w':
                opts.useWholeImageAsOS = true;
                break;
            case 'n':
                opts.doNotSaveRoot = true;
                break;
            case 'S':
                opts.saveFitsSamples = true;
                break;
            case 'r':
                opts.rawMode = true;
                break;
            case 'b':
                opts.useRunningBL = true;
                break;
            case 'R':
                opts.blWRadius = std::strtol(optarg,0,10);
                break;
            case 'B':
                opts.newRunningBL = true;
                break;
            case 'c':
                opts.fitCols = true;
                break;
            case 'g':
                opts.gain = atof(optarg);
                break;
            case 'Q':
                opts.trimQuantile = atof(optarg);
                break;
            case 'H':
                opts.highCharge = true;
                break;
            case 'f':
                opts.dropFirst = true;
                break;
            case 'N':
                opts.dropNoise = true;
                break;
            case 'z':
                if(opts.autoZeroThr || opts.useZeroThr){
                    cerr << red << "\nError, can not set more than one zeroCut!\nWill NOT continue.\n\n"  << normal;
                    return 4;
                }
                else{
                    opts.useZeroThr = true;
                    char * e;
                    errno = 0;
                    opts.zeroCut = std::strtod(optarg, &e);
                    if (*e != '\0' /* didn't consume the entire string */ || errno != 0 /* error, overflow or underflow */ ){
                        cerr << red << "\nError reading zeroCut value.\nWill NOT continue.\n\n"  << normal;
                        return 3;
                    }
                }
                break;
            case 'a':
                if(opts.autoZeroThr || opts.useZeroThr){
                    cerr << red << "\nError, can not set more than one zeroCut!\nWill NOT continue.\n\n"  << normal;
                    return 4;
                }
                else{
                    opts.autoZeroThr = true;
                    char * e;
                    errno = 0;
                    opts.zeroCut = std::strtod(optarg, &e);
                    if (*e != '\0' /* didn't consume the entire string */ || errno != 0 /* error, overflow or underflow */ ){
                        cerr << red << "\nError reading zeroCut value.\nWill NOT continue.\n\n"  << normal;
                        return 3;
                    }
                }
                break;
            case 'v':
                gVerbosity++;
                break;
            case 'q':
                gVerbosity = 0;
                break;
            case 'h':
                return 1;
                break;
            default: /* '?' */
                cerr << red << "For help use -h option.\n\n" << normal;
                return 3;
        }
    }


    inFile="";

    if(argc-optind==0){
        cerr << red << "Error: no input file provided!\n\n" << normal;
        return 1;
    }
    else if(argc-optind>1){
        cerr << red << "Error: more than one input file provided!\n\n" << normal;
        return 1;
    }

    if(opts.useRunningBL){
        if( !(opts.useZeroThr) && !(opts.autoZeroThr) ){
            cout << red << "\nError: if Running-Baseline option \'-b\' is chosen user must also provide a zero-threshold value using \'-z\' or \'-a\' option.\n\n" << normal;
            return -108;
        }
    }

    inFile=argv[optind];
    if(!fileExist(inFile.c_str())){
        cout << red << "\nError reading input file: " << inFile <<"\nThe file doesn't exist!\n\n" << normal;
        return -1;
    }

    if(!outFileFlag){
        cout << endl;
        cout << whiteOnPurple << "Warning: output filename missing." << normal << endl;
        cout << whiteOnViolet << "Will use default output base name:" << normal;
        const int nameStart = std::max((int)(inFile.rfind("/")+1), 0);
        const int nameEnd   = std::max((int)(inFile.rfind(".")+1), 0);
        outFileBaseName = inFile.substr(nameStart, nameEnd-nameStart-1);
        outDir          = "";
        cout << whiteOnViolet << " " << outFile << normal << endl;
    }

    if(outFileFlag){
        const int nameStart = std::max((int)(outFile.rfind("/")+1), 0);
        const int nameEnd   = (outFile.substr(outFile.size()-min(5,(int)(outFile.size()))) == ".root")? outFile.size()-5 : outFile.size();
        outFileBaseName = outFile.substr(nameStart, nameEnd-nameStart);
        outDir          = outFile.substr(0, nameStart);
    }

    /* Handle the cases when output ROOT file will NOT be created */
    if( opts.doNotSaveRoot){
        if( !(opts.saveFits) && !(opts.saveFitsSamples)){
            cerr << red << "\nError: \'-n\' option must be used together with \'-i\' or \'-S\' flag\n" << normal;
            cerr << red << "Will NOT continue\n\n" << normal;
            return -102;
        }
    }

    /* Handle the cases of existing output root and/or fits files */
    // output root file
    if( !(opts.doNotSaveRoot) ) {
        int errcode = checkForExistingFileAndHandle((outDir+outFileBaseName+".root").c_str(), owOutFileFlag, "root");
        if(errcode!=0) return errcode;
    }

    // output mean fits file
    if(opts.saveFits){
        int errcode = checkForExistingFileAndHandle((outDir+"proc_"+outFileBaseName+".fits").c_str(), owOutFileFlag, "mean fits");
        if(errcode!=0) return errcode;

        // output samples fits file
        if(opts.saveFitsSamples){
            errcode = checkForExistingFileAndHandle((outDir+"smpl_"+outFileBaseName+".fits").c_str(), owOutFileFlag, "samples fits");
            if(errcode!=0) return errcode;
        }
    }

    /* Handle the cases when kAutoZeroThrFlag is used */
    if( opts.autoZeroThr){
        cout << yellow << "\nWarning: \'-a\' option, will try to compute the zero threshold. Could fail miserably, use with caution.\n" << normal;
    }

    /* Handle the cases when kUseWholeImageAsOS is used */
    if( opts.useWholeImageAsOS){
        cout << yellow << "\nWarning: \'-w\' option provided, will use whole image to compute zero baseline. Could fail miserably if image is not mostly empty.\n" << normal;
    }
    return 0;
}

int main(int argc, char *argv[])
{
    time_t start,end;
    double dif;
    time (&start);


    string outFileBaseName;
    string outDir;
    string inFile;
    options_t opts;

    basicImageInfo_t bii;
    int returnCode = processCommandLineArgs( argc, argv, inFile, outFileBaseName, opts, outDir, bii);
    if(returnCode!=0){
        if(returnCode == 1) printHelp(argv[0],true);
        if(returnCode == 2) printHelp(argv[0]);
        return returnCode;
    }

    if(gVerbosity){
        cout << bold << "\nWill read the following file:\n" << normal;
        cout << "\t" << inFile << endl;
        cout << bold << "\nWill create the following output file(s):\n";
        if( !(opts.doNotSaveRoot)       ) cout <<"\t" << normal << outDir+outFileBaseName+".root" << endl;
        if(  (opts.saveFits)        ) cout <<"\t" << normal << outDir+"proc_"+outFileBaseName+".fits" << endl;
        if(  (opts.saveFitsSamples) ) cout <<"\t" << normal << outDir+"smpl_"+outFileBaseName+".fits" << endl;
        cout << endl;
    }

    int status = 0;
    status = procSkipperImage       (inFile.c_str(), outFileBaseName, bii, opts, outDir);

    if (status != 0){
        if(status>0) fits_report_error(stderr, status);
        return status;
    }

    /* Report */
    time (&end);
    dif = difftime (end,start);
    if(gVerbosity) cout << green << "\nAll done!\n" << bold << "-> It took me " << dif << " seconds to do it!\n\n" << normal;

    // return status;
    return 0;
}
