//  Copyright 2002 by Natasa Kovacevic

// Compile with: /srutility/gnu/gcc-2.95.2/bin/c++ -O2 -o T1seg T1seg.cpp

#include "analyze_utils.h" // *****NK*****
#include <numeric>
#include <assert.h> // *****NK*****
using namespace nk;

#define MAX_ITER  300  // maximal num of iterations for EM
#define EPSILON 0.0001 // convergence criterion for EM
#define NUM_GAUSSIANS 4
#define MIN_NUM_DATA 20000  // minimal number of non-zero data points in local box

////////////// segmentation values set as global variables /////////////////
unsigned char BGD=0, CSF = 75, GRAY = 150, WHITE = 255;

/////////////////// forward declarations ///////////////////////////////////

void print_usage_message(const char *command);
void parse_command_line(int argc, char **argv,
                        string &in_fname, string &out_fname,
                        double &cutoff, int &histo_start, int &win_size,
                        int *step, int *increment,
                        bool &v, bool &global);
inline bool fit_global_gaussians(unsigned short *Img, int *sz, double cutoff,
				 int histo_start, int win_size,
                                 double *global_weight,
                                 double *global_mean,
                                 double *global_stdev);
void find_borders(const unsigned short *Img, const int *sz, int *B);
inline void output_borders(const int *B);
inline void output_gaussian_parameters(const double *weight,
                                       const double *mean,
                                       const double *stdev);
inline void segment(const unsigned short *Img, unsigned char *Seg,
                    const int *B, const int *sz,
                    const double *mean, const double *stdev);
inline int num_of_nonzero_voxels(unsigned short *Img, const int *B,
                                 const int *sz);
inline void calculate_local_histo(vector<int> &histo,int win_size,
                                  const unsigned short *Img,
                                  const int *LB, const int *sz);
bool fit_gaussians_to_histo(const vector <int> &histo,
                            double *weight,
                            double *mean,
                            double *stdev,
                            int max_iter, double epsilon);

////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
    if(argc < 3 ){
        print_usage_message(argv[0]);
        exit(1);
    }

    std::string in_fname, out_fname;
    int histo_start, win_size = 0;
    double cutoff;
    int step[3], increment[3];
    bool v, global;

    parse_command_line(argc, argv, in_fname, out_fname, cutoff, histo_start, win_size,step, increment, v, global);

    int sz[3];
    // read input file
    info hdr;
    unsigned short *Img = read_analyze<unsigned short>(in_fname, hdr);
    copy(hdr.sz, hdr.sz + 3, sz);
    int volume = sz[0] * sz[1] * sz[2];

    bool success;
    // Fit gaussians to the global histogram
    double global_weight[NUM_GAUSSIANS], global_mean[NUM_GAUSSIANS], global_stdev[NUM_GAUSSIANS];
    success = fit_global_gaussians(Img, sz, cutoff, histo_start, win_size,
                         global_weight, global_mean, global_stdev);
    if(!success) {
        cout << "EM algorithm failed to fit global parameters. Exiting ... "
             << endl;
        exit(1);
    }

    // Craete seg_image of the same size as in_img and fill it with BGD
    unsigned char *Seg = new unsigned char [volume];
    if(!Seg){
        cout << "Memory allocation error. Exiting ... " << endl;
        exit(1);
    }
    fill(Seg, Seg + volume, BGD);

    // find borders of brain
    int GB[6];
    find_borders(Img, sz, GB);

    if( global )  {
        cout << "Performing global segmentation ... " << endl;
        segment(Img, Seg, GB, sz, global_mean, global_stdev) ;
    }
    else {
         cout << "Performing localized segmentations. Local box size ["
              << step[0] << ", " << step[1] << ", " << step[2] << "]"<< endl;

         // set upper bound for all local histos
        int histo_end = (int)ceil(1.2 * global_mean[3]);

        int shift[3];
        for(int t=0; t < 3; t++)  shift[t] = step[t]/3;

        double weight[NUM_GAUSSIANS], mean[NUM_GAUSSIANS], stdev[NUM_GAUSSIANS];

        // LB=borders of local box, CLB=borders of the core of local box
        int LB[6], CLB[6];

        int t, index[3];

        int boxnumber=0;
        for(index[2]=GB[4]; index[2] < GB[5] - shift[2] + 1; index[2] += shift[2])
            for(index[1]=GB[2]; index[1] < GB[3] - shift[1] + 1; index[1] += shift[1])
                for(index[0]=GB[0]; index[0] < GB[1] - shift[0] + 1; index[0] += shift[0]) {

                    boxnumber++;

                    for(t=0; t < 3; ++t) {
                        // lower local box border
                        CLB[2*t] = LB[2*t] = index[t];
                        // upper local box border
                        CLB[2*t+1] = LB[2*t+1] = minval(LB[2*t] + step[t] - 1, GB[2*t+1]);
                    }

                    while( num_of_nonzero_voxels(Img, LB, sz) < MIN_NUM_DATA ) {
                        // extend local box size by increment
                        for(t=0; t < 3; ++t){
                            LB[2*t] = maxval(LB[2*t] - increment[t], GB[2*t]);
                            LB[2*t+1] = minval(LB[2*t+1] + increment[t], GB[2*t+1]);
                        }
                    }

                    // Do local histogram fitting. If neccessary,
                    // increase local box until fitting is successful
                    vector<int> histo(histo_end + 1, 0);
                    success = false;
                    while( !success ) {
                        // calculate local histo
                        calculate_local_histo(histo, win_size, Img, LB, sz);
                        // initilaze local params by global values
                        for(t=0; t < NUM_GAUSSIANS; ++t) {
                            weight[t] = global_weight[t];
                            mean[t] = global_mean[t];
                            stdev[t] = global_stdev[t];
                        }
                        // do EM algorithm to fit NUM_GAUSSIANS gaussans to local histo
                        success = fit_gaussians_to_histo(histo, weight, mean, stdev, MAX_ITER, EPSILON);
                        if(!success) {
                            // extend local box size some more
                            for(t=0; t < 3; ++t){
                                LB[2*t] = maxval(LB[2*t] - increment[t], GB[2*t]);
                                LB[2*t+1] = minval(LB[2*t+1] + increment[t], GB[2*t+1]);
                            }
                        }
                    }

                    // Reduce box to its middle core for all but edge boxes
                    for(t=0; t < 3; ++t){
                        if(CLB[2*t]==GB[2*t]) // inital box was on the lower edge
                            CLB[2*t+1] = GB[2*t] + 2*shift[t];
                        else if(CLB[2*t+1]==GB[2*t+1])// inital box was on the upper edge
                            CLB[2*t] += shift[t];
                        else { // inital box was away from the edges
                            CLB[2*t] += shift[t];
                            CLB[2*t+1] = CLB[2*t] + shift[t];
                        }

                    }
                    // Now do segmentation on the core box
                    segment(Img, Seg, CLB, sz, mean, stdev);

                    if(v){
                        cout << endl;
                        output_borders(CLB);
                        output_gaussian_parameters(weight, mean, stdev);
                    }
                }
        cout << endl;
    }
    hdr.datatype = 2;
    write_analyze(out_fname, Seg, hdr);

    delete [] Img;
    delete [] Seg;
    return 0;
}  // end of main

///////////////////////////////////////////////////////////////////////////

void print_usage_message(const char *command) {
  cout << "\nUsage: " << command << " input output  [options] \n" << endl;
  cout << "where input is brain-extracted T1 image (all non-brain matter voxels set to 0)"<<endl;
  cout << "format:  unsigned short (16bit) Analyze data format (.img/.hdr extensions not needed) \n"<< endl;
  cout << "-------------------------------------------------" << endl;
  cout << "options: \n" << endl;
  cout << "[-win_size] -> window size for histogram smoothing (default 0, i.e. no smoothing)" << endl;
  cout << "[-c] -> cutoff expressed as fraction of total brain volume (default 0.0002)" << endl;
  cout << "[-rc] -> start value for right cutoff (default 0, i.e. tries to find rc automatically). Good value is WM peak.\n" << endl;
  cout << "[-v] -> option to monitor local gaussian parameters\n" << endl;
  cout << "[-x xboxdim] -> set x-dimension of local box (default 48)" << endl;
  cout << "[-y yboxdim] -> set y-dimension of local box (default 48)" << endl;
  cout << "[-z zboxdim] -> set z-dimension of local box (default 36)\n" << endl;
  cout << "[-incx incx] -> set x_increment for local box enlargement (default 2)" << endl;
  cout << "[-incy incy] -> set y_increment for local box enlargement (default 2)" << endl;
  cout << "[-incz incz] -> set z_increment for local box enlargement (default 1)\n" << endl;
  cout << "[-bfgw bgd csf gray white ] -> set segmentation values for background, csf, grey, white "<< endl;
  cout << "                               (defaults: 0 75 150 255)\n" << endl;
  cout << "[-global] -> segmentation based on global gaussians only (no localized gaussian fitting)\n" << endl;
  cout << "-------------------------------------------------" << endl;
  cout << "Example1: T1seg subjectX_T1_masked subject_T1_seg  -x 42 -y 42 -z 42 -incx 1 -incy 1 -incz 1 \n" << endl;
  cout << "Performs segmentation in local boxes of size 42x42x42." << endl;
  cout << "When a local box has less than MIN_NUM_DATA=20000 brain voxels, " << endl;
  cout << "its size is enlarged by increments +/-1 in all three directions\n" << endl;
  cout << "Example2: T1seg subjectX_T1_masked subjectX_T1_seg 256 256 124 -bfgw 0 1 2 3 \n" << endl;
  cout << "Performs segmentation so that output image voxels with value 0 represent background," << endl;
  cout << "voxels with value 1 represent csf, voxels with value 2 represent gray matter, " << endl;
  cout << "and voxels with value 3 represent white matter \n" << endl;
}
///////////////////////////////////////////////////////////////////////
void parse_command_line(int argc, char **argv,
                        std::string &in_fname, std::string &out_fname,
                        double &cutoff, int &histo_start, int &win_size,
                        int *step, int *increment,
                        bool &v, bool &global)
{
    CommandLine cl(argc, argv);

    in_fname = std::string(argv[1]); out_fname = std::string(argv[2]);
    cout << "Input image: " << in_fname << endl;
    cout << "Segmented image: " << out_fname << endl;

    if(cl.find("-bfgw")) {
        BGD = atoi(cl[1]);
        CSF = atoi(cl[2]);
        GRAY = atoi(cl[3]);
        WHITE = atoi(cl[4]);
    }
 
    if(cl.find("-win_size")) {
      win_size  = atoi(cl[1]);
      cout << "Histograms will be smoothed with moving window average with window size " << win_size << endl;
    }

    cutoff = 0.0002;
    if(cl.find("-c")) cutoff = atof(cl[1]);

    histo_start = 0;
    if(cl.find("-rc")) { 
      histo_start = atoi(cl[1]);
      cout << "starting right cutoff search at " << histo_start << endl;
    }

    step[0] = 48; step[1] = 48; step[2] = 36;
    if(cl.find("-x")) step[0] = atoi(cl[1]); step[0] = (step[0]/3)*3;
    if(cl.find("-y")) step[1] = atoi(cl[1]); step[1] = (step[1]/3)*3;
    if(cl.find("-z")) step[2] = atoi(cl[1]); step[2] = (step[2]/3)*3;

    increment[0] = 2; increment[1] = 2; increment[2] = 1;
    if(cl.find("-incx")) increment[0] = atoi(cl[1]);
    if(cl.find("-incy")) increment[1] = atoi(cl[1]);
    if(cl.find("-incz")) increment[2] = atoi(cl[1]);

    v = false;
    if(cl.find("-v")) v = true;
    global = false;
    if(cl.find("-global")) global = true;
}
///////////////////////////////////////////////////////////////////////
void find_borders(const unsigned short *Img, const int *sz, int *B)
{
    int area = sz[0] * sz[1], volume = area * sz[2];
    int i, j, k;

    // find lowest non-zero z-level
    i=0;
    while( Img[i]==0 && i<volume) i++;
    B[4] = i / area;

    // find highest non-zero z-level
    i=volume-1;
    while( Img[i] == (unsigned short)0 && i>=0) i--;
    B[5] = i / area;

    // find lowest non-zero x-level
    bool found=false;
    for(i=0; i < sz[0] && !found; i++)
        for(j=0; j < sz[1] && !found; j++)
            for(k=0; k < sz[2] && !found; k++)
                if( Img[k*area + j*sz[0] + i] ) found=true;
    B[0] = maxval(i-1, 0);

   // find highest non-zero x-level
    found=false;
    for(i=sz[0]-1; i>=0 && !found; i--)
        for(j=0; j < sz[1] && !found; j++)
            for(k=0; k < sz[2] && !found; k++)
                if( Img[k*area + j*sz[0] + i] ) found=true;
    B[1] = minval(i+1, sz[0]-1);

    // find lowest non-zero y-level
    found=false;
    for(j=0; j < sz[1] && !found; j++)
        for(i=0; i < sz[0] && !found; i++)
                  for(k=0; k < sz[2] && !found; k++)
                       if( Img[k*area + j*sz[0] + i] ) found=true;
    B[2] = maxval(j-1, 0);

    // find highest non-zero y-level
    found=false;
    for(j=sz[1]-1; j >=0 && !found; j--)
        for(i=0; i < sz[0] && !found; i++)
                  for(k=0; k < sz[2] && !found; k++)
                       if( Img[k*area + j*sz[0] + i] ) found=true;
    B[3] = minval(j+1, sz[1]-1);

    cout << "Borders :"
         << "  x_range=[" << B[0] << ":" << B[1] << "]"
         << "  y_range=[" << B[2] << ":" << B[3] << "]"
         << "  z_range=[" << B[4] << ":" << B[5] << "]"
         << endl;

    // check that borders in x, y direction are not same as image edges
    if(B[0]==0 || B[1]==sz[0]-1 ||
       B[2]==0 || B[3]==sz[1]-1) {
        cout << "Image borders suspicious. Check your mask." << endl;
        cout << "Algorithm will proceed, results unpredictable. " << endl;
    }
}
///////////////////////////////////////////////////////////////////////////
inline void output_borders(const int *B){
    cout << "["<< B[0] << ":" << B[1] << ", "
               << B[2] << ":" << B[3] <<", "
               << B[4] << ":" << B[5] << "]";
}
///////////////////////////////////////////////////////////////////////////
inline void output_gaussian_parameters(const double *weight,
                                       const double *mean,
                                       const double *stdev)
{
    ostream_iterator<double> out(cout, " ");
    cout << "\nweights: ";
    copy(weight, weight + NUM_GAUSSIANS, out);    cout << endl;
    cout << "means:   ";
    copy(mean, mean + NUM_GAUSSIANS, out);        cout << endl;
    cout << "stdevs: ";
    copy(stdev, stdev + NUM_GAUSSIANS, out);      cout << endl;
}
///////////////////////////////////////////////////////////////////////////
inline void order(double *weight, double *mean, double *stdev)
    // Utility function for ordering gaussian parameters so that
    // means are in ascending order
{
        double E=0.00000000001;
        for(int i=0; i< NUM_GAUSSIANS; i++) {
            for(int j=i+1; j < NUM_GAUSSIANS; j++){
                if( mean[j] < mean[i] - E){
                    swap(mean[j], mean[i]);
                    swap(weight[j], weight[i]);
                    swap(stdev[j], stdev[i]);
                }
            }
        }
}
/////////////////////////////////////////////////////////////////////////////
bool fit_global_gaussians(unsigned short *Img, int *sz, double cutoff, 
			  int histo_start, int win_size,
                          double *global_weight,  double *global_mean,
                          double *global_stdev)
{
    // calculate global image histogram and at the same time
    // calculate total count of brain (non-zero) voxels

    int raw_histo[65536]; // initalize histogram with 0's
    for(int i=0; i < 65536; i++) raw_histo[i] = 0;
    int volume = sz[0] * sz[1] * sz[2];
    unsigned short *ptr = Img;
    for(int t=0; t < volume; ++t){
      ++(raw_histo[*ptr]);
      ++ptr;
    }
    // set histo to 0 at end points
    raw_histo[0] =  raw_histo[65535] = 0 ;
    std::vector<int> histo(raw_histo, raw_histo + 65536 );

    // get total(histo)= sum of all values
    int total= accumulate(histo.begin(), histo.end(), 0);

    // smooth histogram with 10 point half window
    if (win_size > 0) {
      cout << "smoothing histo with window = " << win_size << endl;
      double tmp;
      int cnt;
      for(int w=0; w < histo.size(); w++) {
	tmp=0.0;
	cnt=0;
	for(int i=maxval(0,w-win_size/2); i <= minval(histo.size()-1, w+win_size/2); i++, cnt++) { 
	  tmp += (double)(raw_histo[i]); 
	}      
	if(cnt > 0) histo[w] = round(tmp/cnt); 
	else histo[w] = 0;
      }
    }

    int rc;
    if ( histo_start == 0 )
    {
       // try finding a good starting rc value
	    
	std::cout << "Attempting to find starting rc value... " << std::endl;
       
        // find maximum histogram value
        int max_bin = 0;
        for (int i=0; i < histo.size(); ++i)
        {
  	  if ( histo[i] > histo[max_bin] )
	  {
		max_bin = i;
	  }
      }
    
      // use the intensity value corresponding to the bin containing
      // the specified percentage of the max histogram value for the starting rc
      // -- should land about midway between WM peak and actual rc value
      /////////////////////////////////////////////////
      int auto_rc_bin_val = histo[max_bin] * 0.5;
      ////////////////////////////////////////////////
      
      int auto_rc = 0;
      for (int i=histo.size()-1; i>=0; --i)
      {
	  if ( histo[i] >= auto_rc_bin_val )
	  {
		auto_rc = i;
		break;
	  } 
        }
       rc = auto_rc;
     }
    else
    {
	// use starting rc from command line
	rc = histo_start;	 
    } 

    std::cout << "starting rc value = " << rc << std::endl;
    
    // find histogram's right cutoff point
    while( histo[rc] > cutoff * total ) ++rc;
    
    if (rc == histo_start) {
        cout << "Finding right cutoff point of the global histogram failed."
             << endl;
        exit(1);
    }
    
    std::cout << "rc value = " << rc << std::endl;

    // trim histo from the right
    histo.erase(histo.begin() + rc, histo.end());

    // Now do standard initialization for NUM_GAUSSIANS gaussians fitted
    // to the histo of the entire T1
    global_weight[0] = 0.15;
    global_weight[1] = 0.05;
    global_weight[2] = 0.4;
    global_weight[3] = 0.4;
    global_mean[0] = 0.25 * rc;
    global_mean[1] = 0.35 * rc;
    global_mean[2] = 0.67 * rc;
    global_mean[3] = 0.83 * rc;
    global_stdev[0] = 0.07 * rc;
    global_stdev[1] = 0.03 * rc;
    global_stdev[2] = 0.12 * rc;
    global_stdev[3] = 0.12 * rc;

    // Perform EM algorithm on histo
    bool success = fit_gaussians_to_histo(histo,
					  global_weight, global_mean,
					  global_stdev,
					  MAX_ITER, EPSILON);

    cout << "Global model parameters for the whole brain : ";
    output_gaussian_parameters(global_weight, global_mean, global_stdev);
    cout << endl;

    if(!success || !(global_mean[0] > 15.0 && global_mean[3] < 60000.0)) {
      cout << "Global paremeters found suspicious." << endl;
      cout << "Check whether brain was properly extracted." << endl;
      cout << "Algorithm will proceed but results unpredictable." << endl << endl;
    }
    
    return success;
}
///////////////////////////////////////////////////////////////////////////
inline void segment(const unsigned short *Img, unsigned char *Seg,
                    const int *B, const int *sz, const double *mean,
                    const double *stdev)
    // Perform segmentation within borders defined by borders B
    {
        // csf/gray cutoff: midpoint between mean[0] and mean[2]
        // gray/white cutoff: midpoint between mean[2] and mean[3]
        unsigned short fg_cutoff, gw_cutoff, wh_cutoff;
        fg_cutoff = round((mean[0] + mean[2]) * 0.5);
        gw_cutoff = round((mean[2] + mean[3]) * 0.5);
	int x, y, z, t, area = sz[0]*sz[1];

        for(z=B[4]; z <= B[5]; ++z)
            for(y=B[2]; y <= B[3]; ++y)
                for(x=B[0]; x <= B[1]; ++x) {
                    t = z * area + y * sz[0] + x;
                    if( Img[t] ) {
                        if( Img[t] < fg_cutoff )  Seg[t] = CSF;
                        if( Img[t] >= fg_cutoff &&
                            Img[t] < gw_cutoff )  Seg[t] = GRAY;
                        if( Img[t] >= gw_cutoff ) Seg[t] = WHITE;
		    }
                }
    }
/////////////////////////////////////////////////////////////////////////////
inline int num_of_nonzero_voxels(unsigned short *Img, const int *B, const int *sz)
{
    int count=0;
    int x, y, z, area = sz[0]*sz[1];
    for(z=B[4]; z <= B[5]; ++z)
            for(y=B[2]; y <= B[3]; ++y)
                for(x=B[0]; x <= B[1]; ++x) {
                    if( Img[ z * area + y * sz[0] + x] ) ++count;
                }
    return count;
}
/////////////////////////////////////////////////////////////////////////////
void calculate_local_histo(vector<int> &histo, int win_size,
			   const unsigned short *Img,
                           const int *B, const int *sz)
    // calculate histogram values within box defined by borders B
    // disregard any voxels with intensity > histo.size()
{
    int x, y, z, t, area = sz[0]*sz[1], histo_end = histo.size() - 1;

    for(z=B[4]; z <= B[5]; ++z)
        for(y=B[2]; y <= B[3]; ++y)
            for(x=B[0]; x <= B[1]; ++x) {
                t = z * area + y * sz[0] + x;
                if( Img[t] <= histo_end ) histo[ Img[t] ] += 1;
            }
    histo[0] = 0;
    // smooth local histo
    if (win_size > 0) {
      vector<int> copy_histo(histo.begin(), histo.end());
      double tmp;
      int cnt;
      for(int w=0; w < histo.size(); w++) {
	tmp=0.0;
	cnt=0;
	for(int i=maxval(0,w-win_size/2); i <= minval(histo.size()-1, w+win_size/2); i++, cnt++) { 
	  tmp += (double)(copy_histo[i]); 
	}      
	if(cnt > 0) histo[w] = round(tmp/cnt); 
	else histo[w] = 0;
      }
    }
}
////////////////////////////////////////////////////////////////////////////
inline double gauss(const double &mean, const double &stdev, const int &x) {
    return  exp( -0.5*(x-mean)*(x-mean)/(stdev*stdev) ) /
        (sqrt(M_PI*2.0)*stdev);
}
////////////////////////////////////////////////////////////////////////////
inline void EM_iterate(unsigned int num_of_gaussians,
                double *weight,
                double *mean,
                double *stdev,
                const vector <int> &histo,
                const int &card,
                vector<double> *P,
                double *sumP,
                double *sumPx)
 {
     unsigned int x, i;
    // calculate S and P
    for(i=0; i<num_of_gaussians; ++i) sumP[i]=sumPx[i]=0.0;
    double S;
    // get each P(x) and also  get sum of all P_i(x) and also of all P_i(x)*x
    for(x=0; x < histo.size(); ++x) {
        S = 0.0;
        for(i=0; i<num_of_gaussians; ++i) {
            P[i][x] = weight[i] * gauss(mean[i], stdev[i], x);
            S += P[i][x];
        }
        for(i=0; i<num_of_gaussians; ++i) {
            P[i][x] /= S;
            sumP[i] += histo[x] * P[i][x];
            sumPx[i] += histo[x] * P[i][x] * x;
        }
    }

    // update weights and means
    for(i=0; i<num_of_gaussians; ++i) {
        weight[i] = sumP[i]/card;
        mean[i] = sumPx[i] / sumP[i];
    }

    // update stdevs
    double sum;
    for(i=0; i<num_of_gaussians; i++) {
        sum = 0.0;
        for(x=0; x < histo.size(); x++) {
            sum += histo[x] * P[i][x] * (x - mean[i]) * (x - mean[i]);
        }
        stdev[i] = sqrt( sum / (sumP[i] - 1.0) );
    }
 }
////////////////////////////////////////////////////////////////////////////
inline bool fit_gaussians_to_histo(const vector <int> &histo,
                            double *weight, double *mean, double *stdev,
                            int max_iter, double epsilon)
{
    int num_of_gaussians = NUM_GAUSSIANS;
    vector<double> *P = new vector<double> [num_of_gaussians];
    if(!P) {
        cout << "Memory allocation error " << endl;
        exit(1);
    }

    // do trick to get all P[i] of the same size
    vector<double> Temp(histo.size());
    for(int i=0; i<num_of_gaussians; i++) P[i] = Temp;
    double *sumP = new double[num_of_gaussians];
    assert(sumP);
    double *sumPx = new double[num_of_gaussians];
    assert(sumPx);

    // get cardinality
    int card = 0;
    for(unsigned short x=0; x < histo.size(); ++x) card += histo[x];

    bool success = true;
    double oldweight[NUM_GAUSSIANS], oldmean[NUM_GAUSSIANS], oldstdev[NUM_GAUSSIANS];
    for(int iter=0; iter < max_iter; iter++){

        for(int t=0; t < NUM_GAUSSIANS; ++t){
            oldweight[t]= weight[t];
            oldmean[t] = mean[t];
            oldstdev[t] = stdev[t];
        }

        EM_iterate(num_of_gaussians,
                   weight, mean, stdev,
                   histo, card,
                   P, sumP, sumPx);

        // make sure that  parameters have alowable values
        // if not,then main program will extend the box size and retry
        for(int j=0; j < num_of_gaussians; ++j) {
           if(!(mean[j] > 0.0 && mean[j] < 65535.0 && 
	     stdev[j] > 0.0 && stdev[j] <10000.0 &&
	     weight[j] > 0.0 && weight[j] < 1.0)) {
	    success = false; 
	    break;
	  }
        }
	if(success) {
	  // measure change between two successive iterates
	  double change = 0.0;
	  double temp;
	  for(int j=0; j < NUM_GAUSSIANS; j++) {
            temp = oldweight[j] - weight[j];
            change += temp*temp;
            temp = oldmean[j] - mean[j];
            change += temp*temp;
            temp = oldstdev[j] - stdev[j];
            change += temp*temp;
	  }
	  if(change < epsilon) {iter = max_iter;}
	}

    } // end of iter loop

    delete [] P;
    delete [] sumP;
    delete [] sumPx;

    // Order parameters, in case some means crossed
    if(success && !(mean[0]<=mean[1] && mean[1]<=mean[2] && mean[2]<=mean[3])) order(weight, mean, stdev);
    return success;
}


