/*=========================================================================
  
  Copyright 2002, 2006-2010, Dr. Sandra Black
  Linda C. Campbell Cognitive Neurology Unit
  Sunnybrook Health Sciences Center
  
  This file is part of the Sunnybrook Image Software Processing (SIPS) package

  SIPS is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
 
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.

=========================================================================*/

/////////////////////////////////////////////////////////////
//
// written by Natasa Kovacevic, 2000
//
////////////////////////////////////////////////////////////
// compile with:
// /home/brainlab/gcc-2.95.3/bin/g++ -O2 -o mask_img mask_img.cpp

#include "../../../include_from_cortex/analyze_utils.h"

using namespace nk;
const char *usage = "\

Usage: mask_img <image fileroot> <mask fileroot>
       [-o <outfileroot>] for naming masked image different than default
       [-m <n1 n2 ...>] values to be masked
                       (put this switch always at the end if command line)
       By default, i.e. without -m switch, masking will be done so that
       image voxels that correspond to maximal value in the mask
       (e.g. 1 in a binary mask) are kept, while the rest are set to 0.
       Using -m 2 3 6 will do the masking so that all image voxels that
       correspond to values 2, 3 and 6 in the mask will be set to 0, while
       the rest will remain unchanged.

Example 1: mask_img bald_P bald_HfB
will read maximal value in bald_HfB (say 8) and set to 0 all voxels in bald_P
whose corresponding value in bald_HfB is not 8.

Example 2: mask_img bald_T1 bald_T1_seg -o bald_T1_white -m 0 75 150 255
will set to 0 all voxels in bald_T1 whose corresponding value in bald_T1_seg
is 0, 57, 150 or 255. (e.g. keep white tissue only (value 225 in T1_seg),
everything else masked).

" ;

template <class T_img, class T_mask>
void Mask_image(const string &img_fname, const string &mask_fname,
                vector<int> &mask_val,
                const string &out_fname);

int main(int argc, char ** argv)
{
    if(argc < 3 ){
        cout << usage;
        exit(1);
    }
    
    info img_info;
    string img_fname(argv[1]);
    read_analyze_header(img_fname + ".hdr", img_info);
    
    info mask_info;
    string mask_fname(argv[2]);
    read_analyze_header(mask_fname + ".hdr", mask_info);

    cout << "Image : " << img_fname << endl;
    cout << "Mask  : " << mask_fname << flush << endl;
    
    if( !equal(img_info.sz, img_info.sz+3, mask_info.sz) ) {
        cout << img_fname << " and " << mask_fname
             << " do not have equal dimensions. Exiting... " << endl;
        exit(1);
    }

    if( !equal(img_info.voxsz, img_info.voxsz+3, mask_info.voxsz) )
        cout << "Warning: " << img_fname << " and " << mask_fname
             << "do not have same voxel size. Proceeding..." << endl;

    // find out which values should be masked out
    vector <int> mask_value;
    
    CommandLine cl(argc, argv);
    if( cl.find("-m") ){ 
        for(int j=cl.current() + 1; j < argc; ++j)
            mask_value.push_back(atoi(argv[j]));
    }
    
    string out_fname;
    if(cl.find("-o")) out_fname = cl[1];
    else out_fname = img_fname + "_masked";
    
    // consider all possibilities for data types and
     // choose proper template types
    if(img_info.datatype == 2 && mask_info.datatype == 2)
        Mask_image<unsigned char, unsigned char>(img_fname, mask_fname, mask_value, out_fname);
    
    if(img_info.datatype == 4 && mask_info.datatype == 2)
        Mask_image<short, unsigned char>(img_fname, mask_fname, mask_value, out_fname);
    
    if(img_info.datatype == 2 && mask_info.datatype == 4)
        Mask_image<unsigned char, short>(img_fname, mask_fname, mask_value, out_fname);
    
    if(img_info.datatype == 4 && mask_info.datatype == 4)
        Mask_image<short, short>(img_fname, mask_fname, mask_value, out_fname);
    return 0;
}
// end of main
//////////////////////////////////////////////////////////////////////////

template <class T_img, class T_mask>
void Mask_image(const string &img_fname, const string &mask_fname,
                vector<int> &mask_val,
                const string &out_fname) {
    
    info img_info, mask_info;
    T_img *img = read_analyze<T_img>(img_fname, img_info);
    T_mask *mask = read_analyze<T_mask>(mask_fname, mask_info);
    
    int volume = img_info.sz[0] * img_info.sz[1] * img_info.sz[2];
    T_img *img_ptr;
    T_mask *mask_ptr;
    
    if(mask_val.size() == 0) {
        // calculate maximum  value in the mask
        T_mask M = 0;
        mask_ptr = mask;
        for(int t=0; t < volume; ++t) {
            if(*mask_ptr > M) M = *mask_ptr;
            ++mask_ptr;
        }
        // set up everything else to be masked out
        for(T_mask i = 0; i < M; ++i) mask_val.push_back(i);
    }
    
    cout << "Masking values: ";
    for(int i=0; i < mask_val.size(); ++i)
    cout << mask_val[i] << ", ";
    cout << flush << endl;
    
    img_ptr = img;
    mask_ptr = mask;
    for(int t=0; t < volume; ++t) {
        // if *mask_ptr has value to be masked then set corr img_ptr to 0
        if( find(mask_val.begin(), mask_val.end(),
                 (int)(*mask_ptr)) != mask_val.end() ) *img_ptr = 0;
        ++img_ptr;   ++mask_ptr;
    }
    
    cout << "Writing masked image: " << out_fname << flush << endl;
    
    write_analyze(out_fname, img, img_info);
    
    delete [] img;
    delete [] mask;
}
    
