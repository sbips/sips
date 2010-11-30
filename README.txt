-----------------------------------------------------------------------

SIPS

The Sunnybrook Image Processing Software (SIPS) package consists primarily of
in-house developed programs that can be used to perform various image 
processing tasks such as tissue segmentation (T1seg), semi-automated brain
parcellation (SABRE), white matter hyperintensity and lesion segmentation 
(lesion explorer/flex), and head-from-brain extraction (auto_hfb_template).  

Some external supporting software is also included in the SIPS package, 
the details of which are provided below.

------

If the SIPS package is used as part of published research, please 
reference the following papers where appropriate.

FLEX:
Gibson E., Gao F., Black S.E., Lobaugh N.J.  (2010).  Automatic segmentation of white
matter hyperintensities in the elderly using FLAIR images at 3T. Journal of Magnetic
Resonance Imaging, 31, 1311-22.

SABRE
Dade L.A., Gao F.Q., Kovacevic N., Roy P., Rockel C., O'Toole C.M., Lobaugh N.J., 
Feinstein A., Levine B., Black S.E.  (2004).  Semiautomatic brain region extraction: 
a method of parcellating brain regions from structural magnetic resonance images. 
Neuroimage, 22, 1492-502.

T1seg / auto_hfb
Kovacevic N, Lobaugh NJ, Bronskill MJ, Levine B, Feinstein A, Black SE.
(2002).  A robust method for extraction and automatic segmentation of 
brain images. Neuroimage, 17, 1087-100.  

------

Terms and Conditions governing the use of the SIPS package:

Copyright 2002, 2006-2010, Dr. Sandra Black
Linda C. Campbell Cognitive Neurology Unit
Sunnybrook Health Sciences Center

The SIPS package is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software 
Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

----------------------------------------------------------------------------

External Supporting Software

Source code and binaries for some FSL version 3.3 programs 
(avwmaths, bet, flirt) have been included and are used in the 
registration pipeline.  

Binaries for N3 and supporting programs are also included and are used 
for inhomogeneity correction.  The source code for these binaries can 
be downloaded from:  http://packages.bic.mni.mcgill.ca/tgz/.    

The following licenses and copyrights apply:

-------

FMRIB Software Library, Release 3.3 (c) 2006, The University of Oxford 
(the "Software")

The Software remains the property of the University of Oxford ("the
University").

The Software is distributed "AS IS" under this Licence solely for
non-commercial use in the hope that it will be useful, but in order
that the University as a charitable foundation protects its assets for
the benefit of its educational and research purposes, the University
makes clear that no condition is made or to be implied, nor is any
warranty given or to be implied, as to the accuracy of the Software,
or that it will be suitable for any particular purpose or for use
under any specific conditions. Furthermore, the University disclaims
all responsibility for the use which is made of the Software. It
further disclaims any liability for the outcomes arising from using
the Software.

The Licensee agrees to indemnify the University and hold the
University harmless from and against any and all claims, damages and
liabilities asserted by third parties (including claims for
negligence) which arise directly or indirectly from the use of the
Software or the sale of any products based on the Software.

No part of the Software may be reproduced, modified, transmitted or
transferred in any form or by any means, electronic or mechanical,
without the express permission of the University. The permission of
the University is not required if the said reproduction, modification,
transmission or transference is done without financial return, the
conditions of this Licence are imposed upon the receiver of the
product, and all original and amended source code is included in any
transmitted product. You may be held legally responsible for any
copyright infringement that is caused or encouraged by your failure to
abide by these terms and conditions.

You are not permitted under this Licence to use this Software
commercially. Use for which any financial return is received shall be
defined as commercial use, and includes (1) integration of all or part
of the source code or the Software into a product for sale or license
by or on behalf of Licensee to third parties or (2) use of the
Software or any derivative of it for research with the final aim of
developing software products for sale or license to a third party or
(3) use of the Software or any derivative of it for research with the
final aim of developing non-software products for sale or license to a
third party, or (4) use of the Software to provide any service to an
external organisation for which payment is received. If you are
interested in using the Software commercially, please contact Isis
Innovation Limited ("Isis"), the technology transfer company of the
University, to negotiate a licence. Contact details are:
innovation@isis.ox.ac.uk quoting reference BS/3497.

------

N3-1.10
Copyright 1996, John G. Sled
McConnell Brain Imaging Centre,
Montreal Neurological Institute, McGill University.

------

ebtks-1.6.1
Copyright 1996, Alex P. Zijdenbos,
McConnell Brain Imaging Centre,
Montreal Neurological Institute, McGill University.

------

hdf5-1.6.8
HDF5 (Hierarchical Data Format 5) Software Library and Utilities
Copyright 2006-2008 by The HDF Group (THG).

------

minc-1.5.1
Copyright 1993-2000 Peter Neelin and David MacDonald, McConnell Brain
Imaging Centre, Montreal Neurological Institute, McGill University.

-------

netcdf-3.6.3
Copyright 1993-2004 University Corporation for Atmospheric Research/Unidata

----------------------------------------------------------------------------



