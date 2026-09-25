#include "EDepSimInterface.h"
#include "TH2D.h"

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#ifdef USE_PYTHON3
#include "numpy/arrayobject.h"
#include <cassert>
#else
#include <numpy/ndarrayobject.h>
#endif

#include <iostream>
#include <cmath>

namespace simpledet {
namespace edepsim {

  bool EDepSimInterface::__loaded_numpy = false;

  EDepSimInterface::EDepSimInterface()
    : _readout_image_isgood(false)
  {
    // This class will convert 3D location of energy deposited as ionization
    // into a 2D image from a simple 2D readout.

    // 
    
    // hard code stuff for now, then generalize later to focus

    // when we crop to a specific size image, make a "pre-image" crop that includes the size of
    // ionization diffusion kernel
    padding[0] = 7; // pixels
    padding[1] = 7; // pixels

    // size of the pixels of our fake 2D "readout"
    pixelsize[0] = 0.3; //cm
    pixelsize[1] = 0.3; //cm

    // define the "pre-image" size
    npixels[0] = 3456 + 2*padding[0];
    npixels[1] = 780 + 2*padding[1];
    
    // define our 2D readout plane (the 2D box within which we "capture" our image)
    for (int i=0; i<2; i++) {
      origin[i] = -pixelsize[i]*(npixels[i]/2);
      planelen[i] = npixels[i]*pixelsize[i];
      max_gridpt[i] = origin[i] + planelen[i];
    }
    
    max_step_size = 0.03;
    min_ionization_sigma_width = 0.01;
    distance_to_readout_plane = 128.0;
    max_depth = 256.0;
    num_edep_quantiles = 20;
  }

  // std::vector<int> EDepSimInterface::getPixelIndices( float x, float y ) const
  // {
  //   std::vector<int> pix_indices;
    
  // }
  void EDepSimInterface::clear()
  {
    // clear containers filled
    _pixinfo_map.clear();
    _edep_map.clear();

  }
  
  bool EDepSimInterface::processSegmentHits( const TG4HitSegmentContainer& hit_container )
  {

    _readout_image_isgood = false;
    
    size_t nhits = hit_container.size();

    size_t totpixels = npixels[0]*npixels[1];
    if ( _readout_image.size()!=totpixels ) {
      _readout_image.resize( totpixels, 0.0 );
    }

    /// clear pixel info containers
    clear();

    // zero out image
    std::memset( _readout_image.data(), 0.0, sizeof(float)*totpixels );
    std::cout << "making image with " << _readout_image.size() << " pixels" << std::endl;

    if ( hit_container.size()==0 )
      return _readout_image_isgood;
    
    // get time of first edep. based on how geant4 tracks particles, this should be the first
    // position of energy deposition
    auto const first_hit = hit_container.at(0);
    // we want to fix the start of the shower at:
    // (z=100 cm, y=0 cm)
    std::vector<float> first_edep(4,0);
    std::cout << "first edep pos: (";
    for (int i=0; i<3; i++) {
      first_edep[i]  = (first_hit.GetStart()[i])*0.1; /// mm -> cm
      std::cout << first_edep[i];
      if ( i<2 )
	std::cout << ", ";
    }
    std::cout << ")" << std::endl;
    first_edep[3] = first_hit.GetStart()[3]; // geant4 time in ns

    // shift the readout plane origin in det coords (cm) such that we place the first_edep position in the center of the image
    std::vector<float> offset(3,0); 
    offset[0] = first_edep[0]; // sets the depth
    offset[1] = first_edep[1];
    offset[2] = first_edep[2];

    // we will center the image on the first_edep
    // we adjust the origin of the readout plane
    // note: origin is in 2d plane coordinates (x,y), first_edep is in 3D det coordinates (all in cm)
    _readout_image_origin_cm[0] = first_edep[0]-distance_to_readout_plane;
    _readout_image_origin_cm[1] = origin[1]+first_edep[1]; 
    _readout_image_origin_cm[2] = origin[0]+first_edep[2]; 
  
    for (size_t ihit=0; ihit<nhits; ihit++) {
      
      auto const& hit = hit_container.at(ihit);
      
      // each hit is a line segment
      std::vector<float> start(4);
      std::vector<float> stop(4);
      std::vector<float> stepdir(4);
      for (int i=0; i<3; i++) {
	start[i] = hit.GetStart()[i]*0.1; // mm->cm
	stop[i]  = hit.GetStop()[i]*0.1;  // mm->cm
	stepdir[i] = stop[i]-start[i];
      }
      start[3] = hit.GetStart()[3]; // ns
      stop[3] = hit.GetStop()[3]; // ns
      stepdir[3] = stop[3]-start[3]; // ns
      
      float steplen = 0.0;
      for (int i=0; i<3; i++)
	steplen += stepdir[i]*stepdir[i];
      steplen = sqrt(steplen);
    
      for (int i=0; i<3; i++)
	stepdir[i] /= steplen;

      float Estep = hit.GetEnergyDeposit();

      // we divide the line segment into
      // small pieces
      // we then assign each segment to one or two pixels
      int ndivs = steplen/max_step_size + 1;
      float divlen = steplen/float(ndivs); // change the division size to be smaller, but fit exaclty on line

      for (int istep=0; istep<ndivs; istep++) {
	
	float fstep = (float)istep;
	std::vector<float> pt0(4,0.0);
	std::vector<float> pt1(4,0.0);
	std::vector<float> midpt(4,0.0);
	for (int i=0; i<4; i++) {
	  pt0[i] = start[i] + divlen*fstep*stepdir[i];
	  pt1[i] = start[i] + divlen*(fstep+1.0)*stepdir[i];
	  midpt[i] = 0.5*pt0[i] + 0.5*pt1[i];
	}

	// distribute via approximate gaussian kernel using tranverse diffusion
	// A. Lister, approximate spread:
	// 1.4 mm longitudinal over full 256.0 cm drift
	// 2.5 mm transverse drift over full 256.0 cm drift
	// each pixel is 3.0 mm, so we can limit our gaussian kernel to -3,+3 bins, a 7x7 Gaussian kernel
	// we make the kernel, taking into account the relative offset
	// based on x-position in middle of detector 125.0
	
	// get pos relative to origin, i.e. our plane coordinates
	std::vector<float> planepos(2,0);
	planepos[0] = (midpt[2]-_readout_image_origin_cm[2]); // det z-pos --> plane x pos
	planepos[1] = (midpt[1]-_readout_image_origin_cm[1]); // det y-pos --> plane y pos
	int izpix = planepos[0]/pixelsize[0];
	int iypix = planepos[1]/pixelsize[1];
	std::vector<float> pixcenter(2,0);
	pixcenter[0] = (float(izpix)+0.5)*pixelsize[0];
	pixcenter[1] = (float(iypix)+0.5)*pixelsize[1];
	
	float dz_offset = planepos[0]-pixcenter[0];
	float dy_offset = planepos[1]-pixcenter[1];
	float depth = midpt[0]-_readout_image_origin_cm[0];
	float sigma = (depth)*(0.25/256.0) + min_ionization_sigma_width;

	if ( depth<0.0 || depth>=max_depth )
	  continue;
	
	// std::cout << "------------------------" << std::endl;
	// std::cout << "midpt: " << midpt[2] << " " << midpt[1] << std::endl;
	// std::cout << "pos: " << pos[0] << " " << pos[1] << std::endl;
	// std::cout << "planepos: " << planepos[0] << " " << planepos[1] << std::endl;
	// std::cout << "pixelindex: " << izpix << " " << iypix << std::endl;
	// std::cout << "pixcenter: " << pixcenter[0] << " " << pixcenter[1] << std::endl;
	// std::cout << "pos relative to pixcenter: " << dz_offset << " " << dy_offset << std::endl;
	// std::cout << "sigma: " << sigma << " x=" << midpt[0] << std::endl;

	// simple binning fill
	// int ipix = (iypix)*npixels[0] + (izpix);
	// if ( ipix>=0 && ipix<image.size() ) {
	// 	image[ipix] = Estep/float(ndivs);
	// }

	// apply approximate kernel
	float sum = 0.0;
	float kernel[7][7] = {0.0};
	for (int ikz=-3; ikz<=3; ikz++) {
	  for (int iky=-3; iky<=3; iky++) {
	    // position of kernel pixel relative to center
	    float kpixz = float(ikz)*pixelsize[0];
	    float kpixy = float(iky)*pixelsize[1];
	    float dist2 = (kpixz-dz_offset)*(kpixz-dz_offset) + (kpixy-dy_offset)*(kpixy-dy_offset);
	    float arg = -0.5*dist2/(sigma*sigma);
	    float w = exp(arg);
	    kernel[3+ikz][3+iky] = w;
	    sum += w;
	  }
	}

	// Efrac
	int trackid = hit.GetContributors().at(0);
	float Efrac = Estep/float(ndivs);
	for (int ikz=-3; ikz<=3; ikz++) {
	  for (int iky=-3; iky<=3; iky++) {
	    int ipix = (iypix+iky)*npixels[0] + (izpix+ikz);
	    if ( ipix>=0 && ipix<(int)_readout_image.size() ) {
	      float frac_edep = Efrac*kernel[3+ikz][3+iky]/sum;
	      _readout_image[ipix] += frac_edep;
	      auto it_pix = _pixinfo_map.find( (unsigned long)ipix );
	      if ( it_pix==_pixinfo_map.end() ) {
		// create a place to store info
		PixInfo_t pixinfo;
		pixinfo.pixid = (unsigned long)ipix;
		pixinfo.indices[0] = izpix;
		pixinfo.indices[1] = iypix;
		pixinfo.indices[2] = 0;
		pixinfo.trackid = trackid;
		pixinfo.edep = frac_edep;
		pixinfo.t = pt0[3];
		std::set<EDepInfo_t> edep_v;
		EDepInfo_t edepinfo;
		edepinfo.pixid = (unsigned long)ipix;
		edepinfo.trackid = pixinfo.trackid;
		edepinfo.edep_sum += frac_edep;
		edep_v.insert( edepinfo );
		_edep_map.push_back( edep_v );
		pixinfo.container_index = _edep_map.size()-1;
		_pixinfo_map[ ipix ] = pixinfo;
	      }
	      else {
		// update
		PixInfo_t& pixinfo = it_pix->second;
		std::set<EDepInfo_t>& edep_v = _edep_map.at( pixinfo.container_index );
		pixinfo.edep += frac_edep;
		//pixinfo.t = std::min( pixinfo.t, midpt[3] );
		pixinfo.t = std::min( pixinfo.t, pt0[3] );
		auto it_edep = edep_v.find( EDepInfo_t(trackid) );
		if ( it_edep==edep_v.end() ) {
		  EDepInfo_t edepinfo;
		  edepinfo.pixid = (unsigned long)ipix;
		  edepinfo.trackid = trackid;
		  edepinfo.edep_sum = frac_edep;
		  edep_v.insert( edepinfo );
		}
		else {
		  (*it_edep).edep_sum += frac_edep;
		}
	      }//end of if pixinfo found for pixel id
	    }//end of if pixel inside image
	  }//end of kernel y-dim (y)
	}//end of kernel x-dim (z)
	
      }//loop over ndivisions
      
    }//loop over hits

    _readout_image_isgood = true;
    return _readout_image_isgood;
  }

  // XZ projection: det-x (drift) -> rows, det-z (beam) -> cols, det-y -> diffusion depth.
  // Identical pipeline to processSegmentHits; only the axis assignments differ.
  bool EDepSimInterface::processSegmentHitsXZ( const TG4HitSegmentContainer& hit_container )
  {
    _readout_image_isgood = false;

    size_t nhits = hit_container.size();

    size_t totpixels = npixels[0]*npixels[1];
    if ( _readout_image.size()!=totpixels )
      _readout_image.resize( totpixels, 0.0 );

    clear();
    std::memset( _readout_image.data(), 0.0, sizeof(float)*totpixels );
    std::cout << "making XZ image with " << _readout_image.size() << " pixels" << std::endl;

    if ( hit_container.size()==0 )
      return _readout_image_isgood;

    auto const first_hit = hit_container.at(0);
    std::vector<float> first_edep(4,0);
    for (int i=0; i<3; i++)
      first_edep[i] = (first_hit.GetStart()[i])*0.1;
    first_edep[3] = first_hit.GetStart()[3];

    // XZ: depth axis is det-y, rows are det-x, cols are det-z
    _readout_image_origin_cm[0] = first_edep[1]-distance_to_readout_plane; // y -> depth
    _readout_image_origin_cm[1] = origin[1]+first_edep[0];                  // x -> rows
    _readout_image_origin_cm[2] = origin[0]+first_edep[2];                  // z -> cols

    for (size_t ihit=0; ihit<nhits; ihit++) {

      auto const& hit = hit_container.at(ihit);

      std::vector<float> start(4);
      std::vector<float> stop(4);
      std::vector<float> stepdir(4);
      for (int i=0; i<3; i++) {
	start[i] = hit.GetStart()[i]*0.1;
	stop[i]  = hit.GetStop()[i]*0.1;
	stepdir[i] = stop[i]-start[i];
      }
      start[3] = hit.GetStart()[3];
      stop[3]  = hit.GetStop()[3];
      stepdir[3] = stop[3]-start[3];

      float steplen = 0.0;
      for (int i=0; i<3; i++)
	steplen += stepdir[i]*stepdir[i];
      steplen = sqrt(steplen);

      for (int i=0; i<3; i++)
	stepdir[i] /= steplen;

      float Estep = hit.GetEnergyDeposit();

      int ndivs = steplen/max_step_size + 1;
      float divlen = steplen/float(ndivs);

      for (int istep=0; istep<ndivs; istep++) {

	float fstep = (float)istep;
	std::vector<float> pt0(4,0.0);
	std::vector<float> pt1(4,0.0);
	std::vector<float> midpt(4,0.0);
	for (int i=0; i<4; i++) {
	  pt0[i] = start[i] + divlen*fstep*stepdir[i];
	  pt1[i] = start[i] + divlen*(fstep+1.0)*stepdir[i];
	  midpt[i] = 0.5*pt0[i] + 0.5*pt1[i];
	}

	std::vector<float> planepos(2,0);
	planepos[0] = (midpt[2]-_readout_image_origin_cm[2]); // z -> cols  (same as XY)
	planepos[1] = (midpt[0]-_readout_image_origin_cm[1]); // x -> rows  (was y in XY)
	int izpix = planepos[0]/pixelsize[0];
	int iypix = planepos[1]/pixelsize[1];
	std::vector<float> pixcenter(2,0);
	pixcenter[0] = (float(izpix)+0.5)*pixelsize[0];
	pixcenter[1] = (float(iypix)+0.5)*pixelsize[1];

	float dz_offset = planepos[0]-pixcenter[0];
	float dy_offset = planepos[1]-pixcenter[1];
	float depth = midpt[1]-_readout_image_origin_cm[0]; // y -> depth (was x in XY)
	float sigma = (depth)*(0.25/256.0) + min_ionization_sigma_width;

	if ( depth<0.0 || depth>=max_depth )
	  continue;

	float sum = 0.0;
	float kernel[7][7] = {0.0};
	for (int ikz=-3; ikz<=3; ikz++) {
	  for (int iky=-3; iky<=3; iky++) {
	    float kpixz = float(ikz)*pixelsize[0];
	    float kpixy = float(iky)*pixelsize[1];
	    float dist2 = (kpixz-dz_offset)*(kpixz-dz_offset) + (kpixy-dy_offset)*(kpixy-dy_offset);
	    float arg = -0.5*dist2/(sigma*sigma);
	    float w = exp(arg);
	    kernel[3+ikz][3+iky] = w;
	    sum += w;
	  }
	}

	int trackid = hit.GetContributors().at(0);
	float Efrac = Estep/float(ndivs);
	for (int ikz=-3; ikz<=3; ikz++) {
	  for (int iky=-3; iky<=3; iky++) {
	    int ipix = (iypix+iky)*npixels[0] + (izpix+ikz);
	    if ( ipix>=0 && ipix<(int)_readout_image.size() ) {
	      float frac_edep = Efrac*kernel[3+ikz][3+iky]/sum;
	      _readout_image[ipix] += frac_edep;
	      auto it_pix = _pixinfo_map.find( (unsigned long)ipix );
	      if ( it_pix==_pixinfo_map.end() ) {
		PixInfo_t pixinfo;
		pixinfo.pixid = (unsigned long)ipix;
		pixinfo.indices[0] = izpix;
		pixinfo.indices[1] = iypix;
		pixinfo.indices[2] = 0;
		pixinfo.trackid = trackid;
		pixinfo.edep = frac_edep;
		pixinfo.t = pt0[3];
		std::set<EDepInfo_t> edep_v;
		EDepInfo_t edepinfo;
		edepinfo.pixid = (unsigned long)ipix;
		edepinfo.trackid = pixinfo.trackid;
		edepinfo.edep_sum += frac_edep;
		edep_v.insert( edepinfo );
		_edep_map.push_back( edep_v );
		pixinfo.container_index = _edep_map.size()-1;
		_pixinfo_map[ ipix ] = pixinfo;
	      }
	      else {
		PixInfo_t& pixinfo = it_pix->second;
		std::set<EDepInfo_t>& edep_v = _edep_map.at( pixinfo.container_index );
		pixinfo.edep += frac_edep;
		pixinfo.t = std::min( pixinfo.t, pt0[3] );
		auto it_edep = edep_v.find( EDepInfo_t(trackid) );
		if ( it_edep==edep_v.end() ) {
		  EDepInfo_t edepinfo;
		  edepinfo.pixid = (unsigned long)ipix;
		  edepinfo.trackid = trackid;
		  edepinfo.edep_sum = frac_edep;
		  edep_v.insert( edepinfo );
		}
		else {
		  (*it_edep).edep_sum += frac_edep;
		}
	      }
	    }
	  }
	}

      }//loop over ndivisions
    }//loop over hits

    _readout_image_isgood = true;
    return _readout_image_isgood;
  }


  TH2D* EDepSimInterface::makeWholeDetectorTH2D( const TG4HitSegmentContainer& hit_container )
  {
    
    TH2D* himage = new TH2D("himage", "", npixels[0], origin[0], max_gridpt[0], npixels[1], origin[1], max_gridpt[1] );
    if ( !_readout_image_isgood )
      return himage;
    
    for (int iiz=0; iiz<npixels[0]; iiz++) {
      for (int iiy=0; iiy<npixels[1]; iiy++) {
	int ipix = iiy*npixels[0] + iiz;
	himage->SetBinContent( iiz+1, iiy+1, _readout_image[ipix] );
      }
    }
    
    return himage;

  }

  PyObject* EDepSimInterface::makeNumpyArrayCrop( const TG4HitSegmentContainer& hit_container, int img_pixdim, 
						  int offset_x_pixels, int offset_y_pixels,
						  float threshold, int rand_pix_from_center )
						  
  {
    
    if ( !__loaded_numpy )  {
      import_array1(0);
      __loaded_numpy = true;
    }

    if (!_readout_image_isgood) {
      Py_INCREF(Py_None);
      return Py_None;
    }
    
    // make the full image
    // returns an unrolled 2D image
    int orig_ncols = npixels[0];
    int orig_nrows = npixels[1];
        
    // crop out near the center
    std::vector<float> data( img_pixdim*img_pixdim, 0.0 );
    std::vector<int>   trackidimage( img_pixdim*img_pixdim, 0 );
    //int ipix = (iypix+iky)*npixels[0] + (izpix+ikz);

    std::vector< PixInfo_t > pix_v;
    pix_v.reserve( long(0.5*img_pixdim*img_pixdim) );
    float pix_edep_sum = 0.0;
    
    for (int irow=0; irow<img_pixdim; irow++) {
      // we copy npixels[0] at a given row
      int crop_row_start = orig_nrows/2-img_pixdim/2+irow-offset_y_pixels;
      int crop_col_start = orig_ncols/2-img_pixdim/2-offset_x_pixels;
      // void* memcpy( void* dest, const void* src, std::size_t count );
      // std::cout << "copy row[" << irow << "]" << std::endl;
      // std::cout << "  dest location: " << irow*img_pixdim << std::endl;
      // std::cout << "  src location: " << crop_row_start*npixels[0] + crop_col_start << std::endl;      
      ///std::memcpy( data.data() + irow*img_pixdim,  vdata.data() + crop_row_start*npixels[0] + crop_col_start, sizeof(double)*img_pixdim );
      for (int icol=0; icol<img_pixdim; icol++) {
	
	unsigned long pixid = crop_row_start*orig_ncols + crop_col_start + icol;
	float edep = _readout_image[ pixid ];
	if ( edep<threshold )
	  continue;
	
	data[ icol*img_pixdim + irow ] = edep;
	pix_edep_sum += edep;

	auto it_pix = _pixinfo_map.find( pixid );
	PixInfo_t croppix = it_pix->second;   // copy
	croppix.pixid = icol*img_pixdim+irow; // set to crop pixel id
	pix_v.emplace_back( std::move(croppix) );
	if ( it_pix!=_pixinfo_map.end() ) {
	  auto edep_v = _edep_map.at( it_pix->second.container_index );
	  int max_id = -1;
	  float max_edep = 0;
	  for ( auto& edepinfo : edep_v ) {
	    if ( edepinfo.edep_sum > max_edep ) {
	      max_id = edepinfo.trackid;
	      max_edep = edepinfo.edep_sum;
	    }
	  }
	  if ( max_id>=0 ) {
	    trackidimage[ icol*img_pixdim + irow ] = max_id+1;
	  }
	}//end of if pixel has info
	else {
	  if ( edep>0.0 ) 
	    std::cout << "No pixelinfo for pixid=" << pixid << " (" << icol << "," << irow << ") edep=" << edep << std::endl;
	}
      }//end of col loop
    }//end of row loop

    std::cout << "sort pix_v by time" << std::endl;
    SortByTime_t pix_sort_fn;
    std::sort( pix_v.begin(), pix_v.end(), pix_sort_fn );
    std::vector<float> quantile_points(21,0.0);
    float t0 = pix_v.front().t;
    // quantiles using charge
    float cdf_tot = pix_edep_sum;    
    // quantiles using time
    //cdf_tot = pix_v.back().t-pix_v.front().t;    

    // make quantile intervals  
    float quantile_binwidth  = cdf_tot/float(num_edep_quantiles);
    float next_quantile_threshold = quantile_binwidth;
    int current_quantile = 1;
    float cumulative_edep = 0.;
    quantile_points[0] = pix_v.front().t;
    quantile_points[20] = pix_v.back().t;
    std::vector< int > time_quantile_mask( img_pixdim*img_pixdim, 0);
    std::cout << "quantiles: time boundaries" << std::endl;
    int npix = 0;
    for ( auto& pixinfo : pix_v ) {
      cumulative_edep += pixinfo.edep; // use edep
      //cumulative_edep += pixinfo.t-t0; // using time 
      npix += 1;
      if ( current_quantile<21 && cumulative_edep > next_quantile_threshold ) {
	// pass threshold, move to next edep quantile
	quantile_points[current_quantile] = pixinfo.t;
	std::cout << "  [" << current_quantile << "]: time-bounds=[" << quantile_points[current_quantile-1] << ", " << quantile_points[current_quantile] << "] npixels=" << npix << " CDF(edep | t): " << cumulative_edep/cdf_tot  << std::endl;
	current_quantile += 1;
	next_quantile_threshold += quantile_binwidth;
	npix = 0;
      }
      time_quantile_mask[ pixinfo.pixid ] = current_quantile;
    }

    std::cout << "make arrays" << std::endl;

    int ndims = 2;
    npy_intp* dims = new npy_intp[2];
    dims[0] = img_pixdim;
    dims[1] = img_pixdim;    
    PyArrayObject* array = (PyArrayObject*)PyArray_SimpleNew( ndims, &dims[0], NPY_FLOAT );
    float* np_data = (float*)PyArray_DATA( array );
    memcpy( np_data, data.data(), sizeof(float)*data.size());
    
    PyArrayObject* array_trackid = (PyArrayObject*)PyArray_SimpleNew( ndims, &dims[0], NPY_INT );
    float* np_data_trackid = (float*)PyArray_DATA( array_trackid );
    memcpy( np_data_trackid, trackidimage.data(), sizeof(int)*trackidimage.size());

    PyArrayObject* array_timestep = (PyArrayObject*)PyArray_SimpleNew( ndims, &dims[0], NPY_INT );
    float* np_data_timestep = (float*)PyArray_DATA( array_timestep );
    memcpy( np_data_timestep, time_quantile_mask.data(), sizeof(int)*time_quantile_mask.size());
    
    delete [] dims;

    PyObject *d = PyDict_New();
    PyObject *key_edep     = Py_BuildValue("s", "edep" );
    PyObject *key_trackid  = Py_BuildValue("s", "trackid");
    PyObject *key_timestep = Py_BuildValue("s", "timestepmask");
    PyDict_SetItem( d, key_edep,     (PyObject*)array );
    PyDict_SetItem( d, key_trackid,  (PyObject*)array_trackid );
    PyDict_SetItem( d, key_timestep, (PyObject*)array_timestep );

    return (PyObject*)d;
    
  }

  
}

// PyObject* EDepSimInterface::makeNumpyArrayCrop( const TG4HitSegmentContainer& hit_container, int img_pixdim, 
//                       int offset_x_pixels, int offset_y_pixels,
//                       float threshold, int rand_pix_from_center )
                      
// {
  
//   if ( !__loaded_numpy )  {
//     import_array1(0);
//     __loaded_numpy = true;
//   }

//   if (!_readout_image_isgood) {
//     Py_INCREF(Py_None);
//     return Py_None;
//   }
  
//   // make the full image
//   // returns an unrolled 2D image
//   int orig_ncols = npixels[0];
//   int orig_nrows = npixels[1];
      
//   // crop out near the center
//   std::vector<float> data( img_pixdim*img_pixdim, 0.0 );
//   std::vector<int>   trackidimage( img_pixdim*img_pixdim, 0 );
//   //int ipix = (iypix+iky)*npixels[0] + (izpix+ikz);

//   std::vector< PixInfo_t > pix_v;
//   pix_v.reserve( long(0.5*img_pixdim*img_pixdim) );
//   float pix_edep_sum = 0.0;
  
//   for (int irow=0; irow<img_pixdim; irow++) {
//     // we copy npixels[0] at a given row
//     int crop_row_start = orig_nrows/2-img_pixdim/2+irow-offset_y_pixels;
//     int crop_col_start = orig_ncols/2-img_pixdim/2-offset_x_pixels;
//     // void* memcpy( void* dest, const void* src, std::size_t count );
//     // std::cout << "copy row[" << irow << "]" << std::endl;
//     // std::cout << "  dest location: " << irow*img_pixdim << std::endl;
//     // std::cout << "  src location: " << crop_row_start*npixels[0] + crop_col_start << std::endl;      
//     ///std::memcpy( data.data() + irow*img_pixdim,  vdata.data() + crop_row_start*npixels[0] + crop_col_start, sizeof(double)*img_pixdim );
//     for (int icol=0; icol<img_pixdim; icol++) {
  
//   unsigned long pixid = crop_row_start*orig_ncols + crop_col_start + icol;
//   float edep = _readout_image[ pixid ];
//   if ( edep<threshold )
//     continue;
  
//   data[ icol*img_pixdim + irow ] = edep;
//   pix_edep_sum += edep;

//   auto it_pix = _pixinfo_map.find( pixid );
//   PixInfo_t croppix = it_pix->second;   // copy
//   croppix.pixid = icol*img_pixdim+irow; // set to crop pixel id
//   pix_v.emplace_back( std::move(croppix) );
//   if ( it_pix!=_pixinfo_map.end() ) {
//     auto edep_v = _edep_map.at( it_pix->second.container_index );
//     int max_id = -1;
//     float max_edep = 0;
//     for ( auto& edepinfo : edep_v ) {
//       if ( edepinfo.edep_sum > max_edep ) {
//         max_id = edepinfo.trackid;
//         max_edep = edepinfo.edep_sum;
//       }
//     }
//     if ( max_id>=0 ) {
//       trackidimage[ icol*img_pixdim + irow ] = max_id+1;
//     }
//   }//end of if pixel has info
//   else {
//     if ( edep>0.0 ) 
//       std::cout << "No pixelinfo for pixid=" << pixid << " (" << icol << "," << irow << ") edep=" << edep << std::endl;
//   }
//     }//end of col loop
//   }//end of row loop

//   std::cout << "sort pix_v by time" << std::endl;
//   SortByTime_t pix_sort_fn;
//   std::sort( pix_v.begin(), pix_v.end(), pix_sort_fn );
//   std::vector<float> quantile_points(21,0.0);
//   float t0 = pix_v.front().t;
//   // quantiles using charge
//   float cdf_tot = pix_edep_sum;    
//   // quantiles using time
//   //cdf_tot = pix_v.back().t-pix_v.front().t;    

//   // make quantile intervals  
//   float quantile_binwidth  = cdf_tot/float(num_edep_quantiles);
//   float next_quantile_threshold = quantile_binwidth;
//   int current_quantile = 1;
//   float cumulative_edep = 0.;
//   quantile_points[0] = pix_v.front().t;
//   quantile_points[20] = pix_v.back().t;
//   std::vector< int > time_quantile_mask( img_pixdim*img_pixdim, 0);
//   std::cout << "quantiles: time boundaries" << std::endl;
//   int npix = 0;
//   for ( auto& pixinfo : pix_v ) {
//     cumulative_edep += pixinfo.edep; // use edep
//     //cumulative_edep += pixinfo.t-t0; // using time 
//     npix += 1;
//     if ( current_quantile<21 && cumulative_edep > next_quantile_threshold ) {
//   // pass threshold, move to next edep quantile
//   quantile_points[current_quantile] = pixinfo.t;
//   std::cout << "  [" << current_quantile << "]: time-bounds=[" << quantile_points[current_quantile-1] << ", " << quantile_points[current_quantile] << "] npixels=" << npix << " CDF(edep | t): " << cumulative_edep/cdf_tot  << std::endl;
//   current_quantile += 1;
//   next_quantile_threshold += quantile_binwidth;
//   npix = 0;
//     }
//     time_quantile_mask[ pixinfo.pixid ] = current_quantile;
//   }

//   std::cout << "make arrays" << std::endl;

//   int ndims = 2;
//   npy_intp* dims = new npy_intp[2];
//   dims[0] = img_pixdim;
//   dims[1] = img_pixdim;    
  
//   // Create numpy arrays
//   PyArrayObject* array = (PyArrayObject*)PyArray_SimpleNew( ndims, &dims[0], NPY_FLOAT );
//   if (!array) {
//     delete [] dims;
//     return NULL;
//   }
//   float* np_data = (float*)PyArray_DATA( array );
//   memcpy( np_data, data.data(), sizeof(float)*data.size());
  
//   PyArrayObject* array_trackid = (PyArrayObject*)PyArray_SimpleNew( ndims, &dims[0], NPY_INT );
//   if (!array_trackid) {
//     Py_DECREF(array);
//     delete [] dims;
//     return NULL;
//   }
//   int* np_data_trackid = (int*)PyArray_DATA( array_trackid );  // FIXED: was float*, should be int*
//   memcpy( np_data_trackid, trackidimage.data(), sizeof(int)*trackidimage.size());

//   PyArrayObject* array_timestep = (PyArrayObject*)PyArray_SimpleNew( ndims, &dims[0], NPY_INT );
//   if (!array_timestep) {
//     Py_DECREF(array);
//     Py_DECREF(array_trackid);
//     delete [] dims;
//     return NULL;
//   }
//   int* np_data_timestep = (int*)PyArray_DATA( array_timestep );  // FIXED: was float*, should be int*
//   memcpy( np_data_timestep, time_quantile_mask.data(), sizeof(int)*time_quantile_mask.size());
  
//   delete [] dims;

//   // Create dictionary and add arrays
//   PyObject *d = PyDict_New();
//   if (!d) {
//     Py_DECREF(array);
//     Py_DECREF(array_trackid);
//     Py_DECREF(array_timestep);
//     return NULL;
//   }
  
//   // Create string keys - these create new references
//   PyObject *key_edep     = Py_BuildValue("s", "edep" );
//   PyObject *key_trackid  = Py_BuildValue("s", "trackid");
//   PyObject *key_timestep = Py_BuildValue("s", "timestepmask");
  
//   if (!key_edep || !key_trackid || !key_timestep) {
//     Py_XDECREF(key_edep);
//     Py_XDECREF(key_trackid);
//     Py_XDECREF(key_timestep);
//     Py_DECREF(d);
//     Py_DECREF(array);
//     Py_DECREF(array_trackid);
//     Py_DECREF(array_timestep);
//     return NULL;
//   }
  
//   // PyDict_SetItem steals a reference to the values but not the keys
//   // So we need to DECREF the keys after setting them
//   if (PyDict_SetItem( d, key_edep,     (PyObject*)array ) < 0 ||
//       PyDict_SetItem( d, key_trackid,  (PyObject*)array_trackid ) < 0 ||
//       PyDict_SetItem( d, key_timestep, (PyObject*)array_timestep ) < 0) {
//     Py_DECREF(key_edep);
//     Py_DECREF(key_trackid);
//     Py_DECREF(key_timestep);
//     Py_DECREF(d);
//     Py_DECREF(array);
//     Py_DECREF(array_trackid);
//     Py_DECREF(array_timestep);
//     return NULL;
//   }
  
//   // Clean up key references since PyDict_SetItem doesn't steal them
//   Py_DECREF(key_edep);
//   Py_DECREF(key_trackid);
//   Py_DECREF(key_timestep);
  
//   // Clean up array references since PyDict_SetItem steals them
//   Py_DECREF(array);
//   Py_DECREF(array_trackid);
//   Py_DECREF(array_timestep);

//   return (PyObject*)d;
// }
}
