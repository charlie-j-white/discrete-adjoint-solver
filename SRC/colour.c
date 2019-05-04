/*
 *     written by Charlie Anderson
 *     APR 2019
 */



#include "lodepng.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>





void encodeOneStep(const char* filename, const unsigned char* image, unsigned width, unsigned height) {
  /*Encode the image*/
  unsigned error = lodepng_encode32_file(filename, image, width, height);

  /*if there's an error, display it*/
  if(error) printf("error %u: %s\n", error, lodepng_error_text(error));
}





int colour_(int *nwp, double fluxjac[][*nwp]) {

    printf("           Call from the C program: Start colouring\n");

    int nw = *nwp;
    int i, j;
    double rf, gf, bf;
	
    unsigned width = nw, height = nw;
    unsigned char* image = malloc(width * height * 4 );
    const char* filename = "flux-jacobian.png";

    int nl = 1, R;
    int nt = nw + 2*nl;
    double Cval, jmax, jmin, coltol;






  
//  Note that C[i][j] is same as *(*(C+i)+j) 
    double *C[nt];
    for (i=0; i<nt; i++) 
         C[i] = (double *)malloc(nt * sizeof(double)); 
  
  


//  set boundary offset if needed

    R = nl;





//  initialise shading array with zeros
    for(j = 0; j < nw+2*nl; j++){
    for(i = 0; i < nw+2*nl; i++){
	
	C[i][j] = 0;
	
    }
    }





//  find max value in flux jacobian
    
    jmax = 0.0;
    jmin = 0.0;

    for(j = R; j < nw+nl; j++){
    for(i = R; i < nw+nl; i++){

        if (fluxjac[i-R][j-R] > jmax) {
	    jmax = fluxjac[i-R][j-R];
        }

        if (fluxjac[i-R][j-R] < jmin) {
	    jmin = fluxjac[i-R][j-R];
        }

    }
    }









//  initialise shading array with values from flux jacobian

    coltol = 0.0;

    for(j = R; j < nw+nl; j++){
    for(i = R; i < nw+nl; i++){
	
//	C[i][j] = pow(fabs(fluxjac[i-R][j-R])/jmax,0.05);

	if (fabs(fluxjac[i-R][j-R]) > coltol && fabs(fluxjac[i-R][j-R]) != 0.0) {
	    C[i][j] = 1.0; 
	} else {
	    C[i][j] = 0.0; 
	}
		
    }
    }






//  start the looping; assign desired colour to relevant pixel
    
    for(j = 0; j < height; j++){
    for(i = 0; i < width; i++){
	
	rf = 255, gf = 255, bf = 255;

	Cval = C[i+R][j+R];

	rf = rf*(1-Cval) + Cval*0;
	gf = gf*(1-Cval);
	bf = bf*(1-Cval);


	unsigned red = rf, green = gf, blue = bf;
	image[4 * width * j + 4 * i + 0] = red;
	image[4 * width * j + 4 * i + 1] = green;
	image[4 * width * j + 4 * i + 2] = blue;
	image[4 * width * j + 4 * i + 3] = 255;
	
    }
    }







//  run the formatting function
    
    encodeOneStep(filename, image, width, height);
    free(image);





//  finish off function

    printf("           Call from the C program: Finished colouring\n");

    return 0;


}
