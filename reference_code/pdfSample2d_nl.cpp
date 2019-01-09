
//call using
// [xdraw, ydraw] = pdfSample2d_nl(xcent,ycent,pdf,[ndraws],[randseq]);
//
// the xcent and ycent vectors are the _centers_ of the sampling bins described by the pdf matrix;
// size(pdf,1) = numel(ycent) and size(pdf,2) = numel(xcent). if ndraws is not specified, it defaults to 1.
//
// the optional 5th argument is a set of random draws (if the built-in rand() sequence isn't desired; 
// examples include using quasi-random sequences, seeded sequences, or a sequence with a higher 
// degree of randomization). the randseq is a two-column array, with the first column providing 
// rand() values for the x direction and the second column for the y direction.
// note that if ndraws is greater than the number of columns in randseq, space-time will cease to occur
// in normal physics. to prevent that occurrance, ndraws is set to min(ndraws,rows(randseq)).
// randseq should contain only values from 0 to 1.
//

//the function assumes that xcent and ycent are equally spaced and monotonically increasing 
//(it only uses xcent[0],xcent[1],ycent[0],ycent[1] for spatial interpolation... if you want to be lazy when creating xcent
//and ycent)


#include "mex.h"
#include "math.h"



//borrowed ruthlessly from
//http://www.velocityreviews.com/forums/t611498-implementation-of-drand48-as-given-in-steve-summits-book.html
//

#define PRECISION 2.82e14

double drand48(void)
{
	double x = 0;
	double denom = RAND_MAX + 1.0;
	double need;

	for(need = PRECISION; need > 1; need /= (RAND_MAX + 1.0))
	{
		x += rand()/denom;
		denom *= RAND_MAX + 1.0 ;
	}

	return x;

}


void computeCDF(double *cdf, double *pdf, int startidx, int nx){
	for (int i=startidx; i<nx; i++)
		cdf[i+1] = cdf[i] + pdf[i];
}


void computeCDF2(double *cdf, double *pdf, int nx, int ny){
	
	for (int xidx=0; xidx<nx; xidx++){
		cdf[xidx*(ny+1)] = 0;
		for (int yidx=0; yidx<ny; yidx++)
			cdf[xidx*(ny+1) + yidx + 1] = cdf[xidx*(ny+1) + yidx] + pdf[xidx*ny + yidx];
	}
}



void sample1d_y(double *yvalue, int *xidx, double *pdf, double *y, int nx, int ny, int ndraws, double *randseq)
{

	double *cdf, yinit = y[0], dy = y[1]-y[0], cdfvalue, dyfrac;
	double *normval;
	int cdfidx, coffset;

	//create the column-wise cdfs
	cdf = (double*)mxMalloc(sizeof(double)*nx*(ny+1));
	computeCDF2(cdf, pdf, nx, ny);
	normval = (double*)mxMalloc(sizeof(double)*nx);
	for (int i=0; i<nx; i++)
		normval[i] = cdf[(i+1)*(ny+1)-1];

	//do the y-direction draws
	for (int i=0; i<ndraws; i++){
		//cdfvalue = drand48() * normval[xidx[i]];
		cdfvalue = randseq[i] * normval[xidx[i]];
		coffset = xidx[i]*(ny+1); //where the particular cdf starts
		cdfidx = 0; //TODO: set this to a starting index

		while (cdfvalue>cdf[coffset + cdfidx + 1])
			cdfidx++;

		dyfrac = (cdfvalue - cdf[coffset + cdfidx])/(cdf[coffset+cdfidx+1]-cdf[coffset+cdfidx]);
		yvalue[i] = yinit + dy*((double)cdfidx - 0.5 + dyfrac);
	}

	mxFree(cdf);
	mxFree(normval);
}

//given a 1d pdf, computes the cdf and does the inverse sampling
void sample1d_x(double *xvalue, int *xidx, double *pdf, double *x, int nx, int ndraws, double *randseq)
{
	int startidx=0, npdf, cdfidx;
	double *cdf, cdfvalue, normsum;
	double xinit = x[0], dx = x[1]-x[0], dxfrac;

	//find where the first non-zero entry occurs
	while ( (pdf[startidx]==0.0) && (startidx<(nx-2)) )
		startidx++;

	
	//compute the cdf
	cdf = (double*)mxCalloc(nx+1, sizeof(double));
	computeCDF(cdf, pdf, startidx, nx);
	normsum = cdf[nx]; //the total value in the cdf

	if (normsum == 0.0)
		mexErrMsgTxt("The PDF does not contain any values.");

	//do some draws of the inverse cdf
	for (int i=0; i<ndraws; i++) {
		//cdfvalue = drand48()*normsum; //draw a new cdf value
		cdfvalue = randseq[i]*normsum; 

		//determine which bin of the cdf this value lies within
		cdfidx = startidx;
		while (cdfvalue>cdf[cdfidx+1]) //note: cdf[end]=normsum, so this is
			cdfidx++;

		//using the cdf bin, interpolate the x position that this would have come from
		dxfrac = (cdfvalue - cdf[cdfidx])/pdf[cdfidx]; //pdf[cdfidx] = cdf[cdfidx+1]-cdf[cdfidx]
		xvalue[i] = xinit + dx*((double)cdfidx - 0.5 + dxfrac);
		xidx[i] = cdfidx;
	}
	
	mxFree(cdf);
}



void mexFunction( int nlhs, mxArray *plhs[], int nrhs, mxArray *prhs[])
{
	double *pdf, *x, *y, *suminy, *xdraw, *ydraw, *randseq;
	int ndraws = 1, nx, ny, seqrows, seqcols;
	int *xidx;
	bool randseqAllocated = false;

	if (nrhs<3)
		mexErrMsgTxt("pdfSample2d_nl(x,y,PDF,[ndraws],[randseqs]) requires at least three inputs.");

	x = mxGetPr(prhs[0]);
	y = mxGetPr(prhs[1]);
	pdf = mxGetPr(prhs[2]);
	nx = (int) mxGetN(prhs[2]);
	ny = (int) mxGetM(prhs[2]);


	if (nrhs>3)
		ndraws = (int) mxGetScalar(prhs[3]);
	
	if (nrhs>4) {
		randseq = mxGetPr(prhs[4]);
		seqrows = (int)mxGetM(prhs[4]);
		seqcols = (int)mxGetN(prhs[4]);
		if (seqcols!=2)
			mexErrMsgTxt("pdfSample2d_nl(x,y,PDF,[ndraws],[randseqs]) needs two columns in randseqs.");
		ndraws = (ndraws>seqrows ? seqrows : ndraws);
	}
	else {
		randseq = (double*) mxMalloc(sizeof(double)*ndraws*2);
		seqrows = ndraws;
		randseqAllocated = true;
		for (int i=0; i<ndraws*2; i++)
			randseq[i] = drand48();
	}



	//start by summing the columns of the PDF
	suminy = (double*) mxCalloc(nx, sizeof(double));
	for (int xidx = 0; xidx<nx; xidx++){
		for (int yidx = 0; yidx<ny; yidx++){
			suminy[xidx] += pdf[xidx*ny + yidx];
		}
	}



	//set up some vectors to remember values
	xidx = (int*) mxMalloc(sizeof(int)*ndraws);
	plhs[0] = mxCreateDoubleMatrix(ndraws,1,mxREAL);
	plhs[1] = mxCreateDoubleMatrix(ndraws,1,mxREAL);
	xdraw = mxGetPr(plhs[0]);
	ydraw = mxGetPr(plhs[1]);

	
	//do the x-direction draws, remembering which pdf column the draws came from
	sample1d_x(xdraw, xidx, suminy, x, nx, ndraws, randseq);

	//do the y-direction draws, using the x-direction info
	sample1d_y(ydraw, xidx, pdf, y, nx, ny, ndraws, &(randseq[seqrows]));

	//clean up
	mxFree(suminy);
	mxFree(xidx);
	if (randseqAllocated)
		mxFree(randseq);

	//done.
}