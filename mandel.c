/*
 
    Name: Steven Michael Kost
    ID:   1001388391
 
*/

#include "bitmap.h"
#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <errno.h>
#include <string.h>
#include <pthread.h>
#include <unistd.h>

int iteration_to_color(int i, int max);
int iterations_at_point(double x, double y, int max);
void compute_image(struct bitmap *bm, double xmin, double xmax, double ymin, double ymax, int max, int n, int *size, int *count);
void *threadA(void *);

typedef struct BITS
{
    struct bitmap *BITMAPS;
    double xvalmin;
    double xvalmax;
    double yvalmin;
    double yvalmax;
    int    maxy;
    int    Hei;
    int    Wid;
    int    countA;
    int    *sizeC;
} BITS;

void show_help()
{
	printf("Use: mandel [options]\n");
	printf("Where options are:\n");
	printf("-m <max>    The maximum number of iterations per point. (default=1000)\n");
	printf("-x <coord>  X coordinate of image center point. (default=0)\n");
	printf("-y <coord>  Y coordinate of image center point. (default=0)\n");
	printf("-s <scale>  Scale of the image in Mandlebrot coordinates. (default=4)\n");
	printf("-W <pixels> Width of the image in pixels. (default=500)\n");
	printf("-H <pixels> Height of the image in pixels. (default=500)\n");
    printf("-n <thread> Specified number of threads.\n");
    printf("-o <file>   Set output file. (default=mandel.bmp)\n");
	printf("-h          Show this help text.\n");
	printf("\nSome examples are:\n");
	printf("mandel -x -0.5 -y -0.5 -s 0.2\n");
	printf("mandel -x -.38 -y -.665 -s .05 -m 100\n");
	printf("mandel -x 0.286932 -y 0.014287 -s .0005 -m 1000\n\n");
}

int main(int argc, char *argv[])
{
	char c;

	// These are the default configuration values used
	// if no command line arguments are given.
	const char *outfile = "mandel.bmp";
	double xcenter = 0;
	double ycenter = 0;
	double scale = 4;
	int    image_width = 500;
	int    image_height = 500;
	int    max = 1000;
    int    n = 1;
    int    i = 0;
    int    *size;
    int    *count;

	// For each command line argument given,
	// override the appropriate configuration value.
	while((c = getopt(argc, argv, "x:y:s:W:H:m:n:o:h")) != -1)
    {
		switch(c)
        {
			case 'x':
				xcenter = atof(optarg);
				break;
			case 'y':
				ycenter = atof(optarg);
				break;
			case 's':
				scale = atof(optarg);
				break;
            case 'W':
				image_width = atoi(optarg);
				break;
			case 'H':
				image_height = atoi(optarg);
				break;
			case 'm':
				max = atoi(optarg);
				break;
            case 'n':
                n = atoi(optarg);
                size = malloc((n*sizeof(size))+1);
                count = malloc((n*sizeof(count))+1);
                for(i = 0; i < (n+1); i++)
                {
                    if(i == 0)
                    {
                        size[i] = (image_height/n)*i;
                        count[i] = i;
                    }
                    else if(i == n)
                    {
                        size[i] = image_height;
                        count[i] = i;
                    }
                    else
                    {
                        size[i] = ((image_height/n)*i);
                        count[i] = i;
                    }
                }
                break;
            case 'o':
				outfile = optarg;
				break;
			case 'h':
				show_help();
				exit(1);
				break;
		}
	}

	// Display the configuration of the image.
	printf("mandel: x=%lf y=%lf scale=%lf max=%d outfile=%s\n", xcenter, ycenter, scale, max, outfile);

	// Create a bitmap of the appropriate size.
    struct bitmap *bm = bitmap_create(image_width, image_height);

	// Fill it with a dark blue, for debugging
    bitmap_reset(bm, MAKE_RGBA(0,0,255,0));

	// Compute the Mandelbrot image
    compute_image(bm, xcenter-scale, xcenter+scale, ycenter-scale, ycenter+scale, max, n, size, count);

    // Save the image in the stated file.
	if(!bitmap_save(bm, outfile))
    {
		fprintf(stderr,"mandel: couldn't write to %s: %s\n", outfile, strerror(errno));
		return 1;
	}
    pthread_exit(NULL);
	return 0;
}

/*
Compute an entire Mandelbrot image, writing each point to the given bitmap.
Scale the image to the range (xmin-xmax,ymin-ymax), limiting iterations to "max"
*/
void compute_image(struct bitmap *bm, double xmin, double xmax, double ymin, double ymax, int max, int n, int *size, int *count)
{
    int i, j, k, rc;
    BITS plugh[50];
    int width  = bitmap_width(bm);
	int height = bitmap_height(bm);
    pthread_t  *thread;

    thread = (pthread_t *) malloc(n*sizeof(pthread_t));

	// For every pixel in the image...
    

    if(n <= 1)
    {
        for(j=0; j<height; j++)
        {
            for(i=0; i<width; i++)
            {
                // Determine the point in x,y space for that pixel.
                double x = xmin + i*(xmax-xmin)/width;
                double y = ymin + j*(ymax-ymin)/height;

                // Compute the iterations at that point.
                int iters = iterations_at_point(x, y, max);

                // Set the pixel in the bitmap.
                bitmap_set(bm, i, j, iters);
            }
        }
    }
    else if(n > 1)
    {
        for(k = 0; k < n; k++)
        {
            plugh[k].countA = k;
            plugh[k].BITMAPS = bm;
            plugh[k].xvalmin = xmin;
            plugh[k].xvalmax = xmax;
            plugh[k].yvalmin = ymin;
            plugh[k].yvalmax = ymax;
            plugh[k].maxy    = max;
            plugh[k].Hei     = height;
            plugh[k].Wid     = width;
            plugh[k].sizeC   = size;
            rc = pthread_create(&(thread[k]), NULL, threadA, (void*)&plugh[k]);
        }
        for(j = 0; j < n; j++)
        {
            pthread_join(thread[j], NULL);
        }
    }
}

void *threadA(void *ptr)
{
    int i, j;
    struct BITS *args = (struct BITS *)ptr;
    struct bitmap *bm = args->BITMAPS;
    double xmin = args->xvalmin;
    double xmax = args->xvalmax;
    double ymin = args->yvalmin;
    double ymax = args->yvalmax;
    int    max = args->maxy;
    int    height = args->Hei;
    int    width = args->Wid;
    int    *size = args->sizeC;
    int    count = args->countA;
    int    sizeA = size[count];
    int    sizeB = size[count+1];
    for(j = sizeA; j < sizeB; j++)
    {
        for(i = 0; i < width; i++)
        {
            // Determine the point in x,y space for that pixel.
            double x = xmin + i*(xmax-xmin)/width;
            double y = ymin + j*(ymax-ymin)/height;
            // Compute the iterations at that point.
            int iters = iterations_at_point(x, y, max);

            // Set the pixel in the bitmap.
            bitmap_set(bm, i, j, iters);
        }
    }
    pthread_exit(NULL);
}

/*
Return the number of iterations at point x, y
in the Mandelbrot space, up to a maximum of max.
*/
int iterations_at_point(double x, double y, int max)
{
	double x0 = x;
	double y0 = y;
	int iter = 0;
    while((x*x + y*y <= 4) && iter < max)
    {
		double xt = x*x - y*y + x0;
		double yt = 2*x*y + y0;

		x = xt;
		y = yt;

		iter++;
	}
    return iteration_to_color(iter, max);
}

/*
Convert a iteration number to an RGBA color.
Here, we just scale to gray with a maximum of imax.
Modify this function to make more interesting colors.
*/
int iteration_to_color(int i, int max)
{
    int gray = 255*i/max;
    return MAKE_RGBA(gray, gray, gray, 0);
}

