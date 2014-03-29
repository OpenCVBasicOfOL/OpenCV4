#include <cv.h>
#include <highgui.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI (4*atan(1.0))

double cosines(int x, int y, int m)
{
	double temp;
	temp = (((2 * x) + 1)*y*PI) / (2*m);
	return cos(temp);
}

double coefficient(int x)
{
	if(x == 0)
	{
		return sqrt((double)2)/(double)2;
	}
	else
	{
		return (double)1;
	}
}

void make_DCT(IplImage *orig, IplImage *result, int q[][8])
{
	int u, v, i, j, m, n;
	double temp;
	int DCT[256][256];

	//DCT
	for ( m = 0 ; m < orig->height ; m += 8)
	{
		for ( n = 0 ; n < orig->width ; n += 8)
		{

			for ( u = 0 ; u < 8 ; u ++ )
			{
				for ( v = 0 ; v < 8 ; v ++ )
				{
					temp = 0.0;
					for ( i = 0 ; i < 8 ; i ++ )
					{
						for ( j = 0 ; j < 8 ; j ++ )
						{
							CvScalar itensity = cvGet2D(orig, i+m, j+n);
							temp += cosines(i, u, 8) * cosines(j, v, 8) * (double)itensity.val[0];
						}
					}
					temp *= 2/(sqrt((double)8 * (double)8)) * coefficient(u) * coefficient(v);
					DCT[m + u][n + v] = (int)temp;

				}
			}
		}
	}
	

	//Quantization
	for ( m = 0 ; m < orig->height ; m += 8)
	{
		for ( n = 0 ; n < orig->width ; n += 8)
		{

			for ( u = 0 ; u < 8 ; u ++ )
			{
				for ( v = 0 ; v < 8 ; v ++ )
				{
					DCT[m + u][n + v] = DCT[m + u][n + v] / q[u][v];
				}
			}
		}
	}
	
	//De-Quantization
	for ( m = 0 ; m < orig->height ; m += 8)
	{
		for ( n = 0 ; n < orig->width ; n += 8)
		{

			for ( u = 0 ; u < 8 ; u ++ )
			{
				for ( v = 0 ; v < 8 ; v ++ )
				{
					DCT[m + u][n + v] = DCT[m + u][n + v] * q[u][v];
				}
			}
		}
	}
	

	//Reverse - DCT
	for ( m = 0 ; m < 256 ; m += 8)
	{
		for ( n = 0 ; n < 256 ; n += 8)
		{
			for ( i = 0 ; i < 8 ; i ++ )
			{
				for ( j = 0 ; j < 8 ; j ++ )
				{
					temp = 0.0;
					for ( u = 0 ; u < 8 ; u ++ )
					{
						for ( v = 0 ; v < 8 ; v ++ )
						{
							temp +=  2/(sqrt((double)8 * (double)8)) * coefficient(u) * coefficient(v) * cosines(i, u, 8) * cosines(j, v, 8) * DCT[u + m][v + n];
						}
					}
					
					cvSet2D(result, i + m, j + n, cvScalar((int)temp));

				}
			}
		}
	}
}

//차이를 구하는 함수
void make_err(IplImage *orig, IplImage *result, IplImage *err)
{
	int i, j;
	int temp;
	for ( i = 0 ; i < orig->height ; i ++ )
	{
		for ( j = 0 ; j < orig->width ; j ++ )
		{
			CvScalar itensity = cvGet2D(orig, i, j);
			CvScalar itensity1 = cvGet2D(result, i, j);
			temp = itensity.val[0] - itensity1.val[0];

			if ( temp < 0 )
			{
				temp *= -1;
			}

			cvSet2D(err , i , j, cvScalar(temp));
		}
	}
}


int main()
{
	IplImage *herrington = cvLoadImage("lena256gray.bmp");
	
	IplImage *DCT_Result = cvCreateImage(cvGetSize(herrington), 8, 1);//DCT 및 복원을 거친다
	IplImage *Err = cvCreateImage(cvGetSize(herrington), 8, 1);//차이영상

	IplImage *DCT_15 = cvCreateImage(cvGetSize(herrington), 8, 1);//DCT 및 복원을 거친다(1.5배)
	IplImage *DCT_15_Result = cvCreateImage(cvGetSize(herrington), 8, 1);//DCT 및 복원을 거친다(1.5배)
	IplImage *Err_15 = cvCreateImage(cvGetSize(herrington), 8, 1);//차이영상(1.5배)
	
	int quan[8][8] = {16, 11, 10, 16, 24, 40, 51, 61,
					12, 12, 14, 19, 26, 58, 60, 55,
					14, 13, 16, 24, 40, 57, 69, 56,
					14, 17, 22, 29, 51, 87, 80, 62,
					18, 22, 37, 56, 68, 109, 103, 77,
					24, 35, 55, 64, 81, 104, 113, 92,
					48, 64, 78, 87, 113, 121, 120, 101,
					72, 92, 95, 98, 112, 100, 103, 99};


	int quan_15[8][8];//1.5배 Quantization Matrix

	int i, j;

	for ( i = 0 ; i < 8 ; i ++ )
	{
		for ( j = 0 ; j < 8 ; j ++ )
		{
			quan_15[i][j] = (quan[i][j] * 1.5);
		}
	}
	
	
	make_DCT(herrington, DCT_Result, quan);
	make_err(herrington, DCT_Result, Err);
	make_DCT(herrington, DCT_15, quan_15);
	make_err(herrington, DCT_15, Err_15);
	cvSaveImage("Result_1_1.jpg", DCT_Result);
	cvSaveImage("Result_1_2.jpg", Err);
	cvSaveImage("Result_2_1.jpg", DCT_15);
	cvSaveImage("Result_2_2.jpg", Err_15);

	

	cvReleaseImage(&DCT_Result);
	cvReleaseImage(&Err);
	cvReleaseImage(&DCT_15);
	cvReleaseImage(&DCT_15_Result);
	cvReleaseImage(&Err_15);
	return 0;
}