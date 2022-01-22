#include <iostream>
#include <iomanip>
#include <set>
#include <cmath>
#include<stdio.h>
#include <stdlib.h>
#include <string>
#include <time.h>
#include<math.h>
#include <complex>
#include <bits/stdc++.h>
#include <vector>
#include <fstream>
#include <ctime>
#include <random>
#include <omp.h>
#include <limits>

using namespace std;

const double PI = 3.1415927;
const double M_2PI = 2*PI;
const double eps=1e-12;

typedef std::complex<double> DComplex;

//---------------------------------------------------------------------------
// useful for testing
 inline DComplex polinom_2(DComplex x, double a, double b)
 {
	 //Horner's scheme for x*x + a*x + b
	 return x * (x + a) + b;
 }

//---------------------------------------------------------------------------
// useful for testing
 inline DComplex polinom_3(DComplex x, double a, double b, double c)
 {
	 //Horner's scheme for x*x*x + a*x*x + b*x + c;
	 return x * (x * (x + a) + b) + c;
 }

//---------------------------------------------------------------------------
// useful for testing
 inline DComplex polinom_4(DComplex x, double a, double b, double c, double d)
 {
	 //Horner's scheme for x*x*x*x + a*x*x*x + b*x*x + c*x + d;
	 return x * (x * (x * (x + a) + b) + c) + d;
 }

//---------------------------------------------------------------------------
// solve cubic equation x^3 + a*x^2 + b*x + c
// x - array of size 3
// In case 3 real roots: => x[0], x[1], x[2], return 3
//         2 real roots: x[0], x[1],          return 2
//         1 real root : x[0], x[1] ± i*x[2], return 1
unsigned int solveP3(double *x,double a,double b,double c) {
	double a2 = a*a;
    	double q  = (a2 - 3*b)/9;
	double r  = (a*(2*a2-9*b) + 27*c)/54;
    	double r2 = r*r;
	double q3 = q*q*q;
	double A,B;
    	if(r2<q3)
    	{
    		double t=r/sqrt(q3);
    		if( t<-1) t=-1;
    		if( t> 1) t= 1;
    		t=acos(t);
    		a/=3; q=-2*sqrt(q);
    		x[0]=q*cos(t/3)-a;
    		x[1]=q*cos((t+M_2PI)/3)-a;
    		x[2]=q*cos((t-M_2PI)/3)-a;
    		return 3;
    	}
    	else
    	{
    		A =-pow(fabs(r)+sqrt(r2-q3),1./3);
    		if( r<0 ) A=-A;
    		B = (0==A ? 0 : q/A);

		a/=3;
		x[0] =(A+B)-a;
		x[1] =-0.5*(A+B)-a;
		x[2] = 0.5*sqrt(3.)*(A-B);
		if(fabs(x[2])<eps) { x[2]=x[1]; return 2; }

		return 1;
        }
}

//---------------------------------------------------------------------------
// Solve quartic equation x^4 + a*x^3 + b*x^2 + c*x + d
// (attention - this function returns dynamically allocated array. It has to be released afterwards)
DComplex* solve_quartic(double a, double b, double c, double d)
{
	double a3 = -b;
	double b3 =  a*c -4.*d;
	double c3 = -a*a*d - c*c + 4.*b*d;

	// cubic resolvent
	// y^3 − b*y^2 + (ac−4d)*y − a^2*d−c^2+4*b*d = 0

	double x3[3];
	unsigned int iZeroes = solveP3(x3, a3, b3, c3);

	double q1, q2, p1, p2, D, sqD, y;

	y = x3[0];
	// THE ESSENCE - choosing Y with maximal absolute value !
	if(iZeroes != 1)
	{
		if(fabs(x3[1]) > fabs(y)) y = x3[1];
		if(fabs(x3[2]) > fabs(y)) y = x3[2];
	}

	// h1+h2 = y && h1*h2 = d  <=>  h^2 -y*h + d = 0    (h === q)

	D = y*y - 4*d;
	if(fabs(D) < eps) //in other words - D==0
	{
		q1 = q2 = y * 0.5;
		// g1+g2 = a && g1+g2 = b-y   <=>   g^2 - a*g + b-y = 0    (p === g)
		D = a*a - 4*(b-y);
		if(fabs(D) < eps) //in other words - D==0
			p1 = p2 = a * 0.5;

		else
		{
			sqD = sqrt(D);
			p1 = (a + sqD) * 0.5;
			p2 = (a - sqD) * 0.5;
		}
	}
	else
	{
		sqD = sqrt(D);
		q1 = (y + sqD) * 0.5;
		q2 = (y - sqD) * 0.5;
		// g1+g2 = a && g1*h2 + g2*h1 = c       ( && g === p )  Krammer
		p1 = (a*q1-c)/(q1-q2);
		p2 = (c-a*q2)/(q1-q2);
	}

    DComplex* retval = new DComplex[4];

	// solving quadratic eq. - x^2 + p1*x + q1 = 0
	D = p1*p1 - 4*q1;
	if(D < 0.0)
	{
		retval[0].real( -p1 * 0.5 );
		retval[0].imag( sqrt(-D) * 0.5 );
		retval[1] = std::conj(retval[0]);
	}
	else
	{
		sqD = sqrt(D);
		retval[0].real( (-p1 + sqD) * 0.5 );
		retval[1].real( (-p1 - sqD) * 0.5 );
	}

	// solving quadratic eq. - x^2 + p2*x + q2 = 0
	D = p2*p2 - 4*q2;
	if(D < 0.0)
	{
		retval[2].real( -p2 * 0.5 );
		retval[2].imag( sqrt(-D) * 0.5 );
		retval[3] = std::conj(retval[2]);
	}
	else
	{
		sqD = sqrt(D);
		retval[2].real( (-p2 + sqD) * 0.5 );
		retval[3].real( (-p2 - sqD) * 0.5 );
	}

    return retval;
}

// Different Temperature Profiles for preparing different Synthetic Images
double get_temperature_profile (double pro_type,double T_zero,double rho)
{
    if(pro_type == 1)
    {
        return T_zero * (1 - pow(rho,2));
    }

    if(pro_type == 2)
    {
        return T_zero * pow((1 - pow(rho,2)), 2);
    }

    if(pro_type == 3)
    {
        return T_zero * (1 - pow(rho,8));
    }

    if(pro_type == 4)
    {
        return T_zero * exp(-1 * pow(((rho-0.5)/0.08), 2));
    }

    if(pro_type == 5)
    {
        return T_zero * exp(-1 * pow(((rho)/0.5),2)) * (1-pow(rho,10));
    }

    if(pro_type == 6)
    {
        return T_zero * exp(-1 * pow((rho/0.3),2));
    }

    if(pro_type == 7)
    {
        return T_zero * exp(-1 * pow((((rho-0.7)/0.05)),2));
    }

    if(pro_type == 8)
    {
        return T_zero * exp(-1 * pow((rho/0.25),2));
    }

    else
    {
        return T_zero * exp(-1 * pow((((rho-0.9)/0.05)),2));
    }
}

//Function for finding the coefficients of the quartic equation  
void function_coeffecients(double &x1, double &y1, double &z1, double &x2, double &y2, double &z2, double &a_e, double &ecc, double &r, double &a, double &b, double &c, double &d)
{
				double dx,dy,dz;
                dx = x2 - x1;
                dy = y2 - y1;
                dz = z2 - z1;

                double xc = 0, yc = 0, zc = 0;

				double b_e = a_e*ecc;
                double A,B,C,D,E,F;
                A = (dx*dx*b_e*b_e) + (dy*dy*b_e*b_e) + (dz*dz*a_e*a_e);
                B = (b_e*b_e*((2*x1*dx) - (2*xc*dx) + (2*y1*dy) - (2*yc*dy))) + (a_e*a_e*((2*z1*dz) - (2*zc*dz)));
                C = (b_e*b_e*(((x1 - xc)*(x1 - xc)) + ((y1 - yc)*(y1 - yc)))) + (a_e*a_e*(z1 - zc)*(z1 - zc)) + (b_e*b_e*r*r) - (a_e*a_e*b_e*b_e);
                D = (dx*dx) + (dy*dy);
                E = (2*x1*dx) - (2*xc*dx) + (2*y1*dy) - (2*yc*dy);
                F = ((x1 - xc)*(x1 - xc)) + ((y1 - yc)*(y1 - yc));

				//Variables that form the coefficients of the quartic equation for solving the equation of Line and Torus 
                a = (2*B)/A;
                b = ((B*B) + (2*A*C) - (4*r*r*b_e*b_e*b_e*b_e*D))/(A*A);
                c = ((2*B*C) - (4*r*r*b_e*b_e*b_e*b_e*E))/(A*A);
                d = ((C*C) - (4*r*r*b_e*b_e*b_e*b_e*F))/(A*A);
}

// Function for computing the intensity at each interim points (which would then be added for all such points on the LoS)
void function_interim(double &tt, double &x1, double &y1, double &z1, double &x2, double &y2, double &z2, double &r, double &a_e, double &pixel_intensity, double &t,double &pro_type,double &T_zero)
{
					int end = 1/t;
          			double temp_intensity = 0;    
					for(int cnt=1; cnt <= end; cnt++)
					{
						double a_x = x1 + (x2 - x1)*t*cnt, a_y = y1 + (y2 - y1)*t*cnt, a_z = z1 + (z2 - z1)*t*cnt;
						double theta = atan(a_y/a_x);

						if(a_x<0)
                        {
							theta = theta + PI;
                        }

						double avec_x = a_x - r*cos(theta), avec_y = a_y - r*sin(theta), avec_z = a_z;
						double length = sqrt((avec_x*avec_x) + (avec_y*avec_y) + (avec_z*avec_z));
						double dist = sqrt((a_x*a_x) + (a_y*a_y) + (a_z*a_z));
            			double rho = length/a_e;
						theta = (theta*180)/PI;

						
						if(rho > 1)
						{
							rho = 0.999999;
						}

                        double value_profile_fun = get_temperature_profile(pro_type, T_zero, rho);
                        temp_intensity += value_profile_fun;
					} 
            
					pixel_intensity = temp_intensity;  
}

int main()
{
	//Defining the different variables as per their physical significance
    double pinhole_x, pinhole_y , pinhole_z , camera_centre_x , camera_centre_y , camera_centre_z , area , pixels, r, a_e_or, ecc, num_of_torus, a_e_or_metal, r_metal, pro_type, T_zero;
  
	//Reading the input from the text file
	std::ifstream file_input("input.txt");

	//Reading a total of 15 different inputs from the 'input.txt' file
    string input = "";
	double input_arr[15] = {0};
	int input_counter = 0;

	string input_number = "";
    
	//Fetching variables by space-break from the input file
    while(std::getline(file_input, input))
    {
        for (auto character : input) 
	    {
            if (character == ' ')
            {
                input_arr[input_counter] = stof(input_number);
                input_number = "";
				input_counter ++;
            }

                                
            else
            {
                input_number = input_number + character;
            }
        }
    }

	input_arr[input_counter] = stof(input_number);

	//Variables are assigned to the respective value fetched from the input file (in '.txt' format)
	pinhole_x = input_arr[0];
	pinhole_y = input_arr[1];
	pinhole_z = input_arr[2];
	camera_centre_x = input_arr[3];
	camera_centre_y = input_arr[4];
	camera_centre_z = input_arr[5];
	area = input_arr[6];
	pixels = input_arr[7];
	r = input_arr[8];
	a_e_or = input_arr[9];
	ecc = input_arr[10];
	num_of_torus = input_arr[11];
    r_metal = input_arr[12];
	a_e_or_metal = input_arr[13];
    pro_type = input_arr[14];
    T_zero = input_arr[15];

	//For generating different random numbers every time
	srand(time(NULL));
  
	//Variable for storing the time elapsed of the entire code
   	double final_total_time = 0;

	//Clock timers for marking the start and end of the code
	struct timespec start, end;

    FILE *fp_static;

    fp_static = fopen("time_dynamic.txt", "a+");
    int flops_cnt = 0;

    
    int pixel_array[7] = {1024, 4096, 16384, 65536, 262144, 1048576, 4194304};

        pixels = pixel_array[4]; 

        //Chunk size for scheduling 
    	clock_gettime(CLOCK_MONOTONIC, &start);


        //Dimension of each small pixel is computed from the dimension of the camera screen and the number of pixels required in the synthetic image
        double side = sqrt(area)/sqrt(pixels);
        double dimn = sqrt(pixels);

        double x,y,z;
        double multiplier = (dimn - 1)/2;

        //Arrays for storing the center of each pixel is initialised dynamically (each array would store the respective co-ordinate location as per its name)
        double *x_sp = new double [int(pixels)], *y_sp = new double [int(pixels)], *z_sp = new double [int(pixels)], *intensity_value = new double [int(pixels)];

        //Co-ordinate arrays assigned with values as per the concept of image-flipping
        if(pinhole_x > 0)
        {
            x = camera_centre_x;
            y = camera_centre_y + multiplier*side;
            z = camera_centre_z - multiplier*side;
            
            for (int i = 0; i < (int)dimn ; i++)
            {
                for (int j = 0; j < (int)dimn ; j++)
                {
                        x_sp[int(i*dimn+j)] = x;
                        y_sp[int(i*dimn+j)] = y - j*side;
                        z_sp[int(i*dimn+j)] = z + i*side; 
                }
            }
        }

        else
        {
            x = camera_centre_x;
            y = camera_centre_y - multiplier*side;
            z = camera_centre_z + multiplier*side;

            for (int i = 0; i < (int)dimn ; i++)
            {
                for (int j = 0; j < (int)dimn ; j++)
                {
                        x_sp[int(i*dimn+j)] = x;
                        y_sp[int(i*dimn+j)] = y + j*side;
                        z_sp[int(i*dimn+j)] = z - i*side; 
                }
            }
        }

        //File pointer that would store the intensity of each pixel for plotting the synthetic image
        FILE *fp_intensity;

        //Variable to store the total time elapsed for the code
        final_total_time = 0;
        fp_intensity = fopen("pixels_intensity_txt.txt", "w+");

        //Clock timer started to record the code-execution time
	    clock_gettime(CLOCK_MONOTONIC, &start);

        for (int i = 0; i < (int)pixels; i++)
        {
                    double a, b, c, d, dx, dy, dz;
                    double a_e;

                    double pixel_intensity = 0;
                    double x_i = x_sp[i], y_i = y_sp[i], z_i = z_sp[i];

                    double x1 = pinhole_x , x2 = x_i , y1 = pinhole_y , y2 = y_i , z1 = pinhole_z , z2 = z_i;

                    a_e = a_e_or_metal;
                    function_coeffecients(x1, y1, z1, x2, y2, z2, a_e, ecc, r_metal, a, b, c, d);

                    std::complex<double>* solutions_metal = solve_quartic(a, b, c, d);

                    set<double> st_metal;

                    for(int k = 0; k < 4; k++)
                    {
                        if(fabs(solutions_metal[k].imag()) < eps && !std::isnan(solutions_metal[k].real())) st_metal.insert(solutions_metal[k].real());
                    }

            //If the number of intersection points of an LoS and the metal tokamak is not equal to 4
            if(st_metal.size() != 4)
            {	
                for (int j = 0; j < num_of_torus; j++)
                {
                    x1 = pinhole_x , x2 = x_i , y1 = pinhole_y , y2 = y_i , z1 = pinhole_z , z2 = z_i;

                    a_e = a_e_or - 0.1*j;

                    //Computing the coeeficients of the polynomial for equation solving the parametric equations of line and torus
                    function_coeffecients(x1, y1, z1, x2, y2, z2, a_e, ecc, r, a, b, c, d);
                    
                    //Quartic equation solver based on the coefficients of the polynomial terms
                    std::complex<double>* solutions = solve_quartic(a, b, c, d);

                    //'Set' data structure would store unique intersection points
                    set<double> st;

                    //parameter 't' in the equation of line gets stored in the set
                    for(int k = 0; k < 4; k++)
                    {
                        if(fabs(solutions[k].imag()) < eps && !std::isnan(solutions[k].real())) st.insert(solutions[k].real());
                    }

                    //Case for only 2 intersection points between LoS and the plasma Torus
                    if(st.size() == 2)
                    {
                        double minn = 0;
                        double maxx;

                        for(auto it = st.begin(); it != st.end(); it++)
                        {
                            minn = min(minn,*it);
                            maxx = max(minn,*it);
                        }

                        x2 = x1 + (x_i - x1)*minn;
                        y2 = y1 + (y_i - y1)*minn;
                        z2 = z1 + (z_i - z1)*minn;

                        x1 = x1 + (x_i - x1)*maxx;
                        y1 = y1 + (y_i - y1)*maxx;
                        z1 = z1 + (z_i - z1)*maxx;

                        dx = x2-x1;
                        dy = y2-y1;
                        dz = z2-z1;

                        double sigma = sqrt(dx*dx + dy*dy + dz*dz);
                        double t = 0.001/sigma;
                        double tt = t;

                        function_interim(tt, x1, y1, z1, x2, y2, z2, r, a_e, pixel_intensity, t, pro_type, T_zero);
                    }

                    //Case for 4 intersection points between LoS and the plasma Torus
                    else if(st.size() == 4)
                    {
                        double arr_t[4] = {0};
                        int cnt_i = 0;

                        for(auto it = st.begin(); it != st.end(); it++)
                        {
                            arr_t[cnt_i] = *it;
                            cnt_i++;
                        }

                        sort(arr_t, arr_t + 4);

                        double x2_temp = x2, y2_temp = y2, z2_temp = z2, x1_temp = x1, y1_temp = y1, z1_temp = z1;
                    
                        for(int loop_cnt = 0; loop_cnt <= 2; loop_cnt += 2)
                        {
                            double minn = arr_t[loop_cnt];
                            double maxx = arr_t[loop_cnt + 1];

                            x2 = x1_temp + (x_i - x1_temp)*minn;
                            y2 = y1_temp + (y_i - y1_temp)*minn;
                            z2 = z1_temp + (z_i - z1_temp)*minn;

                            x1 = x1_temp + (x_i - x1_temp)*maxx;
                            y1 = y1_temp + (y_i - y1_temp)*maxx;
                            z1 = z1_temp + (z_i - z1_temp)*maxx;

                            dx = x2-x1;
                            dy = y2-y1;
                            dz = z2-z1;

                            double sigma = sqrt(dx*dx + dy*dy + dz*dz);
                            double t = 0.001/sigma;
                            double tt = t;

                            function_interim(tt, x1, y1, z1, x2, y2, z2, r, a_e, pixel_intensity, t, pro_type, T_zero);
                        }
                    }
                    delete[] solutions;
                }
            }

            //If the number of intersection points of an LoS and the metal tokamak is equal to 4
            else
            {
                for (int j = 0; j < num_of_torus; j++)
                {
                    x1 = pinhole_x , x2 = x_i , y1 = pinhole_y , y2 = y_i , z1 = pinhole_z , z2 = z_i;

                    a_e = a_e_or - 0.1*j;

                    //Computing the coeeficients of the polynomial for equation solving the parametric equations of line and torus
                    function_coeffecients(x1, y1, z1, x2, y2, z2, a_e, ecc, r_metal, a, b, c, d);
                    
                    //Quartic equation solver based on the coefficients of the polynomial terms
                    std::complex<double>* solutions = solve_quartic(a, b, c, d);

                    //'Set' data structure would store unique intersection points
                    set<double> st;

                    //parameter 't' in the equation of line gets stored in the set
                    for(int k = 0; k < 4; k++)
                    {
                        if(fabs(solutions[k].imag()) < eps && !std::isnan(solutions[k].real())) st.insert(solutions[k].real());
                    }

                    //Case for only 2 intersection points between LoS and the plasma Torus
                    if(st.size() == 2)
                    {
                        double minn = 0;
                        double maxx;

                        for(auto it = st.begin(); it != st.end(); it++)
                        {
                            minn = min(minn,*it);
                            maxx = max(minn,*it);
                        }

                        x2 = x1 + (x_i - x1)*minn;
                        y2 = y1 + (y_i - y1)*minn;
                        z2 = z1 + (z_i - z1)*minn;

                        x1 = x1 + (x_i - x1)*maxx;
                        y1 = y1 + (y_i - y1)*maxx;
                        z1 = z1 + (z_i - z1)*maxx;


                        dx = x2-x1;
                        dy = y2-y1;
                        dz = z2-z1;

                        double sigma = sqrt(dx*dx + dy*dy + dz*dz);
                        double t = 0.001/sigma;
                        double tt = t;

                        function_interim(tt, x1, y1, z1, x2, y2, z2, r, a_e, pixel_intensity, t, pro_type, T_zero);
                    }
                    
                    //Case for 4 intersection points between LoS and the plasma Torus
                    else if(st.size() == 4)
                    {
                        double arr_t[4] = {0};
                        int cnt_i = 0;

                        for(auto it = st.begin(); it != st.end(); it++)
                        {
                            arr_t[cnt_i] = *it;
                            cnt_i++;
                        }

                        sort(arr_t, arr_t + 4);

                        double x2_temp = x2, y2_temp = y2, z2_temp = z2, x1_temp = x1, y1_temp = y1, z1_temp = z1;
                    
                        for(int loop_cnt = 2; loop_cnt <= 3; loop_cnt += 2)
                        {
                            double minn = arr_t[loop_cnt];
                            double maxx = arr_t[loop_cnt + 1];

                            x2 = x1_temp + (x_i - x1_temp)*minn;
                            y2 = y1_temp + (y_i - y1_temp)*minn;
                            z2 = z1_temp + (z_i - z1_temp)*minn;

                            x1 = x1_temp + (x_i - x1_temp)*maxx;
                            y1 = y1_temp + (y_i - y1_temp)*maxx;
                            z1 = z1_temp + (z_i - z1_temp)*maxx;


                            dx = x2-x1;
                            dy = y2-y1;
                            dz = z2-z1;

                            double sigma = sqrt(dx*dx + dy*dy + dz*dz);
                            double t = 0.001/sigma;
                            double tt = t;

                            function_interim(tt, x1, y1, z1, x2, y2, z2, r, a_e, pixel_intensity, t, pro_type, T_zero);
                        }
                    }
                    delete[] solutions;
                }
            }
                //Pixel intenisty values for each pixel gets stored in the array		
                intensity_value[i] = pixel_intensity;
        }

        //Storing the pixel intensity values in the file
        for(int i = 0; i < (int)pixels; i++)
        {
            fprintf(fp_intensity,"%lf\n", intensity_value[i]);
        }

        fclose(fp_intensity);

        
        clock_gettime(CLOCK_MONOTONIC, &end); 
    
        double time_taken; 
        time_taken = (end.tv_sec - start.tv_sec) * 1e9; 
        time_taken = (time_taken + (end.tv_nsec - start.tv_nsec)) * 1e-9; 
        
        final_total_time += time_taken;
        
        cout<<"Total Serial Time : "<<final_total_time<<" Pixels : "<<(int)sqrt(pixels)<<"\n";
}
