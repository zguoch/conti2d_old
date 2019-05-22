//
//  Conti3D.h
//  Conti3D
//
//  Created by 郭志馗 on 2017/2/1.
//  Copyright © 2017年 Zhikui Guo. All rights reserved.
//

#ifndef Conti3D_h
#define Conti3D_h

#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
using namespace std;
#include <math.h>
#include "omp.h"
// #include <mkl.h>
// #include "errcheck.inc"
// #include "generatedata.inc"
#include "FFTN.h"
#include "MultiProgressBar.h"

#define PI 3.141592653

#define DWC_CGLS 1
#define DWC_TIKHONOV 2
#define DWC_INTEGRALITERATION 3
#define DWC_LANDWEBER 4

//global variables
//define text color
#define ERROR_COUT "["<<"\033[31mError"<<"\033[0m] "
#define WARN_COUT "["<<"\033[33mWarning"<<"\033[0m] "
#define ERROR_COUT "["<<"\033[31mError"<<"\033[0m] "
#define PURPLE "\033[35m"
#define RED "\033[31m"
#define GREEN "\033[32m"
#define YELLOW "\033[33m"
#define BLUE "\033[34m"
#define COLOR_DEFALUT "\033[0m"
#include <unistd.h>
#define MOVEUP(x) printf("\033[%dA", (x))
// 清除屏幕 
#define CLEAR() printf("\033[2J") 
// 上移光标 
#define MOVEUP(x) printf("\033[%dA", (x)) 
// 下移光标 
#define MOVEDOWN(x) printf("\033[%dB", (x)) 
// 左移光标 
#define MOVELEFT(y) printf("\033[%dD", (y)) 
// 右移光标 
#define MOVERIGHT(y) printf("\033[%dC",(y)) 
// 定位光标 
#define MOVETO(x,y) printf("\033[%d;%dH", (x), (y)) 
// 光标复位 
#define RESET_CURSOR() printf("\033[H") 
// 隐藏光标 
#define HIDE_CURSOR() printf("\033[?25l") 
// 显示光标 
#define SHOW_CURSOR() printf("\033[?25h") 
//反显
#define HIGHT_LIGHT() printf("\033[7m")
#define UN_HIGHT_LIGHT() printf("\033[27m")
////////////////////////////////////
struct GrdHead
{
    int cols, rows;
    double bounds[6];
};
// int Derivative1d_col(double* data,GrdHead grdhead,double *der);
// int Derivative1d_row(double* data,GrdHead grdhead,double *der);
// double* ReadGrd(char* filename,GrdHead& grdhead);
double* ReadGrd(string filename,GrdHead& grdhead,int extNum);
// bool SaveGrd(char* filename, GrdHead grdhead,double* data, bool savexxyz=false,bool isInfo=true);
bool SaveGrd(string filename, GrdHead grdhead,double* data,int extNum, bool savexxyz=false,bool isInfo=true);
void GetKernalMatrix(GrdHead grdhead, double* G, const double rph, int num_thread=4);
void GetKernalMatrix_new(GrdHead grdhead, double* G, const double rph, int num_thread=4);
double GetGij(const int i, const int j, double* firstRow, const GrdHead grdhead);
//1. UWC: plane to plane
int Getkernel_p2p(GrdHead grdhead, double rph, double** kernel, int num_thread);
int Getkernel_p2p_new(GrdHead grdhead, double rph, double* kernel_firstRow, int num_thread);
int Getkernel_p2p_new(GrdHead grdhead, double rph, double** kernel, int num_thread);
void UWC_p2p(string inputfilename,string outputfilename,
    double height1,double height2,int extNum, int num_thread, bool isProgress,bool useOldKernel,
    string filename_exact);
void UWC_p2p_f(string inputfilename,string outputfilename,double height1,double height2,
                int extNum,string filename_exact);
int UWC_p2p_f(double* inputdata,double* result,GrdHead grdhead,double rph);

//2. UWC: uneven surface to plane
// double* GetNz(GrdHead grdhead, double* terrain);
// int GetNxyz(GrdHead grdhead, double* terrain, double* nx,double* ny,double* nz);
int Getkernel_u2p(GrdHead grdhead, double* terrain1, double height2, double** kernel, int num_thread);
int Getkernel_u2p_new(GrdHead grdhead, double* terrain1, double height2,
                    double* nx,double* ny,double* nz, 
                    double** kernel, int num_thread);
int Getkernel_p2s_new(GrdHead grdhead, double h1,double* topo2, double** kernel, int num_thread);
// void UWC_u2p(string inputfilename,string outputfilename,
    // string terrainfile1, double height2, int num_thered, bool isProgress);
void GetPmnij(double* Pmnij,int rows,int cols,double dx,double dy,double rph,double xm,double ym);
void UWC(double* datain, double* dataout, GrdHead grdhead,double** G);
void UWC_Gij(double* b, double* G,double* x, GrdHead grdhead, int num_thread=1);
void UWC_Gji(double* b, double* G,double* x, GrdHead grdhead, int num_thread=1);
void UWC_Gij(double* b,double** G,double* x, int modelnum,int num_thread=1);
//compute: b=(GT*G+lamubda*I)*x == b=M*x
//正则化之后的复杂系数矩阵与向量的乘积
void UWC_G_CGLS_Tik(double* b,double** G,double* x, int modelnum,double lambda2,int num_thread=1);
void UWC_Gji(double* b,double** G,double* x, int modelnum,int num_thread=1);
// 3. DWC: Plane to plane
void DWC_u2p(string inputfilename,string outputfilename,double height1,string filename_topo2);
void DWC_p2p_f(string inputfilename,string outputfilename,double height1,
    double height2,int extNum,double TRP, int num_thread, bool isProgress,bool useOldKernel);
void DWC_p2p(string inputfilename,string outputfilename,double height1,
    double height2,int extNum,double DWC_parameter,int DWC_method, int num_thread, bool isProgress,bool useOldKernel,
    string filename_exact);
//4. UWC: Plane to surface
void UWC_p2s(string inputfilename,string outputfilename,
    double height1,string topoFile,int extNum, int num_thread, bool isProgress,
    bool useOldKernel,string filename_exact);
//5. DWC: surface to plane
void DWC_s2p(string inputfilename,string outputfilename,string topo1,
    double height2,int extNum,double DWC_parameter,int DWC_method, int num_thread, bool isProgress,bool useOldKernel,
    string filename_exact);
//不同方法的向下延拓求解函数
void DWC_Tikhonov_old(double* G_firstRow,double* dataout,double* indata,double TRP,int kmax,double daierta,GrdHead grdhead,int num_thread);
void DWC_p2p_CGLS(double* G_firstRow,double* x, double* b,GrdHead grdhead,int extNum,double delta,int num_thread);
void DWC_s2p_CGLS(double** G,double* x, double* b,GrdHead grdhead,int extNum,double delta,int num_thread);
void DWC_s2p_Tikhonov(double** G,double* x, double* b,GrdHead grdhead,int extNum,double lambda,int num_thread);
//积分迭代发
void DWC_s2p_ItegrationIter(double** G,double* x, double* b,
    GrdHead grdhead,int extNum,int num_thread,string outputfile,double iter_number,double* ExactSolution=NULL);
void DWC_p2p_ItegrationIter(double* G,double* x, double* b,
    GrdHead grdhead,int extNum,int num_thread,string outputfile,double iter_number,double* ExactSolution=NULL);
//landweber迭代法
void DWC_p2p_LandweberIter(double* G,double* x, double* b,
    GrdHead grdhead,int extNum,int num_thread,string outputfile,double iter_number,double* ExactSolution=NULL);
void DWC_s2p_LandweberIter(double** G,double* x, double* b,
    GrdHead grdhead,int extNum,int num_thread,string outputfile,double iter_number,double* ExactSolution=NULL);

double Norm2(double* x,const int num);
//分析文件名和路径
//input: figures/figure_uwc_p2p.ps
//output:figure_uwc_p2p
string Path_GetBaseName(string filepath);
//input: figures/figure_uwc_p2p.ps
//output:ps
string Path_GetExtName(string filepath);
//input: figures/figure_uwc_p2p.ps
//output:figures
string Path_GetPath(string filepath);
//input: figures/figure_uwc_p2p.ps
//output:figures/figure_uwc_p2p
string Path_GetFileName(string filepath);
//求结果的梯度的二范数，看看平滑情况
double Norm2_Gradient(double* result,GrdHead grdhead);
//保存中间结果为时间序列的vtk格式
int SaveGrd2VTK(string outputfile,GrdHead grdhead,double* data);
//start text
static void StartText()
{
    //30:黑  31:红  32:绿  33:黄  34:蓝色  35:紫色  36:深绿
    cout<<"\033[33m";       //开启黄色
    cout << "***************************************************\n";
    cout << "*                 program conti2d                 *\n";
    cout << "*                 ~~~~~~~ ~~~~~~~                 *\n";
    cout << "*                                                 *\n";
    cout << "* Analytical continuation of potential field data *\n";
    cout << "*  - Upward Continuation from plane to plane      *\n";
    cout << "*  - Upward Continuation from plane to surface    *\n";
    cout << "*  - Downward Continuation from plane to plane    *\n";
    cout << "*  - Downward Continuation from surface to plane  *\n";
    cout << "*                                                 *\n";
    cout << "* See head file for input/output details.         *\n";
    cout << "* (c) Zhikui Guo, CUG, Aug 2016, Hangzhou         *\n";
    cout << "*                                                 *\n";
    cout << "***************************************************\n";
    cout << "\n\n";
    cout<<"\033[0m";
                                                                                                                                                                                                                      
}
//static function
static void Text_Axis()
{
    cout<<"   /y"<<endl;
    cout<<"  /"<<endl;
    cout<<" /"<<endl;
    cout<<"/0____________x"<<endl;
    cout<<"|"<<endl;
    cout<<"|"<<endl;
    cout<<"|"<<endl;
    cout<<"|z(-)"<<endl;
}
static void helpINFO()
{
    string version="2.0";
    string author="Zhikui Guo";
    string locus="SIO, Hangzhou";
    unsigned int wordWidth=20;
    // time_t now=time(0);
    // char* now_str=ctime(&now);
    string now_str="Aug 8, 2016";

    //30:黑  31:红  32:绿  33:黄  34:蓝色  35:紫色  36:深绿
    cout<<"===================== conti2d ======================"<<endl;
    cout<<"Analytical continuation of potential field data"<<endl;
    cout<<setw(wordWidth)<<setiosflags(ios::left)<<"Author \033[032m"<<author<<"\033[0m"<<endl;
    cout<<setw(wordWidth)<<setiosflags(ios::left)<<"Locus \033[032m"<<locus<<"\033[0m"<<endl;
    cout<<setw(wordWidth)<<setiosflags(ios::left)<<"Date \033[032m"<<now_str<<"\033[0m"<<endl;
    cout<<setw(wordWidth)<<setiosflags(ios::left)<<"Version \033[032m"<<version<<"\033[0m"<<endl;
    Text_Axis();
    cout<<setw(wordWidth)<<setiosflags(ios::left)<<"Synopsis:\033[031m"
        <<"conti2d "
        <<GREEN<<"input.grd "
        <<PURPLE<<"-G"<<GREEN<<"output.grd "
        <<PURPLE<<"-E"<<GREEN<<"extBoundNum "
        <<PURPLE<<"-H"<<GREEN<<"level_outputdata|.grd "
        <<PURPLE<<"-T"<<GREEN<<"level_inputdata|.grd "
        <<PURPLE<<"-D"<<COLOR_DEFALUT<<"["
        <<GREEN<<"+T"<<COLOR_DEFALUT<<"Tik_Par"<<PURPLE<<"|"
        <<GREEN<<"+I"<<COLOR_DEFALUT<<"It_Par"<<PURPLE<<"|"
        <<GREEN<<"+L"<<COLOR_DEFALUT<<"Landweber_Par"<<PURPLE<<"|"
        <<GREEN<<"+C"<<COLOR_DEFALUT<<"CGLS_Par"<<PURPLE<<""
        <<COLOR_DEFALUT<<"] "
        <<PURPLE<<"-t"<<GREEN<<"threads "
        <<PURPLE<<"-f"<<GREEN<<"[FrequencyDomain] "
        <<COLOR_DEFALUT<<endl;
    cout<<"===================================================="<<endl<<endl;
    // exit(1);
}
static void StartText_artASCII()
{
    cout<<GREEN<<" ________          ________          ________           _________        ___           _______          ________     \n"
    <<"|\\   ____\\        |\\   __  \\        |\\   ___  \\        |\\___   ___\\     |\\  \\         /  ___  \\        |\\   ___ \\    \n"
    <<"\\ \\  \\___|        \\ \\  \\|\\  \\       \\ \\  \\\\ \\  \\       \\|___ \\  \\_|     \\ \\  \\       /__/|_/  /|       \\ \\  \\_|\\ \\   \n"
    <<" \\ \\  \\            \\ \\  \\\\\\  \\       \\ \\  \\\\ \\  \\           \\ \\  \\       \\ \\  \\      |__|//  / /        \\ \\  \\ \\\\ \\  \n"
    <<"  \\ \\  \\____        \\ \\  \\\\\\  \\       \\ \\  \\\\ \\  \\           \\ \\  \\       \\ \\  \\         /  /_/__        \\ \\  \\_\\\\ \\ \n"
    <<"   \\ \\_______\\       \\ \\_______\\       \\ \\__\\\\ \\__\\           \\ \\__\\       \\ \\__\\       |\\________\\       \\ \\_______\\\n"
    <<"    \\|_______|        \\|_______|        \\|__| \\|__|            \\|__|        \\|__|        \\|_______|        \\|_______|\n"
    <<COLOR_DEFALUT<<endl;    
}
//error info
static void OutputErrorinfo(string info)
{
    cout<<"\033[31m";       //开启红色
    cout<<info<<endl;
    cout<<"\033[0m";
    exit(0);
}
//warning info
static void OutputWarninginfo(string info)
{
    cout<<"\033[34m";       //开启蓝色
    cout<<info<<endl;
    cout<<"\033[0m";
}

#endif /* Conti3D_h */
