
#include "Conti2D.h"


bool isNum(string str)
{
    stringstream sin(str);
    double d;
    char c;
    if(!(sin >> d))
        return false;
    if (sin >> c)
        return false;
    return true;
}
void testOMP()
{
    int index=0;
    omp_set_num_threads(3);
    #pragma omp parallel for shared(index)
	for (int i = 0; i < 10; i++)
    {
        
        printf("i = %d, I am Thread %d\n", i, index);
        index++;
    }
		
}
int main(int argc, char* argv[])
{
    struct winsize w;
    ioctl(0, TIOCGWINSZ, &w);
    if(w.ws_col>119)
    {
        StartText_artASCII();
    }else
    {
        StartText();
    }
    helpINFO();
    
    //testOMP();
//********************************************
    string inputfilename;
    string outputfilename;
    string topo1;
    string topo2;
    double height1=0;
    double height2=0;
    bool isDownwardConti=false;
    bool isFrequency=false;
    // double TRP=0.00001;
    int num_inputdata=0;
    bool isplane1=true;  //0:plane to plane; 1:plane to uneven; 2:uneven to plane; 3:uneven to uneven
    bool isplane2=true;
    bool isProgress=true;  //是否显示计算进度，用-p参数控制
    int num_thread=omp_get_max_threads(); //default using all available threads
    bool debug_kernel_old=false;//默认不适用老的核矩阵
    string filename_exact="";//exact solution file name
    double DWC_parameter=-1;//CGLS收敛条件; //正则化参数; //积分迭代发的迭代次数
    int method_DWC=DWC_INTEGRALITERATION;//默认使用积分迭代发下延
    int extNumber=0;//扩边多少个点距
//********************************************
    //option
    for (int i=1; i<argc; i++)
    {
        string s(argv[i]);
        string option_str=s.substr(2,s.length());
        switch (argv[i][0]) {
            case '-':
                switch (argv[i][1])
            {
                    case 'E':
                        if(isNum(option_str))
                        {
                            extNumber=abs(int(atof(option_str.data())));
                            cout<<"extend data "<<extNumber<<" spacing"<<endl;
                        }else
                        {
                            cout<<RED<<"You use -E argument, but don't set extNum"<<COLOR_DEFALUT<<endl;
                        }
                        break;
                    case 'e':
                        filename_exact=option_str;
                        break;
                    case 'D':                   //downward continuation
                        isDownwardConti=true;
                        if(argv[i][2]=='+')
                        {
                            switch(argv[i][3])
                            {
                                case 'T':
                                {
                                    string lambda_str=s.substr(4,s.length());
                                    DWC_parameter=atof(lambda_str.data());
                                    method_DWC=DWC_TIKHONOV;
                                    cout<<"Tiknonv regularization parameter:"<<DWC_parameter<<endl;
                                }

                                    // cout<<rtp<<endl;
                                    break;
                                case 'I':
                                {
                                    string iter_str=s.substr(4,s.length());
                                    DWC_parameter=atof(iter_str.data());
                                    method_DWC=DWC_INTEGRALITERATION;
                                    cout<<"Integral iteration number:"<<DWC_parameter<<endl;
                                }

                                    break;
                                case 'C':
                                {
                                    string delta_str=s.substr(4,s.length());
                                    DWC_parameter=atof(delta_str.data());
                                    method_DWC=DWC_CGLS;
                                    cout<<"CGLS convegence condition:"<<DWC_parameter<<endl;
                                }
                                break;
                                case 'L':
                                {
                                    string iter_str=s.substr(4,s.length());
                                    DWC_parameter=atof(iter_str.data());
                                    method_DWC=DWC_LANDWEBER;
                                    cout<<"Landweber iteration number:"<<DWC_parameter<<endl;
                                }
                                    break;
                            }
                        }
                        break;
                    case 'T':                   //topography input data measured
                        if(isNum(option_str))
                        {
                            height1=atof(option_str.data());
                            isplane1=true;
                        }else
                        {
                            topo1=option_str;
                            isplane1=false;
                        }
                        break;
                    case 't':
                        if(isNum(option_str))
                        {
                            num_thread=atoi(option_str.data());
                            // cout<<"thread number: "<<num_thread<<endl;
                            if(num_thread>omp_get_max_threads())
                            {
                                OutputWarninginfo("the num_thread will be set to max_num_thread.\n");
                                num_thread=omp_get_max_threads();
                            }
                        }else
                        {
                            OutputErrorinfo("-t parmaeter must be a number\n");
                        }
                        break;
                    case 'H':
                        if(isNum(option_str))
                        {
                            height2=atof(option_str.data());
                            isplane2=true;
                        }else
                        {
                            topo2=option_str;
                            isplane2=false;
                        }
                        break;
                    case 'G':
                        outputfilename=option_str;
                        break;
                    case 'p':
                        isProgress=true;
                        break;
                    case 'o':
                        debug_kernel_old=true;
                        break;
                    case 'f':
                        isFrequency=true;
                        cout<<"Frequency domain"<<endl;
                        break;
                    default:
                        break;
                }
                break;
            default:
                inputfilename=argv[i];
                num_inputdata++;
                break;
        }
    }
    switch (num_inputdata)
    {
        case 0:
        {
            string errorinfo="you must input at least 1 grd file as inputdata file!";
            OutputErrorinfo(errorinfo);
            exit(0);
        }
            break;
        case 1:
            break;
        default:
        {
            string errorinfo="you input a wrong input data file!";
            OutputErrorinfo(errorinfo);
            exit(0);
        }
            break;
    }
    
    if (isDownwardConti) {
        if(isplane1==true && isplane2==true)
        {
            double start, end;
            start=omp_get_wtime();
            if(isFrequency)
            {
                DWC_p2p_f(inputfilename,outputfilename,height1,height2,extNumber,
                DWC_parameter,num_thread,isProgress,debug_kernel_old);
            }else
            {
                DWC_p2p(inputfilename,outputfilename,height1,height2,extNumber,
                DWC_parameter,method_DWC,num_thread,isProgress,debug_kernel_old,filename_exact);
            }
            end=omp_get_wtime();
            printf("Use Time:%f second\n",((double)(end - start)));
        }else if (isplane1==true && isplane2==false)
        {
            cout<<"DWC plane to surface: not implemented"<<endl;
            exit(0);
            // DWC_p2s(inputfilename,outputfilename,height1,topo2,num_thread,isProgress,debug_kernel_old,filename_exact);
        }else if (isplane1==false && isplane2==true)
        {
            double start, end;
            start=omp_get_wtime();
            cout<<"DWC surface to plane"<<endl;
            DWC_s2p(inputfilename,outputfilename,topo1,height2,extNumber,
            DWC_parameter,method_DWC,num_thread,isProgress,debug_kernel_old,filename_exact);
            end=omp_get_wtime();
            printf("Use Time:%f second\n",((double)(end - start)));
        }else if (isplane1==false && isplane2==false)
        {
            cout<<"uneven to uneven"<<endl;
        }
        
    }else   //upward conti
    {
        if(isplane1==true && isplane2==true)
        {
            double start, end;
            start=omp_get_wtime();
            if(isFrequency)
            {
                UWC_p2p_f(inputfilename, outputfilename, height1, height2,extNumber,
                filename_exact);
            }else
            {
                UWC_p2p(inputfilename, outputfilename, height1, height2,extNumber,
                num_thread,isProgress,debug_kernel_old,filename_exact);
            }
            end=omp_get_wtime();
            printf("Use Time:%f second\n",((double)(end - start)));
        }else if (isplane1==true && isplane2==false)
        {
            double start, end;
            start=omp_get_wtime();
            //cout<<"plane to uneven"<<endl;
            UWC_p2s(inputfilename,outputfilename,height1,topo2,extNumber,
            num_thread,isProgress,debug_kernel_old,filename_exact);
            end=omp_get_wtime();
            printf("Use Time:%f second\n",((double)(end - start)));
        }else if (isplane1==false && isplane2==true)
        {
            cout<<"UWC surface to plane: under construct"<<endl;
            // UWC_u2p(inputfilename,outputfilename,topo1,height2,num_thread,isProgress);
        }else if (isplane1==false && isplane2==false)
        {
            cout<<"UWC surface to surface: under construct"<<endl;
        }
    }
    return 0;
}
