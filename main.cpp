//g++ main.cpp Lab1.cpp Lab2.cpp Lab3.cpp Lab4.cpp -o labs_vmk
//windows:  labs_vmk.exe
//Linux:    ./labs_vmk

#include "lab_vmk.h"

extern vector<vector<int>> Lab1();
extern void Lab2(vector<vector<int>> f_vec);
extern pair<int, vector<int>> Lab3(vector<vector<int>> f_vec);
extern void Lab4(/*vector<int> func_f, vector<int> Stat_Struct_Coeff*/);

int main(){
    system("chcp 65001 >nul 2>&1");

    cout << "\n\n\tЛабораторная работа №1\n\n" << endl;
    vector<vector<int>> f_vec = Lab1();

    cout << "\n\n\tЛабораторная работа №2\n\n" << endl;
    Lab2(f_vec);

    cout << "\n\n\tЛабораторная работа №3\n\n" << endl;
    pair<int, vector<int>> non_corr_immun_func = Lab3(f_vec);

    cout << "\n\n\tЛабораторная работа №4\n\n" << endl;
    /*f_vec[non_corr_immun_func.first]={0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 
                                        0, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 1, 0, 
                                        0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 0, 
                                        0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0};*/
    Lab4(/*f_vec[non_corr_immun_func.first], non_corr_immun_func.second*/);
}