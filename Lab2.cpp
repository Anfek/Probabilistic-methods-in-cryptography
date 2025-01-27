#include "lab_vmk.h"

#define BAN_DIMENSION_LIMITATION 20

extern int weight(const vector<int> vec);

double predominance_of_zeros_over_ones(const vector<int> vec){ //не static для lab3.1
    double weight_f = weight(vec);
    double rezult = 1.0 - weight_f/pow(2, N_NUMBER-1);
    return rezult;
}

string print_last_bits(unsigned int number, int n) { //не static для lab3.1
    // Проверка корректности входных данных
    if (n <= 0) {
        return "";
    }

    string bits = "";

    // Поочередно проверяем каждый бит, начиная с младшего
    for (int i = 0; i < n; ++i) {
        // Используем битовую операцию AND, чтобы получить значение i-го бита
        bool bit = (number & (1 << i)) > 0;

        // Преобразуем бит в строковый символ ('0' или '1')
        bits += bit ? '1' : '0';
    }
    reverse(bits.begin(), bits.end()); 
    // Возвращаем строку с последними n битами
    return bits;
}


static pair<bool, int> ban(vector<int> vec, const int number_current_f){
    //pair<bool, int> result;
    vector<int> statistic(pow(2, number_current_f), 0);
    for(int i=0; i<vec.size(); i++){    //частота встречаемости каждого значения
        statistic[vec[i]]++;
    }
    //int min = statistic[i];
    for(int i=0; i<statistic.size(); i++){
        if(statistic[i]==0){
            //result = make_pair(true, i);   //false, значение f
            //min = statistic[i];
            return make_pair(true, i);
        }
    }
    return make_pair(false, 0);
    /*
    if(min == 0){
        return make_pair{true, 0};
    }
    else{
        for(int i=0; i<vec.size(); i++){
            if(vec[i] == result.second){
                result = make_pair{false, i};   //false, номер значения f = x1x2...xn
                break;
            }
        }
    }
    return result;
    */
}

static bool highly_equiprobable(vector<int> vec){
    int number_current_f = 0;   //2^(NCF+1) = диапазон значений f
    vector<vector<int>> all_f;
    all_f.push_back(vec);
    
    pair<bool, int> check_ban = ban(all_f[number_current_f], number_current_f+1); //получили есть ли значение f с встречаемостью 0 и если да - значение f
    while(number_current_f <= BAN_DIMENSION_LIMITATION && check_ban.first == false){
        vector<int> buf_new_f(2*all_f[number_current_f].size());
        for(int i=0; i<all_f[number_current_f].size(); i++){
            buf_new_f[2*i]   = (all_f[number_current_f][i]<<1) + vec[(2*i) & ((int)pow(2, N_NUMBER)-1)]; //предыдущее значение f || значение f от последних n переменных(бит)
            buf_new_f[2*i+1] = (all_f[number_current_f][i]<<1) + vec[(2*i+1) & ((int)pow(2, N_NUMBER)-1)];
            //cout << print_last_bits(i, N_NUMBER+number_current_f)       << ": " << print_last_bits(all_f[number_current_f][i], number_current_f + 1) << " \t "; 
            //cout << print_last_bits(2*i, N_NUMBER+number_current_f+1)   << ": " << print_last_bits(buf_new_f[2*i], number_current_f+2) << " \t "; 
            //cout << print_last_bits(2*i+1, N_NUMBER+number_current_f+1) << ": " << print_last_bits(buf_new_f[2*i+1], number_current_f+2) << endl;
        }
        all_f.push_back(buf_new_f);
        number_current_f++;
        check_ban = ban(all_f[number_current_f], number_current_f+1);
    }
    //check_ban = ban(all_f[number_current_f], number_current_f+1);
    /*for (int i=0; i<all_f[number_current_f].size(); i++){
        cout << print_last_bits(i, N_NUMBER+number_current_f) << ": " << print_last_bits(all_f[number_current_f][i], number_current_f+1) << endl; 
    }*/
    check_ban.first ? cout << "запрет = " << print_last_bits(check_ban.second, number_current_f+1)<< endl : cout << "запрет не найден " << endl;
    return !check_ban.first;
}

void Lab2(const vector<vector<int>> f_vec){
    cout << "\n\t\t1. Преобладание нулей над единицами для координатных функций\n" << endl;
    for(int i=0; i<N_NUMBER; i++){
        //cout << "f"+to_string(i+1)+": " << (double)predominance_of_zeros_over_ones(f_vec[i]) << endl;
        printf("f%d: %lf\n", i+1, predominance_of_zeros_over_ones(f_vec[i]));
    }

    cout << "\n\t\t2. Сильно равновероятность для координатных функций" << endl;
    cout << "\t\t\t\tИ" << endl;
    cout << "\t\t3. Запреты для координатных функций\n\n" << endl;
    //f_vec[0] = {0,0,0,0,0,0,0,1};
    //f_vec[1] = {0,0,1,1,0,1,1,0};
    for(int i=0; i<N_NUMBER; i++){
        cout << "f"+to_string(i+1)+": \t";
        highly_equiprobable(f_vec[i]) ? cout << "\tфункция СИЛЬНО РАВНОВЕРОЯТНАЯ\n" << endl : cout << "\tфункция НЕ сильно равновероятная\n" << endl;
    }
}