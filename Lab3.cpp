#include "lab_vmk.h"

extern double predominance_of_zeros_over_ones(vector<int> vec);
extern int weight(vector<int> vec);
extern string print_last_bits(unsigned int number, int n);

void print(vector<double> vec, string str){
    cout << str <<": ";
    for(int i=0; i<vec.size(); i++){
        cout << vec[i] << " ";
    }
    cout << endl;
}

void print_i(vector<int> vec, string str){
    cout << str <<": ";
    for(int i=0; i<vec.size(); i++){
        cout << vec[i] << " ";
    }
    cout << endl;
}

vector<double> func_Furie(vector<int> vec){     //не static для lab4
    vector<double> furie(vec.size()); // Создаем вектор double нужного размера
    for (int i = 0; i < vec.size(); ++i) {
        furie[i] = static_cast<double>(vec[i]); // Преобразуем каждый элемент int в double
    }
    int counter = 1;

    //printf("FURIE\n");
    //print(furie, "vec");
    for(int i=0; i<N_NUMBER; i++){
        vector<double> buf(furie.size());
        for(int j=0; j<vec.size(); j+=2*counter){
            for(int k=j; k<j+counter; k++){
                buf[k] = (furie[k]+furie[k+counter])/2;  
                buf[k+counter] = (furie[k]-furie[k+counter])/2;
            }
        }
        counter*=2;
        furie = buf;
        //print(furie, "furie "+to_string(i));
    }
    return furie;
}

vector<int> int_to_vec_bit(int number, int bit_len){  //не static для lab4
    vector<int> bits(bit_len);
    
    for(int i=0; i<bit_len; ++i){
        bool bit = (number & (1 << i)) > 0;

        bits[bit_len-i-1] = bit ? 1 : 0;
    }
    return bits;
}

pair<bool, int> Correlative_immunity (vector<double> Walsh){ //не static для lab4
    pair<bool, int> result = make_pair(false, 0);

    for(int i=1; i <= N_NUMBER; i++){
        bool check_zero_Walsh = true;
        for(int j=0; j<Walsh.size(); j++){
            if(i == weight(int_to_vec_bit(j, N_NUMBER))){ //если текущий вектор имеет нужный нам вес
                check_zero_Walsh = Walsh[j] == 0.0 ? check_zero_Walsh : false; //проверяем нулевое ли у него W
            }
            if(check_zero_Walsh == false){
                //cout << "W("+print_last_bits(j, N_NUMBER)+") = "<< Walsh[j] << " -> ";
                break;
            }
        }
        
        if(check_zero_Walsh) result = make_pair(check_zero_Walsh, i);
        else break;
    }

    return result;
}

static pair<bool, int> Elasticity (vector<int> vec, pair<bool, int> corr_immun){
    return (corr_immun.first == true && predominance_of_zeros_over_ones(vec) != 0) ? make_pair(false, corr_immun.second) : corr_immun;    
}

static void Spectrum_BF(vector<vector<int>>delta){
    cout<< "x1...x" << N_NUMBER << "\t";
    for(int i=0; i<N_NUMBER; i++){
        cout<< "del(f" << i+1 <<")\t";
    }
    cout<< "\n";

    for(int i=0; i<delta[0].size(); i++){
        cout << print_last_bits(i, N_NUMBER) <<"\t";
        for(int j=0; j<N_NUMBER; j++){
            cout << delta[j][i] <<"\t";
        }
        cout<< "\n";
    }
    cout<< "\n";
}

static string str_vec_poly(vector<int> vec, bool one){
    bool string_emty = true;
    string poly = "";
    for(int i=0; i<vec.size(); i++){
        if(vec[i] == 1){
            if(string_emty) string_emty = false; else poly += " + ";
            poly += "x" + to_string(i+1);
        }
    }

    if(one){
        if(!string_emty) poly += " + ";
        poly += "1";
    }    
    return poly;
}

static vector<string> best_linear_approximation(vector<int> delta){
    int max = *max_element(begin(delta), end(delta));
    int min = *min_element(begin(delta), end(delta));
    max = abs(min) > max ? abs(min) : max;
    if(max == 0) return {""};

    vector<string> result;
    for(int i=0; i<delta.size(); i++){
        if(abs(delta[i]) == max){
            string buf = str_vec_poly(int_to_vec_bit(i, N_NUMBER), delta[i]<0);
            result.push_back(buf);
        }
    }
    return result;
}

static bool check_bent(vector<double> Walsh){
    if(N_NUMBER%2 != 0) return false;
    bool result = true;
    double buf = abs(Walsh[0]);
    for(int i=1; i<Walsh.size(); i++){
        result = buf == abs(Walsh[i]) ? result : false;
        if(!result) return result;
    }
    return result;
}

pair<int, vector<int>> Lab3 (vector<vector<int>> f_vec){
    //f_vec[0] = {0,0,0,1,1,1,1,0};
    vector<vector<double>> Furie(N_NUMBER);
    vector<vector<double>> Walsh_Hadamard(N_NUMBER);
    vector<vector<int>> Stat_Struct_Coeff(N_NUMBER);
    int nunber_non_corr_immun_func = 0;

    for(int i=0; i<N_NUMBER; i++){
        Furie[i] = func_Furie(f_vec[i]);
        Walsh_Hadamard[i].resize(Furie[i].size());
        Stat_Struct_Coeff[i].resize(Furie[i].size());

        for(int j=0; j<Furie[i].size(); j++){
            Walsh_Hadamard[i][j] = j == 0 ? 1 - 2*Furie[i][j] : -2*Furie[i][j];
            Stat_Struct_Coeff[i][j] = pow(2, N_NUMBER-1) * Walsh_Hadamard[i][j];
        }
        //print(Furie[i], "F"+to_string(i));
        //print(Walsh_Hadamard[i], "W"+to_string(i));
        //print_i(Stat_Struct_Coeff[i], "delta"+to_string(i));
    }

    cout << "\n\t\t1. Корреляционная иммунность и эластичность для координатных функций" << endl;
    
    for(int i=0; i<N_NUMBER; i++){
        cout << "\n\nf"+to_string(i+1)+": \t";
        pair<bool, int> corr_immun = Correlative_immunity(Walsh_Hadamard[i]);
        pair<bool, int> elastic = Elasticity(f_vec[i], corr_immun);

        corr_immun.first ? cout << "функция корреляционно иммунна порядка "+to_string(corr_immun.second) : cout << "функция НЕ корреляционно иммунна порядка 1";
        elastic.first ? cout << "\n\tфункция эластична порядка "+to_string(corr_immun.second) : cout << "\n\tфункция НЕ эластична";

        if(corr_immun.first == false){
            nunber_non_corr_immun_func = i;
        }
    }
    
    cout << "\n\n\n\t\t2. Построение спектра булевых функций для координатных функций\n\n" << endl;
    Spectrum_BF(Stat_Struct_Coeff);

    cout << "\n\n\t\t3. Нахождение наилучшего линейного приближения для координатных булеых функций\n\n" << endl;
    for(int i=0; i<N_NUMBER; i++){
        vector<string> lin = best_linear_approximation(Stat_Struct_Coeff[i]);
        cout << "\n\nf"+to_string(i+1)+":";
        for(int j=0; j<lin.size(); j++){
            cout << "\tg_"+to_string(j+1)+"(x) = " << lin[j] <<endl;
        }
    }
    
    cout << "\n\n\t\t4. Является ли координатная функция бент-функцией\n\n" << endl;
    for(int i=0; i<N_NUMBER; i++){
        cout << "f"+to_string(i+1)+":\tданная координатная функция ";
        if(!check_bent(Walsh_Hadamard[i])) cout << "НЕ";
        cout << " является бент-функцией" << endl;
    }

    return make_pair(nunber_non_corr_immun_func, Stat_Struct_Coeff[nunber_non_corr_immun_func]);
}