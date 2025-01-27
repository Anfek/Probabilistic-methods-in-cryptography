#include "lab_vmk.h"


static bool sravn(int buf, vector<int> s_vec){
    for(int i=0; i<s_vec.size(); i++){
        if(buf == s_vec[i]){
            return false;
        }
    }
    return true;
}

vector<int> init_begin_vec_rand(int two_pow_n){ //не static для lab4
    vector<int> s_vec; 
    srand(time(nullptr));
    for(int i=0; i<two_pow_n; i++){
        int buf = rand()%two_pow_n;
        //printf("buf = %d\n", buf);
        while (!sravn(buf, s_vec)){
            buf = rand()%two_pow_n;
            //printf("buf = %d\n", buf);
        }
        s_vec.push_back(buf);
    }

    return s_vec;
}

void print(vector<int> vec, string str){ //не static для lab4
    cout << str <<": ";
    for(int i=0; i<vec.size(); i++){
        cout << vec[i] << " ";
    }
    cout << endl;
}

vector<int> func_f(int numb_bit, vector<int> s_vec){//не static для lab4
    vector<int> f;
    int bit = 1 << numb_bit;
    for(int i=0; i<s_vec.size(); i++){
        int buf = bit & s_vec[i];
        buf >>= numb_bit;
        f.push_back(buf);
    }
    return f;
}

int weight(vector<int> vec){ // не static для lab2.1 и lab4
    int counter=0;
    for (int i=0; i<vec.size(); i++){
        if(vec[i]==1){
            counter++;
        }
    }
    return counter;
}

static vector<int> Jegalkin(vector<int> vec){
    const int vec_size = vec.size();
    int counter = vec.size()/2;
    vector<int> buf(vec.size());
    //print(vec, "step 0");
    for(int i=0; i<N_NUMBER; i++){
        for(int j=0; j<vec_size; j+=2*counter){
            for(int k=j; k<j+counter; k++){
                buf[k] = vec[k];  
                buf[k+counter] = vec[k] ^ vec[k+counter];
            }
        }
        vec = buf;
        //print(vec, "step "+to_string(i+1));
        counter/=2;
    }
    return vec;
}

static pair<string, vector<bool>> Jegalkin_poly(vector<int> vec){
    vector<bool> fict_x(N_NUMBER, true);
    bool string_emty = true;
    string poly = "";
    for(int i=0; i<vec.size(); i++){
        if(vec[i] == 1){
            if(i == 0){
                poly += "1";
                string_emty = false;
            }
            else{
                if(string_emty) string_emty = false; else poly += " + ";
                int bit = 1;
                bool mult_emty = true;
                for(int j=0; j<N_NUMBER; j++){
                    int buf_check_j_bit = bit & i;
                    if(buf_check_j_bit != 0){
                        if(mult_emty) mult_emty= false; else poly += "*";
                        poly += "x" + to_string(j+1);
                        fict_x[j] = false;
                    }
                    bit <<= 1;
                }
                
            }
        }
    }
    return make_pair(poly, fict_x);
}

vector<vector<int>>  Lab1(){
    printf("\t\t1. Генерация случайной подстановки V%d->V%d\n\n", N_NUMBER, N_NUMBER);
    int two_pow_n = pow(2, N_NUMBER);
    vector<int> s_vec;
    s_vec = init_begin_vec_rand(two_pow_n);
    
    print(s_vec, "S");
    cout << endl;
    
    cout << "\t\t2. Вектора значений координатных функций\n" << endl;
    vector<vector<int>> f_vec(N_NUMBER, vector<int>(two_pow_n));
    for(int i=0; i<N_NUMBER; i++){
        f_vec[N_NUMBER - 1 - i] = func_f(i, s_vec);
    }
    
    for(int i=0; i<N_NUMBER; i++){
        print(f_vec[i], "f"+to_string(i+1));
    }
    cout << endl;
    
    cout << "\t\t3. Вес координатных функций\n" << endl;
    for(int i=0; i<N_NUMBER; i++){
        cout << "||f"+to_string(i+1)+"|| = "<< weight(f_vec[i])<<endl;
    }
    cout << endl;
    
    
    //vector<int> vec_test1 = {1,0,1,1,0,0,1,0};
    //vector<int> vec_test2 = {1,1,0,0,0,0,1,1};
    vector<vector<bool>> f_fict_x(N_NUMBER);
    cout << "\t\t4. Многочлен Жегалкина для координатных функций\n" << endl;
    for(int i=0; i<N_NUMBER; i++){
        vector<int> f_vec_J = Jegalkin(f_vec[i]);
        pair<string, vector<bool>> buf_pair = Jegalkin_poly(f_vec_J);
        
        string f_vec_J_poly = buf_pair.first;
        f_fict_x[i] = buf_pair.second;
        
        print(f_vec_J, "\nf"+to_string(i+1));
        cout << "f"+to_string(i+1)+": " << f_vec_J_poly << endl;
    }
    
    cout << "\n\t\t5. Фиктивные переменные для координатных функций\n" << endl;
    for(int i=0; i<N_NUMBER; i++){
        bool no_fict = true;
        cout << "f"+to_string(i+1)+": ";
        for(int j=0; j<N_NUMBER; j++){
            if(f_fict_x[i][j] == true){
                no_fict ? cout << "" : cout << ", ";
                cout << "x"+to_string(j+1);
                no_fict = false;
            }
        }
        no_fict ? cout <<"нет фиктивных переменных\n" : cout <<"\n";
    }
    
    return f_vec;
}