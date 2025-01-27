#include "lab_vmk.h"

//#include <windows.h>
#include <cstdlib>
#include "SNDT.h"
#define MAX_GAMMA_SIZE pow(2, 32)

vector<int> func_f_non_corr_immun;

extern vector<int> init_begin_vec_rand(int max_value);
extern vector<int> func_f(int numb_bit, vector<int> s_vec);
extern vector<double> func_Furie(vector<int> vec);
extern pair<bool, int> Correlative_immunity (vector<double> Walsh);
extern int weight(vector<int> vec);
extern vector<int> int_to_vec_bit(int number, int bit_len);
extern void print(vector<int> vec, string str);                     //Lab1

//вес функции
int weight(bool* vec, unsigned long long int vec_size){ // не static для lab2.1 и lab4
    int counter=0;
    for (int i=0; i<vec_size; i++){
        if(vec[i]==1){
            counter++;
        }
    }
    return counter;
}

//подсчёт коэффициентов стат. структуры для функции
vector<int> gen_stat_struct_coeff(vector<int> vec){
    vector<double> Furie = func_Furie(vec);
    vector<double> Walsh_Hadamard(Furie.size());
    vector<int> Stat_Struct_Coeff(Furie.size());

    for(int j=0; j<Furie.size(); j++){
        Walsh_Hadamard[j] = j == 0 ? 1 - 2*Furie[j] : -2*Furie[j];
        Stat_Struct_Coeff[j] = pow(2, N_NUMBER-1) * Walsh_Hadamard[j];
    }
    return Stat_Struct_Coeff;
}


//функция создания регистров
static vector<vector<int>> Create_Reg(vector<int> register_lengths){
    vector<vector<int>> Registers(register_lengths.size());
    int buf;
    srand(time(nullptr));
    for(int i=0; i<register_lengths.size(); i++){
        uniform_int_distribution<> dist(1,pow(2, register_lengths[i])-1);
        do{
            buf = rand()%(int)(pow(2, register_lengths[i])-1);
        }while(buf < 1);
        Registers[i] = int_to_vec_bit(buf, register_lengths[i]);
    }
    return Registers;
}

//фукция создания ключа из регистров
static vector<int> Create_true_key(const vector<vector<int>> Registers, vector<int> register_lengths){
    vector<int> key;
    for(int i=0; i<register_lengths.size(); i++){
        for(int j=register_lengths[i]-1; j>=0; j--){
            key.push_back(Registers[i][j]);
        }
    }
    return key;
}

//функция генерации нового значения внутри регистра
static int func_reg(vector<int> vec, int reg_size){ 
    int count_begin_calc = vec.size() - reg_size;
    int result = 0;
    if(reg_size == 7 || reg_size == 13 || reg_size == 19){
        result = vec[count_begin_calc] ^ vec[count_begin_calc+1];
    }
    else if(reg_size == 5 || reg_size == 11 || reg_size == 17 || reg_size == 23 || reg_size == 29){
        result = vec[count_begin_calc] ^ vec[count_begin_calc+2];
    }
    else if(reg_size == 9 || reg_size == 16 || reg_size == 25){
        result = vec[count_begin_calc] ^ vec[count_begin_calc+1];
    }
    return result;
}

//функция генерации нового элемента гаммы
static int func_gamma_f(vector<int> x){
    int number = 0;
    for (int i = 0; i < x.size(); i++){
        number << 1;
        number += x[i]==0 ? 0 : 1;
    }
        
    return func_f_non_corr_immun[number];
}

//функция подсчёта НЕ совпадений в векторов-значений xi и f в таблице истинности
static vector<double> truth_table_f(int N){
    vector<double> xn(N, 0);
    int size = pow(2, N);
    for(int i=0; i<size; i++){
        vector<int> i_vec = int_to_vec_bit(i, N);
        int f = func_gamma_f(i_vec);
        for(int j=0; j<N; j++){
            if(i_vec[j] == f){
                xn[j]++;
            }
        }
    }
    for(int i=0; i<N; i++){
        xn[i] = 1 - xn[i]/pow(2, N);
        //printf("x%d = %lf\n", i+1, xn[i]);
    }
    return xn;
}

//генерируем не корреляционно имунную пор. 1 функцию длинной 2^(N+1) 128bit 
//с положительными коэффициентами стат. структуры на векторах веса 1 
//и с частотой НЕ совпадений в таблице истинности != 0.5
static vector<int> gen_func_f(){  
    printf("Подбираем функцию f удовлетворяющую условиям...\n");
    bool check_corr_immun;
    bool check_stat_struct_coeff;
    bool check_p;
    vector<vector<int>> result(N_NUMBER);
    do{
        check_corr_immun = false;
        check_stat_struct_coeff = false;
        check_p = false;
        vector<int> s_vec = init_begin_vec_rand(pow(2, N_NUMBER));
        for(int i=0; i<N_NUMBER && !check_corr_immun; i++){
            result[i] = func_f(i, s_vec);
            //print(result[i], "f"+to_string(i+1));

            vector<double> Furie = func_Furie(result[i]);
            vector<double> Walsh_Hadamard(Furie.size());
            vector<int> Stat_Struct_Coeff(Furie.size());

            for(int j=0; j<Furie.size(); j++){
                Walsh_Hadamard[j] = j == 0 ? 1 - 2*Furie[j] : -2*Furie[j];
                Stat_Struct_Coeff[j] = pow(2, N_NUMBER-1) * Walsh_Hadamard[j];
            }
            pair<bool, int> corr_immun = Correlative_immunity(Walsh_Hadamard);
            check_corr_immun = corr_immun.first;
            if(!check_corr_immun){
                check_stat_struct_coeff = true;
                for(int i=0; i<N_NUMBER; i++){
                    check_stat_struct_coeff = Stat_Struct_Coeff[0b100000 >> i] <= 0 ? false : check_stat_struct_coeff;
                }
                if(check_stat_struct_coeff){
                    func_f_non_corr_immun = result[i];
                    vector<double> p = truth_table_f(N_NUMBER);
                    check_p = true;
                    for (int i = 0; i < N_NUMBER; i++){
                        check_p = p[i] == 0.5 ? false : check_p;
                    }
                    if(check_p){
                        return result[i];
                    }
                }
            }
        }
    }while(check_corr_immun || !check_stat_struct_coeff || !check_p);
}

//функция, создающая гамму
vector<int> create_gamma(int gamma_size, int reg_size, const vector<vector<int>> Registers_const, vector<int> register_lengths){
    vector<vector<int>> Registers = Registers_const;
    vector<int> gamma(gamma_size);
    vector<int> current_vec_x(reg_size);
    for(int i=0; i<gamma_size; i++){
        for(int j=0; j<reg_size; j++){
            current_vec_x[j] = Registers[j][i];
            if(Registers[j].size() < gamma_size){
                vector<int> buf_reg(register_lengths[j]);
                copy_n(Registers[j].end() - register_lengths[j], register_lengths[j], buf_reg.begin()); //ОЧЕНЬ ВАЖНО! НЕ УДАЛЯТЬ! ИНАЧЕ БУДЕТ РАБОТАТЬ  МИЛЛИАРД ЛЕТ!
                //printf("test 4\n");
                Registers[j].push_back(func_reg(buf_reg, register_lengths[j]));
                //printf("test 5\n");
            }            
        }
        //printf("test 6\n");
        gamma[i] = func_gamma_f(current_vec_x);
        //printf("test 7\n");

        /*if((i+1)%10000 == 0){
            cout << "\r\tВычислено уже "<<i+1<<" из "<<gamma_size<<" бит";
        }*/
    }
    return gamma;
}

//функция выбора кандидата в регистры
vector<vector<int>> Search_part_key(double p, double q, int candidate_len, int number_candidate, const vector<vector<int>> Registers, vector<int> register_lengths){
    Standard_Normal_Distribution_Table SNDT = Standard_Normal_Distribution_Table();
    double alpha, beta, T, C;
    double Qantil_n_alpha, Qantil_n_beta;
    alpha = 0.0000; beta = 0.0000;
    cout << "\n\tПодибираем кандидатов в регистр " << number_candidate <<endl;

    do{
        alpha += 0.0000001;
        beta  += 0.0000001;
        Qantil_n_alpha = SNDT.find_x(1 - alpha);
        Qantil_n_beta  = SNDT.find_x(1 - beta);
        printf("\r\tFa(%.9Lf)=%Lf, Fb(%.9Lf)=%Lf", 1 - alpha, Qantil_n_alpha, 1 - beta, Qantil_n_beta);
        T = pow((Qantil_n_alpha*sqrt(p*(1-p)) + Qantil_n_beta*sqrt(q*(1-q))), 2) / pow(q-p, 2);
        printf("\talp=%.9Lf, bet=%.9Lf, X(1-a)=%Lf, X(1-b)=%Lf, p=%lf, q=%.1lf, T=%lf", alpha, beta, Qantil_n_alpha, Qantil_n_beta, p, q, T);
    }while(Qantil_n_alpha==ERROR_CODE || Qantil_n_beta==ERROR_CODE);
    cout << "\n\tT = " << T <<endl;

    double min_C = p < q ? T*p : T*q;
    double max_C = p > q ? T*p : T*q;
    C = Qantil_n_alpha * sqrt(T*p*(1-p)) + T*p;
    if(C < min_C || C > max_C){
        //C = (-Qantil_n_beta) * sqrt(T*p*(1-q)) + T*q;
    }
    cout << "\tC = " << C <<endl;

    vector<vector<int>> result;
    int max_candidate = pow(2, candidate_len);
    int positive_count = 0;

    vector<int> gamma = create_gamma(T, register_lengths.size(), Registers, register_lengths);
    int min_mismatches = (int)T;
    int numb_candidate_min_mismatches = 0;
    for(int i=1; i<max_candidate; i++){
        cout << "\r\t\tпрогресс: "<< i+1 <<"/"<<max_candidate;
        vector<int> buf_res = int_to_vec_bit(i, candidate_len);
        for(int j=candidate_len; j<T; j++){
            vector<int> buf_reg(candidate_len);
            copy_n(buf_res.end() - candidate_len, candidate_len, buf_reg.begin()); //ОЧЕНЬ ВАЖНО! НЕ УДАЛЯТЬ! ИНАЧЕ БУДЕТ РАБОТАТЬ  МИЛЛИАРД ЛЕТ!
            buf_res.push_back(func_reg(buf_reg, candidate_len));
        }
        //vector<int> buf_vec_mismatches(buf_res.size(), 0);
        int mismatches = 0;
        for(int j=0; j<T; j++){
            if(gamma[j] != buf_res[j]){
                mismatches++;
                if(mismatches >= C){
                    break;
                }
            }            
        }
        //int mismatches = weight(buf_vec_mismatches);
        if(mismatches < C){
            vector<int> buf_candidate(candidate_len);
            for(int j=0; j<candidate_len; j++){
                buf_candidate[j] = buf_res[j];
            }
            result.push_back(buf_candidate);
            if(mismatches < min_mismatches){
                numb_candidate_min_mismatches = positive_count;
                min_mismatches = mismatches;
            }
            positive_count++;
            cout << "\r\t\t\t\t\t\t\tположительное решение по: "<< positive_count <<" кандитатy(-ам)";
            //print(buf_candidate, "\n\t\t\t\t\t\t\tcandidate "+to_string(i));
        }
    }
    if(positive_count != 0){
        for (int i = 0; i < result.size(); i++){
            print(result[i], "\n\t\tcandidate "+to_string(i+1)+(i==numb_candidate_min_mismatches ? "(min mismatches)" : ""));
        }
    }
    positive_count == 0 ? cout << "\n\tНеудалось найти кандидатов!\nПерезапустите программу, пожалуйста!" <<endl : cout << "\n\tУспешно!" <<endl;
    
    return result;
}

//функция отсеивания ложных кандитатов путём непоредственной проверки
vector<vector<int>> Search_true_key(int count, vector<vector<int>> keys, vector<vector<vector<int>>> supposed_parts_key, vector<int> vec_numb, int reg_size, vector<int> register_lengths){
    pair<bool, vector<int>> result(false, {});
    if(count == reg_size){
        vector<vector<int>> current_key(reg_size);
        for(int i=0; i<reg_size; i++){
            current_key[i] = supposed_parts_key[i][vec_numb[i]];
        }
        
        keys.push_back(Create_true_key(current_key, register_lengths));
        return keys;
    }

    for (int i = 0; i < supposed_parts_key[count].size() && result.first == false; i++){
        vec_numb[count] = i;
        keys = Search_true_key(count+1, keys, supposed_parts_key, vec_numb, reg_size, register_lengths);
    }
    return keys;
}

//main-функция лабы
void Lab4(/*vector<int> func_f, vector<int> Stat_Struct_Coeff*/){
    //func_f_non_corr_immun = func_f;
    //func_f_non_corr_immun = {0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0}; //функция Миши
    
    func_f_non_corr_immun = gen_func_f();
    //func_f_non_corr_immun = {0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1}; //случайная функция из gen_func_f()
    //func_f_non_corr_immun = {0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 1, 1, 1}; //случайная функция из gen_func_f()
    //func_f_non_corr_immun = {0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1}; //случайная функция из gen_func_f()
    //func_f_non_corr_immun = {0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1}; //случайная функция из gen_func_f()
    //func_f_non_corr_immun = {0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1}; //случайная функция из gen_func_f()
    //func_f_non_corr_immun = {0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 1, 0, 1, 1}; //случайная функция из gen_func_f()
    vector<int> Stat_Struct_Coeff = gen_stat_struct_coeff(func_f_non_corr_immun);
    print(func_f_non_corr_immun, "Функция генерации гаммы");

    vector<int> register_lengths = {7, 11, 13, 17, 19, 23};
    int reg_size = register_lengths.size();
    if(reg_size != N_NUMBER) {
        cout << "Ошибка в размерности. Корректный запуск программы невозможен" << endl;
        return;
    }

    //vector<vector<int>> Registers (reg_size);
    vector<int> key;

    //Registers = {{1,1,1}, {1,1}, {1,1}};
    printf("\n");
    const vector<vector<int>> Registers = Create_Reg(register_lengths);
    for (int i = 0; i < reg_size; i++){
        print(Registers[i], "Регистр "+to_string(i+1));
    }
    
    key = Create_true_key(Registers, register_lengths);
    print(key, "\nИстинный ключ");
    printf("\n");
    
    vector<double> P_f_no_equals_xn(reg_size);     //колличество(частота) НЕ совпадений xn и f = p = 0.5
    P_f_no_equals_xn = truth_table_f(reg_size);

    vector<double> Q_f_no_equals_y(reg_size);      //колличество(частота) совпадений y и f = q1
    for(int i=0; i<reg_size; i++){
        if(Stat_Struct_Coeff[0b100000 >> i] > 0){
            Q_f_no_equals_y[i] = 0.5 + (double)abs(Stat_Struct_Coeff[(0b100000 >> i)])/pow(2, reg_size);
        }
        else{
            Q_f_no_equals_y[i] = 1 - P_f_no_equals_xn[i];
            /*printf("delta(x%d) = 0. Извините, мы не можем обработать этот случай.\n", i+1);
            return;*/
        }
        printf("q%d = %lf\tp%d = %lf\n", i+1, Q_f_no_equals_y[i], i+1, P_f_no_equals_xn[i]);
    }
    
    
    vector<vector<vector<int>>> supposed_parts_key(reg_size);
    for(int i=0; i<reg_size; i++){
        supposed_parts_key[i] = Search_part_key(1-Q_f_no_equals_y[i], 0.5, register_lengths[i], i+1, Registers, register_lengths);
        if(supposed_parts_key.empty()){
            cout << "\nВеликий бог рандома не даровал Вам удачу в этот раз... Чтож, попробуйте снова!" <<endl;
            return;
        }
    }

    vector<vector<int>> keys;
    keys = Search_true_key(0, keys, supposed_parts_key, vector<int>(reg_size), reg_size, register_lengths);

    for(int i=0; i<keys.size(); i++){
        print(keys[i], "Выбранный ключ");
        vector<int> check_true_key(key.size());
        for(int j=0; j<key.size(); j++){
            check_true_key[j] = key[j] ^ keys[i][j];
        }
        cout << "Алгоритм ошибся на " << weight(check_true_key) << " бит" << endl;
    }
    
    return;
}