#include <iostream>
#include <vector>
#include <algorithm>    // std::lower_bound, std::upper_bound, std::sort
#include <functional>   // std::greater
#include <utility>
#include <cmath>       // exp
#include <random>
#include <cfloat>
#include <sstream>
#include <fstream>
#include <regex>
#include <chrono>
#include <thread>


using namespace std;

#define MAX_POP_CAP 3000 // must be even to make things easier on procreating phase
#define POP_LIM 2000
#define IND_PARAMS 50
#define PARENTS_NUM 1000
#define CROSSOVER_PROB 0.8
#define MUT_PROB 0.04
#define MAX_GENE_RANGE 255 // gene values can go from 0 to 255
#define MAX_LOOP 10 // max number of loops to generate expression
#define GEN_LIM 300
#define TRAINING_RATIO 0.8
#define PRESET_RATIO 0.30
#define INIT_MUT 10
#define PRESET1 {232,36,248,112,227,105,23,203,98,92,161,243,231,156,54,28,167,178,114,169,36,16,158,218,47,129,26,238,118,23,231,229,59,198,169,5,44,207,103,139,31,178,25,252,49,169,7,227,179,207,61}
//10_0.1 0.11111 0.115986 -(*(*(8.2)(e))(*(849)(h)))(/(1998.9)(p(*(*(29.9)(k))(*(2.6)(x)))))

#define PRESET2 {196,66,72,122,204,153,190,223,147,233,248,148,196,55,59,99,178,232,149,123,181,247,231,241,26,113,102,171,64,86,137,163,223,120,207,35,96,59,21,78,118,104,135,76,30,48,91,96,212,180}
//40_0.35 0.108539 0.10545 -(*(/(8.8)(+(*(31)(k))(*(2)(g))))(-(p(*(5)(x)))(*(608)(g))))(*(-(~03)(*(998)(x)))(p(*(1.6)(e))))

#define PRESET3 {168,46,252,102,185,135,223,194,31,112,81,43,231,205,142,230,145,199,199,119,47,100,60,243,162,245,198,202,148,222,13,31,179,229,38,237,20,77,207,25,34,43,227,39,249,123,171,5,135,193,61}
//10_0.2 0.100633 0.105957 -(*(p(*(2)(e)))(*(999.7)(h)))(/(1998)(p(*(*(7)(k))(p(*(3.1)(x))))))


double PURGED_PARENTS_RATIO  = 0.1;
int CARDS_CONST = 10;

auto time_begin = std::chrono::high_resolution_clock::now();
auto time_end = std::chrono::high_resolution_clock::now();


vector<vector<int>> population; // number of individuals in a generation and each one is a IND_PARAMS-vector
vector<vector<int>> parents; // individuals who will procreate

random_device rd;  //Will be used to obtain a seed for the random number engine
mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()

ofstream output;

vector<vector<double>> parsedCsv;
string grammar(vector<int> ind);
double evaluate_expression(string expression, vector<double>& values);

double sigmoid(double x){
    return 1/(1+exp(-x));
}

void calculate_cost(string& expr, double &cost, int lower_limit, int upper_limit){
    for(int i = lower_limit; i < upper_limit; i++){
        double v = (parsedCsv[i][8] - sigmoid(0.001*evaluate_expression(expr, parsedCsv[i])));
        cost += v*v;
        if(isnan(cost)) break;
    }
}

long double performance_function(int index){
    double cost = 0;
    string expr = grammar(population[index]);
    int count = 0;
    if (expr.find("x") != string::npos) count++;
    if (expr.find("e") != string::npos) count++;
    if (expr.find("f") != string::npos) count++;
    if (expr.find("g") != string::npos) count++;
    if (expr.find("h") != string::npos) count++;
    if (expr.find("i") != string::npos) count++;
    if (expr.find("j") != string::npos) count++;
    if (expr.find("k") != string::npos) count++;
    if (count < 3) return 0;
    //auto t1 = std::chrono::high_resolution_clock::now();

    double cost1, cost2, cost3, cost4;
    cost1 = cost2 = cost3 = cost4 = 0;
    thread thr1(calculate_cost, ref(expr), ref(cost1), 0, (TRAINING_RATIO*parsedCsv.size())/4);
    thread thr2(calculate_cost, ref(expr), ref(cost2), (TRAINING_RATIO*parsedCsv.size())/4, (TRAINING_RATIO*parsedCsv.size())/2);
    thread thr3(calculate_cost, ref(expr), ref(cost3), (TRAINING_RATIO*parsedCsv.size())/2, (3*(TRAINING_RATIO*parsedCsv.size()))/4);
    thread thr4(calculate_cost, ref(expr), ref(cost4), (3*(TRAINING_RATIO*parsedCsv.size()))/4, (TRAINING_RATIO*parsedCsv.size()));
    thr1.join();
    thr2.join();
    thr3.join();
    thr4.join();
    cost = cost1 + cost2 + cost3 + cost4;
    //auto t2 = std::chrono::high_resolution_clock::now();
    //auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
    //cout << duration << endl;

    cost = sqrt(cost/(TRAINING_RATIO*parsedCsv.size()-1));
    //cout << index << ": " << exp << ": " << cost << endl;
    if(cost == 0)
        return LDBL_MAX;
    if(isnan(cost)) return 0;
    return  exp(CARDS_CONST/cost) - 1;
}

long double performance_function2(int index){
    double cost = 0;
    string expr = grammar(population[index]);
    int count = 0;
    if (expr.find("x") != string::npos) count++;
    if (expr.find("e") != string::npos) count++;
    if (expr.find("f") != string::npos) count++;
    if (expr.find("g") != string::npos) count++;
    if (expr.find("h") != string::npos) count++;
    if (expr.find("i") != string::npos) count++;
    if (expr.find("j") != string::npos) count++;
    if (expr.find("k") != string::npos) count++;
    if (count < 3) return 0;
    //auto t1 = std::chrono::high_resolution_clock::now();

    double cost1, cost2, cost3, cost4;
    cost1 = cost2 = cost3 = cost4 = 0;
    double validation_begin =(TRAINING_RATIO*parsedCsv.size()), validation_size = (1-TRAINING_RATIO)*parsedCsv.size();
    thread thr1(calculate_cost, ref(expr), ref(cost1), validation_begin, validation_begin + validation_size/4);
    thread thr2(calculate_cost, ref(expr), ref(cost2), validation_begin + validation_size/4, validation_begin + validation_size/2);
    thread thr3(calculate_cost, ref(expr), ref(cost3), validation_begin + validation_size/2, validation_begin + (3*validation_size)/4);
    thread thr4(calculate_cost, ref(expr), ref(cost4), validation_begin + (3*validation_size)/4, parsedCsv.size());
    thr1.join();
    thr2.join();
    thr3.join();
    thr4.join();
    cost = cost1 + cost2 + cost3 + cost4;
    //auto t2 = std::chrono::high_resolution_clock::now();
    //auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
    //cout << duration << endl;

    cost = sqrt(cost/(validation_size-1));
    //cout << index << ": " << exp << ": " << cost << endl;
    if(cost == 0)
        return LDBL_MAX;
    if(isnan(cost)) return 0;
    return  exp(CARDS_CONST/cost) - 1;
}

/*double cost_function(int index){
    double cost = 0;
    string exp = grammar(population[index]);
    //auto t1 = std::chrono::high_resolution_clock::now();
    for(int i = 0; i < RAND_SAMPLE; i++){
        double v = ((parsedCsv[i][8]) - (evaluate_expression(exp, parsedCsv[i])));
        cost += v*v;
        if(isnan(cost)) break;
    }
    //auto t2 = std::chrono::high_resolution_clock::now();
    //auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
    //cout << duration << endl;
    cost = sqrt(cost/(RAND_SAMPLE-1));
    //cout << index << ": " << exp << ": " << cost << endl;
    if(isnan(cost)) return DBL_MAX;
    return cost;
}
*/



void parents_mutate(){
    uniform_real_distribution<> r_dis(0.0, 1.0);
    uniform_int_distribution<> i_dis(0, MAX_GENE_RANGE);
    for (int i = PARENTS_NUM*(1-PURGED_PARENTS_RATIO); i < PARENTS_NUM; i++){
        for (int j = 0; j < IND_PARAMS; j++){
            // Use dis to transform the random unsigned int generated by gen into a
            // double in [0, 255]. Each call to i_dis(gen) generates a new random int
            if (r_dis(gen) < MUT_PROB)
                parents[i][j] = i_dis(gen);
        }
    }
}


void select_parents(){
    vector<long double> select_lim; // actual probabilities: [0.1 0.3 0.2 0.4] => [0.1 0.4 0.6 1] (sum of probabilities)
    select_lim.resize(POP_LIM);
    long double total_performance = 0;
    for (int i = 0; i < POP_LIM; i++){
        total_performance += performance_function(i);
        select_lim[i] = (total_performance);
    }
    for (int i = 0; i < POP_LIM; i++)
        select_lim[i] /= (total_performance+0.00000001);

    uniform_real_distribution<> r_dis(0.0, 1.0);
    double random_real;
    vector<long double>::iterator selected;
    for (int i = 0; i < PARENTS_NUM; i++){
        random_real = r_dis(gen);
        selected = lower_bound(select_lim.begin(), select_lim.end(), random_real); // finds first element greater or equal than random_real
        parents[i] = (population[selected-select_lim.begin()]);
    }
    uniform_int_distribution<> i_dis(0, MAX_GENE_RANGE);

    int limit0 = PARENTS_NUM*(1-PURGED_PARENTS_RATIO);
    int limit1 = PARENTS_NUM*(1-PURGED_PARENTS_RATIO) + PARENTS_NUM*(PURGED_PARENTS_RATIO)*(1-PRESET_RATIO);
    int limit2 = PARENTS_NUM*(1-PURGED_PARENTS_RATIO) + PARENTS_NUM*(PURGED_PARENTS_RATIO)*(1-PRESET_RATIO) + (1*PARENTS_NUM*(PURGED_PARENTS_RATIO)*(PRESET_RATIO)/3);
    int limit3 = PARENTS_NUM*(1-PURGED_PARENTS_RATIO) + PARENTS_NUM*(PURGED_PARENTS_RATIO)*(1-PRESET_RATIO) + (2*PARENTS_NUM*(PURGED_PARENTS_RATIO)*(PRESET_RATIO)/3);
    int limit4 = PARENTS_NUM*(1-PURGED_PARENTS_RATIO) + PARENTS_NUM*(PURGED_PARENTS_RATIO)*(1-PRESET_RATIO) + (3*PARENTS_NUM*(PURGED_PARENTS_RATIO)*(PRESET_RATIO)/3);
    for (int i = limit0; i < limit1; i++)
            for (int j = 0; j < IND_PARAMS; j++)
                parents[i][j] = i_dis(gen); // random initialization of the individuals

    for (int i = limit1; i < limit2; i++)
        parents[i] = PRESET1; // preset 1

    for (int i = limit2; i < limit3; i++)
        parents[i] = PRESET2; // preset 2

    for (int i = limit3; i < limit4; i++)
        parents[i] = PRESET3; // preset 3

    for(int i = 0; i < INIT_MUT; i++)   parents_mutate();


}




/*
void select_parents(){
    vector<double> cost_table(POP_LIM, -1);
    random_unique(parsedCsv.begin(), parsedCsv.end(), RAND_SAMPLE);

    uniform_int_distribution<> i_dis(0, POP_LIM-1);
    for (int i = 0; i < PARENTS_NUM; i++){
        int rand_int_1 = i_dis(gen);
        int rand_int_2 = i_dis(gen);
        if(cost_table[rand_int_1] == -1){
            cost_table[rand_int_1] = cost_function(rand_int_1);
        }
        if(cost_table[rand_int_2] == -1){
            cost_table[rand_int_2] = cost_function(rand_int_2);
        }
        if(cost_table[rand_int_1] < cost_table[rand_int_2]){
            parents[i] = population[rand_int_1];
        }
        else {
            parents[i] = population[rand_int_2];
        }
    }
}
*/


void procreate(){
    uniform_int_distribution<> i_dis(0, PARENTS_NUM-1);
    uniform_int_distribution<> i_dis2(0, IND_PARAMS-1);
    uniform_real_distribution<> r_dis(0.0, 1.0);
    int p1_index, p2_index, crossover_index;
    double prob;
    for (int i = 0; i < MAX_POP_CAP/2; i++){
        p1_index = i_dis(gen);
        p2_index = i_dis(gen);
        prob = r_dis(gen);
        population[2*i] = parents[p1_index];
        population[2*i+1] = parents[p2_index];
        if (prob < CROSSOVER_PROB){
            crossover_index = i_dis2(gen);
            for (int j = crossover_index; j < IND_PARAMS; j++){
                population[2*i][j] = parents[p2_index][j];
                population[2*i+1][j] = parents[p1_index][j];
            }
        }
    }
}

void mutate(){
    uniform_real_distribution<> r_dis(0.0, 1.0);
    uniform_int_distribution<> i_dis(0, MAX_GENE_RANGE);
    for (int i = 0; i < MAX_POP_CAP; i++){
        for (int j = 0; j < IND_PARAMS; j++){
            // Use dis to transform the random unsigned int generated by gen into a
            // double in [0, 255]. Each call to i_dis(gen) generates a new random int
            if (r_dis(gen) < MUT_PROB)
                population[i][j] = i_dis(gen);
        }
    }
}


void init_mutate(){
    uniform_real_distribution<> r_dis(0.0, 1.0);
    uniform_int_distribution<> i_dis(0, MAX_GENE_RANGE);
    for (int i = 0; i < POP_LIM; i++){
        for (int j = 0; j < IND_PARAMS; j++){
            // Use dis to transform the random unsigned int generated by gen into a
            // double in [0, 255]. Each call to i_dis(gen) generates a new random int
            if (r_dis(gen) < MUT_PROB)
                population[i][j] = i_dis(gen);
        }
    }
}


void select(){
    vector<pair<double, int> > people_performance;
    people_performance.resize(MAX_POP_CAP);
    long double performance;
    for (int i = 0; i < MAX_POP_CAP; i++){
        performance = performance_function(i);
        people_performance[i] = make_pair(performance, i);
    }
    sort(people_performance.begin(), people_performance.end(), greater<pair<long double, int> >());
    vector<vector<int> > population_copy;
    for (auto item : population[people_performance[0].second])
        output << item << " ";
    output << endl;
    output << grammar(population[people_performance[0].second]) << ": " <<  CARDS_CONST*1.0/log(1+performance_function(people_performance[0].second)) << " " <<
           CARDS_CONST*1.0/log(1+performance_function2(people_performance[0].second)) << endl;
    cout << grammar(population[people_performance[0].second]) << ": " <<  CARDS_CONST*1.0/log(1+performance_function(people_performance[0].second)) << " " <<
           CARDS_CONST*1.0/log(1+performance_function2(people_performance[0].second)) << endl;
    population_copy.resize(POP_LIM);
    for (int i = 0; i < POP_LIM; i++){
        population_copy[i] = population[i];
        if (people_performance[i].second <= i)
            population[i] = population_copy[people_performance[i].second];
        else
            population[i] = population[people_performance[i].second];
    }
}


/*
void select(){
    vector<pair<double, int> > people_cost;
    double cost;
    people_cost.resize(MAX_POP_CAP);
    for (int i = 0; i < MAX_POP_CAP; i++){
        cost = cost_function(i);
        people_cost[i] = make_pair(cost, i);
    }
    sort(people_cost.begin(), people_cost.end());
    vector<vector<int> > population_copy;
    population_copy.resize(POP_LIM);
    cout << grammar(population[people_cost[0].second]) << " :" <<  cost_function(people_cost[0].second);
    for (int i = 0; i < POP_LIM; i++){
        population_copy[i] = population[i];
        if (people_cost[i].second <= i)
            population[i] = population_copy[people_cost[i].second];
        else
            population[i] = population[people_cost[i].second];
    }
}
*/

vector<char> non_t = {'E', 'O', 'V', 'C', 'N', 'n', 'd'};
bool mapp(string &exp, int value){
    int i = 0;
    for(; i < exp.size(); i++){
        if (find(non_t.begin(), non_t.end(),exp[i])!=non_t.end()){
            break; //symbol exp[i] is non_t
        }
    }
    if(i == exp.size()) return false; // only terminals

    char actual_exp = exp[i];
    exp.erase(i, 1);
    switch(actual_exp){
        case('E'):
            switch(value % 4){
                case(0):
                    exp.insert(i, "O(E)(E)");
                    break;
                case(1):
                    exp.insert(i, "p(E)");
                    break;
                case(2):
                    exp.insert(i, "C");
                    break;
                case(3):
                    exp.insert(i, "*(C)(V)");
                    break;
            }
            break;
        case('O'):
            switch(value % 5){
                case(0):
                    exp.insert(i, "+");
                    break;
                case(1):
                    exp.insert(i, "-");
                    break;
                case(2):
                    exp.insert(i, "*");
                    break;
                case(3):
                    exp.insert(i, "/");
                    break;
                case(4):
                    exp.insert(i, "^");
                    break;
            }
            break;

        case('V'):
            switch(value % 8){
                case(0):
                    exp.insert(i, "x");
                    break;
                case(1):
                    exp.insert(i, "e");
                    break;
                case(2):
                    exp.insert(i, "f");
                    break;
                case(3):
                    exp.insert(i, "g");
                    break;
                case(4):
                    exp.insert(i, "h");
                    break;
                case(5):
                    exp.insert(i, "i");
                    break;
                case(6):
                    exp.insert(i, "j");
                    break;
                case(7):
                    exp.insert(i, "k");
                    break;
            }
            break;
        case('C'):
            switch(value % 2){
                case(0):
                    exp.insert(i, "~N");
                    break;
                case(1):
                    exp.insert(i, "N");
                    break;
            }
            break;
        case('N'):
            switch(value % 2){
                case(0):
                    exp.insert(i, "n");
                    break;
                case(1):
                    exp.insert(i, "n.d");
                    break;
            }
            break;
        case('n'):
            switch(value % 2){
                case(0):
                    exp.insert(i, "nd");
                    break;
                case(1):
                    exp.insert(i, "d");
                    break;
            }
            break;
        case('d'):
            exp.insert(i, to_string(value%10));
            break;
        default: cout << "explodiu" << endl;
    }
    return true;
}

string grammar(vector<int> ind){
    string exp = "E";
    int i = 0;
    int j = 0;
    while(mapp(exp,ind[i]) && j < MAX_LOOP){
        i++;
        if(i == (int)ind.size()) {i = 0; j++;};
    };
    if(j == MAX_LOOP){
        return "invalid";
    }
    return exp;

}

double evaluate_expression(string expression, vector<double>& values){
//    cout << "begin " << chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - last_saved).count() << endl;
//    last_saved = std::chrono::high_resolution_clock::now();
    if(expression == "invalid") { return NAN;}
    while(expression[0] == '('){
        expression = expression.substr(1, expression.size() - 2);
    }
    bool stop = false, is_number = true;
    int actual = 0, state = 0;
    if (expression[0] == '~')
        is_number = true;
    while (!stop && is_number){
        if (state == 0){
            if (expression[actual] >= '0' && expression[actual] <= '9'){
                actual++;
            }
            else if (expression[actual] == '.'){
                state = 1;
                actual++;
            }
            else{
                stop = true;
                is_number = false;
            }
        }
        else{
            if (expression[actual] >= '0' && expression[actual] <= '9'){
                actual++;
            }
            else{
                stop = true;
                is_number = false;
            }
        }
        if (actual == expression.size())
            stop = true;
    }
    if (is_number){
        string::size_type found = expression.find('~');
        if (found != string::npos){
            expression.replace(found, 1, "-");
        }
        double value = stod(expression);
        return value;
    }
    if (expression == "x"){
        return values[0];
    }
    if (expression == "e" || expression == "f" || expression == "g" || expression == "h" || expression == "i" || expression == "j" || expression == "k") {
        return values[expression[0] - 'd'];
    }
    if (expression[0] == 'p'){
        int stack = 0, begin = -1, end = -1;
        for (int i = 0; i < expression.size(); i++){
            if (stack == 0 && expression[i] == '(')
                begin = i;
            if (expression[i] == '(')
                stack++;
            if (expression[i] == ')')
                stack--;
            if (stack == 0 && expression[i] == ')')
                end = i;
        }
        double expression1Value;
        expression1Value = evaluate_expression(expression.substr(begin, end-begin+1), values);
        return exp(expression1Value);
    }
    if (expression[0] == '+' || expression[0] == '-' || expression[0] == '*' || expression[0] == '/' || expression[0] == '^'){
        int stack = 0, begin1 = -1, end1 = -1, begin2 = -1, end2 = -1;
        for (int i = 0; i < expression.size(); i++){
            if (stack == 0 && expression[i] == '(' && begin1 == -1)
                begin1 = i;
            if (stack == 0 && expression[i] == '(' && begin1 != -1)
                begin2 = i;
            if (expression[i] == '(')
                stack++;
            if (expression[i] == ')')
                stack--;
            if (stack == 0 && expression[i] == ')' && end1 == -1)
                end1 = i;
            if (stack == 0 && expression[i] == ')' && end1 != -1)
                end2 = i;
        }
        double expression1Value, expression2Value;
        expression1Value = evaluate_expression(expression.substr(begin1, end1-begin1+1), values);
        expression2Value = evaluate_expression(expression.substr(begin2, end2-begin2+1), values);
        if (expression[0] == '+')
            return expression1Value+expression2Value;
        if (expression[0] == '-')
            return expression1Value-expression2Value;
        if (expression[0] == '*')
            return expression1Value*expression2Value;
        if (expression[0] == '/')
            return expression1Value/expression2Value;
        if (expression[0] == '^')
            return pow(expression1Value, expression2Value);
    }
    return 0.0;
}

void readCSV(const string& filename){
    ifstream data(filename);
    string line;
    while(getline(data,line))
    {
        stringstream lineStream(line);
        string cell;
        vector<double> parsedRow;
        int i = 0;
        while(getline(lineStream,cell,','))
        {
            if(i == 0){
                i++;
            }
            else
                parsedRow.push_back(stod(cell));
        }

        parsedCsv.push_back(parsedRow);
    }
}


int main(){


    readCSV("../training.csv");
    population.resize(MAX_POP_CAP);
    parents.resize(PARENTS_NUM);

    for (int i = 0; i < POP_LIM; i++)
        population[i].resize(IND_PARAMS); // pre allocate vector size for speed

    for (double purge_ratio = 0.1; purge_ratio < 0.41; purge_ratio += 0.1) {
        CARDS_CONST = 10;
        PURGED_PARENTS_RATIO = purge_ratio;
        for (int const_it = 0; const_it < 3; const_it++) {
            output = ofstream("results"+to_string(CARDS_CONST)+"_"+to_string(purge_ratio)+".txt", ofstream::out);

            uniform_int_distribution<> i_dis(0, MAX_GENE_RANGE);
            for (int i = 0; i < (1-PRESET_RATIO) * POP_LIM; i++)
                for (int j = 0; j < IND_PARAMS; j++)
                    population[i][j] = i_dis(gen); // random initialization of the individuals

            for (int i = (1-PRESET_RATIO) * POP_LIM; i < (int)(1-2*PRESET_RATIO/3) * POP_LIM; i++)
                population[i] = PRESET1; // preset 1

            for (int i = (int)(1-2*PRESET_RATIO/3) * POP_LIM; i < (int)(1-1*PRESET_RATIO/3) * POP_LIM; i++)
                population[i] = PRESET2; // preset 2

            for (int i = (int)(1-1*PRESET_RATIO/3) * POP_LIM; i < (int)(1-0*PRESET_RATIO/3) * POP_LIM; i++)
                population[i] = PRESET3; // preset 3


            for(int i = 0; i < INIT_MUT; i++)   init_mutate();




            int lim_gen = GEN_LIM;
            while (lim_gen--) {
                select_parents();
                procreate();
                mutate();
                select();
                output << lim_gen << endl;
            }
            output.close();
            cout << CARDS_CONST << "; " << PURGED_PARENTS_RATIO << endl;
            CARDS_CONST *= 2;
        }
    }



    return 0;
}
