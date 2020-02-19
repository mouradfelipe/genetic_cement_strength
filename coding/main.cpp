#include <iostream>
#include <vector>
#include <algorithm>    // std::lower_bound, std::upper_bound, std::sort
#include <functional>   // std::greater
#include <utility>
#include <random>
#include <cfloat>
#include <sstream>
#include <fstream>
#include <regex>
#include <chrono>



using namespace std;

#define MAX_POP_CAP 100*2 // must be even to make things easier on procreating phase
#define POP_LIM 30*2
#define IND_PARAMS 50
#define PARENTS_NUM 20*2
#define CROSSOVER_PROB 1
#define MUT_PROB 1
#define MAX_GENE_RANGE 255 // gene values can go from 0 to 255
#define MAX_LOOP 5 // max number of loops to generate expression
#define RAND_SAMPLE 100

template<class bidiiter>
bidiiter random_unique(bidiiter begin, bidiiter end, size_t num_random) {
    size_t left = std::distance(begin, end);
    while (num_random--) {
        bidiiter r = begin;
        std::advance(r, rand()%left);
        std::swap(*begin, *r);
        ++begin;
        --left;
    }
    return begin;
}


vector<vector<int>> population; // number of individuals in a generation and each one is a IND_PARAMS-vector
vector<vector<int>> parents; // individuals who will procreate

random_device rd;  //Will be used to obtain a seed for the random number engine
mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()


vector<vector<double>> parsedCsv;
string grammar(vector<int> ind);
double evaluate_expression(string expression, vector<double>& values);


/*
double old_performance_function(int index){
    double cost = 0;
    string exp = grammar(population[index]);
    auto t1 = std::chrono::high_resolution_clock::now();
    for(int i = 0; i < RAND_SAMPLE; i++){
        double v = (parsedCsv[i][8] - evaluate_expression(exp, parsedCsv[i]));
        cost += v*v;
        if(isnan(cost)) break;
    }
    auto t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
    //cout << duration << endl;

    cost = sqrt(cost/(parsedCsv.size()-1));
    cout << index << ": " << exp << ": " << cost << endl;
    if(cost == 0)
        return DBL_MAX;
    if(isnan(cost)) return 0;
    return 1.0/cost;
}
*/

double sigmoid(double x){
    return 1/(1+exp(-x));
}

double cost_function(int index){
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

/*
void old_select_parents(){
    vector<double> select_lim; // actual probabilities: [0.1 0.3 0.2 0.4] => [0.1 0.4 0.6 1] (sum of probabilities)
    select_lim.resize(POP_LIM);
    double total_performance = 0;
    random_unique(parsedCsv.begin(), parsedCsv.end(), RAND_SAMPLE);

    for (int i = 0; i < POP_LIM; i++){
        total_performance += performance_function(i);
        select_lim[i] = (total_performance);
    }
    for (int i = 0; i < POP_LIM; i++)
        select_lim[i] /= total_performance;

    uniform_real_distribution<> r_dis(0.0, 1.0);
    double random_real;
    vector<double>::iterator selected;
    for (int i = 0; i < PARENTS_NUM; i++){
        random_real = r_dis(gen);
        selected = lower_bound(select_lim.begin(), select_lim.end(), random_real); // finds first element greater or equal than random_real
        parents[i] = (population[selected-select_lim.begin()]);
    }
}
*/

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

/*
void select(){
    vector<pair<double, int> > people_performance;
    people_performance.resize(MAX_POP_CAP);
    double performance;
    for (int i = 0; i < MAX_POP_CAP; i++){
        performance = performance_function(i);
        people_performance[i] = make_pair(performance, i);
    }
    sort(people_performance.begin(), people_performance.end(), greater<pair<double, int> >());
    vector<vector<int> > population_copy;
    population_copy.resize(POP_LIM);
    for (int i = 0; i < POP_LIM; i++){
        population_copy[i] = population[i];
        if (people_performance[i].second <= i)
            population[i] = population_copy[people_performance[i].second];
        else
            population[i] = population[people_performance[i].second];
    }
}
 */


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



bool mapp(string &exp, int value){
    vector<char> non_t = {'E', 'O', 'P', 'V', 'C', 'N', 'n', 'd'};
    int i = 0;
    for(; i < (int)exp.size(); i++){
        if (find(non_t.begin(), non_t.end(),exp[i])!=non_t.end()){
            break; //symbol exp[i] is non_t
        }
    }
    if(i == (int)exp.size()) return false; // only terminals

    char actual_exp = exp[i];
    exp.erase(i, 1);
    switch(actual_exp){
        case('E'):
            switch(value % 4){
                case(0):
                    exp.insert(i, "O(E)(E)");
                    break;
                case(1):
                    exp.insert(i, "P(E)");
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
        case('P'):
            switch(value % 3){
                case(0):
                    exp.insert(i, "s");
                    break;
                case(1):
                    exp.insert(i, "c");
                    break;
                case(2):
                    exp.insert(i, "l");
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
                    exp.insert(i, "n.n");
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
    if(expression == "invalid") { return 0.0/0;}
    while(expression[0] == '('){
        expression = expression.substr(1, expression.size() - 2);
    }



    regex number("^[~]?([0-9]*[.])?[0-9]+$");
    if (regex_match(expression, number)){
        string::size_type found = expression.find('~');
        if (found != string::npos){
            expression.replace(found, 1, "-");
        }
        double value = stod(expression);
        return value;
    }
    if (expression == "x")
        return values[0];
    if (expression == "e" || expression == "f" || expression == "g" || expression == "h" || expression == "i" || expression == "j" || expression == "k")
        return values[expression[0] - 'd'];
    if (expression[0] == 's' || expression[0] == 'c' || expression[0] == 'l'){
        int stack = 0, begin = -1, end = -1;
        for (int i = 0; i <(int) expression.size(); i++){
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
        if (expression[0] == 's')
            return sin(expression1Value);
        if (expression[0] == 'c')
            return cos(expression1Value);
        if (expression[0] == 'l')
            return log(expression1Value);
    }
    if (expression[0] == '+' || expression[0] == '-' || expression[0] == '*' || expression[0] == '/' || expression[0] == '^'){
        int stack = 0, begin1 = -1, end1 = -1, begin2 = -1, end2 = -1;
        for (int i = 0; i < (int)expression.size(); i++){
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
    return 0;
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

    uniform_int_distribution<> i_dis(0, MAX_GENE_RANGE);
    for (int i = 0; i < POP_LIM; i++)
        for (int j = 0; j < IND_PARAMS; j++)
            population[i][j] = i_dis(gen); // random initialization of the individuals

    int lim_gen = 100;
    while (lim_gen--){
        select_parents();
        procreate();
        mutate();
        select();
        cout << lim_gen << endl;
    }

    vector<pair<double, int> > people_cost;
    people_cost.resize(POP_LIM);
    double cost;
    for (int i = 0; i < POP_LIM; i++){
        cost = cost_function(i);
        people_cost[i] = make_pair(cost, i);
    }
    sort(people_cost.begin(), people_cost.end());
    cout << people_cost[0].second << endl;




    return 0;
}