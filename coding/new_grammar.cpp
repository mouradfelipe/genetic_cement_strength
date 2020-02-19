#include <iostream>
#include <vector>
#include <algorithm>    // std::lower_bound, std::upper_bound, std::sort
#include <functional>   // std::greater
#include <utility>
#include <random>
#include <cmath>
#include <climits>
#include <string>
#include <regex>
#include <fstream>
#include <thread>
#include <cfloat>

using namespace std;


int CARDS_CONST = 256;

bool mapp(string &exp, int value){
    vector<char> non_t = {'E', 'O', 'V', 'C', 'N', 'n', 'd'};
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
    bool end = false;
    string exp = "E";
    int i = 0;
    int j = 0;
    while(mapp(exp,ind[i]) && j < 10){
        i++;
        if(i == ind.size()) {i = 0; j++;};
    };
    if(j == 10){
        return "invalid";
    }
    return exp;

}

long double evaluate_expression(string expression, vector<double>& values){
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
}

vector<vector<double>> parsedCsv;
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



double sigmoid(double x){
    cout << "sigmoid " << x << " " << 1/(1+exp(-x)) << endl;
    return 1/(1+exp(-x));
}

void calculate_cost1(string& expr, long double &cost1){
    for(int i = 0; i < (int)parsedCsv.size()/4; i++){
        long double v = (parsedCsv[i][8] - sigmoid(0.001*evaluate_expression(expr, parsedCsv[i])));
        cost1 += v*v;
        if(isnan(cost1)) break;
    }
}

void calculate_cost2(string& expr, long double &cost2){
    for(int i = parsedCsv.size()/4; i < (int) parsedCsv.size()/2; i++){
        long double v = (parsedCsv[i][8] - sigmoid(0.001*evaluate_expression(expr, parsedCsv[i])));
        cost2 += v*v;
        if(isnan(cost2)) break;
    }
}

void calculate_cost3(string& expr, long double &cost3){
    for(int i = parsedCsv.size()/2; i < ((int) 3*parsedCsv.size())/4; i++){
        long double v = (parsedCsv[i][8] - sigmoid(0.001*evaluate_expression(expr, parsedCsv[i])));
        cost3 += v*v;
        if(isnan(cost3)) break;
    }
}

void calculate_cost4(string& expr, long double &cost4){
    for(int i = (3*parsedCsv.size())/4; i < parsedCsv.size(); i++){
        long double v = (parsedCsv[i][8] - sigmoid(0.001*evaluate_expression(expr, parsedCsv[i])));
        cost4 += v*v;
        if(isnan(cost4)) break;
    }
}



long double performance_function(string expr){
    long double cost = 0;
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

    long double cost1, cost2, cost3, cost4;
    cost1 = cost2 = cost3 = cost4 = 0;
    thread thr1(calculate_cost1, ref(expr), ref(cost1));
    thread thr2(calculate_cost2, ref(expr), ref(cost2));
    thread thr3(calculate_cost3, ref(expr), ref(cost3));
    thread thr4(calculate_cost4, ref(expr), ref(cost4));
    thr1.join();
    thr2.join();
    thr3.join();
    thr4.join();
    cost = cost1 + cost2 + cost3 + cost4;
    //auto t2 = std::chrono::high_resolution_clock::now();
    //auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
    //cout << duration << endl;

    cost = sqrt(cost/(parsedCsv.size()-1));
    cout << cost << endl;
    if(cost == 0)
        return DBL_MAX;
    if(isnan(cost)) return 0;
    return  exp(CARDS_CONST/cost) - 1;
}


int main(){

    readCSV("../training.csv");
    //-(*(*(9)(h))(98))(+(p(/(8)(^(p(1.0))(*(*(6.7)(k))(*(2)(x))))))(-(*(98.8)(e))(+(*(838)(e))(+(+(9)(p(*(~33)(g))))(*(*(28.6)(h))(*(98.8)(e)))))))
    //vector<int> teste = {40, 5, 11, 20, 5, 5999, 57, 102};
    //string exp = grammar(teste);
    string exp =  "-(*(p(*(2)(e)))(*(999.9)(h)))(/(2058)(p(*(*(7)(k))(p(*(3.1)(x))))))";
    long double perf =  performance_function(exp);
    cout << exp << ": " << perf << ": " <<  CARDS_CONST*1.0/log(1+perf) << endl;

    return 0;
};

