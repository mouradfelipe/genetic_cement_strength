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

using namespace std;

bool mapp(string &exp, int value){
    vector<char> non_t = {'E', 'O', 'P', 'V', 'C', 'N', 'n', 'd'};
    int i = 0;
    for(i; i < exp.size(); i++){
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
                    exp.insert(i, "P(E)");
                    break;
                case(2):
                    exp.insert(i, "(C)");
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
        if (expression[0] == 's')
            return sin(expression1Value);
        if (expression[0] == 'c')
            return cos(expression1Value);
        if (expression[0] == 'l')
            return log(expression1Value);
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


int main(){

    readCSV("../training.csv");






    vector<int> teste = {40, 5, 131, 20, 5, 5999, 57, 65, 102};
    string exp = grammar(teste);
    cout << exp << ": " << evaluate_expression(exp, parsedCsv[0]) << endl;

    return 0;
};
