#include <math.h>
#include <stdio.h>
#include <vector>
#include <fstream>
#include <iostream>

using namespace std;

int NumberOfVariables = 4;
int NumberOfStadies = 13;

double Tol(1.e-15);
double Fac(0.8);
double Facmin(0.7);
double Facmax(1.5);

const vector<vector<double>> a = {
    {},
    {1. / 18}, 
    {1. / 48, 1. / 16}, 
    {1. / 32, 0., 3. / 32}, 
    {5. / 16, 0., -75. / 64, 75. / 64}, 
    {3. / 80, 0., 0., 3. / 16, 3. / 20},
    {29443841. / 614563906, 0., 0., 77736538. / 692538347, -28693883. / 1125000000, 23124283. / 1800000000},
    {16016141. / 946692911,0., 0., 61564180. / 158732637, 22789713. / 633445777, 545815736. / 2771057229, -180193667. / 1043307555},
    {39632708. / 573591083, 0., 0., -433636366. / 683701615, -421739975. / 2616292301, 100302831. / 723423059, 790204164. / 839813087, 800635310. / 3783071287},
    {246121993. / 1340847787, 0., 0., -37695042795. / 15268766246, -309121744. / 1061227803, -12992083. / 490766935, 6005943493. / 2108947869, 393006217. / 1396673457, 123872331. / 1001029789},
    {-1028468189. / 846180014, 0., 0., 8478235783. / 508512852, 1311729495. / 1432422823, -10304129995. / 1701304382, -48777925059. / 3047939560, 15336726248. / 1032824649, -45442868181. / 3398467696, 3065993473. / 597172653},
    {185892177. / 718116043, 0., 0., -3185094517. / 667107341, -477755414. / 1098053517, -703635378. / 230739211, 5731566787. / 1027545527, 5232866602. / 850066563, -4093664535. / 808688257, 3962137247. / 1805957418, 65686358. / 487910083},
    {403863854. / 491063109, 0., 0., -5068492393. / 434740067, -411421997. / 543043805, 652783627. / 914296604, 11173962825. / 925320556, -13158990841. / 6184727034, 3936647629. / 1978049680, -160528059. / 685178525, 248638103. / 1413531060, 0.}
};

const vector<double> c = {
    0., 1. / 18, 1. / 12, 1. / 8, 5. / 16, 3. / 8, 59. / 400, 93. / 200, 5490023248. / 9719169821, 13. / 20, 1201146811. / 12990119798, 1., 1.
};

const vector<double> b = {
    14005451. / 335480064, 0., 0., 0., 0., -59238493. / 1068277825, 181606767. / 758867731, 561292985. / 797845732, -1041891430. / 1371343529, 760417239. / 1151165299, 118820643. / 751138087, -528747749. / 2220607170, 1. / 4
};

const vector<double> lowB = {
    13451932. / 455176623, 0., 0., 0., 0., -808719846. / 976000145, 1757004468. / 5645159321, 656045339. / 265891186, -3867574721. / 1518517206, 465885868. / 322736535, 53011238. / 667516719, 2. / 45, 0.
};

double Max(vector<double>& arr) {
    double max = arr[0];
    for (int i = 1; i < arr.size(); ++i) {
        if (arr[i] > max) {
            max = arr[i];
        }
    }
    return max;
}

double FuncX(double t, vector<double> &variables) {
    return variables[1];
}

double FuncY(double t, vector<double> &variables) {
    return - variables[0];
}

double FuncPx(double t, vector<double> &variables) {
    return - variables[0];
}

double FuncPy(double t, vector<double> &variables) {
    return - variables[2] - variables[1];
}

double Func(double t, vector<double> &variables, int var) {
    switch(var) {
    case 0 : 
        return FuncX(t, variables);
    case 1 : 
        return FuncY(t, variables);
    case 2 : 
        return FuncPx(t, variables);
    case 3 :
        return FuncPy(t, variables);
    }

    return -1;
}

double FindLinearCombinationForA(int matrixIndex, int var, vector<vector<double>>& k) {
    double lc = 0.;
    for (int j = 0; j < a[matrixIndex].size(); ++j) {
        lc += a[matrixIndex][j] * k[var][j];
    }
    
    return lc;
}

void UpdateLowVariables(vector<double> &variables, vector<vector<double>> &k, double h) {
    for (int i = 0; i < variables.size(); ++i) {
        double lc = 0.;
        for (int j = 0; j < lowB.size(); ++j) {
            lc += lowB[j] * k[i][j];
        }

        variables[i] += lc * h;
    }
}

void UpdateVariables(vector<double> &variables, vector<vector<double>> &k, double h) {
    for (int i = 0; i < NumberOfVariables; ++i) {
        double lc = 0.;
        for (int j = 0; j < b.size(); ++j) {
            lc += b[j] * k[i][j];
        }

        variables[i] += lc * h;
    }
}

vector<vector<double>> FindK(double t, vector<double> &variables, double h) {
    vector<vector<double>> k(4);
    for (int i = 0; i < NumberOfStadies; ++i) {
        for (int j = 0; j < NumberOfVariables; ++j) {
            double tShift = t + c[i] * h;
            vector<double> variablesWithShift = {
                variables[0] + h * FindLinearCombinationForA(i, 0, k), 
                variables[1] + h * FindLinearCombinationForA(i, 1, k), 
                variables[2] + h * FindLinearCombinationForA(i, 2, k),
                variables[3] + h * FindLinearCombinationForA(i, 3, k)
            };
            k[j].push_back(Func(tShift, variablesWithShift, j));
        }
    }
    UpdateLowVariables(variables, k, h);
    return k;
}

double StepChoice(double t, vector<double> &variables, double& currH) {
    double h = currH;
    for (;;) {
        vector<double> variablesCopy = variables;
        vector<double> wVariables = variables; 
        vector<vector<double>> k;
        k = FindK(t, variablesCopy, h);
        UpdateVariables(wVariables, k, h);
        vector<double> d = {
            fabs(variablesCopy[0] - wVariables[0]),
            fabs(variablesCopy[1] - wVariables[1]),
            fabs(variablesCopy[2] - wVariables[2]),
            fabs(variablesCopy[3] - wVariables[3]),
        };
        double err = Max(d);
        currH = h;
        h *= min(Facmax, max(Facmin, Fac * pow((Tol / err), 1. / 9)));
        if (err <= Tol) {
            variables = variablesCopy;
            return h;
        }
    }
}

int main() {
    ofstream out;
    int draw = 0;
    vector<double> variables = {0., 1., 0., 0.};
    double nextH = 0.01, currH = 0.01;
    out.open("auto_rk.txt");
    for (double T = 0.; T <= 1000000 * 3.1415926535897932; T += currH) {
        currH = nextH;
        if (!(draw % 1)) {
            out << variables[0] - sin(T) << " & " << variables[1] - cos(T) << endl;
        }
        nextH = StepChoice(T, variables, currH);
        draw++;
    }
    cout << draw << endl;
    out.close();
    return 0;
}
