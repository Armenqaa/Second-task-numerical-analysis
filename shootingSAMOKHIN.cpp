#include <math.h>
#include <stdio.h>
#include <utility>
#include <vector>
#include <fstream>
#include <iostream>

using namespace std;

ofstream Out;

const double Period = 10.;
const int NumberOfVariables = 4;
const int NumberOfStadies = 13;

const double Tol(1.e-12);
const double Fac(0.8);
const double Facmin(0.7);
const double Facmax(1.5);

const vector<vector<double>> A = {
    {},
    {1. / 18}, 
    {1. / 48, 1. / 16}, 
    {1. / 32, 0., 3. / 32}, 
    {5. / 16, 0., -75. / 64, 75. / 64}, 
    {3. / 80, 0., 0., 3. / 16, 3. / 20},
    {29443841. / 614563906, 0., 0., 77736538. / 692538347, -28693883. / 1125000000, 23124283. / 1800000000},
    {16016141. / 946692911, 0., 0., 61564180. / 158732637, 22789713. / 633445777, 545815736. / 2771057229, -180193667. / 1043307555},
    {39632708. / 573591083, 0., 0., -433636366. / 683701615, -421739975. / 2616292301, 100302831. / 723423059, 790204164. / 839813087, 800635310. / 3783071287},
    {246121993. / 1340847787, 0., 0., -37695042795. / 15268766246, -309121744. / 1061227803, -12992083. / 490766935, 6005943493. / 2108947869, 393006217. / 1396673457, 123872331. / 1001029789},
    {-1028468189. / 846180014, 0., 0., 8478235783. / 508512852, 1311729495. / 1432422823, -10304129995. / 1701304382, -48777925059. / 3047939560, 15336726248. / 1032824649, -45442868181. / 3398467696, 3065993473. / 597172653},
    {185892177. / 718116043, 0., 0., -3185094517. / 667107341, -477755414. / 1098053517, -703635378. / 230739211, 5731566787. / 1027545527, 5232866602. / 850066563, -4093664535. / 808688257, 3962137247. / 1805957418, 65686358. / 487910083},
    {403863854. / 491063109, 0., 0., -5068492393. / 434740067, -411421997. / 543043805, 652783627. / 914296604, 11173962825. / 925320556, -13158990841. / 6184727034, 3936647629. / 1978049680, -160528059. / 685178525, 248638103. / 1413531060, 0.}
};

const vector<double> C = {
    0., 1. / 18, 1. / 12, 1. / 8, 5. / 16, 3. / 8, 59. / 400, 93. / 200, 5490023248. / 9719169821, 13. / 20, 1201146811. / 12990119798, 1., 1.
};

const vector<double> B = {
    14005451. / 335480064, 0., 0., 0., 0., -59238493. / 1068277825, 181606767. / 758867731, 561292985. / 797845732, -1041891430. / 1371343529, 760417239. / 1151165299, 118820643. / 751138087, -528747749. / 2220607170, 1. / 4
};

const vector<double> LowB = {
    13451932. / 455176623, 0., 0., 0., 0., -808719846. / 976000145, 1757004468. / 5645159321, 656045339. / 265891186, -3867574721. / 1518517206, 465885868. / 322736535, 53011238. / 667516719, 2. / 45, 0.
};

double Module(vector<double> &arr) {
    double res = 0.;
    for (const auto& elem : arr) {
        res += elem * elem;
    }
    
    return sqrt(res);
}

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
    // return variables[1];
    return variables[1];
}

double FuncY(double t, vector<double> &variables) {
    // return fabs(variables[3]) > 1 ? (variables[3] >= 1 ? 1 : -1) : variables[3];
    return variables[3];
}

double FuncPx(double t, vector<double> &variables) {
    return - variables[0];
    // return variables[0];
}

double FuncPy(double t, vector<double> &variables) {
    return - variables[2] - variables[1];
    // return variables[1];
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
    for (int j = 0; j < A[matrixIndex].size(); ++j) {
        lc += A[matrixIndex][j] * k[var][j];
    }
    
    return lc;
}

void UpdateLowVariables(vector<double> &variables, vector<vector<double>> &k, double h) {
    for (int i = 0; i < NumberOfVariables; ++i) {
        double lc = 0.;
        for (int j = 0; j < LowB.size(); ++j) {
            lc += LowB[j] * k[i][j];
        }

        variables[i] += lc * h;
    }
}

void UpdateVariables(vector<double> &variables, vector<vector<double>> &k, double h) {
    for (int i = 0; i < NumberOfVariables; ++i) {
        double lc = 0.;
        for (int j = 0; j < B.size(); ++j) {
            lc += B[j] * k[i][j];
        }

        variables[i] += lc * h;
    }
}

vector<vector<double>> FindK(double t, vector<double> &variables, double h) {
    vector<vector<double>> k(4);
    for (int i = 0; i < NumberOfStadies; ++i) {
        for (int j = 0; j < NumberOfVariables; ++j) {
            double tShift = t + C[i] * h;
            vector<double> variablesWithShift = {
                variables[0] + h * FindLinearCombinationForA(i, 0, k), 
                variables[1] + h * FindLinearCombinationForA(i, 1, k), 
                variables[2] + h * FindLinearCombinationForA(i, 2, k),
                variables[3] + h * FindLinearCombinationForA(i, 3, k)
            };
            k[j].push_back(Func(tShift, variablesWithShift, j));
            /* cout << j << " ";
            for (int i = 0; i < k[j].size(); ++i) {
                cout << k[j][i] << " ";
            } 
            cout << endl; */
        }
    }
    UpdateLowVariables(variables, k, h);
    return k;
}

double StepChoice(double t, vector<double> &variables, double &currH) {
    double h = currH;
    for (;;) {
        vector<double> variablesCopy = variables;
        vector<double> wVariables = variables; 
        vector<vector<double>> k(4);
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
        // cout << t << " ";
        h *= min(Facmax, max(Facmin, Fac * pow((Tol / err), 1. / 9)));
        if (err <= Tol) {
            variables = variablesCopy;
            if (Period - (t + currH) < h) {
                h = Period - (t + currH);
            }
            return h;
        }
    }
}

pair<vector<double>, vector<double>> MainLoopRK(vector<double> variablesInitialConditions = {0., 0., 1., 0.00001}) {
    int draw = 0;
    ofstream out1, out2, out3;
    out1.open("auto_rk1.txt");
    out2.open("auto_rk2.txt");
    out3.open("auto_rk3.txt");
    vector<double> variables = variablesInitialConditions;
    double nextH = 0.01, currH = 0.01;
    double T = 0.;
    for (; 1.e-13 <= Period - T; T += currH) {
        currH = nextH;
        if (!(draw % 1)) {
            Out << T << " " << variables[0] << endl;
            out1 << T << " " << variables[1] << endl;
            out2 << T << " " << variables[2] << endl;
            out3 << T << " " << variables[3] << endl;
        }
        nextH = StepChoice(T, variables, currH);
        draw++;
    }
    Out << T << " " << variables[0] << endl;
    out1 << T << " " << variables[1] << endl;
    out2 << T << " " << variables[2] << endl;
    out3 << T << " " << variables[3] << endl;
    // cout << variables[0] - 0. << endl;
    // cout << variables[1] - 1. << endl;
    out1.close();
    out2.close();
    out3.close();
    // cout << draw << endl;
    return {variablesInitialConditions, variables};
}

vector<vector<double>> FindDerivativeAndInverse(const vector<double> &X, const vector<double> &alpha) {
    // const pair<vector<double>, vector<double>> RKWithShift =  MainLoopRK({0., 0., alpha[0], alpha[1]});
    const pair<vector<double>, vector<double>> RKWithShiftX =  MainLoopRK({0., alpha[0] + 1.e-8, alpha[1], 0.});
    const pair<vector<double>, vector<double>> RKWithShiftY =  MainLoopRK({0., alpha[0], alpha[1] + 1.e-8, 0.});
    // cout << "Compare: " << X[0] - RKWithShiftY.second[0] << " " << X[1] - (RKWithShiftY.second[1] - 1.) << endl;
    // cout << RKWithShiftX.first[0] << " " << RKWithShiftX.first[1] << endl;
    const vector<double> withShiftX = {
        RKWithShiftX.second[0]/* - RKWithShiftX.first[0]*/, 
        RKWithShiftX.second[1] - 1.
    };
    const vector<double> withShiftY = {
        RKWithShiftY.second[0]/* - RKWithShiftY.first[0]*/, 
        RKWithShiftY.second[1] - 1.
    };
    const vector<vector<double>> XDer = {
        {
            (withShiftX[0] - X[0]) / 1.e-8, 
            (withShiftY[0] - X[0]) / 1.e-8
        }, 
        {
            (withShiftX[1] - X[1]) / 1.e-8, 
            (withShiftY[1] - X[1]) / 1.e-8
        }
    };
    double XDerDet = (XDer[0][0] * XDer[1][1] - XDer[0][1] * XDer[1][0]);
    cout << "DerInv: " << XDer[1][1] / XDerDet << " " << - XDer[0][1] / XDerDet << " " << - XDer[1][0] / XDerDet << " " << XDer[0][0] / XDerDet << endl;
    return {
        {
            XDer[1][1] / XDerDet, 
            - XDer[0][1] / XDerDet
        }, 
        {
            - XDer[1][0] / XDerDet, 
            XDer[0][0] / XDerDet
        }
    };
}

void Shooting(const vector<double> &a) {
    int counter = 0;
    bool firstStep = true;
    vector<double> alpha = a, alphaNew = a, X, XPrev;
    vector<vector<double>> XDerivativeInversed = {{0., 0.,}, {0., 0.}};
    double gamma = 1.;

    pair<vector<double>, vector<double>> resultRK = MainLoopRK({0., a[0], a[1], 0.});
    X = {resultRK.second[0]/* - resultRK.first[0]*/, resultRK.second[1] - 1.};

    do {

        if (counter > 30) {
            cout << endl << "-----------------------------------------" << endl;
            cout << "Bad initial conditions or smth went wrong" << endl;
            cout << "-----------------------------------------" << endl;
            break;
        }

        if (!firstStep) {
            // cout << "X: " << X[0] << " " << X[1] << endl;
            alphaNew[0] = alpha[0] - gamma * (XDerivativeInversed[0][0] * X[0] + XDerivativeInversed[0][1] * X[1]);
            alphaNew[1] = alpha[1] - gamma * (XDerivativeInversed[1][0] * X[0] + XDerivativeInversed[1][1] * X[1]);
            pair<vector<double>, vector<double>> resultRK = MainLoopRK({0., alphaNew[0], alphaNew[1], 0.});
            X = {resultRK.second[0]/* - resultRK.first[0]*/, resultRK.second[1] - 1.};
        }
        XDerivativeInversed = FindDerivativeAndInverse(X, alphaNew);
        if (XPrev.size() > 0 && Module(X) > Module(XPrev)) {
            gamma /= 2;
            counter++;
        } else {
            XPrev = X;
            gamma = 1.;
            alpha = alphaNew;
            counter = 0;
        }
        firstStep = false;
        cout << "alpha 1: " << alphaNew[0] << " alpha 2: " << alphaNew[1] << endl;
        cout << "Module: " << Module(X) << " " << Module(XPrev) << " " << gamma << endl;
    } while (Module(X) > 1.e-8);
}

int main() {
    Out.open("auto_rk.txt");
    
    Shooting({4., -3.});
    Out.close();
    
    cout << endl << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
    cout << "Fine." << endl;
    cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
    return 0;
}
