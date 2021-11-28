#include <math.h>
#include <stdio.h>
#include <utility>
#include <vector>
#include <fstream>
#include <iostream>
#include <random>

using namespace std;

ofstream OutX, OutY, OutPX, OutPY, OutT;

const double Period = 20.0;
const int NumberOfVariables = 5;
const int NumberOfStadies = 13;
double NonZeroXDerDet;
int PreviousBranch = -2, CurrentBranch = -2;

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

double Max(vector<double> &arr) {
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
    double toReturn;
    if (fabs(variables[3]) > 1) {
        if (variables[3] > 1) {
            toReturn = 1;
            CurrentBranch = 1;
        } else {
            toReturn = -1;
            CurrentBranch = -1;
        }
    } else {
        toReturn = variables[3];
        CurrentBranch = 0;
    }

    return toReturn;
}

double FuncPx(double t, vector<double> &variables) {
    return - variables[0];
}

double FuncPy(double t, vector<double> &variables) {
    return - variables[2] - variables[1];
}

double FuncB(double t, vector<double> &variables) {
    double u = fabs(variables[3]) > 1 ? (variables[3] > 1 ? 1 : -1) : variables[3];
    return u * u - variables[1] * variables[1] - variables[0] * variables[0];
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
    case 4 :
        return FuncB(t, variables);
    }
    return -1;
}

double FindLinearCombinationForA(int matrixIndex, int var, vector<vector<double>> &k) {
    double lc = 0.;
    for (int j = 0; j < A[matrixIndex].size(); ++j) {
        lc += A[matrixIndex][j] * k[var][j];
    }
    
    return lc;
}

vector<double> ShiftVariables(vector<double> &variables, double h, int stage, vector<vector<double>> &k) {
    vector<double> variablesWithShift(NumberOfVariables);
    for (int var = 0; var < NumberOfVariables; ++var) {
        variablesWithShift[var] = variables[var] + h * FindLinearCombinationForA(stage, var, k);
    }

    return variablesWithShift;
}

void UpdateLowVariables(vector<double> &variables, vector<vector<double>> &k, double h) {
    for (int var = 0; var < NumberOfVariables; ++var) {
        double lc = 0.;
        for (int j = 0; j < LowB.size(); ++j) {
            lc += LowB[j] * k[var][j];
        }

        variables[var] += lc * h;
    }
}

void UpdateVariables(vector<double> &variables, vector<vector<double>> &k, double h) {
    for (int var = 0; var < NumberOfVariables; ++var) {
        double lc = 0.;
        for (int j = 0; j < B.size(); ++j) {
            lc += B[j] * k[var][j];
        }

        variables[var] += lc * h;
    }
}

vector<vector<double>> FindK(double t, vector<double> &variables, double h) {
    vector<vector<double>> k(NumberOfVariables);
    for (int stage = 0; stage < NumberOfStadies; ++stage) {
        for (int var = 0; var < NumberOfVariables; ++var) {
            double tShift = t + C[stage] * h;
            vector<double> variablesWithShift = ShiftVariables(variables, h, stage, k);
            k[var].push_back(Func(tShift, variablesWithShift, var));
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
        vector<vector<double>> k(NumberOfVariables);
        k = FindK(t, variablesCopy, h);
        UpdateVariables(wVariables, k, h);
        vector<double> d(NumberOfVariables - 1); 
        /// NumberOfVariables - 1 because last var is a functional 
        for (int var = 0; var < NumberOfVariables - 1; ++var) {
            d[var] = fabs(variablesCopy[var] - wVariables[var]);
        }
        double err = Max(d);
        currH = h;
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

pair<vector<double>, vector<double>> MainLoopAutoRK(vector<double> variablesInitialConditions = {0., 0., 1., 0.00001}) {
    int draw = 0, previousBranch = -2;
    vector<double> variables = variablesInitialConditions;
    double nextH = 1.e-1, currH = 1.e-1;
    double T = 0.;
    for (; 1.e-13 <= Period - T; T += currH) {
        currH = nextH;
        if (!(draw % 1)) {
            OutX << T << " " << variables[0] << endl;
            OutY << T << " " << variables[1] << " " << T << endl;
            OutPX << T << " " << variables[2] << endl;
            OutPY << T << " " << variables[3] << endl;
            OutT << T << endl;
        }
        nextH = StepChoice(T, variables, currH);
        draw++;
        if (CurrentBranch != PreviousBranch && PreviousBranch != -2) {
            cout.precision(16);
            cout << "Control changed. T:" << T << " h:" << currH << endl;
        }
        PreviousBranch = CurrentBranch;
    }

    OutX << variables[0] << endl;
    OutY << variables[1] << endl;
    OutPX << variables[2] << endl;
    OutPY << variables[3] << endl;
    OutT << T << endl;
    cout << "Functional: " << variables.back() << endl;
    return {variablesInitialConditions, variables};
}

vector<vector<double>> FindDerivativeAndInverse(const vector<double> &X, const vector<double> &alpha) {
    // const pair<vector<double>, vector<double>> RKWithShift =  MainLoopAutoRK({0., 0., alpha[0], alpha[1]});
    const pair<vector<double>, vector<double>> RKWithShiftX =  MainLoopAutoRK({0., 0., alpha[0] + 1.e-8, alpha[1], 0.});
    const pair<vector<double>, vector<double>> RKWithShiftY =  MainLoopAutoRK({0., 0., alpha[0], alpha[1] + 1.e-8, 0.});
    // cout << "Compare: " << X[0] - RKWithShiftY.second[0] << " " << X[1] - RKWithShiftY.second[3] << endl;
    // cout << RKWithShiftX.first[0] << " " << RKWithShiftX.first[1] << endl;
    const vector<double> withShiftX = {
        RKWithShiftX.second[0]/* - RKWithShiftX.first[0]*/, 
        RKWithShiftX.second[3]
    };
    const vector<double> withShiftY = {
        RKWithShiftY.second[0]/* - RKWithShiftY.first[0]*/, 
        RKWithShiftY.second[3]
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

    if (XDerDet != 0) {
        NonZeroXDerDet = XDerDet;
    } else {
        XDerDet = NonZeroXDerDet;
    }

    // cout << "DerInv: " << XDer[1][1] / XDerDet << " " << - XDer[0][1] / XDerDet << " " << - XDer[1][0] / XDerDet << " " << XDer[0][0] / XDerDet << endl;
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

int Shooting(const vector<double> &a) {
    int counterGamma = 0;
    bool firstStep = true;
    vector<double> alpha = a, alphaNew = a, X, XPrev;
    vector<vector<double>> XDerivativeInversed = {{0., 0.,}, {0., 0.}};
    double gamma = 1.;

    pair<vector<double>, vector<double>> resultRK = MainLoopAutoRK({0., 0., a[0], a[1], 0.});
    X = {resultRK.second[0]/* - resultRK.first[0]*/, resultRK.second[3]};
    cout << "X: " << X[0] << " " << X[1] << endl;
    // Out << Module(X) << endl;

    // return;
    int counterRKIterations = 0;
    do {

        if (counterGamma > 30) {
            cout << endl << "-----------------------------------------" << endl;
            cout << "Bad initial conditions or smth went wrong" << endl;
            cout << "-----------------------------------------" << endl;
            break;
        }

        if (!firstStep) {
            cout << "X: " << X[0] << " " << X[1] << endl;
            alphaNew[0] = alpha[0] - gamma * (XDerivativeInversed[0][0] * X[0] + XDerivativeInversed[0][1] * X[1]);
            alphaNew[1] = alpha[1] - gamma * (XDerivativeInversed[1][0] * X[0] + XDerivativeInversed[1][1] * X[1]);
            pair<vector<double>, vector<double>> resultRK = MainLoopAutoRK({0., 0., alphaNew[0], alphaNew[1], 0.});
            X = {resultRK.second[0]/* - resultRK.first[0]*/, resultRK.second[3]};
        }
        XDerivativeInversed = FindDerivativeAndInverse(X, alphaNew);
        if (XPrev.size() > 0 && Module(X) > Module(XPrev)) {
            gamma /= 2;
            counterGamma++;
        } else {
            XPrev = X;
            gamma = 1.;
            alpha = alphaNew;
            counterGamma = 0;
        }
        firstStep = false;
        counterRKIterations++;
        // cout << "alpha 1: " << alphaNew[0] << " alpha 2: " << alphaNew[1] << endl;
        // Out << "Module: " << Module(X) << " " << Module(XPrev) << " " << gamma << endl;
    } while (Module(X) > 1.e-8 && counterRKIterations < 500);

    if (Module(X) < 1.e-8) {
        cout.precision(16);
        cout << "Converge to " << alphaNew[0] << " " << alphaNew[1];
        return 1;
    }

    return 0;
}

void FindSomeInitialConditionToConverge() {
    double leftAlpha0 = -5000, rightAlpha0 = 5000;
    double leftAlpha1 = -5000, rightAlpha1 = 5000;
    random_device rd1, rd2;
    mt19937 gen1(rd1()), gen2(rd2());
    uniform_real_distribution<> dis1(leftAlpha0, rightAlpha0);
    uniform_real_distribution<> dis2(leftAlpha1, rightAlpha1);
    for (auto i = 0; i < 500; ++i) {
        auto alpha0 = dis1(gen1);
        auto alpha1 = dis2(gen2);
        int res = Shooting({alpha0, alpha1});
        if (res) {
            cout << " with initial conditions : " << alpha0 << " " << alpha1 << endl;
        }
    }
}

void OpenFilesToInsert() {
    OutX.open("X.txt");
    OutY.open("Y.txt");
    OutPX.open("PX.txt");
    OutPY.open("PY.txt");
    OutT.open("T.txt");
}

void CloseFiles() {
    OutX.close();
    OutY.close();
    OutPX.close();
    OutPY.close();
    OutX.close();
}

int main() {
    OpenFilesToInsert();
    // FindSomeInitialConditionToConverge();
    // Shooting({-233.4461948624531, -1335.591291777942});
    MainLoopAutoRK({0., 0., -233.4461948624531, -1335.591291777942, 0.});
    CloseFiles();
    return 0;
}
