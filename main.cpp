#include <math.h>
#include <stdio.h>
#include <vector>
#include <fstream>
#include <iostream>

using namespace std;

vector<vector<double>> k = {
    {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}, // x
    {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.} // y
};


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

double tol(1.e-15);
double fac(0.8);
double facmin(0.7);
double facmax(1.5);

using namespace std;

double FuncX(double t, double y, double x) {
    return y;
}

double FuncY(double t, double y, double x) {
    return -x;
}

double FindLinearCombinationForA(int matrixIndex, bool var) {
    double lc = 0.;
    for (int j = 0; j < a[matrixIndex].size(); ++j) {
        lc += a[matrixIndex][j] * k[var][j];
    }
    
    return lc;
}

double FindLinearCombinationForLowB(bool var) {
    double lc = 0.;
    for (int i = 0; i < lowB.size(); ++i) {
        lc += lowB[i] * k[var][i];
    }
    
    return lc;
}

double FindLinearCombinationForB(bool var) {
    double lc = 0.;
    for (int i = 0; i < b.size(); ++i) {
        lc += b[i] * k[var][i];
    }
    
    return lc;
}

void FindK(double t, double &y, double &x, double h) {
    for (int i = 0; i < k[0].size(); ++i) {
        k[0][i] = FuncX(t + c[i] * h, y + h * FindLinearCombinationForA(i, 1), x + h * FindLinearCombinationForA(i, 0));
        k[1][i] = FuncY(t + c[i] * h, y + h * FindLinearCombinationForA(i, 1), x + h * FindLinearCombinationForA(i, 0));
    }
    x += h * FindLinearCombinationForLowB(0);
    y += h * FindLinearCombinationForLowB(1);
}

double StepChoice(double t, double &y, double &x, double& currH) {
    double h = currH;
    for (;;) {
        double locX = x, locY = y, wX = x, wY = y;
        FindK(t, locY, locX, h);
        wX += h * FindLinearCombinationForB(0);
        wY += h * FindLinearCombinationForB(1);
        double d1 = fabs(wX - locX);
        double d2 = fabs(wY - locY);
        double err = max(d1, d2);
        currH = h;
        h *= min(facmax, max(facmin, fac * pow((tol / err), 1. / 9)));
        if (err <= tol) {
            x = locX;
            y = locY;
            return h;
        }
    }
}

int main() {
    ofstream out;
    int draw = 0;
    double nextH = 0.01, x = 0., y = 1., currH = 0.01;
    out.open("auto_rk.txt");
    for (double T = 0.; T <= 100000 * 3.1415926535897932; T += currH) {
        currH = nextH;
        if (!(draw % 400)) {
            out << x - sin(T) << " " << y - cos(T) << endl;
        }
        nextH = StepChoice(T, y, x, currH);
        draw++;
    }
    out.close();
    return 0;
}
