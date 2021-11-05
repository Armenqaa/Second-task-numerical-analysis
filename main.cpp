#include <math.h>
#include <stdio.h>
#include <vector>
#include <fstream>
#include <iostream>

using namespace std;

double kY[] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
double kX[] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};


vector<vector<double>> a = {
    {1. / 18}, 
    {1. / 48, 1. / 16}, 
    {1. / 32, 0., 3. / 32}, 
    {5. / 16, 0., -75. / 64, 75. / 64}, 
    {3. / 80, 0., 0., 3. / 16, 3. / 20},
    {29443841. / 614563906, 0., 0., 77736538. / 692538347, -28693883. / 1125000000, 23124283. / 1800000000},
    {16016141. / 946692911,0., 0., 61564180. / 158732637, 22789713. / 633445777, 545815736. / 2771057229, -180193667. / 1043307555},
    {39632708. / 573591083, 0., 0., -433636366. / 683701615, -421739975. / 2616292301, 100302831. / 723423059, 790204164. / 839813087, 800635310. / 3783071287},
    {246121993. / 1340847787, 0., 0., -37695042795. / 15268766246, -309121744. / 1061227803, -12992083. / 490766935, 6005943493. / 2108947869, 393006217. / 1396673457, 123872331. / 1001029789},
    {-1028468189. / 846180014, 0., 0., 8478235783. / 508512852, 1311729495. / 1432422823, -10304129995. / 1701304382, -48777925059. / 3047939560, 15336726248. / 1032824649, -45442868181. / 3398467696},
    {185892177. / 718116043, 0., 0., -3185094517. / 667107341, -477755414. / 1098053517, -703635378. / 230739211, 5731566787. / 1027545527, 5232866602. / 850066563, -4093664535. / 808688257, 3962137247. / 1805957418, 65686358. / 487910083},
    {403863854. / 491063109, 0., 0., -5068492393. / 434740067, -411421997. / 543043805, 652783627. / 914296604, 11173962825. / 925320556, -13158990841. / 6184727034, 3936647629. / 1978049680, -160528059. / 685178525, 248638103. / 1413531060, 0.}
};

double a21(1. / 18);
double a31(1. / 48); double a32(1. / 16);
double a41(1. / 32); double a42(0.); double a43(3. / 32);
double a51(5. / 16); double a52(0.); double a53(-75. / 64); double a54(75. / 64);
double a61(3. / 80); double a62(0.); double a63(0.); double a64(3. / 16); double a65(3. / 20);
double a71(29443841. / 614563906); double a72(0.); double a73(0.); double a74(77736538. / 692538347); double a75(-28693883. / 1125000000); double a76(23124283. / 1800000000);
double a81(16016141. / 946692911); double a82(0.); double a83(0.); double a84(61564180. / 158732637); double a85(22789713. / 633445777); double a86(545815736. / 2771057229); double a87(-180193667. / 1043307555);
double a91(39632708. / 573591083); double a92(0.); double a93(0.); double a94(-433636366. / 683701615); double a95(-421739975. / 2616292301); double a96(100302831. / 723423059); double a97(790204164. / 839813087); double a98(800635310. / 3783071287);
double a10_1(246121993. / 1340847787); double a10_2(0.); double a10_3(0.); double a10_4(-37695042795. / 15268766246); double a10_5(-309121744. / 1061227803); double a10_6(-12992083. / 490766935); double a10_7(6005943493. / 2108947869); double a10_8(393006217. / 1396673457); double a10_9(123872331. / 1001029789);
double a11_1(-1028468189. / 846180014); double a11_2(0.); double a11_3(0.); double a11_4(8478235783. / 508512852); double a11_5(1311729495. / 1432422823); double a11_6(-10304129995. / 1701304382); double a11_7(-48777925059. / 3047939560); double a11_8(15336726248. / 1032824649); double a11_9(-45442868181. / 3398467696); double a11_10(3065993473. / 597172653);
double a12_1(185892177. / 718116043); double a12_2(0.); double a12_3(0.); double a12_4(-3185094517. / 667107341); double a12_5(-477755414. / 1098053517); double a12_6(-703635378. / 230739211); double a12_7(5731566787. / 1027545527); double a12_8(5232866602. / 850066563); double a12_9(-4093664535. / 808688257); double a12_10(3962137247. / 1805957418); double a12_11(65686358. / 487910083);
double a13_1(403863854. / 491063109); double a13_2(0.); double a13_3(0.); double a13_4(-5068492393. / 434740067); double a13_5(-411421997. / 543043805); double a13_6(652783627. / 914296604); double a13_7(11173962825. / 925320556); double a13_8(-13158990841. / 6184727034); double a13_9(3936647629. / 1978049680); double a13_10(-160528059. / 685178525); double a13_11(248638103. / 1413531060); double a13_12(0.);

double c2(1. / 18);
double c3(1. / 12);
double c4(1. / 8);
double c5(5. / 16);
double c6(3. / 8);
double c7(59. / 400);
double c8(93. / 200);
double c9(5490023248. / 9719169821);
double c10(13. / 20);
double c11(1201146811. / 12990119798);
double c12(1.);
double c13(1.);

double b1(14005451. / 335480064);
double b2(0.);
double b3(0.);
double b4(0.);
double b5(0.);
double b6(-59238493. / 1068277825);
double b7(181606767. / 758867731);
double b8(561292985. / 797845732);
double b9(-1041891430. / 1371343529);
double b10(760417239. / 1151165299);
double b11(118820643. / 751138087);
double b12(-528747749. / 2220607170);
double b13(1. / 4);

double lowB1(13451932. / 455176623);
double lowB2(0.);
double lowB3(0.);
double lowB4(0.);
double lowB5(0.);
double lowB6(-808719846. / 976000145);
double lowB7(1757004468. / 5645159321);
double lowB8(656045339. / 265891186);
double lowB9(-3867574721. / 1518517206);
double lowB10(465885868. / 322736535);
double lowB11(53011238. / 667516719);
double lowB12(2. / 45);
double lowB13(0.);

double tol(1.e-15);
double fac(0.8);
double facmin(0.7);
double facmax(1.5);

using namespace std;

double FuncY(double t, double y, double x) {
    return -x;
}

double FuncX(double t, double y, double x) {
    return y;
}

void FindK(double t, double &y, double &x, double h) {
    for (int i = 0; i < 13; ++i) {
        // FindLinearCombination()
    }
    k1X = FuncX(t, y, x);
    k1Y = FuncY(t, y, x);
    k2X = FuncX(t + c2 * h, y + h * a21 * k1Y, x + h * a21 * k1X);
    k2Y = FuncY(t + c2 * h, y + h * a21 * k1Y, x + h * a21 * k1X);
    k3X = FuncX(t + c3 * h, y + h * (a31 * k1Y + a32 * k2Y), x + h * (a31 * k1X + a32 * k2X));
    k3Y = FuncY(t + c3 * h, y + h * (a31 * k1Y + a32 * k2Y), x + h * (a31 * k1X + a32 * k2X));
    k4X = FuncX(t + c4 * h, y + h * (a41 * k1Y + a42 * k2Y + a43 * k3Y), x + h * (a41 * k1X + a42 * k2X + a43 * k3X));
    k4Y = FuncY(t + c4 * h, y + h * (a41 * k1Y + a42 * k2Y + a43 * k3Y), x + h * (a41 * k1X + a42 * k2X + a43 * k3X));
    k5X = FuncX(t + c5 * h, y + h * (a51 * k1Y + a52 * k2Y + a53 * k3Y + a54 * k4Y), x + h * (a51 * k1X + a52 * k2X + a53 * k3X + a54 * k4X));
    k5Y = FuncY(t + c5 * h, y + h * (a51 * k1Y + a52 * k2Y + a53 * k3Y + a54 * k4Y), x + h * (a51 * k1X + a52 * k2X + a53 * k3X + a54 * k4X));
    k6X = FuncX(t + c6 * h, y + h * (a61 * k1Y + a62 * k2Y + a63 * k3Y + a64 * k4Y + a65 * k5Y), x + h * (a61 * k1X + a62 * k2X + a63 * k3X + a64 * k4X + a65 * k5X));
    k6Y = FuncY(t + c6 * h, y + h * (a61 * k1Y + a62 * k2Y + a63 * k3Y + a64 * k4Y + a65 * k5Y), x + h * (a61 * k1X + a62 * k2X + a63 * k3X + a64 * k4X + a65 * k5X));
    k7X = FuncX(t + c6 * h, y + h * (a61 * k1Y + a62 * k2Y + a63 * k3Y + a64 * k4Y + a65 * k5Y), x + h * (a61 * k1X + a62 * k2X + a63 * k3X + a64 * k4X + a65 * k5X));
    k7Y = FuncY(t + c6 * h, y + h * (a61 * k1Y + a62 * k2Y + a63 * k3Y + a64 * k4Y + a65 * k5Y), x + h * (a61 * k1X + a62 * k2X + a63 * k3X + a64 * k4X + a65 * k5X));
    y = y + h * (lowB1 * k1Y + lowB2 * k2Y + lowB3 * k3Y + lowB4 * k4Y + lowB5 * k5Y + lowB6 * k6Y);
    x = x + h * (lowB1 * k1X + lowB2 * k2X + lowB3 * k3X + lowB4 * k4X + lowB5 * k5X + lowB6 * k6X);
}

double StepChoice(double t, double &y, double &x, double h) {
    // cout << "-----------------------" << endl;
    for (;;) {
        // cout << h << endl;
        double locX = x, locY = y, wX = x, wY = y;
        FindK(t, locY, locX, h);
        double k7X = FuncX(t + c7 * h, y + h * (a71 * k1Y + a72 * k2Y + a73 * k3Y + a74 * k4Y + a75 * k5Y + a76 * k6Y), x + h * (a71 * k1X + a72 * k2X + a73 * k3X + a74 * k4X + a75 * k5X + a76 * k6X));
        double k7Y = FuncY(t + c7 * h, y + h * (a71 * k1Y + a72 * k2Y + a73 * k3Y + a74 * k4Y + a75 * k5Y + a76 * k6Y), x + h * (a71 * k1X + a72 * k2X + a73 * k3X + a74 * k4X + a75 * k5X + a76 * k6X));
        wX = x + h * (b1 * k1X + b2 * k2X + b3 * k3X + b4 * k4X + b5 * k5X + b6 * k6X + b7 * k7X);
        wY = y + h * (b1 * k1Y + b2 * k2Y + b3 * k3Y + b4 * k4Y + b5 * k5Y + b6 * k6Y + b7 * k7Y);
        double d1 = fabs(wX - locX);
        double d2 = fabs(wY - locY);
        double err = max(d1, d2);
        // h *= min(facmax, max(facmin, fac * pow((tol / err), 1. / 6)));
        h *= pow((tol / err), 1. / 5);
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
        if (!(draw % 10000)) {
            // cout << x - sin(T) << " " << y - cos(T) << endl;
        }
        nextH = StepChoice(T, y, x, currH);
        draw++;
    }
    cout << x - sin(100000 * 3.1415926535897932) << " " << y - cos(100000 * 3.1415926535897932) << endl;
    out.close();
}
