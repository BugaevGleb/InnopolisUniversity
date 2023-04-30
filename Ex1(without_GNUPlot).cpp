// Gleb Bugaev, g.bugaev@innopolis.university, CS-02

#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace std;

int main() {

    cout << fixed;
    cout << setprecision(2  );

    int v0, k0, n;
    double alpha1, beta1, alpha2, beta2, timeLimit;
    cin >> v0 >> k0;
    cin >> alpha1 >> beta1 >> alpha2 >> beta2;
    cin >> timeLimit >> n;

    vector<double> times;
    double time = 0;
    double delta = timeLimit / n;
    while(time <= timeLimit) {
        times.push_back(time);
        time += delta;
    }

    vector<double> victims;
    vector<double> killers;
    for(double t: times) {
        victims.push_back((v0 - (alpha2 / beta2)) * cos(sqrt(alpha1 * alpha2) * t) -
        (k0 - (alpha1 / beta1)) * ((sqrt(alpha2) * beta1) / (beta2 * sqrt(alpha1))) *
        sin(sqrt(alpha1 * alpha2) * t) + (alpha2 / beta2));
        killers.push_back((v0 - (alpha2 / beta2)) * ((sqrt(alpha1) * beta2) / (beta1 * sqrt(alpha2)))
        * sin(sqrt(alpha1 * alpha2) * t) + (k0 - (alpha1 / beta1)) * cos(sqrt(alpha1 * alpha2) * t) + (alpha1 / beta1));
    }

    cout << "t:" << endl;
    for(double t: times) {
        cout << t << " ";
    }
    cout << endl;

    cout << "v:" << endl;
    for(double v: victims) {
        cout << v << " ";
    }
    cout << endl;

    cout << "k:" << endl;
    for(double k: killers) {
        cout << k << " ";
    }
    cout << endl;



    return 0;

}
