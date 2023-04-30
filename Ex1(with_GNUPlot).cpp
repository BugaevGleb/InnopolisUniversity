// Gleb Bugaev, g.bugaev@innopolis.university, CS-02

#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

#define GNUPLOT_NAME "C:\\gnuplot\\bin\\gnuplot -persist"

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


    FILE* pipe = _popen(GNUPLOT_NAME, "w");

    if(pipe != nullptr) {


        fprintf(pipe, "set terminal wxt size 2000,2000\n");

        fprintf(pipe, "set multiplot\n");

        fprintf(pipe, "set size 0.5, 0.5\n");
        fprintf(pipe, "set origin 0, 0.5\n");
        fprintf(pipe, "plot [0 : 50] [0 : 150] '-' using 1:2 title 'v(t)' with lines\n");

        for(int i = 0; i < times.size(); ++i) {
            fprintf(pipe, "%f\t%f\n", times[i], victims[i]);
        }

        fprintf(pipe, "e\n");

        fprintf(pipe, "set size 0.5, 0.5\n");
        fprintf(pipe, "set origin 0.5, 0.5\n");
        fprintf(pipe, "plot [0 : 50] [0 : 150] '-' using 1:2 title 'k(t)' with lines\n");

        for(int i = 0; i < times.size(); ++i) {
            fprintf(pipe, "%f\t%f\n", times[i], killers[i]);
        }

        fprintf(pipe, "e\n");

        fprintf(pipe, "set size 0.5, 0.5\n");
        fprintf(pipe, "set origin 0, 0\n");
        fprintf(pipe, "plot [0 : 150] [0 : 150] '-' using 1:2 title 'v(k)' with lines\n");

        for(int i = 0; i < times.size(); ++i) {
            fprintf(pipe, "%f\t%f\n", killers[i], victims[i]);
        }

        fprintf(pipe, "e\n");

        fprintf(pipe, "set size 0.5, 0.5\n");
        fprintf(pipe, "set origin 0.5, 0\n");
        fprintf(pipe, "plot [0 : 150] [0 : 150] '-' using 1:2 title 'k(v)' with lines\n");

        for(int i = 0; i < times.size(); ++i) {
            fprintf(pipe, "%f\t%f\n", victims[i], killers[i]);
        }

        fprintf(pipe, "e\n");

        fprintf(pipe, "unset multiplot\n");

        fflush(pipe);
        _pclose(pipe);

    }


    return 0;

}
