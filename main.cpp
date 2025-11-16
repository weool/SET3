#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <fstream>
#include <algorithm>

using namespace std;

bool isInsideCircle(double cx, double cy, double r, double x, double y) {
    double dx = x - cx;
    double dy = y - cy;
    return dx*dx + dy*dy <= r*r;
}

void monteCarloExperiment(vector<double>& circle1, vector<double>& circle2, vector<double>& circle3, bool narrow, const string& filename) {
    double x_min, x_max, y_min, y_max;

    if (narrow) {
        x_min = max(max(circle1[0] - circle1[2], circle2[0] - circle2[2]), circle3[0] - circle3[2]);
        x_max = min(min(circle1[0] + circle1[2], circle2[0] + circle2[2]), circle3[0] + circle3[2]);
        y_min = max(max(circle1[1] - circle1[2], circle2[1] - circle2[2]), circle3[1] - circle3[2]);
        y_max = min(min(circle1[1] + circle1[2], circle2[1] + circle2[2]), circle3[1] + circle3[2]);
    } else {
        x_min = min(min(circle1[0] - circle1[2], circle2[0] - circle2[2]), circle3[0] - circle3[2]);
        x_max = max(max(circle1[0] + circle1[2], circle2[0] + circle2[2]), circle3[0] + circle3[2]);
        y_min = min(min(circle1[1] - circle1[2], circle2[1] - circle2[2]), circle3[1] - circle3[2]);
        y_max = max(max(circle1[1] + circle1[2], circle2[1] + circle2[2]), circle3[1] + circle3[2]);
    }

    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> dist_x(x_min, x_max);
    uniform_real_distribution<double> dist_y(y_min, y_max);

    ofstream outfile(filename);

    for (int n = 100; n <= 100000; n += 500) {
        int points_inside = 0;

        for (int i = 0; i < n; ++i) {
            double x = dist_x(gen);
            double y = dist_y(gen);

            if (isInsideCircle(circle1[0], circle1[1], circle1[2], x, y) &&
                isInsideCircle(circle2[0], circle2[1], circle2[2], x, y) &&
                isInsideCircle(circle3[0], circle3[1], circle3[2], x, y)) {
                points_inside++;
            }
        }

        double total_area = (x_max - x_min) * (y_max - y_min);
        double area = (double)points_inside / n * total_area;
        outfile << n << " " << area << endl;
    }
    outfile.close();
}

int main() {
    vector<double> circle1 = {1.0, 1.0, 1.0};
    vector<double> circle2 = {1.5, 2.0, sqrt(5)/2};
    vector<double> circle3 = {2.0, 1.5, sqrt(5)/2};

    monteCarloExperiment(circle1, circle2, circle3, true, "data_1.txt");
    monteCarloExperiment(circle1, circle2, circle3, false, "data_2.txt");

    return 0;
}