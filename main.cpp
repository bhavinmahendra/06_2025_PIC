#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include </opt/homebrew/Cellar/boost/1.88.0/include/boost/math/special_functions/lambert_w.hpp>

#include "TAxis.h"
#include "TCanvas.h"
#include "TError.h"
#include "TF1.h"
#include "TGraph.h"
#include "TLatex.h"
#include "TPaveLabel.h"
#include "TPaveStats.h"
#include "TStyle.h"

using namespace std;

const double c = 1.00;
const double M = 1.00;
const double G = 1.00;
const double rs = (2.00 * G * M)/pow(c, 2);
const double rStarMax = 2000;
const double cfl = 0.50;
const double tolVal = 1e-50;

double conversionToNormal(double rStar) {
    double L1, L2, arg = exp(rStar/rs - 1), W;

    if (rStar < -1000) {
        return rs * (1 + exp(rStar/rs - 1));
    }

    else if (rStar > 1000) {
        return r;
    }
    else
        W = boost::math::lambert_w0(arg);

    return rs + rs * W;
}

double finiteDifferenceForward(double space0, double space1, double space2, double drStar) {
    return c * ((-space2) + (4 * space1) - (3 * space0))/(2 * drStar);
}

double finiteDifferenceBackward(double space0, double space1, double space2, double drStar) {
    return - c * ((space0) - (4 * space1) + (3 * space2))/(2 * drStar);
}

double spatialSecondOrder(double space0, double space1, double space2, double drStar) {
    return pow(c/drStar, 2) * (space0 - (2 * space1) + space2);
}

double fFunction(double r, int l) {
    if (r - rs < tolVal)
        return 0;
    return (1 - (rs/r)) * (((l * (l + 1))/pow(r,2)) + (rs/pow(r, 3)));
}

void saveData(const vector<vector<double>>& psi, string details, double dt, double drStar) {
    double timeStamp = 0;
    int timeSteps = psi.size(), spaceSteps = psi[0].size();
    string title_txt = "bin/data_" + to_string(timeSteps) + "_" + to_string(spaceSteps) + "_" + details + ".txt";
    string title_csv = "bin/data_" + to_string(timeSteps) + "_" + to_string(spaceSteps) + "_" + details + ".csv";
    ofstream file_txt(title_txt);
    ofstream file_csv(title_csv);
    file_txt << scientific << setprecision(15);
    file_csv << scientific << setprecision(15);

    for (int i = 0; i < psi.size(); i++) {
        timeStamp = i * dt;
        file_txt << "\"Time = " << timeStamp << "\n";
        for (int j = 0; j < psi[0].size(); j++) {
            file_txt << -rStarMax + (j * drStar) << " " << psi[i][j] << "\n";
            file_csv << psi[i][j] << ",";
        }
        file_txt << "\n";
        file_csv << "\n";
    }
    file_txt.close();
    file_csv.close();
}

void fileToMatrix(string filename, vector<vector<double>>& matrix, int rows, int columns) {
    char* end = nullptr;
    int position = 0;
    string line, delimiter = ",";
    
    ifstream file(filename);
    while (getline(file, line)) {
        vector<double> row;
        position = 0;
        for (int j = 0; j < columns; j++) {
            position = line.find(delimiter);
            row.push_back(strtod(line.substr(0, position).c_str(), &end));
            line.erase(0, position + delimiter.length());
        }
        matrix.push_back(row);
    }
}

void initializeVals(vector<vector<double>>& psi, double dt, double drStar, double norm, double mu, double sigma) {
    double rStar = 0.00;
    for (int j = 0; j < psi[0].size(); j++) {
        rStar = -rStarMax + (j * drStar);
        psi[0][j] = norm * exp(-0.50 * pow((rStar - mu)/sigma, 2));
        psi[1][j] = psi[0][j];
    }
}

void solvePDE(vector<vector<double>>& psi, double dt, double drStar, int l) {
    double last = psi[0].size() - 1;
    vector<double> rStarVal(psi[0].size(), 0.0), fVals(psi[0].size(), 0.0);

    for (int j = 0; j < psi[0].size(); j++) {
        rStarVal[j] = -rStarMax + j * drStar;
        double r = conversionToNormal(rStarVal[j]);
        fVals[j] = fFunction(r, l);
    }

    vector<vector<double>> aux(psi.size(), vector<double>(psi[0].size(), 0.0));
    vector<double> k1(psi[0].size(), 0.00), k1d(psi[0].size(), 0.00), k2(psi[0].size(), 0.00), k2d(psi[0].size(), 0.00), k3(psi[0].size(), 0.00), k3d(psi[0].size(), 0.00), k4(psi[0].size(), 0.00), k4d(psi[0].size(), 0.00);

    for (int i = 0; i < psi.size() - 1; i++) {
        k1[0]  = dt * aux[i][0];
        k1d[0] = dt * finiteDifferenceForward(aux[i][0], aux[i][1], aux[i][2], drStar);
        k2[0]  = dt * (aux[i][0] + 0.5 * k1d[0]);
        k2d[0] = dt * finiteDifferenceForward(aux[i][0] + 0.5 * k1d[0], aux[i][1] + 0.5 * k1d[1], aux[i][2] + 0.5 * k1d[2], drStar);
        k3[0]  = dt * (aux[i][0] + 0.5 * k2d[0]);
        k3d[0] = dt * finiteDifferenceForward(aux[i][0] + 0.5 * k2d[0], aux[i][1] + 0.5 * k2d[1], aux[i][2] + 0.5 * k2d[2], drStar);
        k4[0]  = dt * (aux[i][0] + k3d[0]);
        k4d[0] = dt * finiteDifferenceForward(aux[i][0] + k3d[0], aux[i][1] + k3d[1], aux[i][2] + k3d[2], drStar);
        
        for (int j = 1; j < psi[0].size() - 1; j++) {
            k1[j]  = dt * aux[i][j];
            k1d[j] = dt * (spatialSecondOrder(psi[i][j-1], psi[i][j], psi[i][j+1], drStar) - fVals[j] * psi[i][j]);
        }

        for (int j = 1; j < psi[0].size() - 1; j++) {
            k2[j]  = dt * (aux[i][j] + 0.5 * k1d[j]);
            k2d[j] = dt * (spatialSecondOrder(psi[i][j-1] + 0.5 * k1[j-1],  psi[i][j] + 0.5 * k1[j], psi[i][j+1] + 0.5 * k1[j+1], drStar) - fVals[j] * (psi[i][j] + 0.5 * k1[j]));
        }

        for (int j = 1; j < psi[0].size() - 1; j++) {
            k3[j]  = dt * (aux[i][j] + 0.5 * k2d[j]);
            k3d[j] = dt * (spatialSecondOrder(psi[i][j-1] + 0.5 * k2[j-1], psi[i][j] + 0.5 * k2[j], psi[i][j+1] + 0.5 * k2[j+1], drStar) - fVals[j] * (psi[i][j] + 0.5 * k2[j]));
        }

        for (int j = 1; j < psi[0].size() - 1; j++) {
            k4[j]  = dt * (aux[i][j] + k3d[j]);
            k4d[j] = dt * (spatialSecondOrder(psi[i][j-1] + k3[j-1],  psi[i][j] + k3[j], (psi[i][j+1] + k3[j+1]), drStar) - fVals[j] * ( psi[i][j] + k3[j]));
        }

        k1[last] = dt * aux[i][last];
        k1d[last] = dt * finiteDifferenceBackward(aux[i][last-2], aux[i][last-1], aux[i][last], drStar);
        k2[last] = dt * (aux[i][last] + 0.5 * k1d[last]);
        k2d[last] = dt * finiteDifferenceBackward(aux[i][last-2] + 0.5 * k1d[last-2], aux[i][last-1] + 0.5 * k1d[last-1], aux[i][last] + 0.5 * k1d[last], drStar);
        k3[last] = dt * (aux[i][last] + 0.5 * k2d[last]);
        k3d[last] = dt * finiteDifferenceBackward(aux[i][last-2] + 0.5 * k2d[last-2], aux[i][last-1] + 0.5 * k2d[last-1], aux[i][last] + 0.5 * k2d[last], drStar);
        k4[last] = dt * (aux[i][last] + k3d[last]);
        k4d[last] = dt * finiteDifferenceBackward(aux[i][last-2] + k3d[last-2], aux[i][last-1] + k3d[last-1], aux[i][last] + k3d[last], drStar);

        for (int j = 0; j < psi[0].size(); j++) {
            psi[i+1][j] = psi[i][j] + (k1[j] + 2*k2[j] + 2*k3[j] + k4[j]) / 6.0;
            aux[i+1][j] = aux[i][j] + (k1d[j] + 2*k2d[j] + 2*k3d[j] + k4d[j]) / 6.0;
        }    
    }
}

void convergenceTest(double drStarQuarter, string fileMain, string fileHalf, string fileQuarter, string type, double time, int rows, int columnsQuarter) {
    gErrorIgnoreLevel = kError;
    int position = time/(cfl * drStarQuarter);
    double error1, error2;
    vector<double> e1, e2, rStarVal;
    vector<vector<double>> psiMain, psiHalf, psiQuarter;

    fileToMatrix(fileMain, psiMain, rows, columnsQuarter * 4);
    fileToMatrix(fileHalf, psiHalf, rows, columnsQuarter * 2);
    fileToMatrix(fileQuarter, psiQuarter, rows, columnsQuarter);

    for (int i = 0; i < psiQuarter[0].size(); i++) {
        if (type == "normal") {
            rStarVal.push_back(-rStarMax + i * drStarQuarter);
            error1 = 4 * abs(psiMain[position][i*4] - psiHalf[position][i*2]);
            error2 = abs(psiHalf[position][i*2] - psiQuarter[position][i]);
            e1.push_back(error1);
            e2.push_back(error2);
        }
        else if (type == "log") {
            error1 = 4 * abs(psiMain[position][i*4] - psiHalf[position][i*2]);
            error2 = abs(psiHalf[position][i*2] - psiQuarter[position][i]);
            if (error1 > tolVal && error2 > tolVal) {
                rStarVal.push_back(-rStarMax + i * drStarQuarter);
                e1.push_back(log(4 * abs(psiMain[position][i*4] - psiHalf[position][i*2])));
                e2.push_back(log(abs(psiHalf[position][i*2] - psiQuarter[position][i])));
            }
        }
    }

    string title = "Convergence Test at t = " + to_string(time);
    TCanvas* c = new TCanvas("c", "Canvas", 2500, 1500);
    c->SetGrid();
    TGraph* g1 = new TGraph(rStarVal.size(), rStarVal.data(), e1.data());
    g1->SetTitle(title.c_str());
    if (type == "normal") {
        g1->GetXaxis()-> SetTitle("Space");
        g1->GetYaxis()-> SetTitle("|Errors|");
    }
    else if (type == "log") {
        g1->GetXaxis()-> SetTitle("Space");
        g1->GetYaxis()-> SetTitle("Log |Errors|");
    }
    g1->GetXaxis()-> SetDecimals(kTRUE);
    g1->GetYaxis()-> SetDecimals(kTRUE);
    g1->SetLineColor(kRed);
    g1->SetLineWidth(5);
    g1->Draw("AL");
    TGraph* g2 = new TGraph(rStarVal.size(), rStarVal.data(), e2.data());
    g2->SetLineColor(kBlue);
    g2->SetLineWidth(3);
    g2->Draw("L SAME");

    string output = "convergenceTest/convergence_test_"  + type + "_" + to_string(time) + ".png", command = "open " + output;
    c->SaveAs(output.c_str());
    delete g2;
    delete g1;
    delete c;
    system(command.c_str());
}

string dampedExponential(double peak, string option) {
    if (option == "mode")
        return "[0] + (x-[4])*[1] + log(abs(sin(([2] * (x-[4])) + [3])))";
    else if (option == "modes")
        return "[0] + (x-[4])*[1] + log(abs(sin(([2] * (x-[4])) + [3]))) + [0] + (x-[4])*[5] + log(abs(sin(([6] * (x-[4])) + [3])))";
    else if (option == "tail")
         return "[0]/(x*x) + [1]";
    else
        return "[0]";
}

void quasiNormalModes(string filename, string details, string type, int rows, int columns, double rStar, double dt, double drStar) {
    int position = (rStar + rStarMax)/drStar;
    string command, output, title;
    vector<double> timeVector, waveVector;
    vector<vector<double>> psi;
    fileToMatrix(filename, psi, rows, columns);

    title = "Evolution of the Wave Function at r* = " + to_string(rStar);
    command = "open " + output;

    if (position < 0 || position >= columns) {
        cout << "Position out of limits. No output." << endl;
        return;
    }

    ofstream file("quasinormalmodes/file_" + to_string(rows) + "_" + to_string(columns) + "_" + details + "_" + type + ".csv");
    file << scientific << setprecision(15);
    for (int i = 0; i < psi.size(); i++) {
        if (type == "normal") {
            timeVector.push_back(i * dt);
            waveVector.push_back(psi[i][position]);
            file << timeVector[i] << "," << waveVector[i] << endl;
             
        }
        else if (type == "log") {
            timeVector.push_back(i * dt);
            waveVector.push_back(log(abs(psi[i][position])));
            file << timeVector[i] << "," << waveVector[i] << endl;
        }
    }
    file.close();
}

void quasiNormalModesFit(string filename, string details, string type, int size, int columns, double rStar, double tMin, double tMax, double psiMin, double psiMax, bool scale, string option) {
    gErrorIgnoreLevel = kError;
    double peakTime = 0, peakSpace = 0;
    string command, output, title;
    vector<double> timeVector, waveVector;
    vector<vector<double>> psi;
    fileToMatrix(filename, psi, size, 2);

    for (int i = 0; i < psi.size(); i++) {
        timeVector.push_back(psi[i][0]);
        waveVector.push_back(psi[i][1]);
        if (waveVector[i] > peakSpace) {
            peakTime = timeVector[i];
            peakSpace = waveVector[i];
        }
    }

    title = "Evolution of the Wave Function at r* = " + to_string(rStar);
    output = "quasinormalmodes/file_" + to_string(size) + "_" + to_string(columns) + "_" + details  + "_" + type + ".jpg";
    command = "open " + output;

    TCanvas* c = new TCanvas("c", "Canvas", 2500, 1500);
    TF1* f = new TF1("f", dampedExponential(peakTime, option).c_str(), timeVector[0], timeVector[timeVector.size()-1]);
    TGraph* g = new TGraph(psi.size(), timeVector.data(), waveVector.data());

    c-> SetGrid();
    g-> SetTitle(title.c_str());
    f-> SetNpx(10000);
    if (type == "normal") {
        g-> GetXaxis()-> SetTitle("Time");
        g-> GetYaxis()-> SetTitle("Wave Function Value");
    }
    else if (type == "log") {
        g-> GetXaxis()-> SetTitle("Time");
        g-> GetYaxis()-> SetTitle("Log |Wave Function Value|");
    }

    if (scale) {
        g-> GetXaxis()-> SetRangeUser(tMin, tMax);
        g-> GetYaxis()-> SetRangeUser(psiMin, psiMax);
    }

    if (option == "mode") {
        f-> SetParLimits(0, peakSpace - 0.001, peakSpace + 0.001);
        f-> SetParLimits(3, -M_PI, M_PI);
        f-> SetParLimits(4, peakTime - 0.001, peakTime + 0.001);
        //f-> SetParLimits(1, -0.009, -0.005);
        f-> SetParLimits(2, 0.060, 0.070);
        //f-> SetParameters(peakSpace, -0.007, 0.065, 0.000, peakTime);
    }

    else if (option == "modes"){
        f-> SetParLimits(0, peakSpace - 0.001, peakSpace + 0.001);
        f-> SetParLimits(3, -M_PI, M_PI);
        f-> SetParLimits(4, peakTime - 0.001, peakTime + 0.001);
        f-> SetParLimits(1, -0.009, -0.005);
        f-> SetParLimits(2, 0.060, 0.070);
        f-> SetParLimits(5, -0.009, -0.005);
        f-> SetParLimits(6, 0.060, 0.070);
        f-> SetParameters(peakSpace, -0.007, 0.065, 0.000, peakTime, -0.007, 0.065);
    }

    else if (option == "tail") {
        f-> SetParLimits(0, 0, 1e10);
        f-> SetParLimits(1, -100, 100);
    }

    else
        cout << "None" << endl;

    g-> GetYaxis()-> SetDecimals(kTRUE);
    g-> SetLineColor(kBlack);
    g-> SetLineWidth(5);
    g-> Fit(f, "", "R", peakTime, 605);
    g-> Draw("AL");
    c-> SaveAs(output.c_str());
    delete f;
    delete g;
    delete c;
    system(command.c_str());
}

int main() {
    bool getSample = false, 
         getConvergenceTestData = false, 
         getConvergenceTest = false, 
         getQuasiNormalModes = false,
         getQuasiNormalModesFit = false;

    int l = 4, pointsTime = 10000, pointsSpace = 32000;
    double norm = 100.0, mu = 100.00, sigma = 25.00;
    double dt, drStar;

    if (getSample) {
        vector<vector<double>> psi(pointsTime, vector<double>(pointsSpace, 0.00));
        string details = "l_" + to_string(l) + "_norm_" + to_string(norm) + "_mu_" + to_string(mu) + "_sigma_" + to_string(sigma);

        drStar = (2 * rStarMax)/(psi[0].size() + 1), dt = cfl * drStar;
        initializeVals(psi, dt, drStar, norm, mu, sigma);
        solvePDE(psi, dt, drStar, l);
        saveData(psi, details, dt, drStar);
    }

    if (getConvergenceTestData) {
        vector<vector<double>> psiQuarter(pointsTime, vector<double>(pointsSpace/4, 0.00));
        vector<vector<double>> psiHalf(pointsTime, vector<double>(pointsSpace/2, 0.00));
        vector<vector<double>> psiMain(pointsTime, vector<double>(pointsSpace, 0.00));
        string details = "convergence";

        drStar = (2 * rStarMax)/(psiMain[0].size() + 1), dt = cfl * drStar;
        initializeVals(psiMain, dt, drStar, norm, mu, sigma);
        solvePDE(psiMain, dt, drStar, l);
        saveData(psiMain, details, dt, drStar);
        
        drStar *= 2;
        initializeVals(psiHalf, dt, drStar, norm, mu, sigma);
        solvePDE(psiHalf, dt, drStar, l);
        saveData(psiHalf, details, dt, drStar);
        
        drStar *= 2;
        initializeVals(psiQuarter, dt, drStar, norm, mu, sigma);
        solvePDE(psiQuarter, dt, drStar, l);
        saveData(psiQuarter, details, dt, drStar);
    }

    if (getConvergenceTest) {
        double drStarQuarter = (2 * rStarMax)/((pointsSpace + 1)/4);
        double time = 500.00;
        int rows = pointsTime, columnsQuarter = pointsSpace/4;
        string fileQuarter = "bin/data_" + to_string(pointsTime) + "_" + to_string(pointsSpace/4) + "_convergence.csv", fileHalf = "bin/data_" + to_string(pointsTime) + "_" + to_string(pointsSpace/2) + "_convergence.csv", fileMain = "bin/data_" + to_string(pointsTime) + "_" + to_string(pointsSpace) + "_convergence.csv", type = "log";
        convergenceTest(drStarQuarter, fileMain, fileHalf, fileQuarter, type, time, rows, columnsQuarter);
    }

    if (getQuasiNormalModes) {
        double rStar = 500.00;
        string details = "l_" + to_string(l) + "_norm_" + to_string(norm) + "_mu_" + to_string(mu) + "_sigma_" + to_string(sigma), type = "log", filename = "bin/data_" + to_string(pointsTime) + "_" + to_string(pointsSpace) + "_" + details + ".csv";
        drStar = (2 * rStarMax)/(pointsSpace + 1), dt = cfl * drStar;
        quasiNormalModes(filename, details, type, pointsTime, pointsSpace, rStar, dt, drStar);
    }

    if (getQuasiNormalModesFit) {
        bool scale = true;
        double rStar = 500.00, tMin = 250.00, tMax = 850.00, psiMin = -15.00, psiMax = 15.00;
        string details = "l_" + to_string(l) + "_norm_" + to_string(norm) + "_mu_" + to_string(mu) + "_sigma_" + to_string(sigma), option = "mode", type = "log", filename = "quasinormalmodes/file_" + to_string(pointsTime) + "_" + to_string(pointsSpace) + "_" + details + "_" + type + ".csv";
        quasiNormalModesFit(filename, details, type, pointsTime, pointsSpace, rStar, tMin, tMax, psiMin, psiMax, scale, option);
    }

    return 0;
}
