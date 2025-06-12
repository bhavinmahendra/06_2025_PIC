#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "TAxis.h"
#include "TCanvas.h"
#include "TError.h"
#include "TF1.h"
#include "TGraph.h"
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
const double tolVal = 1e-10;

double conversionToNormal(double rStar) {
    int maxIterations = 1000;
    double r = rStar, rNew, guess1, guess2, rTerm, noneTerm, gFunction, gDerivative;

    rTerm = rs - (2 * rStar) + pow(rs, 2);
    noneTerm = pow(rs, 2) + (rs * rStar) + pow(rStar, 2) + pow(rs, 3);
    guess1 = (rTerm + sqrt(4 * noneTerm))/2;
    guess2 = (rTerm - sqrt(4 * noneTerm))/2;

    if (guess1 > 0)
        r = guess1;
    else if (guess2 > 0)
        r = guess2;

    if (r <= rs)
        r = rs + tolVal;

    for (int i = 0; i < maxIterations; i++) {
        gFunction = r + rs*log(abs(r-rs))-rStar;
        gDerivative = 1 + rs/abs(r-rs);    
        rNew = r - gFunction/gDerivative;
        if (abs(rNew - r) < tolVal)
            return rNew;
        r = rNew;
    }
    cout << "Solution did not converge." << endl;
    return rNew;
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
    int position;
    string line, delimiter = ",";
    
    ifstream file(filename);
    for (int i = 0; i < rows; i++) {
        position = 0;
        getline(file, line);
        for (int j = 0; j < columns; j++) {
            position = line.find(delimiter);
            matrix[i][j] = strtod(line.substr(0, position).c_str(), &end);
            line.erase(0, position + delimiter.length());
        }
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

    for (int i = 0; i < psi.size() - 1; ++i) {
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
    vector<double> e1(columnsQuarter, 0.00), e2(columnsQuarter, 0.00), rStarVal(columnsQuarter, 0.00);
    vector<vector<double>> psiMain(rows, vector<double> (columnsQuarter * 4, 0.00));
    vector<vector<double>> psiHalf(rows, vector<double> (columnsQuarter * 2, 0.00));
    vector<vector<double>> psiQuarter(rows, vector<double> (columnsQuarter, 0.00));

    fileToMatrix(fileMain, psiMain, rows, columnsQuarter * 4);
    fileToMatrix(fileHalf, psiHalf, rows, columnsQuarter * 2);
    fileToMatrix(fileQuarter, psiQuarter, rows, columnsQuarter);

    for (int i = 0; i < columnsQuarter; i++) {
        if (type == "normal") {
            rStarVal[i] = -rStarMax + i * drStarQuarter;
            error1 = 4 * abs(psiMain[position][i*4] - psiHalf[position][i*2]);
            error2 = abs(psiHalf[position][i*2] - psiQuarter[position][i]);
            e1[i] = error1;
            e2[i] = error2;
        }
        else if (type == "log") {
            error1 = 4 * abs(psiMain[position][i*4] - psiHalf[position][i*2]);
            error2 = abs(psiHalf[position][i*2] - psiQuarter[position][i]);
            if (error1 > tolVal && error2 > tolVal) {
                rStarVal[i] = -rStarMax + i * drStarQuarter;
                e1[i] = log(4 * abs(psiMain[position][i*4] - psiHalf[position][i*2]));
                e2[i] = log(abs(psiHalf[position][i*2] - psiQuarter[position][i]));
            }
        }
    }

    string title = "Convergence Test at t = " + to_string(time);
    TCanvas* c = new TCanvas("c", "Canvas", 2500, 1500);
    c->SetGrid();
    TGraph* g1 = new TGraph(columnsQuarter, rStarVal.data(), e1.data());
    g1->SetTitle(title.c_str());
    if (type == "normal") {
        g1->GetXaxis()-> SetTitle("Space");
        g1->GetYaxis()-> SetTitle("Errors");
    }
    else if (type == "log") {
        g1->GetXaxis()-> SetTitle("Space");
        g1->GetYaxis()-> SetTitle("Log |Errors|");
    }
    g1->GetXaxis()-> SetDecimals(kTRUE);
    g1->GetYaxis()-> SetDecimals(kTRUE);
    g1->SetMarkerColor(kRed);
    g1->SetMarkerSize(0.8);
    g1->SetMarkerStyle(8);
    g1->Draw("AP");
    TGraph* g2 = new TGraph(columnsQuarter, rStarVal.data(), e2.data());
    g2->SetMarkerColor(kBlue);
    g2->SetMarkerSize(0.8);
    g2->SetMarkerStyle(8);
    g2->Draw("P SAME");

    string output = "convergenceTest/convergence_test_"  + type + ".png", command = "open " + output;
    c->SaveAs(output.c_str());
    delete g2;
    delete g1;
    delete c;
    system(command.c_str());
}

string dampedExponential(double peak) {
    return "[0]*exp(-(x - " + to_string(peak) + ")/[1])*sin([2]*(x - " + to_string(peak) + ") + [3])";
}

void quasiNormalModes(string filename, string details, string type, int rows, int columns, double rStar, double dt, double drStar) {
    int position = (rStar + rStarMax)/drStar;
    string command, output, title;
    vector<double> timeVector(rows, 0.00), waveVector(rows, 0.00);
    vector<vector<double>> psi(rows, vector<double> (columns, 0.00));
    fileToMatrix(filename, psi, rows, columns);

    title = "Evolution of the Wave Function at r* = " + to_string(rStar);
    output = "quasiNormalModes/file.jpg";
    command = "open " + output;

    if (position < 0 || position >= rows) {
        cout << "Position out of limits. No output." << endl;
        return;
    }

    ofstream file("quasiNormalModes/file_" + details + "_" + type + ".csv");
    for (int i = 0; i < rows; i++) {
        if (type == "normal") {
            timeVector[i] = i * dt;
            waveVector[i] = psi[i][position];
            file << timeVector[i] << "," << waveVector[i] << endl;        
        }
        else if (type == "log") {
            if (i * dt > tolVal) {
                timeVector[i] = log(abs(i * dt));
                waveVector[i] = log(abs(psi[i][position]));
            }
            file << timeVector[i] << "," << waveVector[i] << endl;
        }
    }
    file.close();
}

void quasiNormalModesFit(string filename, string details, string type, int size, double rStar, double tMin, double tMax, double psiMin, double psiMax) {
    gErrorIgnoreLevel = kError;
    double peakTime = 0, peakSpace = 0;
    string command, output, title;
    vector<double> timeVector(size, 0.00), waveVector(size, 0.00);
    vector<vector<double>> psi(size, vector<double> (2, 0.00));
    fileToMatrix(filename, psi, size, 2);

    for (int i = 0; i < size; i++) {
        timeVector[i] = psi[i][0];
        waveVector[i] = psi[i][1];
        if (abs(waveVector[i]) > peakSpace) {
            peakTime = timeVector[i];
            peakSpace = waveVector[i];
        }
    }

    title = "Evolution of the Wave Function at r* = " + to_string(rStar);
    output = "quasiNormalModes/file_" + details  + "_" + type + ".jpg";
    command = "open " + output;

    TCanvas* c = new TCanvas("c", "Canvas", 2500, 1500);
    TF1* f = new TF1("f", dampedExponential(peakTime).c_str(), timeVector[0], timeVector[timeVector.size()-1]);
    TGraph* g = new TGraph(psi.size(), timeVector.data(), waveVector.data());

    c-> SetGrid();
    g-> SetTitle(title.c_str());
    if (type == "normal") {
        g-> GetXaxis()-> SetTitle("Time");
        g-> GetYaxis()-> SetTitle("Wave Function Value");
    }
    else if (type == "log") {
        g-> GetXaxis()-> SetTitle("Log(Time)");
        g-> GetYaxis()-> SetTitle("Log |Wave Function Value|");
    }

    //g-> GetXaxis()-> SetRangeUser(tMin, tMax);
    //g-> GetYaxis()-> SetRangeUser(psiMin, psiMax);

    f-> SetParLimits(0, peakSpace-0.001, peakSpace+0.001);
    f-> SetParLimits(1, 0.00, 500.00);
    f-> SetParLimits(2, 0.00, 1.00);
    f-> SetParLimits(3, -M_PI, M_PI);

    f-> SetParameters(peakSpace, 175.00, (2 * M_PI)/300, 0.00);

    g-> GetYaxis()-> SetDecimals(kTRUE);
    g-> SetLineColor(kBlack);
    g-> SetLineWidth(5);
    g-> Fit(f, "Q", "R", peakTime, 750);
    g-> Draw("AL");
    gStyle-> SetFitFormat(".5g");
    gStyle-> SetOptFit(111);
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

    int l = 4, pointsTime = 5000, pointsSpace = 4000;
    double norm = 100.0, mu = 500.00, sigma = 0.25 * 0.50 * rStarMax;
    double dt, drStar;

    if (getSample) {
        vector<vector<double>> psi(pointsTime, vector<double>(pointsSpace, 0.00));
        string details = "l_" + to_string(l) + "_norm_" + to_string(norm) + "_mu_" + to_string(mu) + "_sigma_" + to_string(sigma);

        drStar = (2 * rStarMax)/(psi[0].size() - 1), dt = cfl * drStar;
        initializeVals(psi, dt, drStar, norm, mu, sigma);
        solvePDE(psi, dt, drStar, l);
        saveData(psi, details, dt, drStar);
    }

    if (getConvergenceTestData) {
        vector<vector<double>> psiQuarter(pointsTime, vector<double>(pointsSpace/4, 0.00));
        vector<vector<double>> psiHalf(pointsTime, vector<double>(pointsSpace/2, 0.00));
        vector<vector<double>> psiMain(pointsTime, vector<double>(pointsSpace, 0.00));
        string details = "";

        drStar = (2 * rStarMax)/(psiMain[0].size()), dt = cfl * drStar;
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
        double drStarQuarter = (2 * rStarMax)/(pointsSpace/4);
        double time = 500.00;
        int rows = pointsTime, columnsQuarter = pointsSpace/4;
        string fileQuarter = "bin/data_5000_1000_.csv", fileHalf = "bin/data_5000_2000_.csv", fileMain = "bin/data_5000_4000_.csv", type = "log";
        convergenceTest(drStarQuarter, fileMain, fileHalf, fileQuarter, type, time, rows, columnsQuarter);
    }

    if (getQuasiNormalModes) {
        double rStar = 50.00, tMin = 6.00, tMax = 10.00, psiMin = -20.00, psiMax = 10.00;
        string details = "l_" + to_string(l) + "_norm_" + to_string(norm) + "_mu_" + to_string(mu) + "_sigma_" + to_string(sigma), type = "normal", filename = "bin/data_5000_4000_" + details + ".csv";
        drStar = (2 * rStarMax)/(pointsSpace - 1), dt = cfl * drStar;
        quasiNormalModes(filename, details, type, pointsTime, pointsSpace, rStar, dt, drStar);
    }

    if (getQuasiNormalModesFit) {
        double rStar = 50.00, tMin = 6.00, tMax = 10.00, psiMin = -20.00, psiMax = 10.00;
        int size = pointsTime;
        string details = "l_" + to_string(l) + "_norm_" + to_string(norm) + "_mu_" + to_string(mu) + "_sigma_" + to_string(sigma), type = "normal", filename = "quasiNormalModes/file_" + details + "_" + type + ".csv";
        quasiNormalModesFit(filename, details, type, size, rStar, tMin, tMax, psiMin, psiMax);
    }

    return 0;
}
