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
const double cfl = 0.5;
const double tolVal = 1e-5;

string dampedExponential(double peak) {
    return "[0]*exp(-(x-" + to_string(peak) + ")/[1])*sin(([2]*(x-" + to_string(peak) + "))+[3])";
}

string format(int time) {
    if (time < 10)
        return "000" + to_string(time);
    else if (time < 100)
        return "00" + to_string(time);
    else if (time < 1000)
        return "0" + to_string(time);
    else
        return to_string(time);
}

double conversionToNormal(double rStar) {
    int maxIterations = 1000;
    double r = rStar, gFunction, gDerivative, rNew;
    if (r <= rs)
        r = rs + tolVal;
    for (int i = 0; i < maxIterations; i++) {
        gFunction = exp((rStar - r)/rs) + rs - r;
        gDerivative = -(1/rs) * exp((rStar - r)/rs) - 1;    
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

void saveData(const vector<vector<double>>& psi, double dt, double drStar) {
    double timeStamp = 0;
    int timeSteps = psi.size(), spaceSteps = psi[0].size();
    string title_txt = "bin/data_" + to_string(timeSteps) + "_" + to_string(spaceSteps) + ".txt";
    string title_csv = "bin/data_" + to_string(timeSteps) + "_" + to_string(spaceSteps) + ".csv";
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

void saveMP4() {
    system("ffmpeg -y -framerate 30 -i spaceEvolution/file.folder/file%04d.jpg -c:v libx264 -pix_fmt yuv420p spaceEvolution/file.mp4 > /dev/null 2>&1");
    system("open spaceEvolution/file.mp4");
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

void spaceEvolution(string filename, int rows, int columns, double dt, double drStar, double psiMax) {
    gErrorIgnoreLevel = kError;
    double maxVal = 0;
    vector<double> spaceArray(columns, 0.00), waveArray(columns, 0.00);
    vector<double> timeVector(rows, 0.00), waveVector(rows, 0.00);
    vector<vector<double>> psi(rows, vector<double> (columns, 0.00));
    string command, output, title;
    fileToMatrix(filename, psi, rows, columns);

    TCanvas* c = new TCanvas("c", "Canvas", 250, 150);
    c-> SetGrid();

    for (int i = 0; i < rows; i++) {
        title = "Space Evolution of the Wave Function at t = " + to_string(i * dt);
        output = "spaceEvolution/file.folder/file" + format(i) + ".jpg";
        command = "open " + output;

        for (int j = 0; j < columns; j++) {
            spaceArray[j] = -rStarMax + j * drStar;
            waveArray[j] = psi[i][j];
            if (abs(waveArray[j]) > maxVal)
                maxVal = abs(waveArray[j]);      
        }

        TGraph* g = new TGraph(columns, spaceArray.data(), waveArray.data());
        g-> SetTitle(title.c_str());
        g-> GetXaxis()-> SetTitle("Space");
        g-> GetYaxis()-> SetTitle("Wave Function Value");
        g-> GetXaxis()-> SetDecimals(kTRUE);
        g-> GetYaxis()-> SetDecimals(kTRUE);
        if (maxVal < psiMax)
            g-> GetYaxis()-> SetRangeUser(-psiMax, psiMax);
        g-> SetMarkerColor(kBlack);
        g-> SetMarkerSize(0.8);
        g-> SetMarkerStyle(8);
        g-> Draw("APL");
        c-> SaveAs(output.c_str());
        delete g;
    }
    delete c;
}

void timeEvolution(string filename, int rows, int columns, double rStar, double dt, double drStar, string type, double tMin, double tMax, double psiMin, double psiMax) {
    gErrorIgnoreLevel = kError;
    int position = (rStar + rStarMax)/drStar;
    double peakTime = 0, peakSpace = 0;
    string command, output, title;
    vector<double> timeVector(rows, 0.00), waveVector(rows, 0.00);
    vector<vector<double>> psi(rows, vector<double> (columns, 0.00));
    fileToMatrix(filename, psi, rows, columns);

    title = "Time Evolution of the Wave Function at r* = " + to_string(rStar);
    output = "timeEvolution/file.jpg";
    command = "open " + output;

    if (position < 0 || position >= rows) {
        cout << "Position out of limits. No output." << endl;
        return;
    }

    ofstream file("timeEvolution/file.csv");
    for (int i = 0; i < rows; i++) {
        if (type == "log") {
            if (i * dt > tolVal) {
                timeVector[i] = log(i * dt);
                waveVector[i] = log(abs(psi[i][position]));
            }
            file << timeVector[i] << "," << waveVector[i] << endl;
        }
        else if (type == "normal") {
            timeVector[i] = i * dt;
            waveVector[i] = psi[i][position];
            file << timeVector[i] << "," << waveVector[i] << endl;        
        }
        if (abs(waveVector[i]) > peakSpace) {
            peakTime = timeVector[i];
            peakSpace = abs(waveVector[i]);
        }
    }
    file.close();

    TCanvas* c = new TCanvas("c", "Canvas", 2500, 1500);
    TF1* f = new TF1("Damped Exponential", "gaus", timeVector[0], timeVector[timeVector.size()-1]);
    TGraph* g = new TGraph(psi.size(), timeVector.data(), waveVector.data());

    c-> SetGrid();
    f-> SetParLimits(0, 0.00, 1.00);
    f-> SetParLimits(1, 0.00, 1.00);
    f-> SetParLimits(2, 0.00, 1.00);
    f-> SetParLimits(3, 0.00, 1.00);
    g-> SetTitle(title.c_str());
    if (type == "log") {
        g-> GetXaxis()-> SetTitle("Log(Time)");
        g-> GetYaxis()-> SetTitle("Log(Wave Function Value)");
    }
    else if (type == "normal") {
        g-> GetXaxis()-> SetTitle("Time");
        g-> GetYaxis()-> SetTitle("Wave Function Value");
    }
    g-> GetXaxis()-> SetRangeUser(tMin, tMax);
    g-> GetYaxis()-> SetRangeUser(psiMin, psiMax);
    g-> GetYaxis()-> SetDecimals(kTRUE);
    g-> SetMarkerColor(kBlack);
    g-> SetMarkerSize(0.1);
    g-> SetMarkerStyle(8);
    g-> Fit(f, "Q");
    g-> Draw("AL");
    gStyle-> SetFitFormat(".5g");
    gStyle-> SetOptFit(111);
    c-> SaveAs(output.c_str());
    delete f;
    delete g;
    delete c;
    system(command.c_str());
}

void fixedSpaceConvergence(double dtQuarter, string fileMain, string fileHalf, string fileQuarter, double rStar, int rowsQuarter, int columns) {
    gErrorIgnoreLevel = kError;
    int position = (rStarMax + rStar)/(dtQuarter/cfl);
    vector<double> e1(rowsQuarter, 0.00), e2(rowsQuarter, 0.00), timeVals(rowsQuarter, 0.00);
    vector<vector<double>> psiMain(rowsQuarter * 4, vector<double> (columns, 0.00));
    vector<vector<double>> psiHalf(rowsQuarter * 2, vector<double> (columns, 0.00));
    vector<vector<double>> psiQuarter(rowsQuarter, vector<double> (columns, 0.00));

    fileToMatrix(fileMain, psiMain, rowsQuarter * 4, columns);
    fileToMatrix(fileHalf, psiHalf, rowsQuarter * 2, columns);
    fileToMatrix(fileQuarter, psiQuarter, rowsQuarter, columns);

    for (int i = 0; i < rowsQuarter; i++) {
        double error1 = 4 * (psiMain[i*4][position] - psiHalf[i*2][position]);
        double error2 = psiHalf[i*2][position] - psiQuarter[i][position];

        timeVals[i] = i * dtQuarter;
        e1[i] = error1;
        e2[i] = error2;
    }

    string title = "Spatial Convergence Test at r* = " + to_string(rStar);
    TCanvas* c = new TCanvas("c", "Canvas", 2500, 1500);
    c->SetGrid();
    TGraph* g1 = new TGraph(rowsQuarter, timeVals.data(), e1.data());
    g1->SetTitle(title.c_str());
    g1->GetXaxis()-> SetTitle("Time");
    g1->GetYaxis()-> SetTitle("Errors");
    g1->GetXaxis()-> SetDecimals(kTRUE);
    g1->GetYaxis()-> SetDecimals(kTRUE);
    g1->SetMarkerColor(kRed);
    g1->SetMarkerSize(0.8);
    g1->SetMarkerStyle(8);
    g1->Draw("AL");
    TGraph* g2 = new TGraph(rowsQuarter, timeVals.data(), e2.data());
    g2->SetMarkerColor(kBlue);
    g2->SetMarkerSize(0.8);
    g2->SetMarkerStyle(8);
    g2->Draw("P SAME");

    string output = "convergenceTest/spatial_convergence_test.png", command = "open " + output;
    c->SaveAs(output.c_str());
    delete g2;
    delete g1;
    delete c;
    system(command.c_str());
}

void fixedTimeConvergence(double drStarQuarter, string fileMain, string fileHalf, string fileQuarter, double time, int rows, int columnsQuarter) {
    gErrorIgnoreLevel = kError;
    int position = time/(cfl * drStarQuarter);
    vector<double> e1(columnsQuarter, 0.00), e2(columnsQuarter, 0.00), rStarVal(columnsQuarter, 0.00);
    vector<vector<double>> psiMain(rows, vector<double> (columnsQuarter * 4, 0.00));
    vector<vector<double>> psiHalf(rows, vector<double> (columnsQuarter * 2, 0.00));
    vector<vector<double>> psiQuarter(rows, vector<double> (columnsQuarter, 0.00));

    fileToMatrix(fileMain, psiMain, rows, columnsQuarter * 4);
    fileToMatrix(fileHalf, psiHalf, rows, columnsQuarter * 2);
    fileToMatrix(fileQuarter, psiQuarter, rows, columnsQuarter);

    for (int i = 0; i < columnsQuarter; i++) {
        double error1 = (4 * abs(psiMain[position][i*4] - psiHalf[position][i*2]));
        double error2 = (abs(psiHalf[position][i*2] - psiQuarter[position][i]));

        rStarVal[i] = -rStarMax + i * drStarQuarter;
        e1[i] = error1;
        e2[i] = error2;
    }

    string title = "Temporal Convergence Test at t = " + to_string(time);
    TCanvas* c = new TCanvas("c", "Canvas", 2500, 1500);
    c->SetGrid();
    TGraph* g1 = new TGraph(columnsQuarter, rStarVal.data(), e1.data());
    g1->SetTitle(title.c_str());
    g1->GetXaxis()-> SetTitle("Space");
    g1->GetYaxis()-> SetTitle("Errors");
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

    string output = "convergenceTest/temporal_convergence_test.png", command = "open " + output;
    c->SaveAs(output.c_str());
    delete g2;
    delete g1;
    delete c;
    system(command.c_str());
}

int main() {
    bool getFixedSpaceConvergenceData = false, 
         getFixedTimeConvergenceData = false, 
         getFixedSpaceConvergenceTest = false, 
         getFixedTimeConvergenceTest = false, 
         getSample = true, 
         getSpaceEvolution = false, 
         getTimeEvolution = false;

    int l = 2, pointsTime = 5000, pointsSpace = 4000;
    double norm = 100.0, mu = rStarMax * 0.50, sigma = 0.05 * 0.50 * rStarMax;

    double dt, drStar;

    if (getSample) {
        vector<vector<double>> psi(pointsTime, vector<double>(pointsSpace, 0.00));

        drStar = (2 * rStarMax)/(psi[0].size() - 1), dt = cfl * drStar;
        initializeVals(psi, dt, drStar, norm, mu, sigma);
        solvePDE(psi, dt, drStar, l);
        saveData(psi, dt, drStar);
    }

    if (getFixedSpaceConvergenceData) {
        vector<vector<double>> psiQuarter(pointsTime/4, vector<double>(pointsSpace, 0.00));
        vector<vector<double>> psiHalf(pointsTime/2, vector<double>(pointsSpace, 0.00));
        vector<vector<double>> psiMain(pointsTime, vector<double>(pointsSpace, 0.00));

        drStar = (2 * rStarMax)/(psiMain[0].size()), dt = cfl * drStar;
        initializeVals(psiMain, dt, drStar, norm, mu, sigma);
        solvePDE(psiMain, dt, drStar, l);
        saveData(psiMain, dt, drStar);
        
        dt *= 2;
        initializeVals(psiHalf, dt, drStar, norm, mu, sigma);
        solvePDE(psiHalf, dt, drStar, l);
        saveData(psiHalf, dt, drStar);
        
        dt *= 2;
        initializeVals(psiQuarter, dt, drStar, norm, mu, sigma);
        solvePDE(psiQuarter, dt, drStar, l);
        saveData(psiQuarter, dt, drStar);
    }

    if (getFixedTimeConvergenceData) {
        vector<vector<double>> psiQuarter(pointsTime, vector<double>(pointsSpace/4, 0.00));
        vector<vector<double>> psiHalf(pointsTime, vector<double>(pointsSpace/2, 0.00));
        vector<vector<double>> psiMain(pointsTime, vector<double>(pointsSpace, 0.00));

        drStar = (2 * rStarMax)/(psiMain[0].size()), dt = cfl * drStar;
        initializeVals(psiMain, dt, drStar, norm, mu, sigma);
        solvePDE(psiMain, dt, drStar, l);
        saveData(psiMain, dt, drStar);
        
        drStar *= 2;
        initializeVals(psiHalf, dt, drStar, norm, mu, sigma);
        solvePDE(psiHalf, dt, drStar, l);
        saveData(psiHalf, dt, drStar);
        
        drStar *= 2;
        initializeVals(psiQuarter, dt, drStar, norm, mu, sigma);
        solvePDE(psiQuarter, dt, drStar, l);
        saveData(psiQuarter, dt, drStar);
    }

    if (getFixedSpaceConvergenceTest) {
        double dtQuarter = cfl * (2 * rStarMax)/(pointsSpace/4);
        double rStar = 500.00;
        int rowsQuarter = pointsTime/4, columns = 4000;
        string fileQuarter = "bin/data_1250_4000.csv", fileHalf = "bin/data_2500_4000.csv", fileMain = "bin/data_5000_4000.csv";
        fixedSpaceConvergence(dtQuarter, fileMain, fileHalf, fileQuarter, rStar, rowsQuarter, columns);
    }

    if (getFixedTimeConvergenceTest) {
        double drStarQuarter = (2 * rStarMax)/(pointsSpace/4);
        double time = 500.00;
        int rows = pointsTime, columnsQuarter = pointsSpace/4;
        string fileQuarter = "bin/data_5000_1000.csv", fileHalf = "bin/data_5000_2000.csv", fileMain = "bin/data_5000_4000.csv";
        fixedTimeConvergence(drStarQuarter, fileMain, fileHalf, fileQuarter, time, rows, columnsQuarter);
    }

    if (getSpaceEvolution) {
        double psiMax = 50.00;
        string filename = "bin/data_5000_4000.csv";
        drStar = (2 * rStarMax)/(pointsSpace - 1), dt = cfl * drStar;
        spaceEvolution(filename, pointsTime, pointsSpace, dt, drStar, psiMax);
        saveMP4();
    }

    if (getTimeEvolution) {
        double rStar = 50.00, tMin = 6.00, tMax = 10.00, psiMin = -20.00, psiMax = 10.00;
        string type = "log", filename = "bin/data_5000_4000.csv";
        drStar = (2 * rStarMax)/(pointsSpace - 1), dt = cfl * drStar;
        timeEvolution(filename, pointsTime, pointsSpace, rStar, dt, drStar, type, 6, 10, -20, 10);
    }

    return 0;
}
