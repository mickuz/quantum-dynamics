#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>

using namespace std;

int N = 101; //number of points
int n = 1; //number of initial state
double w = 3*M_PI*M_PI/2; //resonant angular frequency
double k = 1; //parameter of oscillating field
double dx = 1.0/(N-1); //position step
double dt = 0.0001; //time step

vector<double> positions; //points on [0,1] segment
vector<double> psiR; //real part of wave function
vector<double> psiI; //imaginary part of wave function
vector<double> hamiltonianR; //action of hamiltonian operator on real part of wave function
vector<double> hamiltonianI; //action of hamiltonian operator on imaginary part of wave function

void SetPositions()
{
    for(int i = 0; i < N; i++)
        positions.push_back(i*dx);
}

void SetWaveFunction()
{
    for(auto i : positions)
    {
        psiR.push_back(sqrt(2)*sin(n*i*M_PI));
        psiI.push_back(0);
    }
}

void ComputeHamiltonian(double time)
{
    hamiltonianR.clear();
    hamiltonianI.clear();
    hamiltonianR.push_back(0);
    hamiltonianI.push_back(0);
    for(int i = 1; i < N-1; i++)
    {
        double hR = -0.5*(psiR[i+1]+psiR[i-1]-2*psiR[i])/(dx*dx)+k*(positions[i]-0.5)*psiR[i]*sin(w*time);
        double hI = -0.5*(psiI[i+1]+psiI[i-1]-2*psiI[i])/(dx*dx)+k*(positions[i]-0.5)*psiI[i]*sin(w*time);
        hamiltonianR.push_back(hR);
        hamiltonianI.push_back(hI);
    }
    hamiltonianR.push_back(0);
    hamiltonianI.push_back(0);
}

void UpdateWaveFunction(double time)
{
    for(int i = 0; i < N; i++)
        psiR[i] = psiR[i]+dt/2*hamiltonianI[i];
    ComputeHamiltonian(time+dt/2);
    for(int i = 0; i < N; i++)
        psiI[i] = psiI[i]-dt*hamiltonianR[i];
    ComputeHamiltonian(time+dt);
    for(int i = 0; i < N; i++)
        psiR[i] = psiR[i]+dt/2*hamiltonianI[i];
}

double ComputeNorm()
{
    double norm = 0;
    for(int i = 0; i < N; i++)
        norm += dx*(psiR[i]*psiR[i]+psiI[i]*psiI[i]);
    return norm;
}

double ComputeAvgPosition()
{
    double avgPosition = 0;
    for(int i = 0; i < N; i++)
        avgPosition += dx*positions[i]*(psiR[i]*psiR[i]+psiI[i]*psiI[i]);
    return avgPosition;
}

double ComputeEnergy()
{
    double energy = 0;
    for(int i = 0; i < N; i++)
        energy += dx*(psiR[i]*hamiltonianR[i]+psiI[i]*hamiltonianI[i]);
    return energy;
}

void SaveData(ofstream &outputFile, double time)
{
    outputFile << time << "  " << ComputeNorm() << "  " << ComputeAvgPosition() << "  " << ComputeEnergy() << endl;
}

int main()
{
    double time = 0;
    ofstream outputFile;

    SetPositions();
    SetWaveFunction();
    ComputeHamiltonian(time);

    outputFile.open("avs.dat");
    while(time < 20)
    {
        SaveData(outputFile, time);
        UpdateWaveFunction(time);
        time += dt;
    }
    return 0;
}